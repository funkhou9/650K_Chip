#' # Script: 3-compare_calls
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20161025*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [genotype_analysis](../../genotype_analysis.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#'   - [Filter data](#filter-data)
#'   - [Regression](#regression)
#'   - [Estimate GWBC of animal](#estimate-gwbc-of-animal)
#' 5. [Visualize](#visualize)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/genotype_analysis/scripts")

#' ## Objectives
#'
#' 1. From Yorkshire animals genotyped with the 650K chip, isolate those that
#' were previously genotyped on the Illumina 60K beadchip. Previous 60K data
#' may be loaded from the `SF_PG_Industry` project
#' 2. For each SNP genotyped on both Illumina 60K and Affymetrix 650K platforms,
#' test genotype call consistancy using a linear regression.
#' 3. As seen in `2-PCA.R`, one Landrace animal resembles Yorkshire
#' much more than Landrace. Use common SNPs between Illumina and
#' Affy platforms to estimate the genome-wide breed composition (GWBC)
#' of the animal (*see `SF_PG_Industry`*)
#'

#' ## Install libraries
library(devtools)
library(magrittr)

#' Install [snpTools](https://github.com/funkhou9/snpTools/commit/6603afd1db77fb6a93ece38b4a3eeafc7fbc92f2)
# install_github("funkhou9/snpTools")
library(snpTools)

#' Install [breedTools](https://github.com/funkhou9/breedTools/commit/00b77d774e31b69f885b3ccbe413f7caf92abbbb)
# install_git("/mnt/research/pigsnp/raw_data/SF_working_dir/breed_compos/breedTools")
library(breedTools)

#' ## Load data
#' Load prepped Affymetrix genotyping data from `1-data_prep`
load("../1-data_prep.RData")

#' Load Yorkshire sire genotypes from `SF_PG_Industry`. About 90 animals from
#' this set of ~900 should have also been used for Affymetrix 650K genotyping.
load(paste0("/mnt/research/pigsnp/raw_data/SF_working_dir",
            "/Yorkshire_dataset_cleaning/yorkshireDataForModel.RData"))

#' Load previous work from this project that assembled all physical positions
#' for Affy and Illumina chips
load("/mnt/research/pigsnp/NSR/650K_Chip/snp_inspection/pos_list.RData")

#' Load Sample table from Yorkshire genotyped on plate1 of the
#' Affymetrix assay. The "Sample Name" in this dataset will correspond to
#' registration numbers
plate1_samples <-
   read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                   "PigIA&MIAX_Sample_Table.csv"),
            header = TRUE,
            stringsAsFactors = FALSE)

#' Load Sample table from plate2. Those six of eight Yorkshire on plate2 will
#' have registration numbers in the "MSU.ID.." column.
plate2_samples <-
    read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                    "PigIA&MIAX_Sample_Table.csv"),
             header = TRUE,
             stringsAsFactors = FALSE)

#' Load reference panel of allele frequencies required to compute
#' GWBC.
data("GWBC_ref_B")

#' ## Analysis
#' ### Filter data
#' Obtain registration numbers and "affy array IDs" for all Yorkshire animals genotyped
#' on the Affy platforms.
hd_yorks_reg <- c(plate1_samples$Sample.Name, plate2_samples$MSU.ID..[91:96])
hd_yorks_id <- c(plate1_samples$Best.Array, plate2_samples$Best.Array[91:96])

#' Not all numbers in `hd_yorks_reg` are necessarily in `yorkshireGenoDose`. Find which
#' `hd_yorks_reg` and corresponding `hd_yorks_id` are present in `yorkshireGenoDose`.
idx <- hd_yorks_reg %in% rownames(yorkshireGenoDose)

#' Subset both illumina and affy datasets with animals present in both
geno_60K <- yorkshireGenoDose[hd_yorks_reg[idx], ]
geno_650K <- affy_geno[paste0(hd_yorks_id[idx], ".CEL"), ]

#' Find physical positions that match between Affy and Illumina
#' platforms. Note that there are duplicated positions in the Illumina
#' map (mostly '0:0' positions and some non-'0:0' positions) that
#' have differing marker names. We'll exclude these physical positions.
comm_pos <- intersect(pos_list$Affy650, pos_list$SNP60[!duplicated(pos_list$SNP60)])

#' Find which markers (from each platform) correspond to each common
#' physical position
illum_markers <- marker_list$SNP60[match(comm_pos, pos_list$SNP60)]
affy_markers <- marker_list$Affy650[match(comm_pos, pos_list$Affy650)]

#' Identify markers that are present in both genotyping datasets
idx <- illum_markers %in% colnames(yorkshireGenoDose) &
        affy_markers %in% colnames(affy_geno)

#' Subset SNPs according to common markers present in both datasets.
illum_calls <- geno_60K[, illum_markers[idx]]
affy_calls <- geno_650K[, affy_markers[idx]]

#' Identify markers with at least a 90% call rate in both datasets
idx2 <- apply(illum_calls, 2,
              function(x) {
                  sum(!is.na(x)) / length(x) >= 0.9
              })
idx3 <- apply(affy_calls, 2,
              function(x) {
                  sum(!is.na(x)) / length(x) >= 0.9
              })

#' Identify markers that are unfixed in both datasets
idx4 <- apply(illum_calls, 2,
              function(x) {
                  length(unique(x)) > 1
              })
idx5 <- apply(affy_calls, 2,
              function(x) {
                  length(unique(x)) > 1
              })

#' Keep only SNPs that have sufficient call rate in both datasets and
#' that are unfixed in both datasets.
illum_calls <- illum_calls[, idx2 & idx3 & idx4 & idx5]
affy_calls <- affy_calls[, idx2 & idx3 & idx4 & idx5]

#' ### Regression
#' Apply regression for the ith column of each datasets. When using
#' `mapply()`, the `...` passed must be data.frames if `mapply()` is
#' to iterate over columns.
fits <- mapply(function(y, x) lm(y ~ x),
               as.data.frame(illum_calls),
               as.data.frame(affy_calls))

#' Tabulate marker name, "x coefficient", and standard error
coef_se_r2 <- lapply(fits,
                 function(x) {
                     c(summary(x)$coefficients[2, 1:2], summary(x)$r.squared)
                 })
results <- do.call(rbind, coef_se_r2)

#' ### Estimate GWBC of animal
#' From `2-PCA.R`, I know that animal `a550588-4269754-110716-107_F06.CEL`
#' has clusters much closer with Yorkshire genotyped on the Affy chip
#'
#' Isoloate this "suspect Landrace" animal and compute its GWBC
#' based on the reference panel that has been previously developed
#' for the National Swine Registry and used in `SF_PG_Industry`.
sus_landrace <- affy_geno["a550588-4269754-110716-107_F06.CEL", ]

#' Keep only SNPs present on the Illumina platform and convert to
#' Illumina names
sus_landrace <- sus_landrace[affy_markers[idx]]
names(sus_landrace) <- illum_markers[idx]

#' Keep only SNPs that are used in the reference panel to estimate
#' GWBC
sus_landrace <- sus_landrace[names(sus_landrace) %in% rownames(GWBC_ref_B)]

#' Keep only SNPs that have models in `fits`
sus_landrace <- sus_landrace[names(sus_landrace) %in% names(fits)]

#' Transform each SNP according to model fit provided in `fits`.
trans_geno <- c()
for (i in 1:length(sus_landrace)) {
    # Obtain ith marker name
    marker <- names(sus_landrace)[i]
    # Convert ith marker using model in fits[[marker]]
    trans_geno[i] <-
        fits[[marker]]$coefficients[1] + (fits[[marker]]$coefficients[2] * sus_landrace[i])
}
names(trans_geno) <- names(sus_landrace)


#' Use `breedTools` to estimate GWBC of suspect landrace animal using transformed
#' genotypes
breedTools:::QPsolve(trans_geno, GWBC_ref_B)

#' ## Visualize
#' Joint distribution of coefficient estimates and R squared values
plot(results[, 3],
     results[, 1],
     xlab = "R squared",
     ylab = "Coefficient estimate")
