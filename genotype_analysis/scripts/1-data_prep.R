#' # Script: 1-data-prep
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20161013*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [genotype_analysis](../../genotype_analysis.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#'   - [From Plate one](from-plate-one)
#'   - [From Plate two](from-plate-two)
#' 4. [Analysis](#analysis)
#'   - [Filter samples](#filter-samples)
#'   - [Filter SNPs](#filter-snps)
#'   - [Add breed info](#add-breed-info)
#'   - [Compute and save allele frequencies](#compute-and-save-allele-frequencies)
#' 5. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/genotype_analysis/scripts")

#' ## Objectives
#'
#' 1. Load raw data from Affymetrix containing 650K genotypes for 180 animals
#' of Yorkshire, Landrace, Hampshire, and Duroc breeds.
#' 2. Attach breed information to each sample.
#' 3. Calculate and save allele frequencies for each SNP within each breed
#'

#' ## Install libraries
library(devtools)
library(dplyr)
library(magrittr)
#' load [`snpTools`](https://github.com/funkhou9/snpTools/commit/f3e70629361e023219cea7f5c27f5c56bbd93323)
# devtools::install_github("funkhou9/snpTools")
library(snpTools)
#' Load ['breedTools'](https://github.com/funkhou9/breedTools/commit/fa2d2cf69ef8abd14354aeb80b980b4a6074f28f)
# devtools::install_github("funkhou9/breedTools")
library(breedTools)
library(methods)

#' ## Load data
#'
#' **NOTE:** a decent amount of data must be loaded from each plate, used to:
#'
#' 1. Attach breed information to each sample
#' 2. Keep only PolyHighResolution SNPs
#'
#' ### From plate one
#' Load **genotyping calls** from first 96-well plate
#'
#' > The first 96-well plate contains all Yorkshire animals, however some are from
#' > ISU, which we want to exclude. These samples are identified by their 96-well
#' > position. `check.names = FALSE` is used to preserve original column names,
#' > which contain "-".
#'
plate1_geno <-
    read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                      "cluster/AxiomGT1.calls.txt"),
               header = TRUE,
               stringsAsFactors = FALSE,
               check.names = FALSE)

#' Load **well information**
#'
#' > Contains information to connect the 96-well position with animal ID
#'
plate1_wells <-
   read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                   "PigIA&MIAX_Sample_Table.csv"),
            header = TRUE,
            stringsAsFactors = FALSE)

#' Load **SNP performance dataset**
#'
#' > Contains a list of each probe and how it was classified (PolyHighResolution,
#' > MonoHighResolution, NoMinorHom, CallRate Below Threshold, etc.)
#'
plate1_SNP_performance <-
    read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                      "cluster/filtered/Ps.performance.txt"),
               header = TRUE,
               stringsAsFactors = FALSE)

#' ### From plate two
#' Load **supplemental key** `pt2_IDs_and_breeds.csv` provided by Nancy Raney.
#'
#' > The key contains information necessary to connect *animal ID* (1, 2, 3, 4)
#' > with information about the pig (name, breed, etc.)
#'
plate2_breeds <-
    read.csv("/mnt/research/pigsnp/raw_data/affymetrix_hd/pt2_IDs_and_breeds.csv",
             header = TRUE,
             stringsAsFactors = FALSE)

#' Load **genotyping calls** from the second 96-well plate
#'
#' > Contains Landrace, Hampshire, and Duroc animals identified by the well they
#' > were positioned in in the 96-well plate
#'
plate2_geno <-
    read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                      "cluster_final/AxiomGT1.calls.txt"),
               header = TRUE,
               stringsAsFactors = FALSE,
               check.names = FALSE)

#' Load **well information**
#'
#' > Contains information to connect the 96-well position with animal ID
#'
plate2_wells <-
    read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                    "PigIA&MIAX_Sample_Table.csv"),
             header = TRUE,
             stringsAsFactors = FALSE)

#' Load **PolyHighResolution SNP list**
#'
#' > Contains names of probes that passed all filtering criteria
#'
plate2_SNP_performance <-
 read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                   "cluster_final/SNPolisher/Ps.performance.txt"),
            header = TRUE,
            stringsAsFactors = FALSE)

#' ## Analysis

#' ### Filter samples

#' Several samples from plate 1 failed according to Affy criteria: they did not
#' have a sufficient dish QC value and/or animal-wise calling rate. The animals
#' that failed were all the animals from ISU. No animals failed QC from plate2.

#' From plate 1, ensure that samples from ISU are not included in `plate1_geno`.
#' Samples from ISU were all of the samples that failed QC.
fail_samples <- plate1_wells[plate1_wells$Sample.Status == "Fail", "Best.Array"]
fail_samples <- paste0(fail_samples, ".CEL")
sum(fail_samples %in% colnames(plate1_geno)[-1])

#' Conversely, confirm that all samples that Pass QC are present in `plate1_geno`.
pass_samples <- plate1_wells[plate1_wells$Sample.Status == "Pass", "Best.Array"]
pass_samples <- paste0(pass_samples, ".CEL")
sum(pass_samples %in% colnames(plate1_geno)[-1])

#' ### Filter SNPs

#' Plate1 and Plate2 were processed by Affymetrix seperately. The consequence
#' of this is that PolyHighResolution SNPs (recommended SNPs) were classified
#' for each plate seperately. Plate1 contains all Yorkshire animals, so even
#' if a SNP is not PolyHighResolution on Plate1 (at least two examples of the
#' minor allele), that doesn't mean it isn't interesting, especially if it is
#' PolyHighResolution on Plate2, containing Duroc, Landrace, and Hampshire
#' animals.
#'
#' Find the number of SNPs that are PolyHighResolution in at least one plate,
#' and MonoHighResolution (no minor het or homozygotes), NoMinorHom (no minor
#' homozygotes), or PolyHighResolution in the other plate.

plate_stats <-
  left_join(plate1_SNP_performance,
            plate2_SNP_performance,
            by = "probeset_id") %>%
  filter(ConversionType.x == "PolyHighResolution" | ConversionType.y == "PolyHighResolution") %>%
  filter(ConversionType.x %in% c("MonoHighResolution", "NoMinorHom", "PolyHighResolution")) %>%
  filter(ConversionType.y %in% c("MonoHighResolution", "NoMinorHom", "PolyHighResolution"))

knitr::kable(table(plate_stats$ConversionType.y, plate_stats$ConversionType.x))

#' From both plates, only keep SNPs present on filtered `plate_stats`
plate1_geno <- plate1_geno[plate1_geno$probeset_id %in% plate_stats$probeset_id, ]
dim(plate1_geno)
plate2_geno <- plate2_geno[plate2_geno$probeset_id %in% plate_stats$probeset_id, ]
dim(plate2_geno)

#' Note that for each plate, there are n+1 columns since the first column is used
#' for the probe ID

#' ### Add breed info
#' All animals in `plate1` should be Yorkshire. Conversely, animals in `plate2` are
#' a combination of breeds that need identified. Will use `plate2_breeds` to identify
#' the breed of each animal in `plate2`. Firstly, find which "MSU ID#" correspond with
#' each breed in {Duroc, Hampshire, and Landrace}. Then, find which "Best Array IDs"
#' correspond with each "MSU ID#"

# sapply returns list of length 3 with MSU IDs for each of the 3 breeds.
array_ids <-
    sapply(c("Duroc", "Hampshire", "Landrace"),
           function(x) {
               plate2_breeds$MSU.ID..[plate2_breeds$Breed == x]
           }) %>%
    # lapply returns list of length 3 with Best.Array IDs for each of the 3 breeds.
    lapply(function(x) {
               plate2_wells$Best.Array[plate2_wells$MSU.ID.. %in% x]
           })

#' Isolate Yorkshire "Best Array IDs" from `plate1_wells`
array_ids$Yorkshire <-
    plate1_wells$Best.Array[plate1_wells$Sample.Status == "Pass"]

#' Add Yorkshires from `plate2` using the last 8 elements in `plate2_wells$Best.Array`
array_ids$Yorkshire <- c(array_ids$Yorkshire, plate2_wells$Best.Array[89:96])

#' Combine genotyping datasets into one data.frame `affy_geno`, which will have
#' probe IDs as colnames, animal names as rownames;
#' likewise SNPs in columns, and animals in rows.
#'
#' Transpose both genotyping matrices
plate1_geno <- t(plate1_geno)
plate2_geno <- t(plate2_geno)

#' Add colnames
colnames(plate1_geno) <- unname(plate1_geno[1, ])
colnames(plate2_geno) <- unname(plate2_geno[1, ])

#' Remove the first row of each
plate1_geno <- plate1_geno[-1, ]
plate2_geno <- plate2_geno[-1, ]

#' Merge genotyping datasets with `snpTools::merge_geno`
affy_geno <- snpTools::merge_geno(plate1_geno, plate2_geno)

#' In order to agree with previous code and implementation, `affy_geno` will require
#' numeric coding of genotypes, with missing genotypes represented as `NA` rather
#' than `-1`
storage.mode(affy_geno) <- "numeric"
affy_geno[affy_geno == -1] <- NA

#' How many SNPs contain at least one missing genotype?
sum(apply(affy_geno, 2, function(x) any(is.na(x))))

#' Remove from both `affy_geno` and `array_ids` the two Yorkshire that also have
#' genome sequence (not relevent for this study) and a Landrace animal that
#' was determined to be Yorkshire by PCA analysis and breed composition estimation

#' Yorkshire animals, identified by IDs 348 and 350
offending_animals <- plate2_wells$Best.Array[plate2_wells$MSU.ID.. %in% c(348, 350)]

#' Add suspicious Landrace animal, determined to be "a550588-4269754-110716-107_F06"
offending_animals <- c(offending_animals, "a550588-4269754-110716-107_F06")

#' Remove from `affy_geno`
affy_geno <- affy_geno[!rownames(affy_geno) %in% paste0(offending_animals, ".CEL"), ]

#' Remove from `array_ids`
array_ids$Landrace <-
  array_ids$Landrace[!array_ids$Landrace %in% offending_animals[3]]
array_ids$Yorkshire <-
  array_ids$Yorkshire[!array_ids$Yorkshire %in% offending_animals[1:2]]

#' After removing three animals, remove SNPs with MAF < 0.05
afs <- colMeans(affy_geno, na.rm = TRUE) / 2
affy_geno <- affy_geno[, !(afs > 0.95 | afs < 0.05)]

#' ### Compute and save allele frequencies
freq <-
  breedTools::allele_freq(affy_geno,
                          lapply(array_ids, function(x) paste0(x, ".CEL"))) %>%
  as.data.frame()
freq <- cbind("probeset_id" = rownames(freq), freq)

#' ## Save data
save(affy_geno, array_ids, file = "../1-data_prep.RData")
write.table(freq,
            file = "../allele_frequencies.txt",
            row.names = FALSE,
            quote = FALSE)
