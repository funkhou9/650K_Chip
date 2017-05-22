#' # Script: 2-imputation_accuracy.R
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-05-16*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [imputation](../../imputation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/imputation/scripts")

#' ## Objectives
#' Read imputed genotypes generated in `1-run_fimpute` and calculate correlations
#' (both SNP-wise and animal-wise) between imputed and measured genotypes to
#' assess imputation accuracy.
#'

#' ## Install libraries
#+ message=FALSE
library(devtools)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

#' ## Load data
#' Load SNP info file generated from FImpute. Contains the order of SNPs in
#' FImpute output.
snp_info <- read.table("../output/impute_1/fimpute_run/stat_snp.txt",
                       header = TRUE,
                       skip = 1,
                       stringsAsFactors = FALSE)

#' Load SNPs and physical positions from both Affy and Geneseek platforms
load("../../snp_inspection/pos_list.RData")

#' Load imputed data one at a time. Retain only imputed genotype and convert to
#' matrix. Append imputed genotypes together to create one imputed dataset.
imp_geno <- matrix(0, nrow = 90, ncol = 516813)
for (file in 1:90) {
  geno <- scan(paste0("../output/impute_",
                      file,
                      "/fimpute_run/genotypes_imp.txt"),
               what = character(),
               nlines = 1,
               skip = file)
  g <- unlist(strsplit(geno[3], split = "")) %>%
         as.numeric()
  imp_geno[file, ] <- g
}

#' Load input genotypes prior to masking
load("../../genotype_analysis/1-data_prep.RData")

#' Load and prep map in the same manner that was used previously
map <- read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/",
                       "Axiom_PigHD_v1_Annotation.r1.1.csv"),
                header = TRUE,
                comment.char = "#",
                stringsAsFactors = FALSE)

affy_map <- data.frame("marker" = map$Probe.Set.ID,
                       "chr" = map$Chromosome,
                       "pos" = map$Physical.Position,
                       row.names = "marker",
                       stringsAsFactors = FALSE)

affy_map <- affy_map[affy_map$chr %in% c(as.character(1:18), "X", "Y"), ]

#' ## Analysis
#' Keep only raw genotypes that were used as input to FImpute. When FImpute
#' runs, it re-sorts the order of SNPs according to physical position raw genotypes
#' will likewise need to be re-sorted to match to the order of the imputed
#' genotypes.
raw_geno <- affy_geno[paste0(array_ids$Yorkshire, ".CEL"), ]
raw_geno <- raw_geno[, colnames(raw_geno) %in% rownames(affy_map)]
raw_geno <- raw_geno[, snp_info$SNPID]

#' Obtain SNPs shared by both platforms and SNPs unique to the Affymetrix
matched_pos <- intersect(pos_list$Affy650, pos_list[["GGP-HD"]])
matched_markers <- marker_list$Affy650[pos_list$Affy650 %in% matched_pos]
affy_markers <- marker_list$Affy650[!pos_list$Affy650 %in% matched_pos]

#' Keep only `affy_markers` from both imputed and raw genotyping datasets
idx <- colnames(raw_geno) %in% affy_markers
raw_geno_affy <- raw_geno[, idx]
imp_geno_affy <- imp_geno[, idx]

#' Convert imputed data to 0, 1, 2 allele dosage codes.
imp_geno_affy[imp_geno_affy == 3] <- 1
imp_geno_affy[imp_geno_affy == 4] <- 1

#' Compute animal-wise correlations for all 90 animals. Then snp-wise correlations
#' for all 479282 SNPs
anim_corr <- data.frame("animal" = character(),
                        "correlation" = numeric(),
                        stringsAsFactors = FALSE)
for (i in 1:nrow(raw_geno_affy)) {
  anim_corr[i, ] <- list("animal" = rownames(raw_geno_affy)[i],
                         "correlation" = cor(imp_geno_affy[i, ],
                                             raw_geno_affy[i, ],
                                             use = "complete.obs"))
}

snp_corr <- c()
#+ warning=FALSE
for (i in 1:ncol(raw_geno_affy)) {
  snp_corr[i] <- cor(imp_geno_affy[, i],
                     raw_geno_affy[, i],
                     use = "complete.obs")
}

#' For plotting, append color codes for each SNP. Colors should alternate with
#' alternating chromosomes.
correlations <-
  tibble("position" = pos_list$Affy650[match(colnames(raw_geno_affy), marker_list$Affy650)],
         "corr" = snp_corr) %>%
    separate(col = position,
             into = c("chr", "pos"),
             sep = ":",
             convert = TRUE) %>%
    filter(chr != "X")

correlations$chr <- as.numeric(correlations$chr)

correlations <- arrange(correlations, chr) %>%
                  mutate(color = chr %% 2 == 0) %>%
                  na.omit()

#' ## Visualize
#' Plot SNP-wise correlations and provide table of animal-wise correlations
#+ snp_impute, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
ggplot(correlations, aes(x = seq_along(corr),
                         y = corr,
                         color = color)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_hline(aes(yintercept = mean(corr, na.rm = TRUE))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(breaks = c(-0.25, 0.0, 0.25, 0.5, 0.75, 1.0)) +
  xlab("position") +
  ylab("correlation")

knitr::kable(anim_corr)

#' ## Save data
save(anim_corr, correlations, file = "../output/2-imputation_accuracy.RData")
