#' # Script: 2-PCA
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20161019*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [genotype_analysis](../../genotype_analysis.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/genotype_analysis/scripts")

#' ## Objectives
#'
#' 1. Using genotypes prepared in `1-data_prep.R`, visually inspect the variation
#' observed among genotyping data using eigenvector decomposition of the **G**
#' matrix (genetic relationship matrix, calculated among all 180 animals)
#'

#' ## Install libraries
library(devtools)
library(magrittr)
library(ggplot2)

#' ## Load data
#' Load intermediate data produced from `1-data_prep.R`. Data contains two objects:
#' a matrix of genotypes for all 180 animals, and a list containing animal IDs
#' organized by breed.
load("../1-data_prep.RData")

#' ## Analysis
#' Calculation of the **G** matrix will require a complete dataset with no
#' missing genotype calls. To account for this, I will fill missing genotypes
#' with the mean value for that SNP
#'
#' Genotype call rates should already be high for each SNP
call_rates <- apply(affy_geno, 2, function(x) {
                sum(!is.na(x)) / length(x)
              })
min(call_rates)

#' There should already be no fixed SNPs
fix_snps <- apply(affy_geno, 2, function(x) {
                length(unique(x)) == 1
            })
sum(fix_snps)

#' Fill in missing genotype calls with the mean for that SNP
X <- apply(affy_geno, 2, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            return(x)
     })

#' Compute **G**
X <- scale(X, center = TRUE, scale = TRUE)
G <- X %*% t(X) / ncol(X)

#' Eigendecomposition of **G**
evd <- eigen(G)
pcs <- as.data.frame(evd$vectors)
colnames(pcs) <- paste0(rep("PC", ncol(pcs)), seq(1:ncol(pcs)))

#' Using `array_ids`, provide breed information for each row of `pcs` (or each
#' row of `affy_geno`)
breeds <- c()
# For each rowname in affy_geno
for (i in rownames(affy_geno)) {
    # For each breed in array_ids
    for (j in 1:length(array_ids)) {
        # Does the ith rowname match one of the IDs for the jth breed?
        if (i %in% paste0(array_ids[[j]], ".CEL"))
            breeds[i] <- names(array_ids)[j]
    }
}
# Remove names, which are automatically generated from the for loop and append
# breeds to eigenvectors
breeds <- unname(breeds)
pcs <- cbind(breeds, pcs)

#' ## Visualize
#+ dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
ggplot(pcs, aes(x = PC1, y = PC2, color = breeds)) +
    geom_point(size = 3, alpha = 0.7)
