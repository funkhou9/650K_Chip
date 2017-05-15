#!/usr/bin/env Rscript
#PBS -N 8-process_pofphase_windows
#PBS -l nodes=1:ppn=1,walltime=04:00:00,mem=200Gb
#PBS -j oe
#PBS -o ../output/8-process_pofphase_windows
#PBS -m abe
#PBS -t 1-6
#PBS -M funkhou9@msu.edu

#' # Script: 8-process_pofphase_windows.R
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-05-04*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#' For each pairwise combination of breeds - estimate persistance of phase
#' using a few binning strategies. Windows in this script will match those
#' used in Badke et al to construct Table 3. In such a table, SNP pairs assuming
#' specific distances were used to compute correlation of phase and proportion
#' of phase agreement between breed pairs.
#'

#' ## Install libraries
library(devtools)
library(magrittr)
library(dplyr)

#' ## Load data
#' Obtain list of breed pairs and current job number
pairs <- list(c("Duroc", "Hampshire"),
              c("Duroc", "Landrace"),
              c("Duroc", "Yorkshire"),
              c("Hampshire", "Landrace"),
              c("Hampshire", "Yorkshire"),
              c("Landrace", "Yorkshire"))

job_num <- as.numeric(Sys.getenv("PBS_ARRAYID"))
ith_pair <- pairs[[job_num]]

#' For the current pairwise breed combination, each chromosome's data must be
#' read and combined for both breeds.
ld_data_1 <- data.frame()
ld_data_2 <- data.frame()
for (i in c(paste0("0", 1:9), 10:18)) {
  dat1 <- read.table(paste0("../", ith_pair[1], "_snp1101_", i, "_run/ld_syn_c", i, ".txt"),
                    header = TRUE,
                    stringsAsFactors = FALSE)
  print(paste0("chr", i, ith_pair[1]))
  ld_data_1 <- rbind(ld_data_1, dat1)
  dat2 <- read.table(paste0("../", ith_pair[2], "_snp1101_", i, "_run/ld_syn_c", i, ".txt"),
                     header = TRUE,
                     stringsAsFactors = FALSE)
  print(paste0("chr", i, ith_pair[2]))
  ld_data_2 <- rbind(ld_data_2, dat2)
}

#' ## Analysis
#' Specific windows of various lengths are defined.
groups <- list(c(10000, 50000),
               c(50000, 100000),
               c(900000, 1000000))

#' For each window...
#'
#' 1. Find all SNP pairs in that window.
#' 2. Ensure that SNP pairs are the same between the ith breed pairs.
#' 3. Add sign to each SNP pair r2 value
#' 4. Compute correlation between breed pairs
#' 5. Compute proportion of SNP pairs with opposite sign between breeds
#'

obj <- paste0(paste(ith_pair, collapse = "_"), "_windows")
assign(obj,
       lapply(groups,
              function(x) {
                 sub1 <- ld_data_1[ld_data_1$Dist >= x[1] & ld_data_1$Dist <= x[2], ]
                 sub2 <- ld_data_2[ld_data_2$Dist >= x[1] & ld_data_2$Dist <= x[2], ]
                 rownames(sub1) <- paste0(sub1$X..Chr,
                                          "_",
                                          sub1$SNP1,
                                          "_",
                                          sub1$SNP2,
                                          "_",
                                          sub1$Dist)
                 rownames(sub2) <- paste0(sub2$X..Chr,
                                          "_",
                                          sub2$SNP1,
                                          "_",
                                          sub2$SNP2,
                                          "_",
                                          sub2$Dist)
                 pairs <- intersect(rownames(sub1), rownames(sub2))
                 sub1 <- sub1[pairs, ]
                 sub2 <- sub2[pairs, ]
                 stopifnot(all(sub1[, 1:4] == sub2[, 1:4]))
                 sub1$Sign <- case_when(
                    sub1$Sign == "+" ~ 1,
                    sub1$Sign == "-" ~ -1
                 )
                 sub2$Sign <- case_when(
                    sub2$Sign == "+" ~ 1,
                    sub2$Sign == "-" ~ -1
                 )
                 sub1$r2 <- sub1$Sign * sub1$r2
                 sub2$r2 <- sub2$Sign * sub2$r2
                 oppo <- sum(sub1$Sign != sub2$Sign) / length(sub1$Sign)
                 corr <- cor(sub1$r2, sub2$r2)
                 c("num_snps" = nrow(sub1),
                   "prop_opp_sign" = oppo,
                   "corr" = corr)
              }) %>%
       do.call(rbind, .) %>%
       `rownames<-`(lapply(groups, mean) %>%
                      unlist() %>%
                      as.character()))

#' ## Save data
save(list = obj,
     file = paste0("../", obj, ".RData"))
