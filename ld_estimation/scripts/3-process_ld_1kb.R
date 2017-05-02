#!/usr/bin/env Rscript
#PBS -N 3-process_ld_1kb
#PBS -l nodes=1:ppn=1,walltime=04:00:00,mem=100Gb
#PBS -j oe
#PBS -o ../output/3b-process_ld_1kb
#PBS -m abe
#PBS -t 1-4
#PBS -M funkhou9@msu.edu

#' # Script: 3-process_ld_1kb
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-03-30*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 6. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#' Read LD calculations for each chromosome for each breed. Establish windows
#' of 1kb from 0 to 50Kb with which to count the number of SNP pairs in each
#' window and calculate mean and sd.

#' ## Install libraries
library(devtools)
library(magrittr)
library(dplyr)
library(purrr)

#' ## Load data
#' Map ith job to the ith breed
breeds <- c("Yorkshire", "Landrace", "Hampshire", "Duroc")
job_num <- Sys.getenv("PBS_ARRAYID")
ith_breed <- breeds[as.numeric(job_num)]

#' Each chromosome's data must be read and combined.
ld_data <- data.frame()
for (i in c(paste0("0", 1:9), 10:18)) {
  dat <- read.table(paste0("../", ith_breed, "_snp1101_", i, "_run/ld_syn_c", i, ".txt"),
                    header = TRUE,
                    stringsAsFactors = FALSE)
  print(paste0("chr", i))
  ld_data <- rbind(ld_data, dat)
}

#' ## Analysis
#' 1 KB windows need to be defined, from 0 to 50KB. First establish
#' breakpoints. Breakpoints are then used to establish windows.
breakp <- seq(0, 50000, 1000)
groups <- list()
for (i in 1:50) {
    groups[[i]] <- c(breakp[i], breakp[i + 1])
}

#' For each window, find the number of SNP pairs in that window,
#' calculate mean and sd of r^2. Rename rownames to be the midpoint
#' of each window.
assign(paste0(ith_breed, "_1kb"),
       lapply(groups,
              function(x) {
                 sub <- ld_data[ld_data$Dist >= x[1] & ld_data$Dist <= x[2], ]
                 c("num_snps" = length(sub$r2),
                   "mean" = mean(sub$r2),
                   "sd" = sd(sub$r2))
              }) %>%
       do.call(rbind, .) %>%
       `rownames<-`(lapply(groups, mean) %>%
                      unlist() %>%
                      as.character()))

#' ## Save data
save(list = paste0(ith_breed, "_1kb"),
     file = paste0("../", ith_breed, "_1kb", ".RData"))
