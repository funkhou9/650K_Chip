#' # Script: 2-process_haplotypes
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20161102*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#'
#' 1. Load haplotypes generated from `1-phase_genotypes.R`
#'

#' ## Install libraries
library(devtools)
library(magrittr)

#' Install [snpTools](https://github.com/funkhou9/snpTools/commit/f3e70629361e023219cea7f5c27f5c56bbd93323)
# install_github("funkhou9/snpTools")
library(snpTools)

#' ## Load data
#' Read maker names in the order that they are presented in the haplotype files.
#' This can be done by reading the first line of the haplotype file from any breed.
#' All breeds should have the same marker order.
con <- file("../Duroc_fimpute_run/hap_library.txt", open = "r")
snps <- scan(con, what = "character", nlines = 1)

#' Continue to read haplotypes for each breed. Initially stored as an
#' n x 1 data.frame with a single character column.
duroc <- read.table("../Duroc_fimpute_run/hap_library.txt",
                    skip = 1,
                    colClasses = "character")

landrace <- read.table("../Landrace_fimpute_run/hap_library.txt",
                       skip = 1,
                       colClasses = "character")




#' ## Analysis

#' ## Visualize

#' ## Save data
