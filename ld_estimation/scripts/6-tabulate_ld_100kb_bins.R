#' # Script: 6-tabulate_ld_100kb_bins
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-05-03*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#'

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#'
#' Retrieve mean ld calculations generated from `5-process_ld_100kb` and tabulate
#' results.
#'

#' ## Install libraries
library(devtools)
library(magrittr)
library(tidyr)
library(dplyr)

#' ## Load data
#' Load results from each breed. Results consist of mean LD measurements for
#' all pairwise SNPs that fall in each of the three categories (pairs of SNPs
#' that are .499Mb and .501Mb apart, .999Mb and 1.001Mb apart, and 4.999Mb and
#' 5.001Mb apart)
load("../Yorkshire_100kb.RData")
load("../Duroc_100kb.RData")
load("../Landrace_100kb.RData")
load("../Hampshire_100kb.RData")

#' ## Analysis
#' Assemble matrix with colnames and print table to screen
mean_ld_100kb <- rbind(round(Duroc_100kb[, "mean"], 2),
                       round(Hampshire_100kb[, "mean"], 2),
                       round(Landrace_100kb[, "mean"], 2),
                       round(Yorkshire_100kb[, "mean"], 2))
colnames(mean_ld_100kb) <- c("0.5Mb", "1Mb", "5Mb")
mean_ld_100kb <- cbind(c("Duroc", "Hampshire", "Landrace", "Yorkshire"),
                       mean_ld_100kb)
knitr::kable(mean_ld_100kb)
