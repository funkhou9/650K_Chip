#' # Script: 9-visualize_tabulate_correlations
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-05-11*
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
#' Process results generated in both `7-process_pofphase.R` and
#' `8-process_pofphase_windows.R`. Results are used to plot pairwise correlations
#' between breeds using 1KB windows from 0 to 50KB, and to tabulate correlations
#' using specified windows.
#'

#' ## Install libraries
library(devtools)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)

#' ## Load data
#' Load datasets generated from both `7-process_pofphase.R` and
#' `8-process_pofphase_windows.R`
load("../Duroc_Hampshire_1kb.RData")
load("../Duroc_Landrace_1kb.RData")
load("../Duroc_Yorkshire_1kb.RData")
load("../Landrace_Yorkshire_1kb.RData")
load("../Hampshire_Landrace_1kb.RData")
load("../Hampshire_Yorkshire_1kb.RData")
load("../Duroc_Hampshire_windows.RData")
load("../Duroc_Landrace_windows.RData")
load("../Duroc_Yorkshire_windows.RData")
load("../Landrace_Yorkshire_windows.RData")
load("../Hampshire_Landrace_windows.RData")
load("../Hampshire_Yorkshire_windows.RData")

#' ## Analysis
#' Assemble data from `7-process_pofphase.R` into a format suitable for plotting
#' mean r squared value against distance in KB.
ld_means <- cbind("distance" = as.numeric(rownames(Duroc_Hampshire_1kb)),
                  "Duroc_Hampshire" = Duroc_Hampshire_1kb[, "corr"],
                  "Duroc_Landrace" = Duroc_Landrace_1kb[, "corr"],
                  "Duroc_Yorkshire" = Duroc_Yorkshire_1kb[, "corr"],
                  "Hampshire_Landrace" = Hampshire_Landrace_1kb[, "corr"],
                  "Hampshire_Yorkshire" = Hampshire_Yorkshire_1kb[, "corr"],
                  "Landrace_Yorkshire" = Landrace_Yorkshire_1kb[, "corr"]) %>%
             as.data.frame(stringsAsFactors = FALSE) %>%
             gather(Breed, LD, Duroc_Hampshire:Landrace_Yorkshire)
ld_means$distance <- ld_means$distance / 1000

#+ global_corr, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
ggplot(ld_means, aes(x = distance, y = LD, color = Breed)) +
    geom_line(size = 1, alpha = 0.7) +
    labs(x = "Marker distance in kb",
         y = "Correlation of Phase") +
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13))

#' For data from `8-process_pofphase_windows.R`, tabulate results from each
#' window.
mean_ld_windows <-
  rbind(round(Duroc_Hampshire_windows[, c("prop_opp_sign", "corr")], 3),
        round(Duroc_Landrace_windows[, c("prop_opp_sign", "corr")], 3),
        round(Duroc_Yorkshire_windows[, c("prop_opp_sign", "corr")], 3),
        round(Hampshire_Landrace_windows[, c("prop_opp_sign", "corr")], 3),
        round(Hampshire_Yorkshire_windows[, c("prop_opp_sign", "corr")], 3),
        round(Landrace_Yorkshire_windows[, c("prop_opp_sign", "corr")], 3))

mean_ld_windows <- cbind("Breeds_Compared" = c(rep("Duroc-Hampshire", 3),
                                               rep("Duroc-Landrace", 3),
                                               rep("Duroc-Yorkshire", 3),
                                               rep("Hampshire-Landrace", 3),
                                               rep("Hampshire_Yorkshire", 3),
                                               rep("Landrace-Yorkshire", 3)),
                         "Distance" = rownames(mean_ld_windows), mean_ld_windows)

rownames(mean_ld_windows) <- NULL
mean_ld_windows <- as.data.frame(mean_ld_windows, stringsAsFactors = FALSE)

knitr::kable(mean_ld_windows)
