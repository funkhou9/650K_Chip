#' # Script: 4-visualize_global_ld
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-04-01*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 5. [Visualize](#visualize)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#' Read mean LD estimates from script `3b-process_ld_1kb`.

#' ## Install libraries
library(devtools)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)

#' ## Load data
#' Load datasets from each breed, where mean LD estimates were generated from
#' 1KB distance windows. Merge data into long format.
load("../Yorkshire_1kb.RData")
load("../Duroc_1kb.RData")
load("../Landrace_1kb.RData")
load("../Hampshire_1kb.RData")

ld_means <- cbind("distance" = as.numeric(rownames(Yorkshire_1kb)),
                  "Yorkshire" = Yorkshire_1kb[, "mean"],
                  "Landrace" = Landrace_1kb[, "mean"],
                  "Duroc" = Duroc_1kb[, "mean"],
                  "Hampshire" = Hampshire_1kb[, "mean"]) %>%
             as.data.frame(stringsAsFactors = FALSE) %>%
             gather(Breed, LD, Yorkshire:Hampshire)
ld_means$distance <- ld_means$distance / 1000

#' ## Visualize
#+ global_ld, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
ggplot(ld_means, aes(x = distance, y = LD, color = Breed)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = "Marker distance in kb",
         y = expression(Average~r^2)) +
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13))

#' Print the maximum LD for each breed
ld_means %>%
  group_by(Breed) %>%
  summarize(max(LD))

#' Obtain differences between each breed, for each window, then print the average
#' distance
york <-
  ld_means %>%
    filter(Breed == "Yorkshire") %>%
    select(LD)

duroc <-
  ld_means %>%
    filter(Breed == "Duroc") %>%
    select(LD)

abs(mean(unlist(york - duroc)))
