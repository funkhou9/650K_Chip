#' # Script: 3-process_ld
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
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#' Read LD calculations for each chromosome for each breed. Provide analysis
#' and visualizations and save all data in a table for easy retreival.

#' ## Install libraries
library(devtools)
library(magrittr)
library(dplyr)
library(purrr)

#' ## Load data
#' Each chromosome's data must be read and combined.
ld_data <- data.frame()
for (i in c(paste0("0", 1:9), 10:18)) {
  dat <- read.table(paste0("../Yorkshire_snp1101_", i, "_run/ld_syn_c", i, ".txt"),
                    header = TRUE,
                    stringsAsFactors = FALSE)
  print(paste0("chr", i))
  ld_data <- rbind(ld_data, dat)
}

#' ## Analysis
#' 100 KB windows need to be defined, from 0 to 5.1 MB. First establish
#' breakpoints. Breakpoints are then used to establish windows.
breakp <- seq(0, 50000, 2000)
groups <- list()
for (i in 1:25) {
    groups[[i]] <- c(breakp[i], breakp[i + 1])
}

#' For each window, find
yorkshire_ld_means <-
  lapply(groups,
         function(x) {
            sub <- ld_data[ld_data$Dist >= x[1] & ld_data$Dist <= x[2], ]
            c("num_snps" = length(sub$r2),
              "mean" = mean(sub$r2),
              "sd" = sd(sub$r2))
         }) %>%
  do.call(rbind, .)
rownames(yorkshire_ld_means) <- lapply(groups, mean) %>%
                                  unlist() %>%
                                  as.character()

#' ## Visualize
plot(yorkshire_ld_means$mean)

#' ## Save data
save(yorkshire_ld_means, file = "../yorkshire_ld_means.RData")
