# Script: 3-process_ld

- *Author: Scott Funkhouser*
- *Date: 2017-03-30*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [ld_estimation](../../ld_estimation.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")
```

## Objectives
Read LD calculations for each chromosome for each breed. Provide analysis
and visualizations and save all data in a table for easy retreival.
## Install libraries


```r
library(devtools)
library(magrittr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(purrr)
```

```
## 
## Attaching package: 'purrr'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     contains, order_by
```

```
## The following object is masked from 'package:magrittr':
## 
##     set_names
```

## Load data
Each chromosome's data must be read and combined.


```r
ld_data <- data.frame()
for (i in c(paste0("0", 1:9), 10:18)) {
  dat <- read.table(paste0("../Yorkshire_snp1101_", i, "_run/ld_syn_c", i, ".txt"),
                    header = TRUE,
                    stringsAsFactors = FALSE)
  print(paste0("chr", i))
  ld_data <- rbind(ld_data, dat)
}
```

```
## [1] "chr01"
## [1] "chr02"
## [1] "chr03"
## [1] "chr04"
## [1] "chr05"
## [1] "chr06"
## [1] "chr07"
## [1] "chr08"
## [1] "chr09"
## [1] "chr10"
## [1] "chr11"
## [1] "chr12"
## [1] "chr13"
## [1] "chr14"
## [1] "chr15"
## [1] "chr16"
## [1] "chr17"
## [1] "chr18"
```

## Analysis
100 KB windows need to be defined, from 0 to 5.1 MB. First establish
breakpoints. Breakpoints are then used to establish windows.


```r
breakp <- seq(0, 50000, 2000)
groups <- list()
for (i in 1:25) {
    groups[[i]] <- c(breakp[i], breakp[i + 1])
}
```

For each window, find


```r
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
```

## Visualize


```r
plot(yorkshire_ld_means$mean)
```

```
## Error in yorkshire_ld_means$mean: $ operator is invalid for atomic vectors
```

## Save data


```r
save(yorkshire_ld_means, file = "../yorkshire_ld_means.RData")
```

