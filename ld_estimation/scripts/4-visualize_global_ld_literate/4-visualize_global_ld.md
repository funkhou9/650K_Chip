# Script: 4-visualize_global_ld

- *Author: Scott Funkhouser*
- *Date: 2017-04-01*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [ld_estimation](../../ld_estimation.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
5. [Visualize](#visualize)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")
```

## Objectives
Read mean LD estimates from script `3b-process_ld_1kb`.
## Install libraries


```r
library(devtools)
library(magrittr)
library(ggplot2)
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
library(tidyr)
```

```
## 
## Attaching package: 'tidyr'
```

```
## The following object is masked from 'package:magrittr':
## 
##     extract
```

## Load data
Load datasets from each breed, where mean LD estimates were generated from
1KB distance windows. Merge data into long format.


```r
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
```

## Visualize


```r
ggplot(ld_means, aes(x = distance, y = LD, color = Breed)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(x = "Marker distance in kb",
         y = expression(Average~r^2)) +
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13))
```

![plot of chunk global_ld](figure/global_ld-1.tiff)

Print the maximum LD for each breed


```r
ld_means %>%
  group_by(Breed) %>%
  summarize(max(LD))
```

```
## # A tibble: 4 × 2
##       Breed `max(LD)`
##       <chr>     <dbl>
## 1     Duroc 0.7481810
## 2 Hampshire 0.7370145
## 3  Landrace 0.6702474
## 4 Yorkshire 0.6679377
```

Obtain differences between each breed, for each window, then print the average
distance


```r
york <-
  ld_means %>%
    filter(Breed == "Yorkshire") %>%
    select(LD)

duroc <-
  ld_means %>%
    filter(Breed == "Duroc") %>%
    select(LD)

abs(mean(unlist(york - duroc)))
```

```
## [1] 0.08390722
```

