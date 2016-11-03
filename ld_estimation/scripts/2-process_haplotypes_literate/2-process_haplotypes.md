# Script: 2-process_haplotypes

- *Author: Scott Funkhouser*
- *Date: 20161102*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [ld_estimation](../../ld_estimation.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
6. [Save data](#save-data)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")
```

## Objectives

1. Load haplotypes generated from `1-phase_genotypes.R` and save to disk
as an animal x snp data.frame

## Install libraries


```r
library(devtools)
library(magrittr)
```

## Load data
Read maker names in the order that they are presented in the haplotype files.
This can be done by reading the first line of the haplotype file from any breed.
All breeds should have the same marker order.


```r
con <- file("../Duroc_fimpute_run/hap_library.txt", open = "r")
snps <- scan(con, what = "character", nlines = 1)
```

Continue to read haplotypes for each breed. Initially stored as an
n x 1 data.frame with a single character column.


```r
duroc <- read.table("../Duroc_fimpute_run/hap_library.txt",
                    skip = 1,
                    colClasses = "character")

landrace <- read.table("../Landrace_fimpute_run/hap_library.txt",
                       skip = 1,
                       colClasses = "character")

hampshire <- read.table("../Hampshire_fimpute_run/hap_library.txt",
                        skip = 1,
                        colClasses = "character")

yorkshire <- read.table("../Yorkshire_fimpute_run/hap_library.txt",
                        skip = 1,
                        colClasses = "character")
```

## Analysis
For each row in each of the four haplotype data.frames, convert the haplotypes
from a single character string to a vector using `strsplit()`. Applying
strsplit across each haplotype will result in a list of lists. Remove the inner
list layer with `lapply(., unlist)` then assemble into a data.frame with
`do.call('rbind', .)`


```r
duroc_hap <- lapply(apply(duroc, 1, strsplit, split = ""), unlist) %>%
                do.call('rbind', .)
colnames(duroc_hap) <- snps
storage.mode(duroc_hap) <- "numeric"
```

Repeat process for remaining 3 breeds


```r
landrace_hap <- lapply(apply(landrace, 1, strsplit, split = ""), unlist) %>%
                    do.call('rbind', .)
colnames(landrace_hap) <- snps
storage.mode(landrace_hap) <- "numeric"

hampshire_hap <- lapply(apply(hampshire, 1, strsplit, split = ""), unlist) %>%
                    do.call('rbind', .)
colnames(hampshire_hap) <- snps
storage.mode(hampshire_hap) <- "numeric"

yorkshire_hap <- lapply(apply(yorkshire, 1, strsplit, split = ""), unlist) %>%
                    do.call('rbind', .)
colnames(yorkshire_hap) <- snps
storage.mode(yorkshire_hap) <- "numeric"
```

## Save data


```r
save(duroc_hap,
     landrace_hap,
     hampshire_hap,
     yorkshire_hap,
     file = "../2-process_haplotypes.RData")
```

