# Script: 3-ld_calc

- *Author: Scott Funkhouser*
- *Date: 20161107*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [ld_estimation](../../ld_estimation.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")
```

## Objectives

1. For each chromosome, calculate all pairwise measurements of distance and
save to disk, one binary file per chromosome.]


## Install libraries


```r
library(devtools)
library(magrittr)
```

Install ogbox, an assortment of R tools from open source developer Ogan Mancarci,
containing `latexImg()`, a function to parse a line of latex and return a link
to an image renderable in Github markdown.


```r
# install_github("oganm/ogbox")
library(ogbox)
```

```
## Loading required package: RCurl
```

```
## Loading required package: methods
```

```
## Loading required package: bitops
```

```
## Loading required package: stringr
```

```
## Loading required package: parallel
```

```
## Loading required package: ggplot2
```

```
## Loading required package: tools
```

## Load data
Load haplotypes phased in `1-phase_genotypes` and processed in `2-process_haplotypes`


```r
load("../2-process_haplotypes.RData")
```

Load Affy map containing physical posistions for each SNP


```r
affy_map <-
    read.csv("/mnt/research/pigsnp/raw_data/affymetrix_hd/Axiom_PigHD_v1_Annotation.r1.1.csv",
			 header = TRUE,
			 comment.char = "#")
```

## Analysis
The amount of LD between SNP i and SNP j are calculated using the following:

[1] "![](http://latex.codecogs.com/gif.latex?r%5E2_%7Bij%7D%20%3D%20%5Cfrac%7B(p_%7Bij%7D%20-%20p_ip_j)%5E2%7D%7Bp_i(1%20-%20p_i)p_j(1%20-%20p_j)%7D)"

Correlation statistics will be obtained (mean and variance) for each set of
SNP pairs that fall between breakpoints defined in `breakp`


```r
breakp <- seq(0, 10000000, 100000)
```

Define a list of chromosomes...


```r
chr <- list()
```

Define a list of groups...


```r
groups <- list()
for (i in 1:100) {
    groups[[i]] <- c(breakp[i], breakp[i + 1])
}
```

For each chromosome:

1. calcualte a pairwise distances for all pairs of SNPs on of that chromosome
2. sort pairs into elements of `groups` (ex. groups[1] will correspond to
SNP pairs with a physical distance between 0 and 100kb)
3. append `groups` to `chr`



```r
for (i in c(seq(1:18), "X", "Y")) {
    # identify positions on chromosome i, store as a vector
    pos <- affy_map$Physical.Position[affy_map$Chromosome == i]
    # calculate pairwise distances
    distance <- dist(pos)
    g <- lapply(groups, function(x) {
            which(distance >= x[1] & distance <= x[2])
         })
    chr[[i]] <- g
}
```

```
## Error in which(distance >= x[1] & distance <= x[2]): long vectors not supported yet: ../../src/include/Rinlinedfuns.h:137
```

```r
length(chr)
```

```
## [1] 0
```

```r
length(chr[[1]])
```

```
## Error in chr[[1]]: subscript out of bounds
```
