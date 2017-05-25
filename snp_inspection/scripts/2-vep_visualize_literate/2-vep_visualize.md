# Script: 2-vep_visualize.R

- *Author: Scott Funkhouser*
- *Date: 2017-05-25*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [snp_inspection](../../snp_inspection.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/snp_inspection/scripts")
```

## Objectives
Load VEP results (which were generated with the VEP online application using
an input file prepped from `1-ggp_affy_beadchip_comparison.R`) and plot
the distribution of various consequences using a bar or pie chart.

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

```r
library(readr)
```

## Load data
Load output file from VEP


```r
vep <- read_delim("/mnt/research/pigsnp/raw_data/affymetrix_hd/VEP_ENS_annot.txt",
                  col_names = TRUE,
                  delim = "\t")
```

```
## Parsed with column specification:
## cols(
##   .default = col_character()
## )
```

```
## See spec(...) for full column specifications.
```

## Analysis
Specific data so that only consequences of interest have labels, and remaining
are categorized as "other". First define keywords of interest, then
obtain consequences containing keywords. If consequences don't contain keywords
they are replaced as "other"


```r
keywords <- c("intron_variant",
              "intergenic_variant",
              "upstream_gene_variant",
              "downstream_gene_variant",
              "synonymous_variant",
              "missense_variant",
              "3_prime_UTR_variant",
              "5_prime_UTR_variant")

vep$Consequence <-
  case_when(
    grepl(keywords[1], vep$Consequence) ~ keywords[1],
    grepl(keywords[2], vep$Consequence) ~ keywords[2],
    grepl(keywords[3], vep$Consequence) ~ keywords[3],
    grepl(keywords[4], vep$Consequence) ~ keywords[4],
    grepl(keywords[5], vep$Consequence) ~ keywords[5],
    grepl(keywords[6], vep$Consequence) ~ keywords[6],
    grepl(keywords[7], vep$Consequence) ~ keywords[7],
    grepl(keywords[8], vep$Consequence) ~ keywords[8]
  )

vep$Consequence[is.na(vep$Consequence)] <- "other"

consec <- vep %>%
            group_by(Consequence) %>%
            summarize(value = length(Consequence))
```

## Visualize


```r
ggplot(consec, aes(x = Consequence, y = value)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = Consequence,
                y = value,
                label = value,
                color = "red",
                fontface = "bold"),
            nudge_y = 8000,
            angle = 270) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
```

![plot of chunk consequences](figure/consequences-1.tiff)

