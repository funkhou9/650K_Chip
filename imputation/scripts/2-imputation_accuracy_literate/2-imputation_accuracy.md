# Script: 2-imputation_accuracy.R

- *Author: Scott Funkhouser*
- *Date: 2017-05-16*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [imputation](../../imputation.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/imputation/scripts")
```

## Objectives
Read imputed genotypes generated in `1-run_fimpute` and calculate correlations
(both SNP-wise and animal-wise) between imputed and measured genotypes to
assess imputation accuracy.

## Install libraries


```r
library(devtools)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
```

## Load data
Load SNP info file generated from FImpute. Contains the order of SNPs in
FImpute output.


```r
snp_info <- read.table("../output/impute_1/fimpute_run/stat_snp.txt",
                       header = TRUE,
                       skip = 1,
                       stringsAsFactors = FALSE)
```

Load SNPs and physical positions from both Affy and Geneseek platforms


```r
load("../../snp_inspection/pos_list.RData")
```

Load imputed data one at a time. Retain only imputed genotype and convert to
matrix. Append imputed genotypes together to create one imputed dataset.


```r
imp_geno <- matrix(0, nrow = 90, ncol = 516813)
for (file in 1:90) {
  geno <- scan(paste0("../output/impute_",
                      file,
                      "/fimpute_run/genotypes_imp.txt"),
               what = character(),
               nlines = 1,
               skip = file)
  g <- unlist(strsplit(geno[3], split = "")) %>%
         as.numeric()
  imp_geno[file, ] <- g
}
```

Load input genotypes prior to masking


```r
load("../../genotype_analysis/1-data_prep.RData")
```

Load and prep map in the same manner that was used previously


```r
map <- read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/",
                       "Axiom_PigHD_v1_Annotation.r1.1.csv"),
                header = TRUE,
                comment.char = "#",
                stringsAsFactors = FALSE)

affy_map <- data.frame("marker" = map$Probe.Set.ID,
                       "chr" = map$Chromosome,
                       "pos" = map$Physical.Position,
                       row.names = "marker",
                       stringsAsFactors = FALSE)

affy_map <- affy_map[affy_map$chr %in% c(as.character(1:18), "X", "Y"), ]
```

## Analysis
Keep only raw genotypes that were used as input to FImpute. When FImpute
runs, it re-sorts the order of SNPs according to physical position raw genotypes
will likewise need to be re-sorted to match to the order of the imputed
genotypes.


```r
raw_geno <- affy_geno[paste0(array_ids$Yorkshire, ".CEL"), ]
raw_geno <- raw_geno[, colnames(raw_geno) %in% rownames(affy_map)]
raw_geno <- raw_geno[, snp_info$SNPID]
```

Obtain SNPs shared by both platforms and SNPs unique to the Affymetrix


```r
matched_pos <- intersect(pos_list$Affy650, pos_list[["GGP-HD"]])
matched_markers <- marker_list$Affy650[pos_list$Affy650 %in% matched_pos]
affy_markers <- marker_list$Affy650[!pos_list$Affy650 %in% matched_pos]
```

Keep only `affy_markers` from both imputed and raw genotyping datasets


```r
idx <- colnames(raw_geno) %in% affy_markers
raw_geno_affy <- raw_geno[, idx]
imp_geno_affy <- imp_geno[, idx]
```

Convert imputed data to 0, 1, 2 allele dosage codes.


```r
imp_geno_affy[imp_geno_affy == 3] <- 1
imp_geno_affy[imp_geno_affy == 4] <- 1
```

Compute animal-wise correlations for all 90 animals. Then snp-wise correlations
for all 479282 SNPs


```r
anim_corr <- data.frame("animal" = character(),
                        "R2" = numeric(),
                        stringsAsFactors = FALSE)
for (i in 1:nrow(raw_geno_affy)) {
  anim_corr[i, ] <- list("animal" = rownames(raw_geno_affy)[i],
                         "R2" = cor(imp_geno_affy[i, ],
                                    raw_geno_affy[i, ],
                                    use = "complete.obs")^2)
}

snp_accuracy <- c()
```

```r
for (i in 1:ncol(raw_geno_affy)) {
  snp_accuracy[i] <- cor(imp_geno_affy[, i],
                         raw_geno_affy[, i],
                         use = "complete.obs")^2
}
```

For plotting, append color codes for each SNP. Colors should alternate with
alternating chromosomes. Provide both chromosomal position and scaled position
for multiple ways of showing results.


```r
accuracies <-
  tibble("position" = pos_list$Affy650[match(colnames(raw_geno_affy), marker_list$Affy650)],
         "R2" = snp_accuracy) %>%
    separate(col = position,
             into = c("chr", "pos"),
             sep = ":",
             convert = TRUE) %>%
    filter(chr != "X")

accuracies$chr <- as.numeric(accuracies$chr)
```

```
## Warning: NAs introduced by coercion
```

```r
accuracies <- arrange(accuracies, chr) %>%
                mutate(color = chr %% 2 == 0) %>%
                na.omit()

accuracies <- accuracies %>%
                group_by(chr) %>%
                mutate(scaled_pos = seq_along(pos) / length(pos))
```

## Visualize
Plot SNP-wise accuracies two different ways and provide table of
animal-wise accuracies


```r
ggplot(accuracies) +
  geom_point(aes(x = seq_along(R2),
                 y = R2,
                 color = color),
             size = 1.2,
             alpha = 0.6) +
  geom_hline(aes(yintercept = mean(R2, na.rm = TRUE))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(breaks = c(-0.25, 0.0, 0.25, 0.5, 0.75, 1.0)) +
  xlab("Position") +
  ylab(expression(Imputation~accuracy~R^2))
```

![plot of chunk snp_impute](figure/snp_impute-1.tiff)

```r
ggplot(accuracies, aes(x = scaled_pos, y = R2)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(color = "red") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(breaks = c(-0.25, 0.0, 0.25, 0.5, 0.75, 1.0)) +
  xlab("Scaled position") +
  ylab(expression(Imputation~accuracy~R^2))
```

```
## `geom_smooth()` using method = 'gam'
```

![plot of chunk snp_impute_scaled](figure/snp_impute_scaled-1.tiff)

```r
knitr::kable(anim_corr)
```



|animal                             |        R2|
|:----------------------------------|---------:|
|a550588-4269754-110716-105_A01.CEL | 0.9289151|
|a550588-4269754-110716-105_A04.CEL | 0.9298846|
|a550588-4269754-110716-105_A05.CEL | 0.9621273|
|a550588-4269754-110716-105_A07.CEL | 0.9297618|
|a550588-4269754-110716-105_A09.CEL | 0.9448253|
|a550588-4269754-110716-105_A10.CEL | 0.9412726|
|a550588-4269754-110716-105_A11.CEL | 0.9369074|
|a550588-4269754-110716-105_B01.CEL | 0.9506981|
|a550588-4269754-110716-105_B02.CEL | 0.9333830|
|a550588-4269754-110716-105_B03.CEL | 0.9449916|
|a550588-4269754-110716-105_B06.CEL | 0.9219489|
|a550588-4269754-110716-105_B07.CEL | 0.9402744|
|a550588-4269754-110716-105_B08.CEL | 0.9532266|
|a550588-4269754-110716-105_B10.CEL | 0.9270891|
|a550588-4269754-110716-105_B11.CEL | 0.9296514|
|a550588-4269754-110716-105_B12.CEL | 0.9431122|
|a550588-4269754-110716-105_C01.CEL | 0.9460133|
|a550588-4269754-110716-105_C03.CEL | 0.9448240|
|a550588-4269754-110716-105_C04.CEL | 0.9124763|
|a550588-4269754-110716-105_C05.CEL | 0.9310128|
|a550588-4269754-110716-105_C06.CEL | 0.9468440|
|a550588-4269754-110716-105_C09.CEL | 0.9394991|
|a550588-4269754-110716-105_C10.CEL | 0.9353776|
|a550588-4269754-110716-105_C11.CEL | 0.9515034|
|a550588-4269754-110716-105_C12.CEL | 0.9322890|
|a550588-4269754-110716-105_D02.CEL | 0.9331778|
|a550588-4269754-110716-105_D03.CEL | 0.9596197|
|a550588-4269754-110716-105_D04.CEL | 0.9444521|
|a550588-4269754-110716-105_D05.CEL | 0.8973778|
|a550588-4269754-110716-105_D06.CEL | 0.8938496|
|a550588-4269754-110716-105_D07.CEL | 0.8979647|
|a550588-4269754-110716-105_D08.CEL | 0.9115497|
|a550588-4269754-110716-105_D09.CEL | 0.9020450|
|a550588-4269754-110716-105_D10.CEL | 0.9171192|
|a550588-4269754-110716-105_D11.CEL | 0.9121404|
|a550588-4269754-110716-105_D12.CEL | 0.9048225|
|a550588-4269754-110716-105_E01.CEL | 0.9038128|
|a550588-4269754-110716-105_E02.CEL | 0.9288131|
|a550588-4269754-110716-105_E03.CEL | 0.9220417|
|a550588-4269754-110716-105_E04.CEL | 0.9136065|
|a550588-4269754-110716-105_E05.CEL | 0.9384916|
|a550588-4269754-110716-105_E06.CEL | 0.9259868|
|a550588-4269754-110716-105_E07.CEL | 0.9367155|
|a550588-4269754-110716-105_E08.CEL | 0.9046034|
|a550588-4269754-110716-105_E09.CEL | 0.8617778|
|a550588-4269754-110716-105_E10.CEL | 0.9194756|
|a550588-4269754-110716-105_E11.CEL | 0.8982805|
|a550588-4269754-110716-105_E12.CEL | 0.9335333|
|a550588-4269754-110716-105_F01.CEL | 0.8937393|
|a550588-4269754-110716-105_F02.CEL | 0.9161912|
|a550588-4269754-110716-105_F03.CEL | 0.9167762|
|a550588-4269754-110716-105_F04.CEL | 0.9259619|
|a550588-4269754-110716-105_F05.CEL | 0.8977660|
|a550588-4269754-110716-105_F06.CEL | 0.8877201|
|a550588-4269754-110716-105_F07.CEL | 0.8541539|
|a550588-4269754-110716-105_F08.CEL | 0.8914719|
|a550588-4269754-110716-105_F09.CEL | 0.9088454|
|a550588-4269754-110716-105_F10.CEL | 0.9280576|
|a550588-4269754-110716-105_F11.CEL | 0.8917019|
|a550588-4269754-110716-105_F12.CEL | 0.9113141|
|a550588-4269754-110716-105_G01.CEL | 0.9121279|
|a550588-4269754-110716-105_G02.CEL | 0.9154951|
|a550588-4269754-110716-105_G03.CEL | 0.8837447|
|a550588-4269754-110716-105_G04.CEL | 0.8954372|
|a550588-4269754-110716-105_G05.CEL | 0.9283311|
|a550588-4269754-110716-105_G06.CEL | 0.8947773|
|a550588-4269754-110716-105_G07.CEL | 0.9003670|
|a550588-4269754-110716-105_G08.CEL | 0.9317155|
|a550588-4269754-110716-105_G09.CEL | 0.8869043|
|a550588-4269754-110716-105_G10.CEL | 0.7392466|
|a550588-4269754-110716-105_G11.CEL | 0.9222145|
|a550588-4269754-110716-105_G12.CEL | 0.9071299|
|a550588-4269754-110716-105_H01.CEL | 0.9156858|
|a550588-4269754-110716-105_H02.CEL | 0.9218949|
|a550588-4269754-110716-105_H03.CEL | 0.8936217|
|a550588-4269754-110716-105_H04.CEL | 0.9313570|
|a550588-4269754-110716-105_H05.CEL | 0.9017028|
|a550588-4269754-110716-105_H06.CEL | 0.9521678|
|a550588-4269754-110716-105_H07.CEL | 0.9322312|
|a550588-4269754-110716-105_H08.CEL | 0.9213097|
|a550588-4269754-110716-105_H09.CEL | 0.9283607|
|a550588-4269754-110716-105_H10.CEL | 0.8911468|
|a550588-4269754-110716-105_H11.CEL | 0.9069065|
|a550588-4269754-110716-105_H12.CEL | 0.9301648|
|a550588-4269754-110716-107_H07.CEL | 0.9076942|
|a550588-4269754-110716-107_H08.CEL | 0.9320503|
|a550588-4269754-110716-107_H09.CEL | 0.8955587|
|a550588-4269754-110716-107_H10.CEL | 0.9364907|
|a550588-4269754-110716-107_H11.CEL | 0.9349745|
|a550588-4269754-110716-107_H12.CEL | 0.9394109|

## Save data


```r
save(anim_corr, accuracies, file = "../output/2-imputation_accuracy.RData")
```

