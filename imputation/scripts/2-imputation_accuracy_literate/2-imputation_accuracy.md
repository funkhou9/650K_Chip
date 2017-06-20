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
  xlab("position") +
  ylab("imputation accuracy")
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
  xlab("scaled position") +
  ylab("imputation accuracy")
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
|a550588-4269754-110716-105_A01.CEL | 0.9590743|
|a550588-4269754-110716-105_A04.CEL | 0.9578656|
|a550588-4269754-110716-105_A05.CEL | 0.9806176|
|a550588-4269754-110716-105_A07.CEL | 0.9637168|
|a550588-4269754-110716-105_A09.CEL | 0.9740004|
|a550588-4269754-110716-105_A10.CEL | 0.9697819|
|a550588-4269754-110716-105_A11.CEL | 0.9651517|
|a550588-4269754-110716-105_B01.CEL | 0.9721161|
|a550588-4269754-110716-105_B02.CEL | 0.9669730|
|a550588-4269754-110716-105_B03.CEL | 0.9699430|
|a550588-4269754-110716-105_B06.CEL | 0.9545157|
|a550588-4269754-110716-105_B07.CEL | 0.9672614|
|a550588-4269754-110716-105_B08.CEL | 0.9743092|
|a550588-4269754-110716-105_B10.CEL | 0.9579780|
|a550588-4269754-110716-105_B11.CEL | 0.9610442|
|a550588-4269754-110716-105_B12.CEL | 0.9714463|
|a550588-4269754-110716-105_C01.CEL | 0.9681291|
|a550588-4269754-110716-105_C03.CEL | 0.9702946|
|a550588-4269754-110716-105_C04.CEL | 0.9516580|
|a550588-4269754-110716-105_C05.CEL | 0.9624640|
|a550588-4269754-110716-105_C06.CEL | 0.9743233|
|a550588-4269754-110716-105_C09.CEL | 0.9677266|
|a550588-4269754-110716-105_C10.CEL | 0.9628438|
|a550588-4269754-110716-105_C11.CEL | 0.9784988|
|a550588-4269754-110716-105_C12.CEL | 0.9626181|
|a550588-4269754-110716-105_D02.CEL | 0.9678187|
|a550588-4269754-110716-105_D03.CEL | 0.9798513|
|a550588-4269754-110716-105_D04.CEL | 0.9730636|
|a550588-4269754-110716-105_D05.CEL | 0.9360162|
|a550588-4269754-110716-105_D06.CEL | 0.9352021|
|a550588-4269754-110716-105_D07.CEL | 0.9386487|
|a550588-4269754-110716-105_D08.CEL | 0.9448051|
|a550588-4269754-110716-105_D09.CEL | 0.9428724|
|a550588-4269754-110716-105_D10.CEL | 0.9521092|
|a550588-4269754-110716-105_D11.CEL | 0.9428858|
|a550588-4269754-110716-105_D12.CEL | 0.9379167|
|a550588-4269754-110716-105_E01.CEL | 0.9445176|
|a550588-4269754-110716-105_E02.CEL | 0.9583778|
|a550588-4269754-110716-105_E03.CEL | 0.9552215|
|a550588-4269754-110716-105_E04.CEL | 0.9474559|
|a550588-4269754-110716-105_E05.CEL | 0.9644766|
|a550588-4269754-110716-105_E06.CEL | 0.9603829|
|a550588-4269754-110716-105_E07.CEL | 0.9681327|
|a550588-4269754-110716-105_E08.CEL | 0.9392075|
|a550588-4269754-110716-105_E09.CEL | 0.9102582|
|a550588-4269754-110716-105_E10.CEL | 0.9570978|
|a550588-4269754-110716-105_E11.CEL | 0.9405448|
|a550588-4269754-110716-105_E12.CEL | 0.9604677|
|a550588-4269754-110716-105_F01.CEL | 0.9337175|
|a550588-4269754-110716-105_F02.CEL | 0.9528669|
|a550588-4269754-110716-105_F03.CEL | 0.9526688|
|a550588-4269754-110716-105_F04.CEL | 0.9602465|
|a550588-4269754-110716-105_F05.CEL | 0.9391479|
|a550588-4269754-110716-105_F06.CEL | 0.9339789|
|a550588-4269754-110716-105_F07.CEL | 0.9089079|
|a550588-4269754-110716-105_F08.CEL | 0.9351913|
|a550588-4269754-110716-105_F09.CEL | 0.9511208|
|a550588-4269754-110716-105_F10.CEL | 0.9600205|
|a550588-4269754-110716-105_F11.CEL | 0.9337982|
|a550588-4269754-110716-105_F12.CEL | 0.9464058|
|a550588-4269754-110716-105_G01.CEL | 0.9480213|
|a550588-4269754-110716-105_G02.CEL | 0.9497627|
|a550588-4269754-110716-105_G03.CEL | 0.9306985|
|a550588-4269754-110716-105_G04.CEL | 0.9349626|
|a550588-4269754-110716-105_G05.CEL | 0.9564464|
|a550588-4269754-110716-105_G06.CEL | 0.9308049|
|a550588-4269754-110716-105_G07.CEL | 0.9441886|
|a550588-4269754-110716-105_G08.CEL | 0.9600901|
|a550588-4269754-110716-105_G09.CEL | 0.9305566|
|a550588-4269754-110716-105_G10.CEL | 0.7934482|
|a550588-4269754-110716-105_G11.CEL | 0.9541743|
|a550588-4269754-110716-105_G12.CEL | 0.9391173|
|a550588-4269754-110716-105_H01.CEL | 0.9511263|
|a550588-4269754-110716-105_H02.CEL | 0.9552141|
|a550588-4269754-110716-105_H03.CEL | 0.9381345|
|a550588-4269754-110716-105_H04.CEL | 0.9647660|
|a550588-4269754-110716-105_H05.CEL | 0.9418207|
|a550588-4269754-110716-105_H06.CEL | 0.9729396|
|a550588-4269754-110716-105_H07.CEL | 0.9597190|
|a550588-4269754-110716-105_H08.CEL | 0.9528651|
|a550588-4269754-110716-105_H09.CEL | 0.9588399|
|a550588-4269754-110716-105_H10.CEL | 0.9277263|
|a550588-4269754-110716-105_H11.CEL | 0.9480905|
|a550588-4269754-110716-105_H12.CEL | 0.9620997|
|a550588-4269754-110716-107_H07.CEL | 0.9448157|
|a550588-4269754-110716-107_H08.CEL | 0.9577452|
|a550588-4269754-110716-107_H09.CEL | 0.9358414|
|a550588-4269754-110716-107_H10.CEL | 0.9654534|
|a550588-4269754-110716-107_H11.CEL | 0.9600325|
|a550588-4269754-110716-107_H12.CEL | 0.9623887|

## Save data


```r
save(anim_corr, accuracies, file = "../output/2-imputation_accuracy.RData")
```

