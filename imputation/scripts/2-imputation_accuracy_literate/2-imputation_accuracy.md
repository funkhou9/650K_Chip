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

Load SNPs and physical positions from both Affy and Illumina platforms


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
                        "correlation" = numeric(),
                        stringsAsFactors = FALSE)
for (i in 1:nrow(raw_geno_affy)) {
  anim_corr[i, ] <- list("animal" = rownames(raw_geno_affy)[i],
                         "correlation" = cor(imp_geno_affy[i, ],
                                             raw_geno_affy[i, ],
                                             use = "complete.obs"))
}

snp_corr <- c()
```

```r
for (i in 1:ncol(raw_geno_affy)) {
  snp_corr[i] <- cor(imp_geno_affy[, i],
                     raw_geno_affy[, i],
                     use = "complete.obs")
}
```

For plotting, append color codes for each SNP. Colors should alternate with
alternating chromosomes.


```r
correlations <-
  tibble("position" = pos_list$Affy650[match(colnames(raw_geno_affy), marker_list$Affy650)],
         "corr" = snp_corr) %>%
    separate(col = position,
             into = c("chr", "pos"),
             sep = ":",
             convert = TRUE) %>%
    filter(chr != "X")

correlations$chr <- as.numeric(correlations$chr)
```

```
## Warning: NAs introduced by coercion
```

```r
correlations <- arrange(correlations, chr) %>%
                  mutate(color = chr %% 2 == 0) %>%
                  na.omit()
```

## Visualize
Plot SNP-wise correlations and provide table of animal-wise correlations


```r
ggplot(correlations, aes(x = seq_along(corr),
                         y = corr,
                         color = color)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_hline(aes(yintercept = mean(corr, na.rm = TRUE))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(breaks = c(-0.25, 0.0, 0.25, 0.5, 0.75, 1.0)) +
  xlab("position") +
  ylab("correlation")
```

![plot of chunk snp_impute](figure/snp_impute-1.tiff)

```r
knitr::kable(anim_corr)
```



|animal                             | correlation|
|:----------------------------------|-----------:|
|a550588-4269754-110716-105_A01.CEL |   0.9770497|
|a550588-4269754-110716-105_A04.CEL |   0.9793291|
|a550588-4269754-110716-105_A05.CEL |   0.9896949|
|a550588-4269754-110716-105_A07.CEL |   0.9797389|
|a550588-4269754-110716-105_A09.CEL |   0.9842002|
|a550588-4269754-110716-105_A10.CEL |   0.9822569|
|a550588-4269754-110716-105_A11.CEL |   0.9796473|
|a550588-4269754-110716-105_B01.CEL |   0.9839482|
|a550588-4269754-110716-105_B02.CEL |   0.9829427|
|a550588-4269754-110716-105_B03.CEL |   0.9830469|
|a550588-4269754-110716-105_B06.CEL |   0.9730902|
|a550588-4269754-110716-105_B07.CEL |   0.9803156|
|a550588-4269754-110716-105_B08.CEL |   0.9849708|
|a550588-4269754-110716-105_B10.CEL |   0.9757811|
|a550588-4269754-110716-105_B11.CEL |   0.9770628|
|a550588-4269754-110716-105_B12.CEL |   0.9842776|
|a550588-4269754-110716-105_C01.CEL |   0.9820012|
|a550588-4269754-110716-105_C03.CEL |   0.9838915|
|a550588-4269754-110716-105_C04.CEL |   0.9727716|
|a550588-4269754-110716-105_C05.CEL |   0.9798602|
|a550588-4269754-110716-105_C06.CEL |   0.9835050|
|a550588-4269754-110716-105_C09.CEL |   0.9808143|
|a550588-4269754-110716-105_C10.CEL |   0.9800583|
|a550588-4269754-110716-105_C11.CEL |   0.9860624|
|a550588-4269754-110716-105_C12.CEL |   0.9788666|
|a550588-4269754-110716-105_D02.CEL |   0.9822790|
|a550588-4269754-110716-105_D03.CEL |   0.9877293|
|a550588-4269754-110716-105_D04.CEL |   0.9837055|
|a550588-4269754-110716-105_D05.CEL |   0.9676918|
|a550588-4269754-110716-105_D06.CEL |   0.9618515|
|a550588-4269754-110716-105_D07.CEL |   0.9652827|
|a550588-4269754-110716-105_D08.CEL |   0.9711666|
|a550588-4269754-110716-105_D09.CEL |   0.9688939|
|a550588-4269754-110716-105_D10.CEL |   0.9728114|
|a550588-4269754-110716-105_D11.CEL |   0.9687311|
|a550588-4269754-110716-105_D12.CEL |   0.9650085|
|a550588-4269754-110716-105_E01.CEL |   0.9660048|
|a550588-4269754-110716-105_E02.CEL |   0.9768156|
|a550588-4269754-110716-105_E03.CEL |   0.9762445|
|a550588-4269754-110716-105_E04.CEL |   0.9695411|
|a550588-4269754-110716-105_E05.CEL |   0.9800382|
|a550588-4269754-110716-105_E06.CEL |   0.9763698|
|a550588-4269754-110716-105_E07.CEL |   0.9803544|
|a550588-4269754-110716-105_E08.CEL |   0.9675161|
|a550588-4269754-110716-105_E09.CEL |   0.9491561|
|a550588-4269754-110716-105_E10.CEL |   0.9759595|
|a550588-4269754-110716-105_E11.CEL |   0.9680558|
|a550588-4269754-110716-105_E12.CEL |   0.9781258|
|a550588-4269754-110716-105_F01.CEL |   0.9641445|
|a550588-4269754-110716-105_F02.CEL |   0.9709998|
|a550588-4269754-110716-105_F03.CEL |   0.9741156|
|a550588-4269754-110716-105_F04.CEL |   0.9757727|
|a550588-4269754-110716-105_F05.CEL |   0.9665463|
|a550588-4269754-110716-105_F06.CEL |   0.9606070|
|a550588-4269754-110716-105_F07.CEL |   0.9498672|
|a550588-4269754-110716-105_F08.CEL |   0.9623755|
|a550588-4269754-110716-105_F09.CEL |   0.9707641|
|a550588-4269754-110716-105_F10.CEL |   0.9773186|
|a550588-4269754-110716-105_F11.CEL |   0.9648950|
|a550588-4269754-110716-105_F12.CEL |   0.9698765|
|a550588-4269754-110716-105_G01.CEL |   0.9688186|
|a550588-4269754-110716-105_G02.CEL |   0.9721480|
|a550588-4269754-110716-105_G03.CEL |   0.9604881|
|a550588-4269754-110716-105_G04.CEL |   0.9639117|
|a550588-4269754-110716-105_G05.CEL |   0.9745456|
|a550588-4269754-110716-105_G06.CEL |   0.9611111|
|a550588-4269754-110716-105_G07.CEL |   0.9683700|
|a550588-4269754-110716-105_G08.CEL |   0.9775060|
|a550588-4269754-110716-105_G09.CEL |   0.9617371|
|a550588-4269754-110716-105_G10.CEL |   0.8843701|
|a550588-4269754-110716-105_G11.CEL |   0.9765973|
|a550588-4269754-110716-105_G12.CEL |   0.9664030|
|a550588-4269754-110716-105_H01.CEL |   0.9725181|
|a550588-4269754-110716-105_H02.CEL |   0.9755509|
|a550588-4269754-110716-105_H03.CEL |   0.9642962|
|a550588-4269754-110716-105_H04.CEL |   0.9811046|
|a550588-4269754-110716-105_H05.CEL |   0.9670846|
|a550588-4269754-110716-105_H06.CEL |   0.9845612|
|a550588-4269754-110716-105_H07.CEL |   0.9771563|
|a550588-4269754-110716-105_H08.CEL |   0.9741999|
|a550588-4269754-110716-105_H09.CEL |   0.9775362|
|a550588-4269754-110716-105_H10.CEL |   0.9625826|
|a550588-4269754-110716-105_H11.CEL |   0.9705227|
|a550588-4269754-110716-105_H12.CEL |   0.9768178|
|a550588-4269754-110716-107_H07.CEL |   0.9687114|
|a550588-4269754-110716-107_H08.CEL |   0.9773565|
|a550588-4269754-110716-107_H09.CEL |   0.9641925|
|a550588-4269754-110716-107_H10.CEL |   0.9786741|
|a550588-4269754-110716-107_H11.CEL |   0.9792834|
|a550588-4269754-110716-107_H12.CEL |   0.9817243|

## Save data


```r
save(anim_corr, correlations, file = "../output/2-imputation_accuracy.RData")
```

