# Script: 3-compare_calls

- *Author: Scott Funkhouser*
- *Date: 20161025*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [genotype_analysis](../../genotype_analysis.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
  - [Filter data](#filter-data)
  - [Regression](#regression)
  - [Estimate GWBC of animal](#estimate-gwbc-of-animal)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/genotype_analysis/scripts")
```

## Objectives

1. From Yorkshire animals genotyped with the 650K chip, isolate those that
were previously genotyped on the Illumina 60K beadchip. Previous 60K data
may be loaded from the `SF_PG_Industry` project
2. For each SNP genotyped on both Illumina 60K and Affymetrix 650K platforms,
test genotype call consistancy using a linear regression.
3. As seen in `2-PCA.R`, one Landrace animal resembles Yorkshire
much more than Landrace. Use common SNPs between Illumina and
Affy platforms to estimate the genome-wide breed composition (GWBC)
of the animal (*see `SF_PG_Industry`*)

## Install libraries


```r
library(devtools)
library(magrittr)
library(ggplot2)
```

Install [snpTools](https://github.com/funkhou9/snpTools/commit/6603afd1db77fb6a93ece38b4a3eeafc7fbc92f2)


```r
# install_github("funkhou9/snpTools")
library(snpTools)
```

Install [breedTools](https://github.com/funkhou9/breedTools/commit/00b77d774e31b69f885b3ccbe413f7caf92abbbb)


```r
# install_git("/mnt/research/pigsnp/raw_data/SF_working_dir/breed_compos/breedTools")
library(breedTools)
```

## Load data
Load prepped Affymetrix genotyping data from `1-data_prep`


```r
load("../1-data_prep.RData")
```

Load Yorkshire sire genotypes from `SF_PG_Industry`. About 90 animals from
this set of ~900 should have also been used for Affymetrix 650K genotyping.


```r
load(paste0("/mnt/research/pigsnp/raw_data/SF_working_dir",
            "/Yorkshire_dataset_cleaning/yorkshireDataForModel.RData"))
```

Load previous work from this project that assembled all physical positions
for Affy and Illumina chips


```r
load("/mnt/research/pigsnp/NSR/650K_Chip/snp_inspection/pos_list.RData")
```

Load Sample table from Yorkshire genotyped on plate1 of the
Affymetrix assay. The "Sample Name" in this dataset will correspond to
registration numbers


```r
plate1_samples <-
   read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                   "PigIA&MIAX_Sample_Table.csv"),
            header = TRUE,
            stringsAsFactors = FALSE)
```

Load Sample table from plate2. Those six of eight Yorkshire on plate2 will
have registration numbers in the "MSU.ID.." column.


```r
plate2_samples <-
    read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                    "PigIA&MIAX_Sample_Table.csv"),
             header = TRUE,
             stringsAsFactors = FALSE)
```

Load reference panel of allele frequencies required to compute
GWBC.


```r
data("GWBC_ref_B")
```

## Analysis
### Filter data
Obtain registration numbers and "affy array IDs" for all Yorkshire animals genotyped
on the Affy platforms.


```r
hd_yorks_reg <- c(plate1_samples$Sample.Name, plate2_samples$MSU.ID..[91:96])
hd_yorks_id <- c(plate1_samples$Best.Array, plate2_samples$Best.Array[91:96])
```

Not all numbers in `hd_yorks_reg` are necessarily in `yorkshireGenoDose`. Find which
`hd_yorks_reg` and corresponding `hd_yorks_id` are present in `yorkshireGenoDose`.


```r
idx <- hd_yorks_reg %in% rownames(yorkshireGenoDose)
```

Subset both illumina and affy datasets with animals present in both


```r
geno_60K <- yorkshireGenoDose[hd_yorks_reg[idx], ]
geno_650K <- affy_geno[paste0(hd_yorks_id[idx], ".CEL"), ]
```

Find physical positions that match between Affy and Illumina
platforms. Note that there are duplicated positions in the Illumina
map (mostly '0:0' positions and some non-'0:0' positions) that
have differing marker names. We'll exclude these physical positions.


```r
comm_pos <- intersect(pos_list$Affy650, pos_list$SNP60[!duplicated(pos_list$SNP60)])
```

Find which markers (from each platform) correspond to each common
physical position


```r
illum_markers <- marker_list$SNP60[match(comm_pos, pos_list$SNP60)]
affy_markers <- marker_list$Affy650[match(comm_pos, pos_list$Affy650)]
```

Identify markers that are present in both genotyping datasets


```r
idx <- illum_markers %in% colnames(yorkshireGenoDose) &
        affy_markers %in% colnames(affy_geno)
```

Subset SNPs according to common markers present in both datasets.


```r
illum_calls <- geno_60K[, illum_markers[idx]]
affy_calls <- geno_650K[, affy_markers[idx]]
```

Identify markers with at least a 90% call rate in both datasets


```r
idx2 <- apply(illum_calls, 2,
              function(x) {
                  sum(!is.na(x)) / length(x) >= 0.9
              })
idx3 <- apply(affy_calls, 2,
              function(x) {
                  sum(!is.na(x)) / length(x) >= 0.9
              })
```

Identify markers that are unfixed in both datasets (For a single SNP, at
two animals must have differing conclusive genotype calls)


```r
idx4 <- apply(illum_calls, 2,
              function(x) {
                  length(unique(x[!is.na(x)])) > 1
              })
idx5 <- apply(affy_calls, 2,
              function(x) {
                  length(unique(x[!is.na(x)])) > 1
              })
```

Keep only SNPs that have sufficient call rate in both datasets and
that are unfixed in both datasets.


```r
illum_calls <- illum_calls[, idx2 & idx3 & idx4 & idx5]
affy_calls <- affy_calls[, idx2 & idx3 & idx4 & idx5]
```

### Regression
Apply regression for the ith column of each datasets. When using
`mapply()`, the `...` passed must be data.frames if `mapply()` is
to iterate over columns.


```r
fits <- mapply(function(y, x) lm(y ~ x),
               as.data.frame(illum_calls),
               as.data.frame(affy_calls))
```

Tabulate marker name, "x coefficient", and standard error.
Warnings of "essentially perfect fits". Ignore these to avoid
thousands of lines of output.


```r
suppressWarnings({int_slope_r2 <- lapply(fits,
                     function(x) {
                         c("intercept" = summary(x)$coefficients[1, 1],
                           "slope" = summary(x)$coefficients[2, 1],
                           "R2" = summary(x)$r.squared)
                     })})
results <- as.data.frame(do.call(rbind, int_slope_r2))
```

How many SNPs out of 30,941 in `results` have a slope of exactly 1 and R2 of
exactly 1?


```r
nrow(results[results$slope == 1 & results$R2 == 1, ])
```

```
## [1] 2080
```

Add 4th column, a logical vector that is `TRUE` if the SNP is currently used
in GWBC calculation.


```r
results <- cbind(results, "GWBC" = rownames(results) %in% rownames(GWBC_ref_B))
```

Plot joint distribution of coefficient estimates (slopes) and R squared values across
all models. Highlight interesting points that deserve further investigation.


```r
ggplot(results, aes(x = R2, y = slope, color = GWBC)) +
    geom_point(size = 2.5, alpha = 0.4) +
    geom_point(aes(x = 1, y = 0.5), size = 10, shape = 1, color = "red") +
    geom_point(aes(x = 1, y = -1), size = 10, shape = 1, color = "darkgreen") +
    geom_point(aes(x = 0, y = 0), size = 10, shape = 1, color = "blue")
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)

The data suggests that for those SNPs used to compute GWBC, no transformation
of Affy coding to Illumina coding is needed.
Any models where the estimated slope is above 0.5 may indicate that
both Affy and Illumina platforms are counting the same allele, but genotyping
errors (or the possibility that the platforms are counting slightly different
loci) may prevent the estimated coefficients to be 1.0.
From each of the **red**, **blue**, and **green** points highlighted in the
above plot, provide a contingency table of genotype calls.

For the **red** points (points with roughly 0.5 slope and near perfect R2):


```r
results[round(results$R2, 1) == 1 & round(results$slope, 1) == 0.5, ]
```

```
##                 intercept     slope        R2  GWBC
## ALGA0002747  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0004678  9.855879e-01 0.5044219 0.9808203 FALSE
## H3GA0002791  3.947957e-16 0.5000000 1.0000000 FALSE
## MARC0082076 -3.242965e-16 0.5000000 1.0000000 FALSE
## MARC0073315  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0006547  9.680851e-01 0.5106383 0.9612015  TRUE
## ALGA0007262  1.000000e+00 0.5000000 1.0000000 FALSE
## MARC0033050  1.022114e+00 0.4913755 0.9731162 FALSE
## ASGA0006216  1.000000e+00 0.5000000 1.0000000  TRUE
## H3GA0005011  5.198663e-03 0.5061270 0.9785122 FALSE
## ALGA0010917  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0008132  4.640371e-03 0.5127610 0.9685486 FALSE
## ASGA0105436  6.767926e-16 0.5000000 1.0000000 FALSE
## M1GA0004043  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0108553  9.757306e-01 0.5096582 0.9729839 FALSE
## ALGA0027315  1.000000e+00 0.5000000 1.0000000 FALSE
## MARC0082004  1.550983e-16 0.5000000 1.0000000 FALSE
## CASI0005793 -7.049923e-17 0.5000000 1.0000000 FALSE
## ASGA0095665  9.798735e-01 0.5031627 0.9676206 FALSE
## ALGA0108257  1.000000e+00 0.5000000 1.0000000 FALSE
## MARC0014928  1.000000e+00 0.5000000 1.0000000 FALSE
## MARC0057634  4.424779e-03 0.5176991 0.9614412 FALSE
## M1GA0011382  1.000000e+00 0.5000000 1.0000000 FALSE
## MARC0031511  1.000000e+00 0.5000000 1.0000000 FALSE
## H3GA0027585  9.680851e-01 0.5106383 0.9612015 FALSE
## H3GA0028220  1.014410e+00 0.4956434 0.9800222 FALSE
## ALGA0058439  9.603960e-01 0.5148515 0.9561528 FALSE
## ALGA0058520  1.000000e+00 0.5000000 1.0000000  TRUE
## ALGA0059917  0.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0088514  3.524962e-16 0.5000000 1.0000000 FALSE
## H3GA0056386  9.658754e-01 0.5118694 0.9597552 FALSE
## MARC0039964  1.000000e+00 0.5000000 1.0000000  TRUE
## ASGA0056220  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0094334  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0071323  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0059376  3.877458e-17 0.5000000 1.0000000 FALSE
## MARC0063244  5.685980e-17 0.5000000 1.0000000 FALSE
## MARC0022985  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0062927 -4.229954e-17 0.5000000 1.0000000 FALSE
## ASGA0063418  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0063866  1.141227e-02 0.5092725 0.9619591 FALSE
## ALGA0078141  9.855879e-01 0.5044219 0.9808203 FALSE
## DRGA0013948  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0078662  9.698349e-01 0.5128059 0.9686334 FALSE
## MARC0090899  9.778859e-01 0.5086245 0.9748636 FALSE
## DRGA0014120 -4.344678e-03 0.4822592 0.9559067 FALSE
## INRA0045033  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0065768  1.010683e+00 0.4978915 0.9830165 FALSE
## ALGA0081091  1.038576e-02 0.5118694 0.9597552 FALSE
## MARC0074458  9.525826e-01 0.5215919 0.9562518 FALSE
## MARC0084778  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0087956  9.897592e-01 0.5017991 0.9835262 FALSE
## ASGA0071410  1.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0096571  9.528302e-01 0.5188679 0.9512579 FALSE
## MARC0075925  9.825493e-01 0.5060698 0.9784016 FALSE
## ASGA0089601  4.640371e-03 0.5127610 0.9685486 FALSE
## M1GA0021726  9.730988e-01 0.5111226 0.9711330 FALSE
## ALGA0095183  9.730988e-01 0.5111226 0.9711330 FALSE
## MARC0044488  9.778859e-01 0.5086245 0.9748636 FALSE
## ASGA0078296  1.000000e+00 0.5000000 1.0000000 FALSE
## ALGA0097979  1.000000e+00 0.5000000 1.0000000 FALSE
```

Examples of genotypes for this set of markers:


```r
table(fits$ALGA0002747$model)
```

```
##    x
## y    0  2
##   1  3  0
##   2  0 59
```

```r
table(fits$ALGA0004678$model)
```

```
##    x
## y    0  1  2
##   1 17  1  0
##   2  0  0 44
```

```r
table(fits$H3GA0002791$model)
```

```
##    x
## y    0  2
##   0 41  0
##   1  0 21
```

```r
table(fits$ASGA0065768$model)
```

```
##    x
## y    0  1  2
##   1 23  0  0
##   2  0  1 38
```

For the **blue** points (points that have near-zero R2 values)


```r
results[results$R2 < 0.05, ]
```

```
##             intercept       slope           R2  GWBC
## ASGA0096844 0.1080617  0.06232133 0.0121087782 FALSE
## ALGA0028021 1.7493809 -0.01188707 0.0004052411 FALSE
## DRGA0008687 0.4154786  0.04480652 0.0032858113 FALSE
## ALGA0081437 0.6379310  0.36206897 0.0071093486 FALSE
```

Examples of genotypes for this set of markers:


```r
table(fits$ASGA0096844$model)
```

```
##    x
## y    0  1  2
##   0 11 25 15
##   1  0  8  3
```

```r
table(fits$ALGA0028021$model)
```

```
##    x
## y    0  1  2
##   1  5  0 11
##   2  5 19 20
```

```r
table(fits$DRGA0008687$model)
```

```
##    x
## y    0  1  2
##   0 15 15  6
##   1  7 12  5
##   2  1  1  0
```

```r
table(fits$ALGA0081437$model)
```

```
##    x
## y    1  2
##   0  1 10
##   1  0 17
##   2  1 31
```

For the **green** points (points that have near -1 slope)


```r
results[results$slope < -0.75, ]
```

```
##             intercept      slope        R2  GWBC
## DRGA0000277  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0000883  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0003027  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0003043  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0003901  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0001119  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0001168  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0004732  2.000000 -1.0000000 1.0000000 FALSE
## MARC0007197  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0005395  2.000000 -1.0000000 1.0000000 FALSE
## MARC0016307  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0004323  2.000000 -1.0000000 1.0000000 FALSE
## MARC0026602  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0002681  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0003063  2.000000 -1.0000000 1.0000000 FALSE
## MARC0068306  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0001249  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0007818  2.000000 -1.0000000 1.0000000 FALSE
## MARC0072618  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0008458  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0002269  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0026765  1.959870 -0.9669197 0.9349338 FALSE
## ALGA0015543  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0003491  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0111164  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0103397  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0018491  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0109278  2.000000 -1.0000000 1.0000000 FALSE
## MARC0095377  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0085168  2.000000 -1.0000000 1.0000000 FALSE
## MARC0009481  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0017696  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0015231  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0020071  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0010208  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0010237  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0010289  1.972355 -0.9909623 0.9681099 FALSE
## ALGA0020644  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0023691  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0019889  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0020481  2.000000 -1.0000000 1.0000000 FALSE
## MARC0027566  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0021414  2.000000 -1.0000000 1.0000000 FALSE
## MARC0001008  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0025164  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0026243  2.000000 -1.0000000 1.0000000 FALSE
## DIAS0000686  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0008164  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0053970  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0083116  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0053959  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0103781  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0001502  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0103459  1.710396 -0.8366337 0.8984342 FALSE
## MARC0035152  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0123175  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0038840  1.934492 -0.8533868 0.8320956 FALSE
## MARC0030001  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0039771  1.952998 -0.8978930 0.9077240 FALSE
## ASGA0032442  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0040191  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0021066  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0034372  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0022157  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0022462  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0010933  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0091667  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0119079  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0085363  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0118900  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0051700  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0043130  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0043752  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0043938  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0044381  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0044392  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0118468  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0047043  2.000000 -1.0000000 1.0000000 FALSE
## CAHM0000131  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0117353  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0062449  2.000000 -1.0000000 1.0000000 FALSE
## CASI0000302  2.004425 -0.8783186 0.8673942 FALSE
## SIRI0000605  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0025413  1.936567 -0.9589552 0.9460710 FALSE
## ASGA0057107  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0069929  2.000000 -1.0000000 1.0000000 FALSE
## MARC0028090  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0057447  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0036438  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0090575  2.000000 -1.0000000 1.0000000 FALSE
## MARC0034409  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0091358  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0102775  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0072378  2.000000 -1.0000000 1.0000000 FALSE
## MARC0082696  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0037597  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0097399  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0078246  2.000000 -1.0000000 1.0000000 FALSE
## MARC0002437  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0041247  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0041845  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0066140  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0066204  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0081371  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0081570  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0042303  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0083674  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0083680  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0044392  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0085733  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0001013  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0086180  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0070239  2.000000 -1.0000000 1.0000000 FALSE
## MARC0089468  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0015447  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0020479  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0095426  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0000715  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0001388  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0046417  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0092775  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0096878  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0112265  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0092994  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099483  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099504  2.000000 -1.0000000 1.0000000 FALSE
## CASI0006594  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0090010  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0023684  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0081175  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099766  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099903  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0051929  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0100024  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0100228  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0100258  2.000000 -1.0000000 1.0000000 FALSE
```

Examples of genotypes for this set of markers:


```r
table(fits$DRGA0000277$model)
```

```
##    x
## y    0  1  2
##   0  0  0 50
##   1  0  9  0
##   2  1  0  0
```

```r
table(fits$M1GA0026765$model)
```

```
##    x
## y    0  1  2
##   0  0  1 22
##   1  0 27  1
##   2 10  0  0
```

```r
table(fits$CASI0000302$model)
```

```
##    x
## y    0  1  2
##   0  0  0  1
##   1  0  3  0
##   2 56  1  0
```

```r
table(fits$M1GA0025413$model)
```

```
##    x
## y    0  1  2
##   0  0  0 12
##   1  2 20  0
##   2 26  0  0
```

Plot joint distribution of slopes and intercepts fitted across all models.
Highlight interesting points that deserve further investigation.


```r
ggplot(results, aes(x = intercept, y = slope, color = GWBC)) +
    geom_point(size = 2.5, alpha = 0.4) +
    geom_point(aes(x = 1, y = 0.5), size = 10, shape = 1, color = "red") +
    geom_point(aes(x = 0, y = 0.5), size = 10, shape = 1, color = "darkgreen") +
    geom_point(aes(x = 2, y = -1), size = 10, shape = 1, color = "blue")
```

![plot of chunk unnamed-chunk-32](figure/unnamed-chunk-32-1.png)

From each of the **red**, **blue**, and **green** points highlighted in the
above plot, provide a contingency table of genotype calls.

For the **red** points(slope near 0.5 and intercept near 1.0):


```r
results[round(results$intercept, 1) == 1 & round(results$slope, 1) == 0.5, ]
```

```
##             intercept     slope        R2  GWBC
## ALGA0002747 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0004678 0.9855879 0.5044219 0.9808203 FALSE
## MARC0073315 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0006547 0.9680851 0.5106383 0.9612015  TRUE
## ALGA0007262 1.0000000 0.5000000 1.0000000 FALSE
## MARC0033050 1.0221141 0.4913755 0.9731162 FALSE
## ASGA0006216 1.0000000 0.5000000 1.0000000  TRUE
## ALGA0010917 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0116887 0.9535759 0.5106650 0.9362192 FALSE
## M1GA0004043 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0108553 0.9757306 0.5096582 0.9729839 FALSE
## ALGA0027315 1.0000000 0.5000000 1.0000000 FALSE
## ASGA0095665 0.9798735 0.5031627 0.9676206 FALSE
## ALGA0108257 1.0000000 0.5000000 1.0000000 FALSE
## MARC0014928 1.0000000 0.5000000 1.0000000 FALSE
## M1GA0011382 1.0000000 0.5000000 1.0000000 FALSE
## MARC0031511 1.0000000 0.5000000 1.0000000 FALSE
## H3GA0027585 0.9680851 0.5106383 0.9612015 FALSE
## H3GA0028220 1.0144102 0.4956434 0.9800222 FALSE
## ALGA0058439 0.9603960 0.5148515 0.9561528 FALSE
## ALGA0058520 1.0000000 0.5000000 1.0000000  TRUE
## H3GA0056386 0.9658754 0.5118694 0.9597552 FALSE
## MARC0039964 1.0000000 0.5000000 1.0000000  TRUE
## ASGA0056220 1.0000000 0.5000000 1.0000000 FALSE
## ASGA0094334 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0071323 1.0000000 0.5000000 1.0000000 FALSE
## MARC0066980 1.0000000 0.5000000 0.9194805  TRUE
## MARC0022985 1.0000000 0.5000000 1.0000000 FALSE
## ASGA0063418 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0078141 0.9855879 0.5044219 0.9808203 FALSE
## DRGA0013948 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0078662 0.9698349 0.5128059 0.9686334 FALSE
## MARC0090899 0.9778859 0.5086245 0.9748636 FALSE
## INRA0045033 1.0000000 0.5000000 1.0000000 FALSE
## ASGA0065768 1.0106832 0.4978915 0.9830165 FALSE
## M1GA0019916 0.9531416 0.5154420 0.9449769 FALSE
## ALGA0122296 0.9790301 0.5080821 0.9276288  TRUE
## MARC0074458 0.9525826 0.5215919 0.9562518 FALSE
## MARC0084778 1.0000000 0.5000000 1.0000000 FALSE
## ALGA0087956 0.9897592 0.5017991 0.9835262 FALSE
## ASGA0071410 1.0000000 0.5000000 1.0000000 FALSE
## ASGA0096571 0.9528302 0.5188679 0.9512579 FALSE
## MARC0075925 0.9825493 0.5060698 0.9784016 FALSE
## M1GA0021726 0.9730988 0.5111226 0.9711330 FALSE
## ALGA0095183 0.9730988 0.5111226 0.9711330 FALSE
## MARC0044488 0.9778859 0.5086245 0.9748636 FALSE
## ASGA0078296 1.0000000 0.5000000 1.0000000 FALSE
## ASGA0095831 1.0000000 0.5000000 0.9077381  TRUE
## ALGA0097979 1.0000000 0.5000000 1.0000000 FALSE
```

Examples of genotypes for this set of markers:


```r
table(fits$ALGA0002747$model)
```

```
##    x
## y    0  2
##   1  3  0
##   2  0 59
```

```r
table(fits$ALGA0004678$model)
```

```
##    x
## y    0  1  2
##   1 17  1  0
##   2  0  0 44
```

```r
table(fits$MARC0073315$model)
```

```
##    x
## y    0  2
##   1  1  0
##   2  0 61
```

```r
table(fits$ASGA0006216$model)
```

```
##    x
## y    0  2
##   1 14  0
##   2  0 48
```

For the **green** points(slope near 0.5 and intercept near 0.0):


```r
results[round(results$intercept, 1) == 0 & round(results$slope, 1) == 0.5, ]
```

```
##                 intercept     slope        R2  GWBC
## ASGA0003370  1.639344e-02 0.5491803 0.8585519 FALSE
## H3GA0002791  3.947957e-16 0.5000000 1.0000000 FALSE
## MARC0082076 -3.242965e-16 0.5000000 1.0000000 FALSE
## H3GA0005011  5.198663e-03 0.5061270 0.9785122 FALSE
## ASGA0008132  4.640371e-03 0.5127610 0.9685486 FALSE
## ASGA0105436  6.767926e-16 0.5000000 1.0000000 FALSE
## MARC0082004  1.550983e-16 0.5000000 1.0000000 FALSE
## CASI0005793 -7.049923e-17 0.5000000 1.0000000 FALSE
## MARC0006040  4.827330e-03 0.5414036 0.9178680 FALSE
## MARC0057634  4.424779e-03 0.5176991 0.9614412 FALSE
## ALGA0059917  0.000000e+00 0.5000000 1.0000000 FALSE
## ASGA0088514  3.524962e-16 0.5000000 1.0000000 FALSE
## ASGA0099478  1.629914e-02 0.5381911 0.9027030 FALSE
## ASGA0059376  3.877458e-17 0.5000000 1.0000000 FALSE
## MARC0063244  5.685980e-17 0.5000000 1.0000000 FALSE
## ASGA0062927 -4.229954e-17 0.5000000 1.0000000 FALSE
## ASGA0063866  1.141227e-02 0.5092725 0.9619591 FALSE
## DRGA0014120 -4.344678e-03 0.4822592 0.9559067 FALSE
## ALGA0081091  1.038576e-02 0.5118694 0.9597552 FALSE
## INRA0061306  5.240747e-03 0.5335735 0.9236889 FALSE
## ASGA0089601  4.640371e-03 0.5127610 0.9685486 FALSE
```

Examples of genotypes for this set of markers:


```r
table(fits$ASGA0003370$model)
```

```
##    x
## y    0  1  2
##   0 49  1  0
##   1  0  5  7
```

```r
table(fits$H3GA0002791$model)
```

```
##    x
## y    0  2
##   0 41  0
##   1  0 21
```

```r
table(fits$MARC0082076$model)
```

```
##    x
## y    0  2
##   0 55  0
##   1  0  7
```

```r
table(fits$H3GA0005011$model)
```

```
##    x
## y    0  1  2
##   0 47  0  0
##   1  0  1 14
```

For the **blue** points(slope near 0.5 and intercept near 0.0):


```r
results[round(results$intercept, 1) == 2 & round(results$slope, 1) == -1, ]
```

```
##             intercept      slope        R2  GWBC
## DRGA0000277  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0000883  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0003027  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0003043  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0003901  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0001119  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0001168  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0004732  2.000000 -1.0000000 1.0000000 FALSE
## MARC0007197  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0005395  2.000000 -1.0000000 1.0000000 FALSE
## MARC0016307  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0004323  2.000000 -1.0000000 1.0000000 FALSE
## MARC0026602  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0002681  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0003063  2.000000 -1.0000000 1.0000000 FALSE
## MARC0068306  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0001249  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0007818  2.000000 -1.0000000 1.0000000 FALSE
## MARC0072618  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0008458  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0002269  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0026765  1.959870 -0.9669197 0.9349338 FALSE
## ALGA0015543  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0003491  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0111164  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0103397  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0018491  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0109278  2.000000 -1.0000000 1.0000000 FALSE
## MARC0095377  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0085168  2.000000 -1.0000000 1.0000000 FALSE
## MARC0009481  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0017696  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0015231  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0020071  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0010208  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0010237  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0010289  1.972355 -0.9909623 0.9681099 FALSE
## ALGA0020644  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0023691  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0019889  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0020481  2.000000 -1.0000000 1.0000000 FALSE
## MARC0027566  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0021414  2.000000 -1.0000000 1.0000000 FALSE
## MARC0001008  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0025164  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0026243  2.000000 -1.0000000 1.0000000 FALSE
## DIAS0000686  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0008164  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0053970  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0083116  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0053959  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0103781  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0001502  2.000000 -1.0000000 1.0000000 FALSE
## MARC0035152  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0123175  2.000000 -1.0000000 1.0000000 FALSE
## MARC0030001  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0032442  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0040191  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0021066  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0034372  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0022157  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0022462  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0010933  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0091667  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0119079  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0085363  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0118900  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0051700  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0043130  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0043752  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0043938  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0044381  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0044392  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0118468  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0047043  2.000000 -1.0000000 1.0000000 FALSE
## CAHM0000131  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0117353  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0062449  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0000605  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0057107  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0069929  2.000000 -1.0000000 1.0000000 FALSE
## MARC0028090  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0057447  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0036438  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0090575  2.000000 -1.0000000 1.0000000 FALSE
## MARC0034409  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0091358  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0102775  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0072378  2.000000 -1.0000000 1.0000000 FALSE
## MARC0082696  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0037597  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0097399  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0078246  2.000000 -1.0000000 1.0000000 FALSE
## MARC0002437  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0041247  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0041845  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0066140  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0066204  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0081371  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0081570  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0042303  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0083674  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0083680  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0044392  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0085733  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0001013  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0086180  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0070239  2.000000 -1.0000000 1.0000000 FALSE
## MARC0089468  2.000000 -1.0000000 1.0000000 FALSE
## DRGA0015447  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0020479  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0095426  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0000715  2.000000 -1.0000000 1.0000000 FALSE
## SIRI0001388  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0046417  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0092775  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0096878  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0112265  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0092994  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099483  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099504  2.000000 -1.0000000 1.0000000 FALSE
## CASI0006594  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0090010  2.000000 -1.0000000 1.0000000 FALSE
## M1GA0023684  2.000000 -1.0000000 1.0000000 FALSE
## ASGA0081175  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099766  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0099903  2.000000 -1.0000000 1.0000000 FALSE
## H3GA0051929  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0100024  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0100228  2.000000 -1.0000000 1.0000000 FALSE
## ALGA0100258  2.000000 -1.0000000 1.0000000 FALSE
```

Examples of genotypes for this set of markers:


```r
table(fits$DRGA0000277$model)
```

```
##    x
## y    0  1  2
##   0  0  0 50
##   1  0  9  0
##   2  1  0  0
```

```r
table(fits$M1GA0000883$model)
```

```
##    x
## y    0  1  2
##   0  0  0 41
##   1  0 18  0
##   2  3  0  0
```

```r
table(fits$ASGA0003027$model)
```

```
##    x
## y    0  1  2
##   0  0  0 16
##   1  0 30  0
##   2 16  0  0
```

```r
table(fits$ASGA0003043$model)
```

```
##    x
## y    0  1  2
##   0  0  0 16
##   1  0 29  0
##   2 17  0  0
```

### Estimate GWBC of animal
From `2-PCA.R`, I know that animal `a550588-4269754-110716-107_F06.CEL`
clusters much closer with Yorkshire genotyped on the Affy chip

Isoloate this "suspect Landrace" animal and compute its GWBC
based on the reference panel that has been previously developed
for the National Swine Registry and used in `SF_PG_Industry`.


```r
sus_landrace <- affy_geno["a550588-4269754-110716-107_F06.CEL", ]
```

Keep only SNPs present on the Illumina platform and convert to
Illumina names


```r
sus_landrace <- sus_landrace[affy_markers[idx]]
names(sus_landrace) <- illum_markers[idx]
```

Keep only SNPs that are used in the reference panel to estimate
GWBC


```r
sus_landrace <- sus_landrace[names(sus_landrace) %in% rownames(GWBC_ref_B)]
```

Use `breedTools` to estimate GWBC of suspect landrace animal using transformed
genotypes


```r
breedTools:::QPsolve(sus_landrace, GWBC_ref_B)
```

```
##         Duroc     Hampshire      Landrace     Yorkshire            R2 
##  0.000000e+00 -5.374697e-17  2.369919e-18  1.000000e+00  3.214907e-01
```

