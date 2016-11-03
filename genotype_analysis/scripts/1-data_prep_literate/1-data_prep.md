# Script: 1-data-prep

- *Author: Scott Funkhouser*
- *Date: 20161013*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [genotype_analysis](../../genotype_analysis.md)*

## Table of Contents

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
  - [From Plate one](from-plate-one)
  - [From Plate two](from-plate-two)
4. [Analysis](#analysis)
  - [Filter samples](#filter-samples)
  - [Filter SNPs](#filter-snps)
  - [Add breed info](#add-breed-info)
5. [Save data](#save-data)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/genotype_analysis/scripts")
```

## Objectives

1. Load raw data from Affymetrix containing 650K genotypes for 180 animals
of Yorkshire, Landrace, Hampshire, and Duroc breeds.
2. Attach breed information to each sample.

## Install libraries


```r
library(devtools)
library(magrittr)
```

load [`snpTools`](https://github.com/funkhou9/snpTools/commit/6603afd1db77fb6a93ece38b4a3eeafc7fbc92f2)


```r
library(snpTools)
```

## Load data

**NOTE:** a decent amount of data must be loaded from each plate, used to:

1. Attach breed information to each sample
2. Keep only PolyHighResolution SNPs

### From plate one
Load **genotyping calls** from first 96-well plate

> The first 96-well plate contains all Yorkshire animals, however some are from
> ISU, which we want to exclude. These samples are identified by their 96-well
> position. `check.names = FALSE` is used to preserve original column names,
> which contain "-".



```r
plate1_geno <-
    read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                      "cluster/AxiomGT1.calls.txt"),
               header = TRUE,
               stringsAsFactors = FALSE,
               check.names = FALSE)
```

Load **well information**

> Contains information to connect the 96-well position with animal ID



```r
plate1_wells <-
   read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                   "PigIA&MIAX_Sample_Table.csv"),
            header = TRUE,
            stringsAsFactors = FALSE)
```

Load **PolyHighResolution SNP list**

> Contains names of probes that passed all filtering criteria



```r
plate1_goodSNPs <-
    read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX-selected/",
                      "cluster/filtered/PolyHighResolution.ps"),
               header = TRUE,
               stringsAsFactors = FALSE)
```

### From plate two
Load **supplemental key** `pt2_IDs_and_breeds.csv` provided by Nancy Raney.

> The key contains information necessary to connect *animal ID* (1, 2, 3, 4)
> with information about the pig (name, breed, etc.)



```r
plate2_breeds <-
    read.csv("/mnt/research/pigsnp/raw_data/affymetrix_hd/pt2_IDs_and_breeds.csv",
             header = TRUE,
             stringsAsFactors = FALSE)
```

Load **genotyping calls** from the second 96-well plate

> Contains Landrace, Hampshire, and Duroc animals identified by the well they
> were positioned in in the 96-well plate



```r
plate2_geno <-
    read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                      "cluster_final/AxiomGT1.calls.txt"),
               header = TRUE,
               stringsAsFactors = FALSE,
               check.names = FALSE)
```

Load **well information**

> Contains information to connect the 96-well position with animal ID



```r
plate2_wells <-
    read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                    "PigIA&MIAX_Sample_Table.csv"),
             header = TRUE,
             stringsAsFactors = FALSE)
```

Load **PolyHighResolution SNP list**

> Contains names of probes that passed all filtering criteria



```r
plate2_goodSNPs <-
 read.table(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/PigIA&MIAX_pt2/",
                   "cluster_final/SNPolisher/PolyHighResolution.ps"),
            header = TRUE,
            stringsAsFactors = FALSE)
```

## Analysis
### Filter samples
From plate 1, ensure that samples from ISU are not included in `plate1_geno`.
Samples from ISU were all of the samples that failed QC.


```r
fail_samples <- plate1_wells[plate1_wells$Sample.Status == "Fail", "Best.Array"]
fail_samples <- paste0(fail_samples, ".CEL")
sum(fail_samples %in% colnames(plate1_geno)[-1])
```

```
## [1] 0
```

Conversely, confirm that all samples that Pass QC are present in `plate1_geno`.


```r
pass_samples <- plate1_wells[plate1_wells$Sample.Status == "Pass", "Best.Array"]
pass_samples <- paste0(pass_samples, ".CEL")
sum(pass_samples %in% colnames(plate1_geno)[-1])
```

```
## [1] 84
```

### Filter SNPs
Use `plate1_goodSNPs` to only keep PolyHighResolution SNPs from `plate1_geno`


```r
plate1_geno <- plate1_geno[plate1_geno$probeset_id %in% plate1_goodSNPs$probeset_id, ]
dim(plate1_geno)
```

```
## [1] 514176     85
```

Likewise use `plate2_goodSNPs` to only keep PolyHighResolution SNPs from `plate2_geno`


```r
plate2_geno <- plate2_geno[plate2_geno$probeset_id %in% plate2_goodSNPs$probeset_id, ]
dim(plate2_geno)
```

```
## [1] 532909     97
```

Note that for each plate, there are n+1 columns since the first column is used
for the probe ID
### Add breed info
All animals in `plate1` should be Yorkshire. Conversely, animals in `plate2` are
a combination of breeds that need identified. Will use `plate2_breeds` to identify
the breed of each animal in `plate2`. Firstly, find which "MSU ID#" correspond with
each breed in {Duroc, Hampshire, and Landrace}. Then, find which "Best Array IDs"
correspond with each "MSU ID#"


```r
# sapply returns list of length 3 with MSU IDs for each of the 3 breeds.
array_ids <-
    sapply(c("Duroc", "Hampshire", "Landrace"),
           function(x) {
               plate2_breeds$MSU.ID..[plate2_breeds$Breed == x]
           }) %>%
    # lapply returns list of length 3 with Best.Array IDs for each of the 3 breeds.
    lapply(function(x) {
               plate2_wells$Best.Array[plate2_wells$MSU.ID.. %in% x]
           })
```

Isolate Yorkshire "Best Array IDs" from `plate1_wells`


```r
array_ids$Yorkshire <-
    plate1_wells$Best.Array[plate1_wells$Sample.Status == "Pass"]
```

Add Yorkshires from `plate2` using the last 8 elements in `plate2_wells$Best.Array`


```r
array_ids$Yorkshire <- c(array_ids$Yorkshire, plate2_wells$Best.Array[89:96])
```

Combine genotyping datasets into one data.frame `affy_geno`, which will have
probe IDs as colnames, animal names as rownames;
likewise SNPs in columns, and animals in rows.

Transpose both genotyping matrices


```r
plate1_geno <- t(plate1_geno)
plate2_geno <- t(plate2_geno)
```

Add colnames


```r
colnames(plate1_geno) <- unname(plate1_geno[1, ])
colnames(plate2_geno) <- unname(plate2_geno[1, ])
```

Remove the first row of each


```r
plate1_geno <- plate1_geno[-1, ]
plate2_geno <- plate2_geno[-1, ]
```

Merge genotyping datasets with `snpTools::merge_geno`


```r
affy_geno <- snpTools::merge_geno(plate1_geno, plate2_geno)
```

In order to agree with previous code and implementation, `affy_geno` will require
numeric coding of genotypes, with missing genotypes represented as `NA` rather
than `-1`


```r
storage.mode(affy_geno) <- "numeric"
affy_geno[affy_geno == -1] <- NA
```

How many SNPs contain at least one missing genotype?


```r
sum(apply(affy_geno, 2, function(x) any(is.na(x))))
```

```
## [1] 91620
```

## Save data


```r
save(affy_geno, array_ids, file = "../1-data_prep.RData")
```

