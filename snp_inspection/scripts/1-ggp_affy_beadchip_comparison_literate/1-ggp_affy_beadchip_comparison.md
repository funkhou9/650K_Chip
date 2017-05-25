# Script: 1-ggp_affy_beadchip_comparison.R

- *Author: Scott Funkhouser*
- *Date: 20160504*
- *Project: [650K_Chip](../../../README.md)*
- *Sub Folder: [snp_inspection](../../snp_inspection.md)*

## Table of Contents

1. [Background](#background)
2. [Objectives](#objectives)
3. [Install libraries](#install-libraries)
4. [Load data](#load-data)
5. [Analysis](#analysis)
  - [Objective One](#objective-one)
  - [Objective Two](#objective-two)
    - [Save physical positions](#save-physical-positions)
  - [Objective Three](#objective-three)
	   - [Write VEP input to disk](#write-vep-input-to-disk)
6. [Conclusion](#conclusion)


```r
setwd("/mnt/research/pigsnp/NSR/650K_Chip/snp_inspection/scripts")
```

## Background
The Illumina Porcine BeadChip60 is being replaced by the GGP-HD, which contains
a similar number of SNPs (60K) that aren't all the same as the ones on the BeadChip60. A subset of
BeadChip60 SNPs "KIT SNPs" are used in the package [`breedTools`](https://github.com/funkhou9/breedTools)
to estimate KIT-based breed probabilities. It was discovered that not all of these SNPs are present
on the GGP-HD, instead other SNPs are present on that chip in similar positions ("GGP-HD KIT SNPs").
In order to use both animals genotyped on the older BeadChip60 and GGP-HD as reference animals in
KIT-based breed probabilities, we'd need to impute SNPs between platforms. To determine imputation
accuracy, we would need to have animals that are genotyped on all "KIT SNPs" and "GGP-HD KIT SNPs".

Separately from the need to characterize the Affymetrix 650K chip for the purposes of Kit-based
breed probability work, we wish to characterize the SNPs present on the affy by comparing to those
on the GGP-HD, GGP-LD, and Porcine SNP60 by position.
## Objectives

1. Determine if "KIT SNPs" and "GGP-HD KIT SNPs" are all present on the Affy 650K chip.
2. Determine overlap of SNPs between GGP-HD, GGP-LD, SNP60, and Affy650 using physical position.
3. Use [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)
to annotate Affy650K Chip. Save VEP results to where raw 650K data is stored:
`/mnt/research/pigsnp/raw_data/affymetrix_hd/`

## Install libraries


```r
library(devtools)
library(magrittr)
library(methods)
library(VennDiagram)
```

```
## Loading required package: grid
```

```
## Loading required package: futile.logger
```

## Load data
Paste in KIT snp information as a data.frame. Information comes from GGP-LD/SNP60 map.


```r
kit_snps <- data.frame("snp" = c("ALGA0123881",
					 			 "MARC0034580",
					 			 "ALGA0102731",
					 			 "ALGA0115258",
					 			 "ALGA0047798",
					 			 "ALGA0047807",
					 			 "ALGA0047809"),
					   "chr" = rep(8, 7),
					   "pos" = c(43068687,
								 43425758,
								 43462399,
								 43651639,
								 43730377,
								 43878303,
								 43916646))
```

Read in Affy 650K SNP map


```r
affy_map <- read.csv("/mnt/research/pigsnp/raw_data/affymetrix_hd/Axiom_PigHD_v1_Annotation.r1.1.csv",
					 header = TRUE,
					 comment.char = "#")
```

Read in GGP-HD SNP map


```r
ggphd_map <- read.table("/mnt/research/pigsnp/raw_data/80K_SNP_Map.txt", skip = 1)
```

Read in GGP-LD SNP map, variable name `org_LowD_chip`


```r
load("/mnt/research/pigsnp/raw_data/SF_working_dir/LowD_chip_maps.RData")
```

Read in Porcine SNP60 map


```r
snp60_map <- read.table("/mnt/research/pigsnp/raw_data/SF_working_dir/map_updated_commercial.txt")
```

## Analysis
### Objective One
View "KIT SNPs" within affymetrix 650 panel. Are all 7 SNPs present? Answer: yes but with completely
different names.


```r
affy_map[affy_map$Chromosome %in% kit_snps$chr & affy_map$Physical.Position %in% kit_snps$pos, 1:5]
```

```
##       Probe.Set.ID    Affy.SNP.ID dbSNP.RS.ID Chromosome Physical.Position
## 26100 AX-116691049 Affx-114940851         ---          8          43068687
## 26102 AX-116691051 Affx-115060428         ---          8          43425758
## 26103 AX-116345663 Affx-115223271         ---          8          43462399
## 26104 AX-116345697 Affx-114629526         ---          8          43651639
## 26105 AX-116345721 Affx-115245040         ---          8          43730377
## 26107 AX-116691053 Affx-115203820         ---          8          43878303
## 26108 AX-116345774 Affx-114753388         ---          8          43916646
```

What SNPs are present within the KIT SNP range (chr8:43068687-43916646) on the GGP-HD? Answer:
there are 8 SNPs in this region, 3 of which are part of the original BeadChip60 (ALGA0047809,
ALGA0102731, and ALGA0123881)


```r
ggp_kit <- ggphd_map[ggphd_map[, 3] == 8 &
				(ggphd_map[, 4] >= 43068687 & ggphd_map[, 4] <= 43916646), ]

ggp_kit
```

```
##          V1                 V2 V3       V4     V5    V6  V7  V8 V9
## 6069   6069        ALGA0047809  8 43916646 0.8160 [T/C] BOT TOP  0
## 11474 11474        ALGA0102731  8 43462399 0.3624 [A/G] TOP BOT  0
## 14878 14878        ALGA0123881  8 43068687 0.8646 [T/C] BOT BOT  0
## 64267 64279 WU_10.2_8_43277756  8 43277756 0.8870 [T/C] BOT TOP  0
## 64269 64281 WU_10.2_8_43364095  8 43364095 0.7347 [T/C] BOT BOT  0
## 64270 64282 WU_10.2_8_43545395  8 43545395 0.5629 [A/G] TOP BOT  0
## 64272 64284 WU_10.2_8_43641500  8 43641500 0.8158 [T/C] BOT TOP  0
## 64273 64285 WU_10.2_8_43818637  8 43818637 0.0000 [A/G] TOP TOP  0
```

Are these 8 "GGP KIT" SNPs on the affymetrix 650K? Answer: 7 of the 8 GGP KIT SNPs are present on
the Affy 650K


```r
affy_map[affy_map$Chromosome %in% ggp_kit[, 3] & affy_map$Physical.Position %in% ggp_kit[, 4], 1:5]
```

```
##        Probe.Set.ID    Affy.SNP.ID dbSNP.RS.ID Chromosome
## 26100  AX-116691049 Affx-114940851         ---          8
## 26103  AX-116345663 Affx-115223271         ---          8
## 26108  AX-116345774 Affx-114753388         ---          8
## 95509  AX-116736314 Affx-114850107         ---          8
## 95513  AX-116345643 Affx-114992210         ---          8
## 161054 AX-116345747 Affx-114935317         ---          8
## 555768 AX-116345695 Affx-114958854         ---          8
##        Physical.Position
## 26100           43068687
## 26103           43462399
## 26108           43916646
## 95509           43277756
## 95513           43364095
## 161054          43818637
## 555768          43641500
```

### Objective Two
To compare physical positions across platforms, first combine chromosome and position of each SNP
into a single character (formatted `chr:pos`) which will provide an easier identifyer to work with.


```r
snp60_pos <- paste(snp60_map$chr,
	 			   snp60_map$pos,
	 			   sep = ":")
affy_pos <- paste(affy_map$Chromosome,
				  affy_map$Physical.Position,
				  sep = ":")
ggphd_pos <- paste(ggphd_map$V3,
				   ggphd_map$V4,
				   sep = ":")
ggpld_pos <- paste(org_LowD_chip$chr,
				   org_LowD_chip$pos,
				   sep = ":")

pos_list <- list("SNP60" = snp60_pos,
				 "Affy650" = affy_pos,
				 "GGP-HD" = ggphd_pos,
				 "GGP-LD" = ggpld_pos)
```

Save marker names for each platform as well. Note that `ggphd_map` does
not contain marker names.


```r
marker_list <- list("SNP60" = rownames(snp60_map),
					"Affy650" = as.character(affy_map$Probe.Set.ID),
					"GGP-LD" = rownames(org_LowD_chip))
```

#### Save physical positions


```r
save(pos_list, marker_list, file = "../pos_list.RData")
```

Plot venn diagram to visualize overlap between all 4 platforms


```r
venn <- venn.diagram(pos_list,
					 filename = NULL,
					 fill = c("red", "blue", "green", "yellow"),
					 alpha = 0.6,
					 cat.fontfamily = "Helvetica",
					 fontfamily = "Helvetica",
					 cat.cex = 2,
					 cex = 1.8)
plot.new()
grid.draw(venn)
```

![plot of chunk mapcompare](figure/mapcompare-1.tiff)

### Objective Three
Manipulate data from `affy_map` to prepare input for Variant Effect Predictor


```r
affy_vep <- data.frame(affy_map$Chromosome,
					   affy_map$Physical.Position,
					   affy_map$Physical.Position,
					   paste(affy_map$Allele.A, affy_map$Allele.B, sep = '/'),
					   affy_map$Strand,
					   affy_map$Affy.SNP.ID)
colnames(affy_vep) <- NULL
```

#### Write VEP input to disk


```r
write.table(affy_vep,
			file = "../affy.vep",
			quote = FALSE,
			sep = '\t',
			col.names = FALSE,
			row.names = FALSE)
```

> Not shown: the resulting file written to disk is used as input for VEP. The assembly used is
> Sscrofa10.2, transcript database used is 'Ensemble transcripts', all identifyers and extra
> options selected, SIFT provides prediction and score, conservation assessed using BLOSUM62

## Conclusion
1. The Affymetrix 650K Swine SNP Chip contains 7 SNPs within chr8:43068687-43916646 (the same SNPs
used by [`breedTools`](https://github.com/funkhou9/breedTools)). It also contains 7 SNPs in this
region that are present on the GGP-HD Chip. Therefore, individuals genotyped on the Affymetrix
650K chip can be used to assess accuracy of imputing "KIT SNPs" from "GGP-HD KIT SNPS" and vice versa.
2. The 4 platforms (GGP-HD, GGP-LD, SNP60, Affy650) contain a decent amount of overlap across the
genome. See [Venn diagram above](#objective-two).
