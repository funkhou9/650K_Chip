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
on the GGP-HDv2, GGP-LD, and Porcine SNP60 by position.
## Objectives

1. Determine if "KIT SNPs" and "GGP-HD KIT SNPs" are all present on the Affy 650K chip.
2. Determine overlap of SNPs between GGP-HD, GGP-LD, SNP60, and Affy650 using physical position.
3. Use [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)
to annotate Affy650 Chip. Save VEP results to where raw 650K data is stored:
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
					           comment.char = "#",
					           stringsAsFactors = FALSE)
```

Read in GGP-HDv2 SNP map


```r
ggphd_map <- read.table(paste0("/mnt/research/pigsnp/raw_data/ggp_hdv2/",
                               "GGP_Porcine_HD_Public_E_StrandReport_FDT.txt"),
												header = TRUE,
												stringsAsFactors = FALSE)
```

Read in GGP-LD SNP map, variable name `org_LowD_chip`


```r
load("/mnt/research/pigsnp/raw_data/SF_working_dir/LowD_chip_maps.RData")
```

Read in Porcine SNP60 map


```r
snp60_map <- read.table("/mnt/research/pigsnp/raw_data/SF_working_dir/map_updated_commercial.txt",
                        stringsAsFactors = FALSE)
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

What SNPs are present within the KIT SNP range (chr8:43068687-43916646) on the GGP-HD?


```r
ggp_kit <- ggphd_map[ggphd_map[, 5] == 8 &
				(ggphd_map[, 6] >= 43068687 & ggphd_map[, 6] <= 43916646), ]

ggp_kit
```

```
##       Index           SNP_Name                              Ilmn_ID Build
## 541     541     Affx-114629526      Affx-114629526-1_B_F_2353936408  10.2
## 5195   5195        ALGA0047809         ALGA0047809-0_B_R_2348767067  10.2
## 11709 11709        ALGA0123881 ALGA0123881_IlmnDup-0_B_F_2354426544  10.2
## 47809 47809 WU_10.2_8_43277756  WU_10.2_8_43277756-0_B_R_2348821553  10.2
## 47811 47811 WU_10.2_8_43364095  WU_10.2_8_43364095-0_B_F_2348821555  10.2
## 47813 47813 WU_10.2_8_43641500  WU_10.2_8_43641500-0_B_R_2348821559  10.2
##       Chr    Coord
## 541     8 43651639
## 5195    8 43916646
## 11709   8 43068687
## 47809   8 43277756
## 47811   8 43364095
## 47813   8 43641500
##                                                                                                                                                                                                         Forward_Seq
## 541   AGGTAGAGGGAAGAGGAGACAAATAGAGAAATGAAGGGAAGTAGAGAAGAGAGAAGAGAACACATTACCACCAGGAGGGCTCTGTGTGCCAGACCAAGAT[T/G]TGTGCTTTGTGCTAATATATTCTTGCAGCTTTAATATTTATTCATTCATTCATTGCTTTTTAGGGCCACGCCTGTAGCATATGGACGTTCCCAGGCTAGG
## 5195                                                                                  TTTACTGTAATCCTTGTATTATCAATAGTCTTCCTTCAAGTATTTTCCTCCTAACTCAGA[A/G]GGCTTATCTCTGAGTCAGACCTCAGCTGAGCCTGGTATGTTCTACTACTTTAGCTCCGGA
## 11709                                                                                 TCATTGGCCAATCCGGCTTCAGTCTTCACACTTAGCACGAAAAAGGTGGGGCTTGCATTG[T/C]GATTACACSTGTAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
## 47809                                                                                 GCAAAAAGACAGAAAAAAAAAAAAAAAACAACTCTGTTCTATTGGAAACTCAGGGCACCT[A/G]TCCTTTGGTGTGGTCAAAGACAGTGTCTGTCCTCTCACCCACTCTAGGAGGCAACAGCAA
## 47811                                                                                 GCAGAACTCTGGTTTTGTTAAAGTGGCAGTGCATATAGCCCAGGCGATGATCACCATTGA[T/C]CTAAGCCAGACATGTCAATCTCCTTCCTCTTTACCAGTGATTAATTATGGGCATCTGACC
## 47813                                                                                 TCCTGTGATATAGTTATTATGAAGTTCATTTAGCCAGATGAAGACYGAGACCTAGAGAAG[A/G]GGTTGAATCCATTTCTATAGGTCACTTAGCAACTAGGTGGTTGATAAGGAATCAAACTTA
##       Forward_Allele1 Forward_Allele2
## 541                 T               G
## 5195                A               G
## 11709               T               C
## 47809               A               G
## 47811               T               C
## 47813               A               G
##                                                                                                                                                                                                          Design_Seq
## 541   AGGTAGAGGGAAGAGGAGACAAATAGAGAAATGAAGGGAAGTAGAGAAGAGAGAAGAGAACACATTACCACCAGGAGGGCTCTGTGTGCCAGACCAAGAT[T/G]TGTGCTTTGTGCTAATATATTCTTGCAGCTTTAATATTTATTCATTCATTCATTGCTTTTTAGGGCCACGCCTGTAGCATATGGACGTTCCCAGGCTAGG
## 5195                                                                                  TCCGGAGCTAAAGTAGTAGAACATACCAGGCTCAGCTGAGGTCTGACTCAGAGATAAGCC[T/C]TCTGAGTTAGGAGGAAAATACTTGAAGGAAGACTATTGATAATACAAGGATTACAGTAAA
## 11709                                                                                 TCATTGGCCAATCCGGCTTCAGTCTTCACACTTAGCACGAAAAAGGTGGGGCTTGCATTG[T/C]GATTACACSTGTAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
## 47809                                                                                 TTGCTGTTGCCTCCTAGAGTGGGTGAGAGGACAGACACTGTCTTTGACCACACCAAAGGA[T/C]AGGTGCCCTGAGTTTCCAATAGAACAGAGTTGTTTTTTTTTTTTTTTTCTGTCTTTTTGC
## 47811                                                                                 GCAGAACTCTGGTTTTGTTAAAGTGGCAGTGCATATAGCCCAGGCGATGATCACCATTGA[T/C]CTAAGCCAGACATGTCAATCTCCTTCCTCTTTACCAGTGATTAATTATGGGCATCTGACC
## 47813                                                                                 TAAGTTTGATTCCTTATCAACCACCTAGTTGCTAAGTGACCTATAGAAATGGATTCAACC[T/C]CTTCTCTAGGTCTCRGTCTTCATCTGGCTAAATGAACTTCATAATAACTATATCACAGGA
##       Design_Allele1 Design_Allele2
## 541                T              G
## 5195               T              C
## 11709              T              C
## 47809              T              C
## 47811              T              C
## 47813              T              C
##                                                                                                                                                                                                             Top_Seq
## 541   CCTAGCCTGGGAACGTCCATATGCTACAGGCGTGGCCCTAAAAAGCAATGAATGAATGAATAAATATTAAAGCTGCAAGAATATATTAGCACAAAGCACA[A/C]ATCTTGGTCTGGCACACAGAGCCCTCCTGGTGGTAATGTGTTCTCTTCTCTCTTCTCTACTTCCCTTCATTTCTCTATTTGTCTCCTCTTCCCTCTACCT
## 5195                                                                                  TTTACTGTAATCCTTGTATTATCAATAGTCTTCCTTCAAGTATTTTCCTCCTAACTCAGA[A/G]GGCTTATCTCTGAGTCAGACCTCAGCTGAGCCTGGTATGTTCTACTACTTTAGCTCCGGA
## 11709                                                                                 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTACASGTGTAATC[A/G]CAATGCAAGCCCCACCTTTTTCGTGCTAAGTGTGAAGACTGAAGCCGGATTGGCCAATGA
## 47809                                                                                 GCAAAAAGACAGAAAAAAAAAAAAAAAACAACTCTGTTCTATTGGAAACTCAGGGCACCT[A/G]TCCTTTGGTGTGGTCAAAGACAGTGTCTGTCCTCTCACCCACTCTAGGAGGCAACAGCAA
## 47811                                                                                 GGTCAGATGCCCATAATTAATCACTGGTAAAGAGGAAGGAGATTGACATGTCTGGCTTAG[A/G]TCAATGGTGATCATCGCCTGGGCTATATGCACTGCCACTTTAACAAAACCAGAGTTCTGC
## 47813                                                                                 TCCTGTGATATAGTTATTATGAAGTTCATTTAGCCAGATGAAGACYGAGACCTAGAGAAG[A/G]GGTTGAATCCATTTCTATAGGTCACTTAGCAACTAGGTGGTTGATAAGGAATCAAACTTA
##       Top_AlleleA Top_AlleleB
## 541             A           C
## 5195            A           G
## 11709           A           G
## 47809           A           G
## 47811           A           G
## 47813           A           G
```

Are these "GGP KIT" SNPs on the affymetrix 650K?


```r
affy_map[affy_map$Chromosome %in% ggp_kit[, 5] & affy_map$Physical.Position %in% ggp_kit[, 6], 1:5]
```

```
##        Probe.Set.ID    Affy.SNP.ID dbSNP.RS.ID Chromosome
## 26100  AX-116691049 Affx-114940851         ---          8
## 26104  AX-116345697 Affx-114629526         ---          8
## 26108  AX-116345774 Affx-114753388         ---          8
## 95509  AX-116736314 Affx-114850107         ---          8
## 95513  AX-116345643 Affx-114992210         ---          8
## 555768 AX-116345695 Affx-114958854         ---          8
##        Physical.Position
## 26100           43068687
## 26104           43651639
## 26108           43916646
## 95509           43277756
## 95513           43364095
## 555768          43641500
```

### Objective Two
To compare physical positions across platforms, first combine chromosome and position of each SNP
into a single character (formatted `chr:pos`) which will provide an easier identifyer to work with.
Remove SNPs with ambiguous positions (chr=0, pos=0)


```r
snp60_pos <- paste(snp60_map$chr,
	 			           snp60_map$pos,
	 			           sep = ":")
idx <- snp60_pos != "0:0"
snp60_pos <- snp60_pos[idx]
snp60_markers <- rownames(snp60_map)[idx]

affy_pos <- paste(affy_map$Chromosome,
				          affy_map$Physical.Position,
				          sep = ":")
idx <- affy_pos != "0:0"
affy_pos <- affy_pos[idx]
affy_markers <- as.character(affy_map$Probe.Set.ID)[idx]

ggphd_pos <- paste(ggphd_map$Chr,
				           ggphd_map$Coord,
				           sep = ":")
ggphd_pos <- ggphd_pos[ggphd_pos != "0:0"]

ggpld_pos <- paste(org_LowD_chip$chr,
				           org_LowD_chip$pos,
				           sep = ":")
idx <- ggpld_pos != "0:0"
ggpld_pos <- ggpld_pos[idx]
ggpld_markers <- rownames(org_LowD_chip)[idx]

pos_list <- list("SNP60" = snp60_pos,
				         "Affy650" = affy_pos,
				         "GGP-HD" = ggphd_pos,
				         "GGP-LD" = ggpld_pos)
```

Save marker names for each platform as well.


```r
marker_list <- list("SNP60" = snp60_markers,
					          "Affy650" = affy_markers,
					          "GGP-LD" = ggpld_markers)
```

#### Save physical positions


```r
save(pos_list, marker_list, file = "../pos_list.RData")
```

Plot venn diagram to visualize overlap between all 4 platforms


```r
venn <- venn.diagram(x = list("SNP60" = pos_list[[1]],
                              "PigHD" = pos_list[[2]],
										          "GGP-HDv2" = pos_list[[3]],
										          "GGP-LD" = pos_list[[4]]),
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
