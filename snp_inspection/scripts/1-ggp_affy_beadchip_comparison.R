#' # Script: 1-ggp_affy_beadchip_comparison.R
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20160504*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [snp_inspection](../../snp_inspection.md)*
#'
#' ## Table of Contents
#'
#' 1. [Background](#background)
#' 2. [Objectives](#objectives)
#' 3. [Install libraries](#install-libraries)
#' 4. [Load data](#load-data)
#' 5. [Analysis](#analysis)
#'   - [Objective One](#objective-one)
#'   - [Objective Two](#objective-two)
#'     - [Save physical positions](#save-physical-positions)
#'   - [Objective Three](#objective-three)
#' 	   - [Write VEP input to disk](#write-vep-input-to-disk)
#' 6. [Conclusion](#conclusion)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/snp_inspection/scripts")

#' ## Background
#' The Illumina Porcine BeadChip60 is being replaced by the GGP-HD, which contains
#' a similar number of SNPs (60K) that aren't all the same as the ones on the BeadChip60. A subset of
#' BeadChip60 SNPs "KIT SNPs" are used in the package [`breedTools`](https://github.com/funkhou9/breedTools)
#' to estimate KIT-based breed probabilities. It was discovered that not all of these SNPs are present
#' on the GGP-HD, instead other SNPs are present on that chip in similar positions ("GGP-HD KIT SNPs").
#' In order to use both animals genotyped on the older BeadChip60 and GGP-HD as reference animals in
#' KIT-based breed probabilities, we'd need to impute SNPs between platforms. To determine imputation
#' accuracy, we would need to have animals that are genotyped on all "KIT SNPs" and "GGP-HD KIT SNPs".
#'
#' Separately from the need to characterize the Affymetrix 650K chip for the purposes of Kit-based
#' breed probability work, we wish to characterize the SNPs present on the affy by comparing to those
#' on the GGP-HD, GGP-LD, and Porcine SNP60 by position.

#' ## Objectives
#'
#' 1. Determine if "KIT SNPs" and "GGP-HD KIT SNPs" are all present on the Affy 650K chip.
#' 2. Determine overlap of SNPs between GGP-HD, GGP-LD, SNP60, and Affy650 using physical position.
#' 3. Use [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)
#' to annotate Affy650 Chip. Save VEP results to where raw 650K data is stored:
#' `/mnt/research/pigsnp/raw_data/affymetrix_hd/`
#'

#' ## Install libraries
library(devtools)
library(magrittr)
library(methods)
library(VennDiagram)

#' ## Load data
#' Paste in KIT snp information as a data.frame. Information comes from GGP-LD/SNP60 map.
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

#' Read in Affy 650K SNP map
affy_map <- read.csv("/mnt/research/pigsnp/raw_data/affymetrix_hd/Axiom_PigHD_v1_Annotation.r1.1.csv",
					 header = TRUE,
					 comment.char = "#")

#' Read in GGP-HD SNP map
ggphd_map <- read.table("/mnt/research/pigsnp/raw_data/80K_SNP_Map.txt", skip = 1)

#' Read in GGP-LD SNP map, variable name `org_LowD_chip`
load("/mnt/research/pigsnp/raw_data/SF_working_dir/LowD_chip_maps.RData")

#' Read in Porcine SNP60 map
snp60_map <- read.table("/mnt/research/pigsnp/raw_data/SF_working_dir/map_updated_commercial.txt")

#' ## Analysis
#' ### Objective One
#' View "KIT SNPs" within affymetrix 650 panel. Are all 7 SNPs present? Answer: yes but with completely
#' different names.
affy_map[affy_map$Chromosome %in% kit_snps$chr & affy_map$Physical.Position %in% kit_snps$pos, 1:5]

#' What SNPs are present within the KIT SNP range (chr8:43068687-43916646) on the GGP-HD? Answer:
#' there are 8 SNPs in this region, 3 of which are part of the original BeadChip60 (ALGA0047809,
#' ALGA0102731, and ALGA0123881)
ggp_kit <- ggphd_map[ggphd_map[, 3] == 8 &
				(ggphd_map[, 4] >= 43068687 & ggphd_map[, 4] <= 43916646), ]

ggp_kit

#' Are these 8 "GGP KIT" SNPs on the affymetrix 650K? Answer: 7 of the 8 GGP KIT SNPs are present on
#' the Affy 650K
affy_map[affy_map$Chromosome %in% ggp_kit[, 3] & affy_map$Physical.Position %in% ggp_kit[, 4], 1:5]

#' ### Objective Two
#' To compare physical positions across platforms, first combine chromosome and position of each SNP
#' into a single character (formatted `chr:pos`) which will provide an easier identifyer to work with.
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

#' Save marker names for each platform as well. Note that `ggphd_map` does
#' not contain marker names.
marker_list <- list("SNP60" = rownames(snp60_map),
					"Affy650" = as.character(affy_map$Probe.Set.ID),
					"GGP-LD" = rownames(org_LowD_chip))

#' #### Save physical positions
save(pos_list, marker_list, file = "../pos_list.RData")

#' Plot venn diagram to visualize overlap between all 4 platforms
#+ mapcompare, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
venn <- venn.diagram(x = list("SNP60" = pos_list[[1]],
                              "PigHD" = pos_list[[2]],
										          "GGP-HD" = pos_list[[3]],
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

#' ### Objective Three
#' Manipulate data from `affy_map` to prepare input for Variant Effect Predictor
affy_vep <- data.frame(affy_map$Chromosome,
					   affy_map$Physical.Position,
					   affy_map$Physical.Position,
					   paste(affy_map$Allele.A, affy_map$Allele.B, sep = '/'),
					   affy_map$Strand,
					   affy_map$Affy.SNP.ID)
colnames(affy_vep) <- NULL

#' #### Write VEP input to disk
write.table(affy_vep,
			file = "../affy.vep",
			quote = FALSE,
			sep = '\t',
			col.names = FALSE,
			row.names = FALSE)

#' > Not shown: the resulting file written to disk is used as input for VEP. The assembly used is
#' > Sscrofa10.2, transcript database used is 'Ensemble transcripts', all identifyers and extra
#' > options selected, SIFT provides prediction and score, conservation assessed using BLOSUM62
#'

#' ## Conclusion
#' 1. The Affymetrix 650K Swine SNP Chip contains 7 SNPs within chr8:43068687-43916646 (the same SNPs
#' used by [`breedTools`](https://github.com/funkhou9/breedTools)). It also contains 7 SNPs in this
#' region that are present on the GGP-HD Chip. Therefore, individuals genotyped on the Affymetrix
#' 650K chip can be used to assess accuracy of imputing "KIT SNPs" from "GGP-HD KIT SNPS" and vice versa.
#' 2. The 4 platforms (GGP-HD, GGP-LD, SNP60, Affy650) contain a decent amount of overlap across the
#' genome. See [Venn diagram above](#objective-two).
