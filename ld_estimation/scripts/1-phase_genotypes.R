#' # Script: 1-phase_genotypes
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20161101*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#'
#' 1. Phase genotypes within each breed {Landrace, Yorkshire, Duroc and Landrace}
#' and save the result.
#'

#' ## Install libraries
library(devtools)
library(magrittr)

#' Install [snpTools](https://github.com/funkhou9/snpTools/commit/f3e70629361e023219cea7f5c27f5c56bbd93323)
# install_github("funkhou9/snpTools")
library(snpTools)

#' Install [breedTools](https://github.com/funkhou9/breedTools/commit/00b77d774e31b69f885b3ccbe413f7caf92abbbb)
# install_git("/mnt/research/pigsnp/raw_data/SF_working_dir/breed_compos/breedTools")
library(breedTools)

#' ## Load data
#' Load Affymetrix genotypes from `1-data_prep.R`
load("../../genotype_analysis/1-data_prep.RData")

#' Load SNP map provided by Affymetrix
map <- read.csv(paste0("/mnt/research/pigsnp/raw_data/affymetrix_hd/",
                       "Axiom_PigHD_v1_Annotation.r1.1.csv"),
                header = TRUE,
                comment.char = "#",
                stringsAsFactors = FALSE)

#' ## Analysis
#' Prep `map` to be used as map argument for `fimpute_run()`
affy_map <- data.frame("marker" = map$Probe.Set.ID,
                       "chr" = map$Chromosome,
                       "pos" = map$Physical.Position,
                       row.names = "marker",
                       stringsAsFactors = FALSE)

#' FImpute only accepts sample names that are 30 characters or less. Sample names
#' are identical for the first 23 characters, so we can modify sample names
#' keeping only the unique part of each name.
rownames(affy_geno) <- substr(rownames(affy_geno), 24, 30)

#' Prep `groups` argument passed to `fimpute_run()` so that IDs of each group
#' match newly modified IDs in `affy_geno`
array_ids <- lapply(array_ids, function(x) substr(x, 24, 30))

#' Remove markers on unassigned contigs
affy_map <- affy_map[affy_map$chr %in% c(as.character(1:18), "X", "Y"), ]

#' Phase genotypes with `fimpute_run()`. Stdout will be printed to screen and
#' without returning anything, `fimpute_run()` will save haplotypes to disk
#' in `output_folder`
snpTools::fimpute_run(geno = affy_geno,
                      map = affy_map,
                      groups = array_ids,
                      path = "/mnt/research/pigsnp/bin/FImpute_Linux",
                      output_folder = "..")
