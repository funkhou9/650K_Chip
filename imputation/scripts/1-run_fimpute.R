#!/usr/bin/env Rscript
#PBS -N 1-run_fimpute
#PBS -l nodes=1:ppn=1,walltime=00:30:00,mem=10Gb
#PBS -j oe
#PBS -o ../output/log/1-run_fimpute
#PBS -m abe
#PBS -t 1-90
#PBS -M funkhou9@msu.edu

#' # Script: 1-run_fimpute
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-05-15*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [imputation](../../imputation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/imputation/scripts")

#' ## Objectives
#' For each of the 90 animals genotyped on the Affymetrix 650K Chip, mask all
#' SNPs except those also present on the Geneseek Genomic Profiler HD, then run
#' FImpute to obtain imputed genotypes.
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

#' Process PBS ARRAY ID
job_num <- as.numeric(Sys.getenv("PBS_ARRAYID"))

#' ## Load data
#' Load Yorkshire 650K Genotypes
load("../../genotype_analysis/1-data_prep.RData")

#' Load SNPs and physical positions from both Affy and Illumina platforms
load("../../snp_inspection/pos_list.RData")

#' ## Analysis
#' Filter Yorkshire animals (remove animals of other breeds) and mask the
#' 650K-specific genotypes for the animal corresponding to the current job
#' number
yorks <- affy_geno[paste0(array_ids$Yorkshire, ".CEL"), ]
matched_pos <- intersect(pos_list$Affy650, pos_list[["GGP-HD"]])
matched_markers <- marker_list$Affy650[pos_list$Affy650 %in% matched_pos]
yorks[job_num, !(colnames(yorks) %in% matched_markers)] <- NA

#' Prep additional inputs for FImpute, including Affymetrix map
#' Load SNP map provided by Affymetrix
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

#' FImpute only accepts sample names that are 30 characters or less. Sample names
#' are identical for the first 23 characters, so we can modify sample names
#' keeping only the unique part of each name.
rownames(yorks) <- substr(rownames(yorks), 24, 30)

#' Remove markers on unassigned contigs
affy_map <- affy_map[affy_map$chr %in% c(as.character(1:18), "X", "Y"), ]

#' Prep output file and run FImpute. Note that fimpute_run() should be executed
#' inside the dedicated output file since fimpute_run() will generate several
#' input files in the working directory that will be used for FImpute software.
out_dir <- paste0("../output/impute_", job_num)
dir.create(out_dir)
setwd(out_dir)
snpTools::fimpute_run(geno = yorks,
                      map = affy_map,
                      path = "/mnt/research/pigsnp/bin/FImpute_Linux")
