#!/usr/bin/env Rscript
#PBS -N 2d-landrace_ld_calculation
#PBS -l nodes=1:ppn=1,walltime=00:30:00,mem=10Gb
#PBS -j oe
#PBS -o ../output/2d-landrace_ld_calculation
#PBS -m abe
#PBS -t 1-20
#PBS -M funkhou9@msu.edu

#' # Script: 2d-landrace_ld_calculation
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 2017-03-23*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Analysis](#analysis)
#'   - [Prep input files](#prep-input-files)
#'   - [Write control file](#write-control-file)
#'   - [Submit executable](#submit-executable)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#' Using software provided by The University of Guelph, calculate genome-wide
#' LD using phased genotypes from `1-phase_genotypes`. A control file will
#' be written then used in the software.
#'

#' ## Install libraries
library(devtools)
library(magrittr)

# #' ## Analysis
# #' ### Prep input files
# #' FImpute output will need to be modified to serve as input for snp1101.
# #' Namely, snp1101 cannot load/process the full dataset. Genotypes will
# #' need to be broken into seperate chromosomes for each breed. Use
# #' `snp_info.txt` to subset genotypes.
job_num <- Sys.getenv("PBS_ARRAYID")

geno <- read.table("../Landrace_fimpute_run/genotypes_imp.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   colClasses = rep("character", 3))
geno <- geno[, c(1, 3)]

snps <- read.table("../Landrace_fimpute_run/snp_info.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE)
snps_subset <- snps[snps$Chr == job_num, -4]
sub <- snps[snps$Chr == job_num, "chip_1"]

geno_subset <- vapply(geno[, 2],
                      function(x) {
                        substr(x, sub[1], sub[length(sub)])
                      },
                      FUN.VALUE = character(1)) %>%
                 cbind("ID" = geno[, 1],
                       "Calls..." = .)

write.table(snps_subset,
            file = paste0("../Landrace_fimpute_run/snp_info_",
                          job_num,
                          ".txt"),
            sep = '\t',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(geno_subset,
            file = paste0("../Landrace_fimpute_run/genotypes_imp_",
                          job_num,
                          ".txt"),
            sep = '\t',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

#' ### Write control file
if (as.numeric(job_num) < 10) {
  out_num <- paste0("0", job_num)
} else {
  out_num <- job_num
}

sink(paste0("../ld_Landrace_", job_num, ".ctr"))
cat(
  paste0(
'title
  "LD Landrace";

gfile
  "../Landrace_fimpute_run/genotypes_imp_', job_num, '.txt"
  skip 1;

mapfile
  "../Landrace_fimpute_run/snp_info_', job_num, '.txt"
  skip 1;

ld pair_wise
  hap
  maf_range 0.01 0.5
  dist_range 1 5100000
  plot mean;

nthread
  10;

output_folder
  "../Landrace_snp1101_', out_num, '_run";\n'
  )
)
sink()

#' ### Submit executable
system(paste0("echo ../ld_Landrace_", job_num, ".ctr | ../src/snp1101_Scott170203/snp1101"))

#' Remove control file, which is no longer needed since it is stored automatically
#' in the output file
system(paste0("rm ../ld_Landrace_", job_num, ".ctr"))
