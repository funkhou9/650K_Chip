#' # Script: 3-ld_calc
#'
#' - *Author: Scott Funkhouser*
#' - *Date: 20161107*
#' - *Project: [650K_Chip](../../../README.md)*
#' - *Sub Folder: [ld_estimation](../../ld_estimation.md)*
#'
#' ## Table of Contents
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)

setwd("/mnt/research/pigsnp/NSR/650K_Chip/ld_estimation/scripts")

#' ## Objectives
#'
#' 1. For each chromosome, calculate all pairwise measurements of distance and
#' save to disk, one binary file per chromosome.]
#'
#'

#' ## Install libraries
library(devtools)
library(magrittr)

#' Install ogbox, an assortment of R tools from open source developer Ogan Mancarci,
#' containing `latexImg()`, a function to parse a line of latex and return a link
#' to an image renderable in Github markdown.
# install_github("oganm/ogbox")
library(ogbox)

#' ## Load data
#' Load haplotypes phased in `1-phase_genotypes` and processed in `2-process_haplotypes`
load("../2-process_haplotypes.RData")

#' Load Affy map containing physical posistions for each SNP
affy_map <-
    read.csv("/mnt/research/pigsnp/raw_data/affymetrix_hd/Axiom_PigHD_v1_Annotation.r1.1.csv",
			 header = TRUE,
			 comment.char = "#")

#' ## Analysis
#' The amount of LD between SNP i and SNP j are calculated using the following:
#+ echo=FALSE, results='asis'
ogbox::latexImg("r^2_{ij} = \\frac{(p_{ij} - p_ip_j)^2}{p_i(1 - p_i)p_j(1 - p_j)}")

#' Correlation statistics will be obtained (mean and variance) for each set of
#' SNP pairs that fall between breakpoints defined in `breakp`
breakp <- seq(0, 10000000, 100000)

#' Define a list of chromosomes...
chr <- list()

#' Define a list of groups...
groups <- list()
for (i in 1:100) {
    groups[[i]] <- c(breakp[i], breakp[i + 1])
}

#' Define a compiled function `which_in_range()` that will iterate through a long vector
#' finding elements that fall between `start` and `end` integers.

# Will push_back work with a NumericVector object? How do I append?
Rcpp::cppFunction('NumericVector which_in_range(NumericVector x, double start, double end) {
                    int n = x.size();
                    NumericVector out;
                    for (int i = 0; i < n; i++) {
                        if (x[i] >= start && x[i] <= end)
                            out.push_back(x[i]);
                    }
                    return out;
}')

#' Define a compiled function that will repeatedly apply `which_in_range()` for a set of `start`
#' and `end` values, storing the results each time in a list. If the input of a compiled function
#' is a list, then C++ needs to know the type of elements in the list

# How to access elements in a List object??
Rcpp::cppFunction('List pairwise_dist(List g) {
                    if (!g.inherits("numeric")) {
                        stop("Input list needs to be numeric");
                    }
                    int n = g.size();
                    NumericVector x;
                    for (int i = 0; i < n; i++) {
                        x.push_back(which_in_range(x, g[[i]][1], g[[i]][2]))
                    }
                    return
}')

#' For each chromosome:
#'
#' 1. calcualte a pairwise distances for all pairs of SNPs on of that chromosome
#' 2. sort pairs into elements of `groups` (ex. groups[1] will correspond to
#' SNP pairs with a physical distance between 0 and 100kb)
#' 3. append `groups` to `chr`
#'

for (i in c(seq(1:18), "X", "Y")) {
    # identify positions on chromosome i, store as a vector
    pos <- affy_map$Physical.Position[affy_map$Chromosome == i]
    # calculate pairwise distances
    distance <- dist(pos)
    g <- lapply(groups, function(x) {
            which(distance >= x[1] & distance <= x[2])
         })
    chr[[i]] <- g
}

length(chr)
length(chr[[1]])
