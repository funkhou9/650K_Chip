[back](../README.md)

# ld_estimation

> Scripts to phase genotypes within each breed, used to estimate
> high-resolution linkage-disequilibrium across the genome

[**1-phase_genotypes.R**](./scripts/1-phase_genotypes_literate/1-phase_genotypes.md)
Using wrappers in `snpTools`, phase genotypes with FImpute software for all chromosomes
within each breed and store the result for downstream LD calculations.

[**2a-2d_ld_calculation**](./scripts/2a-yorkshire_ld_calculation.R)
A set of 4 scripts (one for each breed), each one parameterized by PBS array
to estimate LD of each chromosome separately. The yorkshire script is linked
above.

[**3-process_ld_1kb.R**](./scripts/3-process_ld_1kb.R)
Once LD is estimated for all pairwise combinations of SNPs, whose distances are
less than 50KB (scripts 2a-2d), mean and sd LD estimates are calculated for 1KB
windows from 0KB to 50KB. This script is parameterized by PBS Array and has no
literate document.

[**4-visualize_global_ld.R**](./scripts/4-visualize_global_ld_literate/4-visualize_global_ld.md)
Plotting estimates obtained in `3-process_ld_1kb.R` to visualize global patterns
of short range (0-50KB) LD.
