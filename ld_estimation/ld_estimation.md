[back](../README.md)

# ld_estimation

> Scripts to phase genotypes within each breed, used to estimate
> high-resolution linkage-disequilibrium across the genome

[**1-phase_genotypes.R**](./scripts/1-phase_genotypes_literate/1-phase_genotypes.md)
Using wrappers in `snpTools`, phase genotypes with FImpute software for all chromosomes
within each breed and store the result for downstream LD calculations.

[**2-process_haplotypes.R**](./scripts/2-process_haplotypes_literate/2-process_haplotypes.md)
Loads haplotypes phased in the previous script and saves results as a data.frame
to disk
