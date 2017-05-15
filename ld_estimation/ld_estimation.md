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
less than 5Mb (scripts 2a-2d), mean and sd LD estimates are calculated for 1KB
windows from 0KB to 50KB. This script is parameterized by PBS Array and has no
literate document.

[**4-visualize_global_ld.R**](./scripts/4-visualize_global_ld_literate/4-visualize_global_ld.md)
Plotting estimates obtained in `3-process_ld_1kb.R` to visualize global patterns
of short range (0-50KB) LD.

[**5-process_ld_100kb.R**](./scripts/5-process_ld_100kb.R)
New 100Kb windows around distances of 0.5Mb, 1Mb, and 5Mb are used for mean
LD calculations.

[**6-tabulate_ld_100kb_bins.R**](./scripts/6-tabulate_ld_100kb_bins_literate/6-tabulate_ld_100kb_bins.md)
Tabulate estimates obtained in `5-process_ld_100kb.R`.

[**7-process_pofphase.R**](./scripts/7-process_pofphase.R)
PBS Array parameterized script to calculate correlation of phase for breed
pairs using 1KB windows from 0-50KB.

[**8-process_pofphase_windows.R**](./scripts/8-process_pofphase_windows.R)
PBS Array parameterized script to calculate correlation of phase for breed
pairs using specified windows of various sizes.

[**9-visualize_tabulate_correlations.R**](./scripts/9-visualize_tabulate_correlations_literate/9-visualize_tabulate_correlations.md)
Process results obtained in `7-process_pofphase.R` and `8-process_pofphase_windows.R`
to create tables and plots for publication.
