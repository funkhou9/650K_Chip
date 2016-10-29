[back](../README.md)

# genotype_analysis

> Scripts to load 650K data from all animals of each breed, perform any necessary
> data cleaning, inspect genotypes with exploratory plots, and compare genotype
> calls to previous calls made on the 60K SNP chip.

[**1-data_prep.R**](./scripts/1-data_prep_literate/1-data_prep.md)
As the name suggests, loading and cleaning of raw 650K genotyping data provided
by Affymetrix

[**2-PCA.R**](./scripts/2-PCA_literate/2-PCA.md)
Inspect genotype variation between using PCA

[**3-compare_calls.R**](./scripts/3-compare_calls_literate/3-compare_calls.md)
Most Yorkshire animals genotyped on the 650K Affy chip were previously
genotyped on the Illumina SNP60 Beadchip. Isolate those animals and analyze
differences in genotype calls for each SNP in common with the 60K and 650K
chip
