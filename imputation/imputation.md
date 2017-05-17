[back](../README.md)

# imputation

> Scripts to impute SNPs present on the 650K Chip, using SNPs present on the
> 60K Beadchip, and estimate imputation accuracy.

[**1-run_fimpute.R**](./scripts/1-run_fimpute.R)
PBS Array parameterized script to individually mask Affymetrix-specific
SNPs on one Yorkshire animal at a time and impute those SNPs using all
remaining Yorkshire animals.

[**2-imputation_accuracy.R**](./scripts/2-imputation_accuracy_literate/2-imputation_accuracy.md)
Inspect SNP-wise and animal-wise imputation accuracy using data generated
in `1-run_fimpute.R`
