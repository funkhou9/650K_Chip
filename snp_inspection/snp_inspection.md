[back](../README.md)

# snp_inspection

> Inspect the SNPs present on 650K chip. Some SNPs may be especially useful if found across
> several different platforms, including the GGP-LD, GGP-HD, and SNP60 Beadchip

[**1-ggp_affy_beadchip_comparison.R**](./scripts/1-ggp_affy_beadchip_comparison_literate/1-ggp_affy_beadchip_comparison.md) Comparison of SNPs present within *chr8:43068687-43916646* for use in
"KIT-based breed probability" estimation as implemented in 
[`breedTools`](https://github.com/funkhou9/breedTools). Additionally, using physical positions,
a comparison is conducted between the GGP-LD, GGP-HD, Porcine SNP60, and Affy650K. Finally,
this script saves data to disk `affy.vep`, which will be used for Variant Effect Predictor to
annotate SNPs present on Affy650K chip. Annotation will be saved in the same location as raw
Affy650K data: `/mnt/research/pigsnp/raw_data/affymetrix_hd/`
