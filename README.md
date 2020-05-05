# DOT
DOT - Decorrelation by Orthogonal Transformation, a gene-set analysis by combining decorrelated association statistics

# Introduction
This package provides DOT function that combines GWAS summary statistics of P non-independed variants (i.e., SNPs fall in a gene or a chromosome window) to achieve higher power (i.e., lower p-values) than a Bonferroni style variant scan and the TQ (test by quadratic form) summary statistics .

The DOT 

# Usage
An user supplies the GWAS report (e.g., p-values and beta coefficients) for the P variants of interest, and an estimated LD (linkage disquilibrium) between the variants (a PxP correlation matrix), the package gives the DOT summary statistics as the sum of re-weighted P test statistics, followed by a single p-value for all variants involved.

# Examples

(TODO)
