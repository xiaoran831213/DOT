# DOT
DOT - Decorrelation by Orthogonal Transformation, a gene-set analysis by combining decorrelated association statistics

# Introduction
This package provides DOT function that combines GWAS summary statistics of P non-independed variants (i.e., SNPs fall in a gene or a chromosome window) to achieve higher power (i.e., lower p-values) than a Bonferroni style variant scan and the TQ (test by quadratic form) summary statistics .

The DOT 

# Usage
An user inputs a GWAS report (e.g., p-values and beta coefficients) of M genetic variants of interest, and an estimated LD (linkage disquilibrium) between these M variants (a MxM correlation matrix), the package produces DOT summary statistics as the sum of M re-weighted test statistics, and a single p-value for H0: none of the M variants is effective.

# Examples

(TODO)
