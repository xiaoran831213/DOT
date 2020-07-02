# DOT
Decorrelation by Orthogonal Transformation (DOT) on genetic association test statistics.

# Introduction
The DOT package provides functions to decorrelate a set of genetic association test statistics (i.e., Z-score of SNP(s) in a gene) via De-correlation by Orthogonal Transformation (DOT). One could jointly analyze the genetic effect of all variants involved by combining the decorrelated statistics; since the DOT only relies on the original association test statistics and the LD of the genotype variants, it bypasses the logistical issue of accessing genotype and phenotype.

# Usage
An user inputs a GWAS report (e.g., p-values and beta coefficients) of M genetic variants of interest, and an estimated LD (linkage disequilibrium) between these M variants (a MxM correlation matrix), the package produces DOT summary statistics as the sum of M re-weighted test statistics, and a single p-value for H0: none of the M variants is effective.

# Examples

(TODO)
