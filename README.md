# DOT
Decorrelation by Orthogonal Transformation (DOT) on genetic association test statistics.

# Overview
DOT  means  Decorrelation  by  Orthogonal Transformation  (DOT).   The  package
provides function to  decorrelate genetic association test  statistics.  One can
treat the output as independent normal variables.  Because DOT only requires the
association test  statistics and  their estimated  correlation, it  bypasses the
logistical issue of accessing genotype and phenotype.

See  published  work  "[DOT:  Gene-set  analysis  by  combining  decorre-  lated
association statis- tics](dot)" for more details.

[dot]:https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819&rev=2

# Usage

Users inputs  a series  of genetic  association test  statistics of  M variants,
usually the p-values and estimated effects  taken from a GWAS, and the estimated
LD (linkage disequilibrium) among the variants, the package outputs docorrelated
statistics which can combine into a joint statistics to test the overall genetic
effect hypothesis (H_0: none of the M variants is effective).

Here is an example,
```{r}
library(dotgen)
gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))

## estimate the correlation among association test statistics
sgm <- css(gno, cvr)

## get the result of genetic association analysis (p-values and effects)
res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))

## recover association test statistics (z-scores)
stt <- with(res, zsc(P, BETA))

## decorrelate z-scores by DOT
rpt <- dot(stt, sgm)
print(rpt$X)                            # decorrelated statistics
print(rpt$W)                            # orthogonal transformation

## sum of squares of decorrelated statistics is a chi-square
ssq <- sum(rpt$X^2)
pvl <- 1 - pchisq(ssq, df=length(stt))
print(ssq)                              # sum of square = 35.76306
print(pvl)                              # chisq p-value =  0.001132132
```
