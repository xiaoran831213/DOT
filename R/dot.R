#' DOT: Decorrelation by Orthogonal Transformation
#'
#' DOT decorrelated raw  test statistics and summarize a  chi-square value which
#' usually provides  higher power  than raw statistics  after Bonferroni  or FDR
#' correction.
#'
#' DOT does not need the original  genotype that to generate the raw statistics,
#' and only requires the correlation of variants.
#'
#' @param Z raw test statistics (Z-scores)
#' @param C correlation of test statistics
#' @param ... absorb additional parameters
#' @return the chi-square value, p-value, and other reports.
dot <- function(Z, C=NULL, ...)
{
    if(is.null(C))
        C <- diag(length(Z))

    W <- nsp(C)   # DOT weights
    X <- W %*% Z  # decorrelate
    S <- sum(X^2) # sum square
    
    ## P-values
    P <- 1 - pchisq(S, df=length(Z), ncp=0)
    
    list(mtd='dot', Z=Z, C=C, W=W, X=X, ssq=S, pvl=P)
}

#' Get test statistics from p-values and effect sizes
#'
#' Given the p-values (`p`) and effect  sizes (`beta`) of a double tailed test,
#' recover the test statistics `z = Phi^{-1}(p / 2) sign(beta)`.
#'
#' @param P vector of p-values
#' @param BETA vector of effect sizes
zvl <- function(P, BETA)
{
    qnorm(P / 2) * -(sign(BETA))
}
