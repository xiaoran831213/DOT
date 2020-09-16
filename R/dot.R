#' Decorrelation by Orthogonal Transformation (DOT)
#'
#' \code{\link{dot}}  decorrelates   genetic  association  test   statistics  by
#' special symmetric orthogonal transformation.
#'
#' @details
#' Genetic  association studies  typically provide  per-variant test  statistics
#' that can be  converted to asymptotically normal, signed  Z-scores. Once those
#' Z-scores are transformed to independent random variables, various methods can
#' be applied to combine them and obtain SNP-set overall association.
#' 
#' \code{\link{dot}} uses  per-variant genetic association test  statistics and
#' the correlation among them to decorrelate Z-scores.
#'
#' To estimate the  correlation among genetic association  test statistics, use
#' \code{\link{cst}}.    If  P-values   and   estimated   effects  (i.e,   beta
#' coefficients) are given instead of test statistics, \code{\link{zsc}} can be
#' used to recover the test statistics (i.e., Z-scores).
#'
#' XT:
#' A number statistics that combines de-correlated P-values are made available,
#' see \code{\link{dot_sst}} for detailed.
#'
#' For details about DOT, see the reference below.
#'
#' @references
#' \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819}{
#' Vsevolozhskaya, O. A., Shi, M., Hu, F., & Zaykin, D. V. (2020). DOT: Gene-set
#' analysis by combining decorrelated association statistics. PLOS Computational
#' Biology, 16(4), e1007819.}
#'
#' @param Z vector of association test statistics (i.e., Z-scores)
#' @param C matrix of correlation among the association test statistics, as
#'     obtained by \code{\link{cst}}
#' @param ... additional parameters
#'
#' @return
#' a  list  containing the  association  statistics  \code{Z}, its  decorrelated
#' counterpart  \code{X}, and  the  orthogonal  transformation matrix  \code{W},
#' where \code{X == WZ}.
#' 
#' @seealso \code{\link{cst}}, \code{\link{zsc}}, \code{\link{dot_sst}}
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## estimate the correlation among association test statistics
#' sgm <- cst(gno, cvr)$C
#'
#' ## get the result of genetic association analysis (P-values and effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover Z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## decorrelate Z-scores by DOT
#' result <- dot(stt, sgm)
#' print(result$X)          # decorrelated statistics
#' print(result$W)          # orthogonal transformation
#'
#' ## sum of squares of decorrelated statistics is a chi-square
#' ssq <- sum(result$X^2)
#' pvl <- 1 - pchisq(ssq, df=length(stt))
#'
#' print(ssq)            # sum of square = 35.76306
#' print(pvl)            # chisq P-value =  0.001132132
#' @export
dot <- function(Z, C=NULL, ...)
{
    if(is.null(C))
        C <- diag(length(Z))
    W <- nsp(C)                         # orthogonal transformation
    X <- W %*% Z                        # decorrelated statistics
    list(Z=Z, W=W, X=X)
}

#' Calculate Z-scores from P-values and estimated effects
#'
#' \code{\link{zsc}} recovers  Z-scores from  P-values and  corresponding effect
#' directions (or beta coefficients) reported by a genetic association analysis.
#'
#' @details
#'
#' For any  genetic variant,  its two-sided  P-value (\eqn{p})  and the  sign of
#' estimated effect (\eqn{\beta}) is used to recover the Z-score (\eqn{z}), that
#' is, \eqn{z = sign(\beta) \Phi^{-1}(p/2)}.
#'
#' @seealso \code{\link{dot}}
#'
#' @param P vector of P-values.
#' @param BETA vector of effect directions or beta coefficients.
#' @return A vector of Z-scores.
#' @examples
#' ## result of per-variant analysis (P-values and estimated effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover Z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## checking
#' stopifnot(all.equal(pnorm(abs(stt), lower.tail = FALSE) * 2, res$P))
#' 
#' @export
zsc <- function(P, BETA)
{
    qnorm(P / 2) * (sign(BETA))
}
