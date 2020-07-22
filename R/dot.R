#' Decorrelation by Orthogonal Transformation
#'
#' \code{\link{dot}}  decorrelates   genetic  association  test   statistics  by
#' symmetric orthogonal transformation.
#'
#' @details 
#' The genetic association test statistics  (\code{Z}) are typically the variant
#' Z-scores  reported   by  a  genome   wide  association  study   (GWAS).   The
#' decorrelated statistics can be seen  as independent normals.  One could thus,
#' for  example, test  for the  overall genetic  effect by  treating the  sum of
#' squares of decorrelated statistics as a chi-square.
#'
#' \code{\link{dot}} uses  existing genetic association test  statistics and the
#' correlation among these statistics, it does not require original genotype and
#' phenotype that gaven rise to the said statistics.
#'
#' To estimate  the correlation among  genetic association test  statistics, use
#' \code{\link{cst}}.   If  the  p-values   and  estimated  effects  (i.e,  beta
#' coefficients) are given instead of  test statistics, use \code{\link{zsc}} to
#' recover the test statistics (i.e., Z-scores).
#'
#' For details  of DOT and  its impact on  statistical power, see  the reference
#' below.
#'
#' @references
#' \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819}{DOT:
#' Gene-set analysis by combining decorrelated association statistics}
#'
#' @param Z vector of association test statistics (i.e., Z-scores)
#' @param C matrix of correlation among the association test statistics, as
#'     obtained by \code{\link{cst}}
#' @param ... additional parameters
#' @return
#' a  list  containing the  vector  of  association statistics  (\code{Z}),  its
#' decorrelated counterpart (\code{X}), and the orthogonal transformation matrix
#' (\code{W}), where \eqn{X = WZ}.
#' 
#' @seealso \code{\link{cst}}, \code{\link{zsc}}
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## estimate the correlation among association test statistics
#' sgm <- cst(gno, cvr)
#'
#' ## get the result of genetic association analysis (p-values and effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## decorrelate z-scores by DOT
#' rpt <- dot(stt, sgm)
#' print(rpt$X)          # decorrelated statistics
#' print(rpt$W)          # orthogonal transformation
#'
#' ## sum of squares of decorrelated statistics is a chi-square
#' ssq <- sum(rpt$X^2)
#' pvl <- 1 - pchisq(ssq, df=length(stt))
#'
#' print(ssq)            # sum of square = 35.76306
#' print(pvl)            # chisq p-value =  0.001132132
#' @export
dot <- function(Z, C=NULL, ...)
{
    if(is.null(C))
        C <- diag(length(Z))
    W <- nsp(C)                         # orthogonal transformation
    X <- W %*% Z                        # decorrelated statistics
    list(Z=Z, W=W, X=X)
}

#' Z-scores from p-values and estimated effects
#'
#' \code{\link{zsc}} recovers  Z-scores from the p-values  and estimated effects
#' reported by a genetic association analysis.
#'
#' @details
#' In reality, for any genetic variant, it is the p-value (\eqn{p}) and the sign
#' of  estimated  effect  (\eqn{\beta})  being   used  to  recover  the  Z-score
#' (\eqn{z}), that is, \eqn{z = sign(\beta) \Phi^{-1}(p/2)}.
#'
#' @seealso \code{\link{dot}}
#'
#' @param P vector of p-values
#' @param BETA vector of effects (or the signs of effects)
#' @return vector of Z-scores
#' @examples
#' ## result of per-variant analysis (p-values and estimated effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## checking
#' stopifnot(all.equal(pnorm(abs(stt), lower.tail = FALSE) * 2, res$P))
#' 
#' @export
zsc <- function(P, BETA)
{
    stats::qnorm(P / 2) * -(sign(BETA))
}
