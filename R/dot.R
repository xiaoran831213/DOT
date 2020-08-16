#' Decorrelation by Orthogonal Transformation
#'
#' \code{\link{dot}}  decorrelates   genetic  association  test   statistics  by
#' symmetric orthogonal transformation.
#'
#' @details
#' Genetic  association studies  typically provide  per-variant test  statistics
#' that can be  converted to asymptotically normal, signed  Z-scores. Once those
#' Z-scores are transformed to independent random variables, various methods can
#' be applied to combined them and obtain SNP-set overall association.
#' 
#' \code{\link{dot}} uses  existing genetic association test  statistics and the
#' correlation among these statistics to decorrelate Z-scores.
#'
#' To estimate  the correlation among  genetic association test  statistics, use
#' \code{\link{cst}}.    If   p-values   and  estimated   effects   (i.e,   beta
#' coefficients) are given instead of  test statistics, use \code{\link{zsc}} to
#' recover the test statistics (i.e., Z-scores).
#'
#' For details about DOT, see the reference below.
#'
#' @references
#' \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819}{DOT:
#' Gene-set analysis by combining decorrelated association statistics}
#'
#' @param Z vector of association test statistics (i.e., Z-scores)
#' @param C matrix of correlation among the association test statistics, as
#'     obtained by \code{\link{cst}}
#' @param ... additional parameters
#'
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
#' \code{\link{zsc}} recovers  Z-scores from  p-values and  corresponding effect
#' directions (or beta coefficients) reported by a genetic association analysis.
#'
#' @details
#'
#' For any  genetic variant,  its two-sided  p-value (\eqn{p})  and the  sign of
#' estimated effect (\eqn{\beta}) is used to recover the Z-score (\eqn{z}), that
#' is, \eqn{z = sign(\beta) \Phi^{-1}(p/2)}.
#'
#' @seealso \code{\link{dot}}
#'
#' @param P vector of p-values
#' @param BETA vector of effect directions or beta coefficients
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
    stats::qnorm(P / 2) * (sign(BETA))
}


#' De-correlation followed by Chi-square Test
#'
#' Summarize one  p-value from a  series of genetic association  test statistics
#' (i.e., Z-scores, \code{Z}) and their correlation (i.e., LD-values, \code{C}).
#'
#' @details
#'
#' These functions first call \code{\link{dot}(Z, C)} to decorrelate the genetic
#' association test statistics,  then use various techniques  (one per function)
#' to  combine  the  decorrelated  statistics  or the  p-values  into  a  single
#' statistics and one  p-value, which essentially tests the joint  null: none of
#' the variants are associated with the phenotype of interests.
#'
#' The original  association test statistics are  typically normally distributed
#' Z-scores reported by a GWAS.
#'
#' @param Z vector of association test statistics (i.e., Z-scores)
#' @param C matrix of correlation among the association test statistics, as
#'     obtained by \code{\link{cst}}
#' @param  k combine  the  \code{k}  most  significant test  statistics  (i.e.,
#' \code{k} smallist p-values)
#' @param ... additional parameters
#'
#' @return
#' A list with decorrelated association  test statistics \code{X} and orthogonal
#' transformation \code{W}  such that  \eqn{X =  WZ}, the  summarized statistics
#' \code{Y} and p-value \code{P}.
#'
#' @examples
#' ## get the test statistics and pre-calculated LD matrix
#' stt <- readRDS(system.file("extdata", 'art_zsc.rds', package="dotgen"))
#' sgm <- readRDS(system.file("extdata", 'art_ldm.rds', package="dotgen"))
#'
#' ## run de-correlated chi-square test
#' rpt <- dot_chisq(stt, sgm)
#'
#' print(rpt$Y)  # 37.2854
#' print(rpt$P)  #  0.0003736988
dot_chisq <- function(Z, C, ...)
{
    within(dot(Z, C, ...),
    {
        Y <- sum(X^2)
        P <- 1 - pchisq(Y, df=length(Z))
    })
}

#' @describeIn dot_chisq De-correlation followed by Augmented Rank Trucated (ART) Test
#'
#' For details about DOT, see the reference below.
#'
#' @references
#' \href{https://www.frontiersin.org/articles/10.3389/fgene.2019.01051}{Detecting
#' Weak Signals by Combining Small P-Values in Genetic Association Studies}
#'
#' @examples
#' ## get the test statistics and pre-calculated LD matrix
#' stt <- readRDS(system.file("extdata", 'art_zsc.rds', package="dotgen"))
#' sgm <- readRDS(system.file("extdata", 'art_ldm.rds', package="dotgen"))
#'
#' ## run de-correlated chi-square test
#' rpt <- dot_art(stt, sgm)
#'
#' print(rpt$Y)  # 22.50976
#' print(rpt$S)  #  4.484002
#' print(rpt$P)  #  0.0006704994
dot_art <- function(Z, C, k, ...)
{
    L <- length(Z)                      # number of statistics
    if(missing(k))                      # k = L / 2 (default)
        k <- round(L * .5)

    within(dot(Z, C, ...),
    {
        ## de-correlated two-tail p-values, sorted
        P <- sort(1 - pchisq(X^2, df=1))
        
        ## summarized statistics: shape parameter
        S <- (k - 1) * (digamma(L + 1) - digamma(k))

        ## summarized statistics: summation
        Y <- 0
        Y <- Y + (k - 1) * log(P[k])
        Y <- Y - sum(log(P[1:(k - 1)]))
        Y <- Y + qgamma(1 - pbeta(P[k], k, L - k + 1), S)

        ## summarized statistics: p-value
        P <- 1 - pgamma(Y, k + S - 1)
    })
}


#' @describeIn dot_chisq
#' De-correlation followed by Rank Trucated Product (RTP) Test
#'
#' @examples
#' ## get the test statistics and pre-calculated LD matrix
#' stt <- readRDS(system.file("extdata", 'art_zsc.rds', package="dotgen"))
#' sgm <- readRDS(system.file("extdata", 'art_ldm.rds', package="dotgen"))
#'
#' ## run de-correlated chi-square test
#' rpt <- dot_rtp(stt, sgm)
#'
#' print(rpt$Y)  # 22.6757
#' print(rpt$P)  #  0.0007275518
dot_rtp <- function(Z, C, k, ...)
{
    L <- length(Z)                      # number of statistics
    if(missing(k))                      # k = L / 2 (default)
        k <- round(L * .5)

    within(dot(Z, C, ...),
    {
        ## de-correlated two-tail p-values, sorted
        P <- sort(1 - pchisq(X^2, df=1))

        Y <- sum(-log(P[1:k]))

        P <- integrate(function(x, y, m, n)
        {
            1-pgamma(log(qbeta(x, m + 1, n - m)) * m + y, m)
        },
        0, 1, Y, k, L)$va
    })
}
