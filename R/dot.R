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
#' A number statistics that combines de-correlated P-values are made available,
#' see \code{\link{dot_sst}} for detailed.
#'
#' Tow or more genetic variants with correlation higher than \code{1 - tol.cor}
#' are  considered colinear,  and only  one them  is retained  to bring  the LD
#' matrix  closer to  full-rank.   The default  tolerence  of high  correlation
#' (\code{tol.cor})is \code{sqrt(.Machine$double.eps)}
#'
#' Eigenvalues  smaller than  \code{tol.egv} times  the largest  eigenvalue are
#' treated as  non-positve and truncated  to make the  orthogonal tranformation
#' \code{W} possible.  The default tolerence of small eigenvalue \code{tol.egv}
#' is \code{sqrt(.Machine$double.eps)}.
#' 
#' For details about DOT, see
#' the reference below.
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
#' @param tol.cor tolerence for correlation
#' @param tol.egv tolerence for eigenvalues
#' @param ... additional parameters
#'
#' @return
#' a  list,
#' \itemize{
#' \item{\code{Z}:}{association test statistics, original.}
#' \item{\code{X}:}{association test statistics, de-correlated.}
#' \item{\code{W}:}{orthogonal transformation such that \code{X = WZ}.}
#' \item{\code{M}:}{effective number of variants after breaking colinarity.}
#' \item{\code{L}:}{effective number of eigenvalues after zeroing the negatives.}
#' }
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
dot <- function(Z, C, tol.cor=NULL, tol.egv=NULL, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    if(is.null(tol.egv))
        tol.egv <- sqrt(.Machine$double.eps)

    ## trim colinear variants
    m <- dvt(C, tol.cor)
    M <- sum(M)                      # effective number of variants
    C <- C[m, m]
    Z <- Z[m]

    ## get orthogonal transformation
    d <- nsp(C, eps=tol.egv, ...)
    W <- d$W                         # orthogonal transformation
    L <- d$L                         # effective number of eigenvalues

    X <- W %*% Z                     # decorrelated statistics

    list(Z=Z, W=W, X=X, L=L, M=M)
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


#' Correlation among association test statistics
#'
#' Calculates the correlation among genetic association test statistics.
#'
#' @details
#' When no covariates are present in per-variant association analyses, that is,
#' \code{x==NULL},  correlation  among  test  statistics is  the  same  as  the
#' correlation among variants, \code{cor(g)}.
#'
#' With  covariates, correlation  among  test  statistics is  not  the same  as
#' \code{cor(g)}. In this case, \code{\link{cst}} takes the generalized inverse
#' of the entire correlation matrix, \code{corr(cbind(g, x))}, and then inverts
#' back only the submtarix containing genotype variables, \code{g}.
#'
#' Missed genotype calls are fill  with 'soft' (i.e., continuous) allele dosage
#' values  in the  interval from  0 to  2, imputed  using information  from all
#' non-missing entries in the genotype matrix.
#'
#' By default (\code{imp=0}), missed calls for each variant are filled with the
#' average of non-missing entires; setting \code{imp=1} uses advance techniques
#' based on expected allele count of one variant conditioned on other variants.
#'
#' The advance imputation (\code{imp=1})  improves statistical power in certian
#' cases, but  we highly suggest  the user to  rerun the association  test with
#' imputed genotype to avoid type I error.
#'
#' @param  g matrix of  genotype, one row per  sample, one column  per variant,
#'     missings allowed.
#' @param x matrix of covariates, one row per sample, no missing allowed.
#' @param imp  imputation methods,  0=naive averge,  1=conditional expectation
#'     (def=0).
#'
#' @return a list containing
#' \itemize{
#' \item{G}{imputed genotype matrix, no missing entreis}
#' \item{C}{the correlation matrix \code{cor{G}}}
#' }
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## correlation among association statistics, covariates involved
#' res <- cst(gno, cvr)
#' print(res$C[1:4, 1:4])
#'
#' ## with 2% missed calls
#' g02 <- readRDS(system.file("extdata", 'rs208294_g02.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#' res <- cst(g02, cvr, imp=0)
#' print(res$C[1:4, 1:4])
#' 
#' @export
cst <- function(g, x=NULL, imp=0)
{
    ## impute missed calls
    if(imp == 1)
        imp <- imp.cxp
    else
        imp <- imp.avg
    g <- imp(g)
    
    if(is.null(x))
        r <- stats::cor(g)              # no covariate, use full cor
    else
    {
        r <- stats::cor(cbind(x, g))    # full cor
        i <- seq(ncol(x))               # index of covariates
        r <- stats::cov2cor(scp(r, i))  # cond cor
    }
    list(G=g, C=r)
}
