#' a package to decorrelate genetic association test statistics
#'
#' The package  decorrelate genetic  association test statistics  via orthogonal
#' transformation.
#'
#' The core  function \code{\link{dot}()} decorrelates genetic  association test
#' statistics which typically are Z-scores  of genetic variants (i.e., SNP) came
#' from a genome wide association  study (GWAS). The decorrelated statistics can
#' be analyzed as independent normal variables, for example, to test the overall
#' effect  of all  SNP(s) in  a gene,  one could  use the  chi-square statistics
#' formed by the sum of squares of decorrelated Z-scores.
#'
#' \code{\link{dot}()} requires the correlation  of association test statistics,
#' but does not require the study genotype and phenotype that gave rise to these
#' test statistics, therefore bypassing the logistic issues of data sharing.
#'
#' The  package  also  provides  function \code{\link{css}()}  to  estimate  the
#' correlation among association test statistics  from the genotype, which takes
#' care of missing values and covariates.
#'
#' For details of  the underlying decorrelation algorithm and its  impact on the
#' power of  association test, see the  reference below.  For an  example use of
#' the package, please type \code{?dot}.
#'
#' @references
#'     \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819&rev=2}{DOT:
#'     Gene-set analysis by combining decorrelated association statistics}
#' 
#' @author Olga A. Vsevolozhskaya, \email{ovs222@@g.uky.edu}
#' @author Dmitri V. Zaykin, \email{zaykin@@gmail.com}
#' @author Xiaoran Tong, \email{tongxia1@msu.edu}
#' 
#' @docType package
#' @name dotgen
NULL

#' Decorrelation by Orthogonal Transformation
#'
#' Decorrelate  association  test  statistics   (i.e.,  Z-scores)  by  symmetric
#' orthogonal transformation
#'
#' \code{dot()} takes  association test statistics (\eqn{\strong{z}})  and their
#' correlation   (\eqn{\strong{C}})    to   produce    decorrelated   statistics
#' (\eqn{\strong{x}}), that  is, \eqn{\strong{x} = \strong{Wz}};  the orthogonal
#' transformation \eqn{\strong{W}}  is the  inverted square root  of correlation
#' \eqn{\strong{C}}, that is, \eqn{(\strong{W W'})^{-1} = \strong{C}}.
#'
#' To estimate \eqn{\strong{C}} which is  the correlation among association test
#' statistics, use function \code{\link{css}()}.
#'
#' If the p-values and effects (i.e, beta coefficients) are available instead of
#' the test statistics, recover the Z-scores with \code{\link{zsc}()} first.
#'
#' @param Z vector of association test statistics (i.e., Z-scores)
#' @param  C matrix  of correlation  among the  association test  statistics, as
#'     obtained by \code{\link{css}()}
#' @param ... additional parameters
#'
#' @return a  list containing decorrelated statistics  (\code{X}) and orthogonal
#'     transformation matrix (\code{W})
#' 
#' @seealso \code{\link{css}}, \code{\link{zsc}}
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## estimate the correlation among association test statistics
#' sgm <- css(gno, cvr)
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
#'
#' @author Olga A. Vsevolozhskaya, \email{ovs222@@g.uky.edu}
#' @author Dmitri V. Zaykin, \email{zaykin@@gmail.com}
#' @author Xiaoran Tong, \email{tongxia1@msu.edu}
#' @references
#'     \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819&rev=2}{DOT:
#'     Gene-set analysis by combining decorrelated association statistics}
#' @export
dot <- function(Z, C=NULL, ...)
{
    if(is.null(C))
        C <- diag(length(Z))
    W <- nsp(C)                         # orthogonal transformation
    X <- W %*% Z                        # decorrelated statistics
    list(W=W, X=X)
}

#' Z-scores from p-values and estimated effects
#'
#' \code{zsc} recovers  Z-scores from  the p-values and  estimated effects  of a
#' genetic association analysis.
#'
#' In reality,  for each genetic  variant, it is  the p-value (\eqn{p})  and the
#' sign  of estimated  effect (\eqn{\beta})  being used  to recover  the Z-score
#' (\eqn{z}), that is, \eqn{z = sign(\beta) \Phi^{-1}(p/2)}.
#'
#' @seealso \code{\link{dot}}
#'
#' @param P vector of p-values
#' @param BETA vector of effects, or their signs
#' @return vector of Z-scores
#'
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

