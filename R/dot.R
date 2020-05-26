#' dotgen: a package to combine genetic association tests
#'
#' Summarize genetic association test statistics via Decorrelation by Orthogonal
#' Transformation (DOT).
#'
#' The main function \code{\link{dot}}  decorrelate \code{M} genetic association
#' test statistics  (i.e., z-scores of  \code{M} SNPs in  a gene) and  gives one
#' chi-square statistics of  \code{M} degree of freedom,  which usually provides
#' superior statistical power  than the largest of  \code{M} original statistics
#' after Bonferroni or FDR correction.
#'
#' DOT requires  only the \code{M}  test statistics and the  correlation between
#' involvded genetic variants, not  the genotype and phenotype, therefore it
#' bypasses the logistic and privacy concern of sharing the data.
#'
#' The  variants correlation  may come  from the  undisclosed study  genotype or
#' publicly available genotype such as the 1000 genome project.  The correlation
#' should  also  consider the  prior  information  of covariants.   The  package
#' provides  helpers to  calculate  the correlation  of  genotype variants  with
#' missing values, conditioning on given covariants.
#'
#' @author Olga A. Vsevolozhskaya, \email{ovs222@@g.uky.edu}
#' @author Dmitri V. Zaykin, \email{zaykin@@gmail.com}
#' @references
#'     \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819&rev=2}{DOT:
#'     Gene-set analysis by combining decorrelated association statistics}
#' 
#' @docType package
#' @name dotgen
#'
#' @importFrom plinkFile readBED
#' @importFrom plinkFile readBSM
#' @importFrom plinkFile saveBSM
NULL

#' Decorrelation by Orthogonal Transformation
#'
#' \code{dot} summarizes  genetic association  test statistics  (i.e., z-scores)
#' into a chi-square statistics.
#'
#' Given \code{M} genetic association test statistics (\eqn{\strong{z}}) and the
#' correlation (\eqn{\Sigma})  of the  \code{M} underlying  variants, \code{dot}
#' first de-correlate the test statistics via Orthogonal Transformation
#'
#' \deqn{\strong{x} = \Sigma^{-\frac{1}{2}} \strong{z},}
#'
#' the summarized statistics  is the sum of square  \eqn{\sigma^2 = \sum_{m=1}^M
#' x_m^2}, which follows a chi-square of \code{M} degree of freedom.
#' 
#' In case  the p-values and  effects (i.e, beta)  are available instead  of the
#' test statistics, recover z-scores with \code{\link{zsc}} first.
#'
#' Genotype  owners can  read  and calcuate  and  calculate variant  correlation
#' matrix with \code{\link{cog}}.
#'
#'
#' @seealso \code{\link{zsc}}
#' @seealso \code{\link{cog}}
#'
#' @param Z vector of variant statistics (z-scores)
#' @param C matrix of variant correlation
#' @param ... additional parameters
#' @return the chi-square statistics, p-value, and other information.
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## conditional variate correlation
#' sgm <- cvc(gno, cvr)
#'
#' ## result of per-variant test (p-values and effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## DOT  summarized chi-square  statistics based  on z-scores  and conditional
#' ## variant correlation
#' rpt <- dot(stt, sgm)
#'
#' print(rpt$ssq)  # chi-square = 35.76306
#' print(rpt$pvl)  #    p-value =  0.001132132
#'
#' @export
dot <- function(Z, C=NULL, ...)
{
    if(is.null(C))
        C <- diag(length(Z))

    W <- nsp(C)   # DOT weights
    X <- W %*% Z  # decorrelate
    S <- sum(X^2) # chi-square
    
    ## P-values
    P <- 1 - pchisq(S, df=length(Z), ncp=0)
    list(mtd='dot', Z=Z, C=C, W=W, X=X, ssq=S, pvl=P)
}

#' z-scores from p-values and effects
#'
#' \code{zsc}  recovers  z-scores  from   the  p-values  (\eqn{p})  and  effects
#' (\eqn{beta})  of   a  genetic   associatio  analysis,   that  is,   \eqn{z  =
#' \Phi^{-1}(\frac{p}{2})  sign(\beta)}.  The  z-scores can  then be  summarized
#' into a chi-square statistics by function \code{\link{dot}}.
#'
#' @seealso \code{\link{dot}}
#'
#' @param P vector of p-values
#' @param BETA vector of effects, or their signs
#' @return vector of Z-scores
#'
#' @examples
#' ## result of per-variant test (p-values and effects)
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
    qnorm(P / 2) * -(sign(BETA))
}
