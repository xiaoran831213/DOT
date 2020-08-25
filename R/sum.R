#' Decorrelation followed by Summarized Statistical Test
#'
#' Summarize a  series of genetic  association test statistics  (i.e., Z-scores,
#' \strong{Z}) given their correlation (i.e., LD-values, \strong{C}).
#'
#' @details
#'
#' The original association statistics are normally Z-scores reported by a GWAS.
#'
#' These functions first call \code{\link{dot}(Z, C)} to decorrelate the genetic
#' association test statistics, allowing various of tools to combine independent
#' statistics or corresponding  p-values into one statistics  and p-value, which
#' essentially tests  the joint null: none  of the variants are  associated with
#' the phenotype of interests.
#'
#' The  two  rank  truncated test  \code{\link{dot_art}},  \code{\link{dot_rtp}}
#' require an additional  parameter \code{k} specifying the  number of smallest
#' (decorrelated) p-values  to combine,  which by default  retains half  of the
#' p-values.  The adaptive \code{\link{dot_arta}} on the other hand, decides the
#' appropriate \code{k} automatically.
#'
#' The  property  of  combined statistics  via  \link{dot_art},  \link{dot_rtp},
#' \link{dot_arta} and \link{dot_fisher} are detailed in the reference below.
#'
#' @param Z vector of association test statistics (i.e., Z-scores).
#' @param  C matrix of  correlation among the  test statistics, as  obtained by
#'     \code{\link{cst}}.
#' @param  k  controls  the  number of  statistics  (i.e.,  \code{k}  smallest,
#'     decorrelated p-values) to combine.
#' @param ... additional parameters
#'
#' @return typically a list, with
#' \itemize{
#' \item{\code{X}:} {decorrelated  association statistics}
#' \item{\code{W}:} {orthogonal transformation, such that \code{X == WZ}}
#' \item{\code{Y}:} {the summarized statistics}
#' \item{\code{P}:} {the p-value corresponding to \code{Y}}
#' }
#'
#' @references
#' \href{https://www.frontiersin.org/articles/10.3389/fgene.2019.01051}{
#' Vsevolozhskaya, O.   A., Hu, F.,  & Zaykin,  D.  V. (2019).   Detecting weak
#' signals  by  combining  small   P-values  in  genetic  association  studies.
#' Frontiers in genetics, 10, 1051.}
#'
#' @name dot_sst
#' @seealso \code{\link{dot}}
#' 
#' @examples
#' ## get the test statistics and pre-calculated LD matrix
#' stt <- readRDS(system.file("extdata", 'art_zsc.rds', package="dotgen"))
#' sgm <- readRDS(system.file("extdata", 'art_ldm.rds', package="dotgen"))
NULL


#' @describeIn dot_sst Decorrelation followed by Chi-square Test
#'
#' @examples
#'
#' ## decorrelated chi-square test
#' rpt <- dot_chisq(stt, sgm)
#' print(rpt$Y)  # 37.2854
#' print(rpt$P)  #  0.0003736988
#' @export
dot_chisq <- function(Z, C, ...)
{
    ret <- dot(Z, C, ...)
    Y <- sum(ret$X^2)
    P <- 1 - stats::pchisq(Y, df=length(Z))
    c(list(P=P, Y=Y), ret)
}

#' @describeIn dot_sst
#'
#' Decorrelation followed by Augmented Rank Trucated (ART) Test 
#'
#' @examples
#'
#' ## decorrelated augmented rank trucated (ART) test
#' rpt <- dot_art(stt, sgm, k=6)
#' print(rpt$Y)  # 22.50976
#' print(rpt$P)  #  0.0006704994
#' @export
dot_art <- function(Z, C, k, ...)
{
    L <- length(Z)                      # number of statistics
    if(missing(k))                      # k = L / 2 (default)
        k <- round(L * .5)

    ret <- dot(Z, C, ...)

    ## decorrelated two-tail p-values, sorted
    P <- sort(1 - stats::pchisq(ret$X^2, df=1))
        
    ## summarized statistics: shape parameter
    S <- (k - 1) * (digamma(L + 1) - digamma(k))

    ## summarized statistics: summation
    Y <- 0
    Y <- Y + (k - 1) * log(P[k])
    Y <- Y - sum(log(P[1:(k - 1)]))
    Y <- Y + stats::qgamma(1 - stats::pbeta(P[k], k, L - k + 1), S)

    ## summarized statistics: p-value
    P <- 1 - stats::pgamma(Y, k + S - 1)

    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelation followed by Augmented Rank Truncated Adaptive (ARTA) Test
#'
#' \code{k} limit number of decorrelated p-values to combine; the actual number
#' is decided adaptively.
#'
#' @return
#' For Augmented Rank Truncated Adaptive (ARTA) Test,
#' \itemize{
#' \item{k:} {gives the number of smallest, decorrelated p-values combined}
#' }
#'
#' @param pts number of samples for MVN integration
#' 
#' @examples
#'
#' ## decorrelated augmented rank truncated adaptive (ARTA) test
#' rpt <- dot_arta(stt, sgm, k=6, pts=1e7)
#' print(rpt$Y)  # -1.738662
#' print(rpt$k)  #  5
#' print(rpt$P)  #  0.003174 (varies)
#' @export
dot_arta <- function(Z, C, k, pts=1e6, ...)
{
    L <- length(Z) # number of statistics
    if(missing(k)) # k = L (default)
        k <- L

    wgt <- rep(1, k)
    ret <- dot(Z, C)

    P <- sort(1 - stats::pchisq(ret$X^2, df=1))

    . <- ((1 - P) / c(1, 1 - P[-L]))^(L:1)
    p <- (1 - .)[1:k]

    sumQ <- cumsum(stats::qnorm(p) * wgt)

    pSg <- matrix(cumsum(wgt[1:k]^2), k, k)
    pSg[lower.tri(pSg)] <- t(pSg)[lower.tri(pSg)]
    pCr <- stats::cov2cor(pSg)

    sQ <- sumQ / sqrt(diag(pSg))

    Y <- max(sQ)
    ## P <- pmvnorm(lower = rep(-Inf,k), upper = rep(max(sQ), k), sigma = pCr)[1]
    P <- pmvn(rep(max(Y), k), v=pCr, pts=pts)
    c(list(P=P, k=which.min(sQ), Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelation followed  by Rank Truncated Product (RTP) Test
#'
#' @examples
#'
#' ## decorrelated Rank Truncated Product (RTP)
#' rpt <- dot_rtp(stt, sgm, k=6)
#' print(rpt$Y)  # 22.6757
#' print(rpt$P)  #  0.0007275518
#' @export
dot_rtp <- function(Z, C, k, ...)
{
    L <- length(Z)                      # number of statistics
    if(missing(k))                      # k = L / 2 (default)
        k <- round(L * .5)

    ret <- dot(Z, C, ...)

    ## decorrelated two-tail p-values, sorted
    P <- sort(1 - stats::pchisq(ret$X^2, df=1))
    Y <- sum(-log(P[1:k]))
    P <- stats::integrate(function(x, y, m, n)
    {
        1 - stats::pgamma(log(stats::qbeta(x, m + 1, n - m)) * m + y, m)
    },
    0, 1, Y, k, L)$va
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelation followed by Fisher's Combined p-value Test
#'
#' @examples
#'
#' ## decorrelated Fisher's combined p-value chi-square test
#' rpt <- dot_fisher(stt, sgm)
#' print(rpt$Y)  # 58.44147
#' print(rpt$P)  #  0.0002706851
#' @export
dot_fisher <- function(Z, C, ...)
{
    L <- length(Z)                      # number of statistics
    ret <- dot(Z, C, ...)
    
    ## decorrelated two-tail p-values, sorted
    P <- sort(1 - stats::pchisq(ret$X^2, df=1))
    Y <- -2 * sum(log(P))
    P <- 1 - stats::pchisq(Y, 2 * L)
    c(list(P=P, Y=Y), ret)
}
