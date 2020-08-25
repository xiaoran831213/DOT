#' De-correlation followed by Summarized Statistical Test
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
#' The two rank trancated tools \code{\link{dot_art}} and \code{\link{dot_rtp}},
#' require an addnitional  parameter \code{k} specifying the  number of smallest
#' (de-correlated) p-values  to combine,  which by default  retains half  of the
#' p-values.
#'
#' The property  of combined statistics via  \link{dot_art}, \link{dot_rtp}, and
#' \link{dot_fisher} are detailed in the reference below.
#'
#' @param Z vector of association test statistics (i.e., Z-scores).
#' @param  C matrix of  correlation among the  test statistics, as  obtained by
#'     \code{\link{cst}}.
#' @param k controls how many top statistics (i.e., \code{k} smallist p-values)
#'     to be combined.
#' @param ... additional parameters
#'
#' @return
#' A list  with decorrelated  association statistics \strong{X},  the orthogonal
#' transformation \strong{W} such that  \strong{X} = \strong{WZ}, the summarized
#' statistics Y, and the corresponding p-value P.
#'
#' @references
#' \href{https://www.frontiersin.org/articles/10.3389/fgene.2019.01051}{
#' Vsevolozhskaya, O.   A., Hu, F.,  & Zaykin,  D.  V. (2019).   Detecting weak
#' signals  by  combining  small   P-values  in  genetic  association  studies.
#' Frontiers in genetics, 10, 1051.}
#'
#' @name dot_sst
#' @export
#'
#' @examples
#' ## get the test statistics and pre-calculated LD matrix
#' stt <- readRDS(system.file("extdata", 'art_zsc.rds', package="dotgen"))
#' sgm <- readRDS(system.file("extdata", 'art_ldm.rds', package="dotgen"))
NULL


#' @describeIn dot_sst De-correlation followed by Chi-square Test
#'
#' @examples
#'
#' ## de-correlated chi-square test
#' rpt <- dot_chisq(stt, sgm)
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

#' @describeIn dot_sst
#'
#' De-correlation followed by Augmented Rank Trucated (ART) Test 
#'
#' @examples
#'
#' ## de-correlated augmented rank trucated (ART) test
#' rpt <- dot_art(stt, sgm)
#' print(rpt$Y)  # 22.50976
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


#' @describeIn dot_sst
#'
#' De-correlation followed by Augmented Rank Trucated Adaptive (ARTA) Test
#'
#' The paramter  \code{k} limite the number  of top statistics to  combine. The
#' actual number of statistics included is decided adaptively.
#'
#' @return
#' the actual number (\code{k}) of statistics combined
#'
#' @param pts number of samples for MVN integration
#' 
#' @examples
#'
#' ## de-correlated augmented rank trucated (ART) test
#' rpt <- dot_arta(stt, sgm, round(length(stt)/2), pts=1e7)
#' print(rpt$Y)  # -1.738662
#' print(rpt$k)  #  5
#' print(rpt$P)  #  0.003174 (varies)
dot_arta <- function(Z, C, k, pts=1e6, ...)
{
    L <- length(Z) # number of statistics
    if(missing(k)) # k = L (default)
        k <- L

    wgt <- rep(1, k)
    ret <- dot(Z, C)

    P <- sort(1 - pchisq(ret$X^2, df=1))

    . <- ((1 - P) / c(1, 1 - P[-L]))^(L:1)
    p <- (1 - .)[1:k]

    sumQ <- cumsum(qnorm(p) * wgt)

    pSg <- matrix(cumsum(wgt[1:k]^2), k, k)
    pSg[lower.tri(pSg)] <- t(pSg)[lower.tri(pSg)]
    pCr <- cov2cor(pSg)

    sQ <- sumQ / sqrt(diag(pSg))

    Y <- max(sQ)
    ## P <- pmvnorm(lower = rep(-Inf,k), upper = rep(max(sQ), k), sigma = pCr)[1]
    P <- pmvn(rep(max(Y), k), v=pCr, pts=pts)
    c(list(P=P, k=which.min(sQ), Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' De-correlation followed  by Rank Trucated Product (RTP) Test
#'
#' @examples
#'
#' ## de-correlated Rank Trucated Product (RTP)
#' rpt <- dot_rtp(stt, sgm)
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


#' @describeIn dot_sst
#'
#' De-correlation followed by Fisher's Combined p-value Test
#'
#' @examples
#'
#' ## de-correlated Fisher's combined p-value chi-square test
#' rpt <- dot_fisher(stt, sgm)
#' print(rpt$Y)  # 58.44147
#' print(rpt$P)  #  0.0002706851
dot_fisher <- function(Z, C, ...)
{
    L <- length(Z)                      # number of statistics
    within(dot(Z, C, ...),
    {
        ## de-correlated two-tail p-values, sorted
        P <- sort(1 - pchisq(X^2, df=1))

        Y <- -2 * sum(log(P))

        P <- 1 - pchisq(Y, 2 * L)
    })
}
