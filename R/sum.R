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
#' Reference (a) details how to combine decorrelated test statistics or p-values
#' via  \code{\link{dot_art}},  \code{\link{dot_rtp}} and  \link{dot_arta};  for
#' \code{\link{dot_rtm}} and \code{\link{dot_fisher}}, see reference (b).
#'
#' @param Z vector of association test statistics (i.e., Z-scores).
#' @param  C matrix  of correlation  among the test  statistics, as  obtained by
#'     \code{\link{cst}}.
#' @param k consider \code{k} smallest (decorrelated) p-values.
#'
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
#' (a) \href{https://www.frontiersin.org/articles/10.3389/fgene.2019.01051}{
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
#' result <- dot_chisq(stt, sgm)
#' print(result$Y)  # 37.2854
#' print(result$P)  #  0.0003736988
#' @export
dot_chisq <- function(Z, C, ...)
{
    ret <- dot(Z, C, ...)
    Y <- sum(ret$X^2)
    P <- 1 - pchisq(Y, df=length(Z))
    c(list(P=P, Y=Y), ret)
}

#' @describeIn dot_sst
#'
#' Decorrelated Augmented Rank Trucated (ART) Test 
#'
#' @examples
#'
#' ## decorrelated augmented rank trucated (ART) test
#' result <- dot_art(stt, sgm, k=6)
#' print(result$Y)  # 22.50976
#' print(result$P)  #  0.0006704994
#' @export
dot_art <- function(Z, C, k=NULL, ...)
{
    L <- length(Z)                      # number of statistics
    if(is.null(k))                      # k = L / 2 (default)
        k <- round(L * .5)
    ret <- dot(Z, C, ...)

    ## decorrelated two-tail p-values, sorted
    P <- sort(1 - pchisq(ret$X^2, df=1))
    if(any(P == 0))
        return(c(list(P=0, Y=1, k=k, ret)))

    ## summarized statistics: shape parameter
    S <- (k - 1) * (digamma(L + 1) - digamma(k))

    ## summarized statistics: summation
    Y <- 0
    Y <- Y + (k - 1) * log(P[k])
    Y <- Y - sum(log(P[1:(k - 1)]))
    Y <- Y + qgamma(1 - pbeta(P[k], k, L - k + 1), S)

    ## summarized statistics: p-value
    P <- 1 - pgamma(Y, k + S - 1)
    
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelated Augmented Rank Truncated Adaptive (ARTA) Test
#'
#' @return
#' for Augmented Rank Truncated Adaptive (ARTA) Test,
#' \itemize{
#' \item{k:} {the number of decorrelated p-values adaptively picked}}
#'
#' @param w weight assigned to cumulated statistics, default to 1.
#' 
#' @examples
#'
#' ## decorrelated augmented rank truncated adaptive (ARTA) test
#' result <- dot_arta(stt, sgm, k=6)
#' print(result$Y)  # -1.738662
#' print(result$k)  #  5 smallest p-value retained
#' print(result$P)  #  0.003165 (varies)
#' @export
dot_arta <- function(Z, C, k=NULL, w=NULL, ...)
{
    L <- length(Z)                      # number of statistics
    if(is.null(k))                      # k = L (default)
        k <- L
    if(is.null(w))
        w <- rep(1, k)
    k <- min(k, L)                      # k <= L
    w <- w[k]                           # first k weights

    ret <- dot(Z, C)
    P <- sort(1 - pchisq(ret$X^2, df=1))

    ## z_i = [(1 - p_i) / (1 - p_{i-1})]^(L - i + 1), i = 1 .. L, p_0 = 1
    z <- ((1 - P) / c(1, 1 - P[-L]))^(L:1)

    ## cumulate the top-ranking k statistics, apply weights
    p <- (1 - z)[1:k]

    sumQ <- cumsum(qnorm(p) * w)

    ## L(1) %*% diag(w^2) %*% U(1)
    pSg <- matrix(cumsum(w^2), k, k)
    pSg[lower.tri(pSg)] <- t(pSg)[lower.tri(pSg)]
    pCr <- cov2cor(pSg)

    sQ <- sumQ / sqrt(diag(pSg))

    Y <- max(sQ)
    P <- mvtnorm::pmvnorm(lower=rep(-Inf,k), upper=rep(max(sQ), k), sigma=pCr)[1]
    c(list(P=P, k=which.min(sQ), Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelated Rank Truncated Product (RTP) Test
#'
#' @examples
#'
#' ## decorrelated Rank Truncated Product (RTP)
#' result <- dot_rtp(stt, sgm, k=6)
#' print(result$Y)  # 22.6757
#' print(result$P)  #  0.0007275518
#' @export
dot_rtp <- function(Z, C, k=NULL, ...)
{
    L <- length(Z)                      # number of statistics
    if(is.null(k))                      # k = L / 2 (default)
        k <- round(L * .5)
    k <- min(L, k)

    ## decorrelated two-tail p-values, sorted
    ret <- dot(Z, C, ...)
    P <- sort(1 - pchisq(ret$X^2, df=1))

    Y <- sum(-log(P[1:k]))
    P <- integrate(function(x, y, m, n)
    {
        1 - pgamma(log(qbeta(x, m + 1, n - m)) * m + y, m)
    },
    0, 1, Y, k, L)$va
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelated Fisher's Combined p-value Test
#'
#' @examples
#'
#' ## decorrelated Fisher's combined p-value chi-square test
#' result <- dot_fisher(stt, sgm)
#' print(result$Y)  # 58.44147
#' print(result$P)  #  0.0002706851
#' @export
dot_fisher <- function(Z, C, ...)
{
    L <- length(Z)                      # number of statistics
    ret <- dot(Z, C, ...)               # decorrelated two-tail p-values
    P <- 1 - pchisq(ret$X^2, 1)
    
    Y <- -2 * sum(log(P))               # Fisher statistics
    P <- 1 - pchisq(Y, 2 * L)           # summarized p-value
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' Decorrelated Truncated Product Method (TPM).
#'
#' @examples
#'
#' ## decorrelated Truncated Product Method (TPM)
#' result <- dot_tpm(stt, sgm, tau=0.05)
#' print(result$Y)  #  1.510581e-08
#' print(result$k)  #  6 p-values <= tau
#' print(result$P)  #  0.0007954961
#'
#' @param tau combine (decorrelated) p-values no large than tau (def=0.05).
#' @return
#' for Truncated Product Method (TPM),
#' \itemize{
#' \item{k:} {the number of decorrelated p-values \eqn{\le} \code{tau}}}
#'
#' @references
#' (b) \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.0042}{Zaykin,
#' D. V., Zhivotovsky, L.  A., Westfall, P. H., & Weir,  B. S. (2002). Truncated
#' product  method for  combining P‐values.  Genetic Epidemiology:  The Official
#' Publication  of  the  International   Genetic  Epidemiology  Society,  22(2),
#' 170-185.}
#'
#' @export
dot_tpm <- function(Z, C, tau=0.05, ...)
{
    L <- length(Z)                      # number of statistics
    if(is.null(tau))                    # k = L / 2 (default)
        tau <- 0.05

    ret <- dot(Z, C, ...)               # decorrelated two-tail p-values
    P <- 1 - pchisq(ret$X^2, 1)

    k <- sum(P <= tau)
    if(k == 0)
        return(c(list(P=1, Y=1, k=k, ret)))
    if(any(P == 0))
        return(c(list(P=0, Y=1, k=k, ret)))


    lw <- sum(log(P[P <= tau]))         # log(w), w = prod(p_i | p_i<=tau)
    
    ## sequence: k = 1 .. L
    K <- seq(L)

    ## log Pr(k) = log(L choose k) + log((1 - tau)^(L - k), k = 1..L
    lpk <- lchoose(L, K) + (L - K) * log(1 - tau)

    ## log Pr(W <= w | k)
    lpw <- double(L)
    klt <- K * log(tau)                 # k ln(tau), k = 1..L
    for(k in K[lw <= klt])              # case 1: ln(w) <= k ln(tau)
    {
        S <- seq(0, k - 1)
        ## sum_{s=0..k-1}(k log(tau) - log(w))^s / s!
        . <- S * log(klt[k] - lw) - lfactorial(S)
        . <- exp(.)
        lpw[k] <- lw + log(sum(.))
    }
    lpw[lw > klt] <- klt[lw > klt]      # case 2: ln(w)  > k ln(tau)
    
    ## sum Pr(k) Pr(W <= w|k), k = 1 .. L
    P <- sum(exp(lpk + lpw))
    Y <- exp(lw)
    
    c(list(P=P, Y=Y, k=k, ret))
}
