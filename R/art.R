# RTP via Equ. 2
RTP <- function(Z, k, L)
{
    ## Z <- sum(-log(P[1:k]))
    integrate(function(x,y,m,n) 1-pgamma(log(qbeta(x,m+1,n-m))*m+y,m),0,1,Z,k,L)$va
}

ART.A <- function(P, k, L)
{
    ## "ART-A"
    wgt <- rep(1,k)
    z <- P
    z[1] <- ( 1 - P[1] )^L
    for(j in 2:L) z[j] <- ((1-P[j]) / (1-P[j-1]))^((L-(j-1)))
    p <- (1-z)[1:k]
    k = length(p)
    sumZ <- rep(0, k)
    y <- qnorm(p)
    z <- y
    gz <- z[1] * wgt[1]
    sumZ[1] <- gz
    for(i in 2:k) {
        gz <- p[ 1 : i ]
        for(j in 1:i) gz[j] <- z[j] * wgt[j]
        sumZ[i] <- sum(gz)
    }
    Lo = diag(k); Lo[lower.tri(Lo)] <- 1
    pSg <- Lo %*% diag(wgt[1:k]^2) %*% t(Lo)
    pCr <- cov2cor(pSg)
    sZ <- sumZ
    for(i in 1:k) {
        sZ[i] <- sumZ[i] / sqrt(diag(pSg))[i]
    }
    ppZ <- pmvnorm(lower = rep(-Inf,k), upper = rep(max(sZ), k), sigma = pCr)[1]
    c(ppZ, which.min(sZ))
}

ART <- function(lW, Pk, k, L)
{
   d = (k-1)*(digamma(L+1) - digamma(k))
   ak = (k-1)*log(Pk) - lW + qgamma(1-pbeta(Pk, k, L-k+1), shape=d)
   1 - pgamma(ak, shape=k+d-1)
}

test <- function()
{
    ## read in summary statistics and LD matrix
    stt <- readRDS("inst/extdata/art_res.rds")
    sgm <- readRDS("inst/extdata/art_ldm.rds")
    ## the number of statistics
    Z <- stt$Z
    L <- length(Z)
    ## transform Z-scores to chi-squares
    y <- Z^2
    ## calculate TQ-statistic
    yy <- y %*% rep(1, L)
    ## eigen decomposition of the LD matrix
    ee <- eigen(sgm); eivec <- ee$vectors; eigva <- ee$values
    pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
    ## calculate decorreated statistics
    x <- (Z %*% pc)^2
    k <- round(0.5*L)
                                        #k <- L
    px <- 1-pchisq(x, df=1)
    P <- sort(px)
    P.rtp <- RTP(sum(-log(P[1:k])), k, L)
    P.art <- ART(sum(log(P[1:(k-1)])), P[k], k, L)
    cat(P.art, P.rtp, "\n")
    ## P.arta <- ART.A(P, k, L)
    ## cat(P.art, P.rtp, P.arta[1], P.arta[2], "\n")
}

try1 <- function(N=1e3)
{
    zs <- readRDS("inst/extdata/art_zsc.rds")
    ld <- readRDS("inst/extdata/art_ldm.rds")

    z2 <- zs^2
    l2 <- ld^2

    x2 <- dot(z2, l2)$X
    p2 <- 1 - pchisq(x2, df=1)

    L <- length(zs)
    k <- round(L / 2)
}
