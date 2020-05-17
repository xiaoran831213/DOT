## Imputation by random sample
imp.rnd <- function(g)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- sample(x[-i], length(i))
        x
    })
}

## remove colinear variable
rmv.col <- function(x, t=1.0)
{
    r <- abs(cor(x))
    r[upper.tri(r, TRUE)] <- 0
    b <- which(r >= t, arr.ind=TRUE)
    x[, apply(r, 2, function(.) all(. < t))]
}

## get genotype matrix
get.gmx <- function(N=200, P=N/10, MAF=0.01, gno="i21", bat="003", nlr=1.0, ...)
{
    if(is.null(bat))
        bat <- sub(".bed$", "", dir(gno, "bed$", ful=TRUE))
    gmx <- readBED(file.path(gno, bat), quiet=FALSE)
    
    ## MAF threshold
    alf <- colMeans(gmx, TRUE) / 2.0
    gmx[, alf > 0.5] <- 2 - gmx[, alf > 0.5]
    maf <- pmin(alf, 1 - alf)
    gmx <- gmx[, maf > MAF] 

    ## N x P
    i <- sample.int(nrow(gmx), N)               # N
    j <- seq(sample(ncol(gmx) - P, 1) + 1, l=P) # P
    gmx <- gmx[i, j]                            # N x P
    gmx <- imp.rnd(gmx)                         # no NA

    ## maintain non-linearity
    gmx <- rmv.col(gmx, nlr)
    gmx
}

## get phenotype
get.phe <- function(gmx, frq=.2, hsq=.5)
{
    N <- nrow(gmx)
    P <- ncol(gmx)

    if(hsq > 0 && frq > 0)
    {
        Q <- max(1, P * frq)
        cmx <- gmx[, sample(P, Q), drop=FALSE]
        ## sgm <- tcrossprod(scale(cmx)) / Q
        ## gvt <- mvn(1, 0, sgm)   # genetic effect
        b <- rnorm(Q)
        gvt <- cmx %*% b
        
        gvr <- var(gvt)            # genetic variance
        evr <- gvr / hsq - gvr     # noise variance
    }
    else
    {
        gvt <- 0.0 # no genotype effect
        evr <- 1.0 # pure noise
        print("noise")
    }
    evt <- rnorm(N, 0, evr) # noise vector
    as.vector(scale(gvt + evt))
}

tp1.err <- function(gmx, tms=1e3, int=TRUE)
{
    n <- nrow(gmx)
    r <- t(replicate(tms,
    {
        y <- as.vector(scale(rnorm(n)))
        get.gwa(y, gmx, int)[, 3]
    }))
    r
}

## get GWAS
get.gwa <- function(y, gmx, int=TRUE)
{
    if(int)
        r <- t(apply(gmx, 2, function(x) summary(lm(y ~ x + 1))$coef[2, ]))
    else
        r <- t(apply(gmx, 2, function(x) summary(lm(y ~ x - 1))$coef[1, ]))
    colnames(r) <- c("BETA", "SE", "Z", "P")
    rownames(r) <- rownames(gmx)
    .d(r)
}

## get.gwa <- function(y, x)
## {
##     n <- nrow(x)
##     xx <- colSums(x^2)
##     xy <- colSums(x*y)
##     b <- xy / xx
##     e <- y - sweep(x, 2, b, `*`)
##     ee <- colSums(e^2)
##     s <- sqrt(ee / ((n - 1) * xx))
##     z <- b / s
##     p <- pt(z, n-1, lower.tail = TRUE) * 2
##     .d(BETA=b, SE=s, Z=z, P=p)
## }
## get multi-regression analysis
get.mra <- function(y, gmx, int=TRUE)
{
    if(int)
        ret <- anova(lm(y ~ gmx + 1))[1, c(1, 3, 4, 5)]
    else
        ret <- anova(lm(y ~ gmx - 1))[1, c(1, 3, 4, 5)]
    names(ret) <- c('df', 'ssq', 'F', 'P')
    ret
}

## get simulation data
get.sim <- function(N=500, P=N/20, hsq=.1, frq=1, MAF=0.01, ...)
{
    gmx <- get.gmx(N, P, MAF, ...)
    y <- get.phe(gmx, frq, hsq)
    gwa <- get.gwa(y, gmx)
    mra <- get.mra(y, gmx)

    list(gmx=gmx, y=y, gwa=gwa, mra=mra)
}
