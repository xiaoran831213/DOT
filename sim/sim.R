## a reference, test of the sum quadratic terms
tsq <- function(Z, C, ...)
{
    S <- sum(Z^2) # sum square

    ## chi-square mixture
    L <- eigen(C, TRUE, TRUE)$value
    P <- imhof(S, L, delta=rep(0, length(Z)))$Qq
    list(mtd='tsq', Z=Z, C=C, lmd=L, ssq=S, pvl=P)
}

## False Discovery Rate
fdr <- function(P, ...) min(p.adjust(P, 'fdr'))

## bonferroni correction
bon <- function(P, ...) min(p.adjust(P, 'bon'))


sim <- function(N=1e3, P=N/20, hsq=.1, MAF=0.01, frq=1, ...)
{
    arg <- get.arg()
    set.seed(arg$seed)
    dat <- get.sim(N, P, frq=frq, hsq=hsq, MAF=MAF, ...)
    set.seed(NULL)
    gwa <- dat$gwa
    mra <- dat$mra$P
    gmx <- dat$gmx
    ngv <- ncol(gmx)
    gcv <- cor(gmx)
    
    dot <- .d(mtd='dot', pvl=dot(gwa$Z, gcv)$pvl)
    tsq <- .d(mtd='tsq', pvl=tsq(gwa$Z, gcv)$pvl)
    mra <- .d(mtd='mra', pvl=mra)
    fdr <- .d(mtd='fdr', pvl=fdr(gwa$P))
    bon <- .d(mtd='bon', pvl=bon(gwa$P))

    rpt <- rbind(dot, tsq, mra, fdr, bon)
    cbind(arg, ngv=ngv, rpt)
}


gwa <- function(x, y)
{
    rowSums(y * x) / rowSums(x * x)
}
