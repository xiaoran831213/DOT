# TPM: https://doi.org/10.1002/gepi.0042

# use Tpm.exact for relatively small L, for large L use simulations instead
tpm.ext <- function(p, tau)
{
    lw <- k <- 0
    ltau <- log(tau)
    lp <- log(p)

    ## get log(w), w = prod_i p_i^I(p_i < tau)
    for(i in 1 : (L <- length(p)))
    {
        if(lp[i] <= ltau)
        {
            lw <- lw + lp[i]
            k <- k + 1
        }
    }
    if(k <= 0) return (1)

    w <- exp(lw)
    r <- 0
    for(k in 1:L)
    {
        r1 <- lfactorial(L) - lfactorial(L-k) - lfactorial(k)
        r2 <-  (L-k)*log(1-tau)

        innersum <- 0
        if (lw <= k*ltau)
        {
            for(s in 0 : (k-1))
            {
                a <- k*ltau - lw
                innersum <- innersum + exp(s*log(a) - lfactorial(s))
            }
            innersum <- w * innersum;
        }
        else
        {
            innersum = tau^k
        }
        r <- r + exp(r1 + r2) * innersum
    }
    return (r)
}

# TPM via simulations for large L
tpm.sim <- function(p, tau, pts = 1e5, rs=NULL) {
    set.seed(rs)
    L <- length(p)
    w0len <- length(p[p <= tau])
    w0 <- ifelse(w0len > 0, sum(-log(p[p <= tau])), 0)
    cnt <- 0
    for(i in 1:pts) {
        x <- runif(L)
        wlen <- length(p[p <= tau])
        wi <- ifelse(wlen > 0, sum(-log(x[x <= tau])), 0)
        if(wi > w0) cnt <- cnt+1
    }
    set.seed(NULL)
    cnt / pts
}

new.ext <- function(P, tau, ...)
{
    L <- length(P)                      # number of statistics
    lw <- sum(log(P[P <= tau]))         # log(w), w = prod(p_i | p_i<=tau)
    if(lw == 0)
        return(1)

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
    sum(exp(lpk + lpw))
}

new.sim <- function(P, tau, pts=1e4, rs=NULL, ...)
{
    set.seed(rs)
    L <- length(P)
    lw <- sum(log(P[P <= tau]))

    x <- replicate(pts,
    {
        p <- runif(L)
        sum(log(p[p <= tau]))
    })
    set.seed(NULL)
    mean(x < lw)
}

## check
tpm_check <- function(n=5, L=1000, tau=0.3, rs=NULL, pts=1e4)
{
    cat(sprintf("%10s  %10s  %10s  %10s\n",
                "old.ext", "old.sim", "new.ext", "new.sim"))
    for(i in seq(n))
    {
        p <- rbeta(L, 0.9, 1)
        r1 <- tpm.ext(p, tau)                 # older
        r2 <- tpm.sim(p, tau, pts=pts, rs=rs) # older
        r3 <- new.ext(p, tau)                 # newer
        r4 <- new.sim(p, tau, pts=pts, rs=rs) # newer

        cat(sprintf("%10.8f  %10.8f  %10.8f  %10.8f\n", r1, r2, r3, r4))
    }
    set.seed(NULL)
}
