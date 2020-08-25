pmvn <- function(x, b=NULL, m=NULL, v=NULL, pts=1e5, tol=NULL)
{
    if(is.vector(x)) # quantile
        x <- t(x)
    n <- nrow(x)
    p <- ncol(x)

    ## bounds
    if(is.null(b))
    {
        a <- matrix(-Inf, n, p)
        b <- x
    }
    else
        a <- x
    
    if(is.null(m)) # mean
        m <- 0
    m <- rep(m, p)

    if(is.null(v)) # covariance
        v <- 1
    if(is.vector(v))
        v <- diag(v, p)

    ## take samples
    if(is.null(tol)) # sqrt cov
        tol <- sqrt(.Machine$double.eps)
    z <- with(svd(v),
    {
        d[d < d[1] * tol] <- 0
        u %*% (sqrt(d) * t(v))
    })
    y <- z %*% matrix(stats::rnorm(p * pts), p, pts) + m
    
    r <- double(n)
    for(i in seq(n))
        r[i] <- mean(colSums(a[i, ] < y & y < b[i, ]) == p)
    r
}

ts1 <- function()
{
    x <- c(3, 3, 3)
    b <- -x
    m <- c(0, 0, 0)
    v <- matrix(c(1.0, 0.5, 0.0,
                  0.5, 1.0, 0.5,
                  0.0, 0.5, 1.0), 3, 3)
    pmvn(x, b, m=m, v=v, pts=1e5)
}
