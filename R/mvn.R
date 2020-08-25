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
    y <- z %*% matrix(rnorm(p * pts), p, pts) + m
    
    r <- double(n)
    for(i in seq(n))
        r[i] <- mean(colSums(a[i, ] < y & y < b[i, ]) == p)
    r
}

mpdf <- function(x, m=NULL, v=NULL)
{
    if(is.vector(x)) # points
        x <- t(x)
    n <- nrow(x)
    p <- ncol(x)

    if(is.null(m)) # mean
        m <- 0
    m <- rep(m, len=p)
    x <- sweep(x, 2, m)
    
    if(is.null(v)) # covariance
        v <- 1
    if(is.vector(v))
        v <- diag(v, p)

    u <- chol(v)
    a <- chol2inv(u)
    r <- -.5 * rowSums(x %*% a * x) - sum(log(diag(u))) - .5 * log(2*pi)
    exp(r)
}

wmvn <- function(x, b=NULL, m=NULL, v=NULL, pts=1e3)
{
    if(is.vector(x))                    # quantile
        x <- t(x)
    n <- nrow(x)
    p <- ncol(x)

    a <- 0 * x
    b <- x
    d <- (b - a) / pts

    if(is.null(m))                      # mean
        m <- 0
    m <- rep(m, p)

    if(is.null(v))                      # covariance
        v <- 1
    if(is.vector(v))
        v <- diag(v, p)

    psum <- 0 * d[, 1]
    for(i in seq(0, pts))
    {
        x <- a + i * d
        f <- mpdf(x, m, v)
        psum <- psum + f
    }
    psum / pts
}

mcdf <- function(x, mu=NULL, sigma=diag(ncol(x)), a=NULL)
{
    if(is.vector(x))
        x <- t(x)
    ## Monte-Carlo confidence factor for the standard error: 99 %
    gamma = 2.5;
    ## Tolerance
    err_eps = 1e-3;

    ## Dimension
    l <- nrow(x)
    q <- nrow(sigma)
    cases <- nrow(x)

    ## Default value for mu
    if (is.null(mu))
        mu <- 0
    mu <- rep(mu, len=q)
    
    ## move x to center
    ## x <- scale(x, mu, FALSE)
    x <- sweep(x, 2, mu)

    ## lower bound
    if(is.null(a))
        a <- matrix(-Inf, cases, q)

    ## Cholesky lower Tri
    u <- t(chol(sigma))
    
    ## Number of integral transformations
    n <- 1

    p <- rep(0, l)
    varsum <- rep(0, l)

    err <- rep(err_eps, l)
    
    ## Apply crude Monte-Carlo estimation
    while(any(err >= err_eps))
    {
        ## Sample from q-1 dimensional unit hypercube
        w <- matrix(runif(cases * (q - 1)), cases, q - 1)
        
        ## Transformation of the multivariate normal integral
        d <- pnorm(a[, 1] / u[1, 1])
        e <- pnorm(x[, 1] / u[1, 1])
        f <- e - d
        y <- matrix(0, cases, q - 1)

        for(i in seq(q - 1))
        {
            y[, i] <- qnorm(d + w[, i] * (e - d))
            d <- pnorm((a[, i+1] - sum(u[i+1, 1:i] * y[, 1:i])) / u[i+1, i+1])
            e <- pnorm((x[, i+1] - sum(u[i+1, 1:i] * y[, 1:i])) / u[i+1, i+1])
            f <- (e - d) * f;
        }
        n <- n + 1

        ## Estimate standard error
        varsum <- varsum + (n - 1) * ((f - p)^2) / n
        err <- gamma * sqrt (varsum / (n * (n - 1)))
        p <- p + (f - p) / n
    }
    list(p=p, err=err)
}
