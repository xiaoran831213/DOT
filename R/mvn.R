pmvn <- function(x, m=NULL, v=NULL, pts=1e5, tol=NULL)
{
    if(is.vector(x))                    # quantile
        x <- t(x)
    n <- nrow(x)
    p <- ncol(x)

    if(is.null(m))                      # mean
        m <- 0
    m <- rep(m, p)

    if(is.null(v))                      # covariance
        v <- 1
    if(is.vector(v))
        v <- diag(v, p)

    if(is.null(tol))                    # sqrt cov
        tol <- sqrt(.Machine$double.eps)
    z <- with(svd(v),
    {
        i <- d > d[1] * tol
        d <- d[i]
        u <- u[i, i]
        v <- v[i, i]
        u %*% (sqrt(d) * t(v))
    })

    y <- z %*% matrix(rnorm(p * pts), p, pts)

    r <- apply(x, 1, function(.)
    {
        mean(colSums(y < .) == p)
    })
    r
}


mcdf <- function(x, mu=NULL, sigma=diag(ncol(x)), a=NULL)
{
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

    ## Cholesky Upper Tri
    u <- chol(sigma)

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
        dvev <- pnorm(cbind(a[, 1] / u[1, 1], x[, 1] / u[1, 1]))
        dv <- dvev[, 1]
        ev <- dvev[, 2]
        fv <- ev - dv
        y <- matrix(0, cases, q - 1)

        for(i in seq(q - 1))
        {
            y[, i] <- qnorm(dv + w[, i] * (ev - dv))

            dvev <- pnorm(cbind(
            (a[, i + 1, drop=FALSE] - u[i + 1, 1:i, drop=FALSE] * y[, 1:i, drop=FALSE]) / u[i + 1, i + 1],
            (x[, i + 1, drop=FALSE] - u[i + 1, 1:i, drop=FALSE] * y[, 1:i, drop=FALSE]) / u[i + 1, i + 1])
            )

            dv <- dvev[, 1];
            ev <- dvev[, 2];
            fv <- (ev - dv) * fv;
        }
        n <- n + 1

        ## Estimate standard error
        varsum <- varsum + (n - 1) * ((fv - p)^2) / n
        err <- gamma * sqrt (varsum / (n * (n - 1)))
        p <- p + (fv - p) / n
    }
    list(p=p, err=err)
}

ts1 <- function()
{
    x <- matrix(c(1, 2), 1, 2)
    mu <- c(0.5, 1.5)
    sigma <- matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2)
    ret <- mcdf(x, mu, sigma)
    ret
}
