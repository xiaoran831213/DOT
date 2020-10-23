library(MASS)
library(Matrix)
options(width=180)

main <- function()
{
    oo <- 2e6
    L <- 5 # number of statistics
    Le <- 3 # number of eigenvalues to keep
                                        #
    ## generate MVN correlation matrix
    R <- matrix(1, L, L)
    R[lower.tri(R)] <- sort(2*rbeta(L*(L-1)/2, 1, 1) - 1)
    R <- (R * lower.tri(R)) + t(R * lower.tri(R)); diag(R) <- 1
    R <- as.matrix(nearPD(cov2cor(R), corr=TRUE, posd.tol = 1e-8)$mat)
    u <- 1.5 *(2*(rbeta(L, 1, 1) - 0.5)) # add jiggle
    R <- cov2cor(R + u %*% t(u))
    LD <- R[lower.tri(R)]
    summary(LD)
    ## sample MVN statistics
    Y <- mvrnorm(n = oo, mu = rep(0,L), Sigma=R)
    ##
    ee <- eigen(R)

    ## drop L-Le positive eigenvalues
    eigva <- ee$values # [1 : Le]
    eivec <- ee$vectors # [, 1 : Le]
    sD <- diag(sqrt(1/eigva))
    H <- eivec %*% sD %*% t(eivec)
    X <- (Y %*% H)
    Cr <- cor(X)
    ## when Le < L, X remain correlated
    summary(Cr[lower.tri(Cr)])

    ## scale X to xx so that xx have unit variance
    ## H %*% R %*% H <-- cov(X) == cov(Y %*% H)
    dg <- diag(1/sqrt(diag(H %*% R %*% H)))
    xx <- X %*% dg # diag(cov(xx)) == 1
    ## sum squared and scaled values
    xx2 <- xx^2 %*% rep(1, L)
    ## compute theoretical eigenvalues of cor(xx) and set negative or small values to zero
    w <- eigen(cov2cor(H %*% R %*% H))$values # <-- eigen(cor(xx))$values
    w[w < 1e-14] <- 0
    ##

    ## empirical cor(xx^2) should approach 'xxc'
    ## when the mean vector of Y ~ MVN is zero
    xxc <- cov2cor(H %*% R %*% H)^2
    xxs <- xxc[lower.tri(xxc)]
    ## Compute theoretical variance of the sum of squared and scaled values (xx2)
    ( v <- 2*sum(w^2) )
    2*sum(diag(xxc)) + 4*sum(xxs) # check that it is equal to v
    var(xx2) # empirical value, should approach 2*sum(w^2)

    ## theoretical and empirical mean
    (m <- sum(w)); mean(xx2)

    ## if G ~ gamma(shape=sh, scale=sc), then mean(G) = sh*sc, var(G) = ash*sc^2
    ## given mean and var, solve to get 'sh' and 'sc'
    ( sh <- m^2 / v ) # gamma shape
    ( sc <- v / m ) # gamma scale
    ## get approximate P-values
    p <- 1 - pgamma(xx2, shape = sh, scale = sc)
    ## check that approximately p ~ U(0,1)
    hist(p[1 : (oo/5),])
    ecdf(p)(0.05)

    ## check that xx^2 follow 1 df chisquare marginally
    p <- 1 - pchisq(xx^2, df=1)
    hist(p[1 : (oo/5),])

    ## Now convert P-values to 2 df chi-squares (Fisher's method)
    x2 <- qchisq(1-p, df=2)

    ## the following covariances should be close
    cov(xx^2); cov(x2/sqrt(2)); cov(x2)/2

    ## theoretical matrix 'xxc' should be close to empirical cor(xx^2)
    xxc; cor(xx^2)

    ## sum the 2 df chisquares and compute approximate P-values
    x2sum <- x2 %*% rep(1,L)
    ## p2 and p3 are equivalent, and also equivalent to
    ## 1-pgamma(x2sum/sc, shape = 2*sh, scale = 1)
    p2 <- 1 - pgamma(x2sum, shape = 2*sh, scale = sc)
    p3 <- 1 - pchisq(x2sum / (0.5*sc), df = 4*sh)
    head(p2); head(p3)
    hist(p2[1 : (oo/5),])

    ## compute combined P via sum -log(P-values) via gamma or chisquare
    xlog.sum <- -log(p) %*% rep(1,L)
    ## Compare: under independence, i.e., when P~U(0,1) and iid, 
    ## then shape = L and scale = 1,
    ## or equivalently, sh = L/2 and sc = 2
    p.log.1 <- 1 - pgamma(xlog.sum, shape = 2*sh, scale = sc/2)
    p.log.2 <- 1 - pchisq(xlog.sum / (0.25*sc), df = 4*sh)
    head(p.log.1); head(p.log.2); head(p3)
    ##
    ## the same:
    ## xlog.sum <- ( -log(p) / (0.5*sc) ) %*% rep(1,L)
    ## p.log.3 <- 1 - pgamma(xlog.sum, shape = 2*sh, scale = 1)
    ## head(p.log.3)
}
