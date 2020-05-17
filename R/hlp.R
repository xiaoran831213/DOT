#' Multivariate Normal Samples
#'
#' A simplified copy of \code{MASS::mvrnorm()}.
#'
#' @param N number of samples to be drawn
#' @param  M vector  of mean,  to be rotated  to the  length of  sample features
#'     \code{D}
#' @param V square matrix of covariance, with \code{D} dimensions
#' @param drop TRUE to drop matrix of single sample to vector
#' 
#' @return matrix of N row samples and D column features, whose
#' covariance is \code{V}
mvn <- function (N=1, M=0, V=NULL, drop=TRUE)
{
    ## default V is for demonstration
    if(is.null(V))
        V = rbind(c(1, .5, .25), c(.5, 1, .5), c(.25, .5, 1))

    ## dimensionality
    D <- nrow(V)

    ## mean vector
    M <- drop(rep(M, length=D))

    ## eigen decomposition
    e <- eigen(V, symmetric = TRUE)
    d <- e$values
    U <- e$vectors
    s <- sqrt(pmax(d, 0))          # square root of V

    ## random values
    X <- matrix(rnorm(D * N), D, N)
    y <- t(M + U %*% (s * X))

    if(drop)
        y <- drop(y)
    y
}


#' Schur complement of matrix
#'     
#' X   = [A  B]
#'       [B' C]
#'
#' X/C = A - BC^{-1}B'
#' @param X the whole matrix
#' @param C mask or index of block C in matrix X
#' @return the Schur complement of block C of matrix X
scp <- function(X, C)
{
    ## C can be names and masks, turn them into index
    i <- seq(ncol(X))
    names(i) <- colnames(X)
    i <- unname(i[C])

    ## A <- X[-C, -C]
    ## B <- X[-C, +C]
    ## B'<- X[+C, -C]
    ## C <- X[+C, +C]

    ## A - B C^{-1} B'
    X[-i, -i] - X[-i, +i] %*% solve(X[+i, +i], X[+i, -i])
}


#' Negative square root of positive definite matrix
#'
#' @param X the positive definite matrix
#' @return Z, such that `(Z Z')^{-1} = X`
nsp <- function(X)
{
    eps <- sqrt(.Machine$double.eps)
    ## eigen decomposition
    . <- svd(X)
    u <- .$u
    v <- .$v
    d <- .$d

    ## positive eigen values
    . <- d > d[1] * eps
    d <- d[  .]
    u <- u[, .]
    v <- v[, .]

    ## square root
    d <- sqrt(1/d)
    Y <- u %*% (d * t(v))       # U diag(d) V'
    Y <- 0.5 * (Y + t(Y))
    Y
}


#' Correlation of genotype variants
#'
#' Given two variants `x_i` and `x_j`,  the allele frequency `p_i` and `p_j`, and
#' the haplotype frequency `p_{ij}`, their correlation is
#'
#' `r_{ij} = \frac{p_{ij} - p_i p_j}{\sqrt{p_i(1 - p_i)p_j(1 - p_j)}}`
#'
#' @param x matrix of genotype
cog <- function(x)
{
    PA_ <- colMeans(x) / 2
    P_B <- PA_
    PAB <- (crossprod(x > 0) + crossprod(x > 1)) / (2 * nrow(x))
    DAB <- PAB - outer(PA_, P_B) # LD
    VA_ <- PA_ * (1 - PA_)       # VA
    V_B <- P_B * (1 - P_B)       # VB
    DAB / sqrt(outer(VA_, V_B))  # cor(AB)
}
