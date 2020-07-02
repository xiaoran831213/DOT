#' Schur complement
#'
#' Given a full matrix \code{X} and a target block \code{C} within,
#' 
#' X = [A  B]
#'     [B' C],
#'
#' calculate the Schur complement \eqn{X/C = A - BC^{-1}B'}
#' 
#' @param X the whole matrix
#' @param C mask or index of the target block
#' @return the Schur complement of block C of matrix X
#' @noRd
scp <- function(X, C)
{
    ## C can be names, numbers, or Boolean masks, turn them into index
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
#' Given a positive  semi-definite (PSD) matrix \code{X},  \code{nsp} gives the
#' negative square root of \code{X}, that is, \eqn{(R R')^{-1} = X}.
#'
#' For  matrix \code{X}  that  is  not PSD,  \code{nsp}  truncates \code{X}  and
#' operates on the PSD part of \code{X}.
#' 
#' @param X supposedly positive definite matrix
#' @return negative square root of \code{X}
#' @noRd
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


#' Impute missing genotype calls
#'
#' Fill missing values considering the correlated variants
#'
#' The imputation can be seen as regressing  each of the \code{M} variant on the
#' rest of  \code{M - 1} variants,  and taking the average  prediction of missed
#' calls from \code{M - 1} linear models.
#'
#' @param g genotype matrix in allele dosage format
#' @return the same matrix with imputed dosage values
#' @noRd
imp <- function(g)
{
    P <- ncol(g)
    N <- nrow(g)
    
    ## NA mask, and pairwise completion count
    Z <- is.na(g)
    C <- N - crossprod(Z)
    a <- colMeans(g, na.rm=TRUE)
    
    ## center and scale X_j, for j = 1 ... P (the SNPs)
    x <- sweep(g, 2, a, `-`)
    x[Z] <- 0
    
    ## d_ij =  x_i' x_j *  (n_ij / c_ij),  where n_ij =  N, and c_ij  counts the
    ## complete pairs;
    ## d_ij resembles a distance measure.
    D <- crossprod(x) * (N / C)

    ## simple regression, use x_i to product x_j
    ## hat{x_j}(x_i) = b_ij x_i, needs regression coefficient b_ij
    ## b_ij = (x_i' x_j) / (x_i' x_i) = d_ij / d_ii, i,j = 1 ... P
    B <- D / diag(D)
    
    ## for a sample h, the ith SNP
    ## h_i = sum_{j=1, j!=i} x_j b_ji / (non-na x_j count)
    diag(B) <- 0
    x <- x %*% B / (P - 1)
    x <- sweep(x, 2, a, `+`)
    
    ## break imputed values into discrete dosage
    h <- matrix(0L, nrow(x), ncol(x))
    for(j in seq(ncol(g)))
    {
        h[x[, j] > max(x[g[, j] == 0, j], -Inf, na.rm=TRUE), j] <- 1L
        h[x[, j] > max(x[g[, j] == 1, j], -Inf, na.rm=TRUE), j] <- 2L
    }
    g[Z] <- h[Z]
    g
}


#' Correlation of association statistics
#'
#' Calculate the correlation among genetic association test statistics
#'
#' When no  covariate was present  in the association  analysis (\code{x=NULL}),
#' correlation among the  test statistics equals the genotype  LD or correlation
#' \code{cor(g)}.
#'
#' With covariates, \code{css()} gives Schur complement of the genotype block in
#' cor(cbind(g, x)) -- the correlation among all predictors.
#'
#' @param g matrix of genotype, one row per sample, one column per variant;
#' @param x matrix of covariates, one row per sample, one column per variable.
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## correlation among association statistics, covariates involved
#' sgm <- css(gno, cvr)
#' print(sgm[1:4, 1:4])
#' @export
css <- function(g, x=NULL)
{
    g <- imp(g) # impute missed calls
    
    if(is.null(x))
        r <- cor(g) # no covariate, use full cor
    else
    {
        r <- cor(cbind(x, g))   # full cor
        i <- seq(ncol(x))       # index of covariates
        r <- cov2cor(scp(r, i)) # cond cor
    }
    r
}
