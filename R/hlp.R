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


#' Schur complement
#'
#' Given a full matrix \code{X} and a target block \code{C} within,
#' 
#' X   = [A  B]
#'       [B' C],
#'
#' \code{scp} calculate the Schur complement \eqn{X/C = A - BC^{-1}B'}
#' 
#' @param X the whole matrix
#' @param C mask or index of the target block
#' @return the Schur complement of block C of matrix X
scp <- function(X, C)
{
    ## C can be names, numbers, or bollean masks, turn them into index
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
#' Given a poitive  (simi) definite (PSE) matrix \code{X},  \code{nsp} gives the
#' negative squre root of \code{X}, that is, \eqn{(R R')^{-1} = X}.
#'
#' For  matrix \code{X}  that  is  not PSD,  \code{nsp}  truncates \code{X}  and
#' operates on the PSD part of \code{X}.
#' 
#' @param X supposedly positive definite matrix
#' @return negative square root of \code{X}
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


#' LD-correlation of genotype variants
#'
#' variant correlation based on Linkage  Disequilibrium (LD).
#' 
#' Given  two genetic  variants \eqn{x_i}  and \eqn{x_j},  the allele  frequency
#' \eqn{p_i}  and \eqn{p_j},  and  the haplotype  frequency \eqn{p_{ij}},  their
#' LD-correlation is
#'
#' \deqn{r_{ij} = \frac{p_{ij} - p_i p_j}{\sqrt{p_i(1 - p_i)p_j(1 - p_j)}}}
#'
#' \code{ldc} is related to Linkage  Disequilibrium (LD), that is, the numerator
#' \eqn{p_{ij} - p_i  p_j} is the LD between variant  \code{i} and \code{j}. The
#' LD-correlation is different from the typical 
#'
#' @param x genotype matrix of N row samples and M column variants.
#' @export
ldc <- function(x)
{
    PA_ <- colMeans(x) / 2
    P_B <- PA_
    PAB <- (crossprod(x > 0) + crossprod(x > 1)) / (2 * nrow(x))
    DAB <- PAB - outer(PA_, P_B) # LD
    VA_ <- PA_ * (1 - PA_)       # VA
    V_B <- P_B * (1 - P_B)       # VB
    DAB / sqrt(outer(VA_, V_B))  # cor(AB)
}


#' Impute missing genotype calls
#'
#' Fill missing values considering the correlated variants
#'
#' The imputation can be seen as regressing  each of the \code{M} variant on the
#' rest of  \code{M - 1} variants,  and taking the average  prediction of missed
#' calls from \code{M - 1} linear models.
#'
#' @param g genotype matrix in allele dosage formate
#' @return the same matrix with imputed dosage values
#' @export
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

    ## simple regression, use x_i to preduct x_j
    ## hat{x_j}(x_i) = b_ij x_i, needs regresion coeficient b_ij
    ## b_ij = (x_i' x_j) / (x_i' x_i) = d_ij / d_ii, i,j = 1 ... P
    B <- D / diag(D)
    
    ## for a sample h, the ith SNP
    ## h_i = sum_{j=1, j!=i} x_j b_ji / (non-na x_j count)
    diag(B) <- 0
    x <- x %*% B / (P - 1)
    x <- sweep(x, 2, a, `+`)

    ## break imputed values into descrete dosage
    h <- matrix(0L, nrow(x), ncol(x))
    for(j in seq(ncol(g)))
    {
        h[x[, j] > max(x[g[, j] == 0, j], na.rm=TRUE), j] <- 1L
        h[x[, j] > max(x[g[, j] == 1, j], na.rm=TRUE), j] <- 2L
    }
    g[Z] <- h[Z]
    g
}


#' Conditional variant correlation
#'
#' calculate the correlation between  genotype variants (\code{g}) conditioned on
#' covariates (\code{x}).
#'
#' With no covariate (\code{x=NULL}),  the conditioned variant correlation (cvc)
#' is identical  to the raw  variant correlation \code{cor(g)};  with covariate,
#' the conditional
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## conditional variate correlation
#' sgm <- cvc(gno, cvr)
#'
#' ## result of per-variant test (p-values and effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## DOT  summarized chi-square  statistics based  on z-scores  and conditional
#' ## variant correlation
#' rpt <- dot(stt, sgm)
#'
#' print(rpt$ssq)  # chi-square = 35.76306
#' print(rpt$pvl)  #    p-value =  0.001132132
cvc <- function(g, x=NULL)
{
    g <- imp(g) # impute missed calls
    
    if(is.null(x))
        r <- cor(g) # no covariate, cond cor = full cor
    else
    {
        r <- cor(cbind(x, g))   # full cor
        i <- seq(ncol(x))       # index of covariates
        r <- cov2cor(scp(r, i)) # cond cor
    }
    r
}
