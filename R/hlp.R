#' Decolinearity by Variant Trimming 
#'
#' Retain only one of highly colinear variants.
#'
#' Decolinearity by variant trimming bring the correlation closer to full-rank,
#' Low rank LD correlation matrix causes inflated type I error for some summary
#' statistics.
#'
#' Tow or  more variants with  correlation higher  than \code{1 -  tol.cor} are
#' considered colienar, and  \code{dvt} marks all but one of  them for removal.
#' The    default   tolerence    of    high   correlation    (\code{tol.cor})is
#' \code{.Machine$double.eps}
#' 
#' @param C matrix of LD-correlation;
#' @param tol.cor tolerence of high correlation
#'
#' @return boolean mask of variants to retain.
#' @noRd
dvt <- function(C, tol.cor=NULL, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    
    C <- abs(C)
    C[upper.tri(C, TRUE)] <- 0
    apply(C, 2, function(.) all(. < 1 - tol.cor))
}


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


#' Negative square root of a positive definite matrix
#'
#' Given a positive  semi-definite (PSD) matrix \code{X},  \code{nsp} gives the
#' negative square root of \code{X}, that is, \eqn{(R R')^{-1} = X}.
#'
#' For matrix \code{X}  that is not PSD, \code{nsp} trims  \code{X} and retains
#' \code{L} only positive eigenvalues.
#' 
#' @param X supposedly positive definite matrix
#' @return a list
#' \itemize{
#' \item{\code{W}: }{negative square root  of \code{X}}
#' \item{\code{L}: }{the effective number of positve eigenvalues}
#' }
#' @noRd
nsp <- function(X, L=NULL, eps=NULL, ...)
{
    if(is.null(eps))
        eps <- sqrt(.Machine$double.eps)
    ## eigen decomposition
    . <- svd(X)
    u <- .$u
    v <- .$v
    d <- .$d
    
    ## positive eigen values
    . <- d > d[1] * eps
    if(!all(.))
    {
        d <- d[  .]
        u <- u[, .]
        v <- v[, .]
    }
    L <- length(d)
    
    ## square root
    d <- sqrt(1/d)
    Y <- u %*% (d * t(v))       # U diag(d) V'
    Y <- 0.5 * (Y + t(Y))
    
    list(W=Y, L=L)
}
