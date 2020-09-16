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


#' Impute by Average Variant
#'
#' @param g genotype matrix in allele dosage format
#' @return matrix with imputed dosage values
imp.avg <- function(g)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- mean(x[-i])
        x
    })
}

#' Impute by Conditional Expectation
#'
#' Soft imputation of missing values considering the correlation among variants
#'
#' @details
#' The imputation can be seen as predicting  each of the \code{M} variant by the
#' rest \code{M  - 1} variants,  and taking  the weighted average  prediction of
#' missed calls.
#'
#' @param g genotype matrix in allele dosage format
#' @return matrix with imputed dosage values
#' @noRd
imp.cxp <- function(g)
{
    N <- nrow(g); M <- ncol(g); Z <- is.na(g); C <- !Z

    .zs <- function(g)
    {
        ## column center; setting NA to 0 after opt out incomplete pairs.
        x <- as.matrix(scale(g, TRUE, FALSE))
        x[Z] <- 0

        ## sum_i(x_ir * x_is * i_rs), sum of product between of x_r and x_s,
        ## complete pairs only
        xy <- crossprod(x)

        ## sum_i(x_ir * x_ir * i_rs), sum of square of x_r in the complete
        ## pairs of x_r and x_s
        xx <- crossprod(x^2, C)
        
        ## sum_i(x_ir * x_is * i_rs) / sum_i(x_ir  * x_ir * i_rs) = xy / xx_rs,
        ## the regression coeficient
        rc <- xy / xx

        ## sum_i(x_ir * i_rs) / sum_i(i_rs), mean  of x_r in the complete pairs
        ## of x_r and x_s
        rs <- crossprod(x, C) # the sum
        nn <- crossprod(C)    # the non-NA count
        mu <- rs / nn         # the mean

        ## sum_i(x_is * x_is * i_rs) - 2 * rc sum_i(x_ir * x_is * i_rs) + rc^2 *
        ## sum_i(x_ir * x_ir * i_rs), squared residual
        e2 <- t(xx) - xy * xy / xx

        ## denominator
        d2 <- xx - 2 * rs * mu + mu^2
        
        ## squared standard error
        s2 <- e2 / d2 / (nn - 2)
        diag(s2) <- (N - diag(nn)) / ((diag(nn) - 2) * (diag(nn) - 1))
        ## z-scores
        s2[s2 <= 0 ] <- min(s2[s2 > 0]) # avoid s2=0 caused by perfect fit
        rc / sqrt(s2)
    }
    
    p <- g
    x <- array(c(g == 0 & C, g == 1 & C, g == 2 & C), c(N, M, 3)) * 1
    
    ## g <- imp.mod(g)
    c_0 <- Inf # consistancy
    w <- 1

    ## complete pairs, for each genotype {0, 1, 2}
    n <- array(apply(x, 3, crossprod, C), c(M, M, 3))

    ## weights
    w <- .zs(p)^2
    w <- w / mean(diag(w))

    ## transition
    y <- array(0, c(N, M, 3))
    for(i in 1:3)
    {
        n_i <- crossprod(x[, , i], C) # pairwise complete count for x == 0, 1, or 2
        r_i <- 1/ n_i
        r_i[is.infinite(r_i)] <- 0
        for(j in 1:3)
        {
            y[, , j] <- y[, , j] + x[, , i] %*% (w * crossprod(x[, , i], x[, , j]) * r_i)
        }
    }

    ## balance the contribution of predictors.
    y <- y / rowSums(C)
    y <- y / array(rowSums(y, dims=2), c(N, M, 3))
    
    ## imputation
    p <- 0 * y[, , 1] + 1* y[, , 2] + 2 * y[, , 3]

    g[Z] <- p[Z]
    g
}


#' Correlation among association test statistics
#'
#' Calculates the correlation among genetic association test statistics.
#'
#' @details
#' When no covariates are present in per-variant association analyses, that is,
#' \code{x==NULL},  correlation  among  test  statistics is  the  same  as  the
#' correlation among variants, \code{cor(g)}.
#'
#' With  covariates, correlation  among  test  statistics is  not  the same  as
#' \code{cor(g)}. In this case, \code{\link{cst}} takes the generalized inverse
#' of the entire correlation matrix, \code{corr(cbind(g, x))}, and then inverts
#' back only the submtarix containing genotype variables, \code{g}.
#'
#' Missed genotype calls are fill  with 'soft' (i.e., continuous) allele dosage
#' values  in the  interval from  0 to  2, imputed  using information  from all
#' non-missing entries in the genotype matrix.
#'
#' By default (\code{imp=0}), missed calls for each variant are filled with the
#' average of non-missing entires; setting \code{imp=1} uses advance techniques
#' based on expected allele count of one variant conditioned on other variants.
#'
#' The advance imputation (\code{imp=1})  improves statistical power in certian
#' cases, but  we highly suggest  the user to  rerun the association  test with
#' imputed genotype to avoid type I error.
#'
#' @param  g matrix of  genotype, one row per  sample, one column  per variant,
#'     missings allewed.
#' @param x matrix of covariates, one row per sample, no missing allowed.
#' @param imp  imputation methods,  0=naive averge,  1=conditional expectation
#'     (def=0).
#'
#' @return a list containing
#' \itemize{
#' \item{G}{imputed genotype matrix, no missing entreis}
#' \item{C}{the correlation matrix \code{cor{G}}}
#' }
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## correlation among association statistics, covariates involved
#' res <- cst(gno, cvr)
#' print(res$C[1:4, 1:4])
#'
#' ## with 2% missed calls
#' g02 <- readRDS(system.file("extdata", 'rs208294_g02.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#' res <- cst(g02, cvr, imp=0)
#' print(res$C[1:4, 1:4])
#' 
#' @export
cst <- function(g, x=NULL, imp=0)
{
    ## impute missed calls
    if(imp == 1)
        imp <- imp.cxp
    else
        imp <- imp.avg
    g <- imp(g)
    
    if(is.null(x))
        r <- stats::cor(g)              # no covariate, use full cor
    else
    {
        r <- stats::cor(cbind(x, g))    # full cor
        i <- seq(ncol(x))               # index of covariates
        r <- stats::cov2cor(scp(r, i))  # cond cor
    }
    list(G=g, C=r)
}
