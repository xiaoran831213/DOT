#' Impute missing genotype
#'
#' Fill missed genotype entries with guesses based on observed entries.
#'
#' The most  naive approach to  fill the  missings in a  variant is to  use the
#' average  of observed  entries  (\link{\code{imp_avg}}). Besides  simplicity,
#' imputation by  average has  the advantage  of approximating  the correlation
#' among test  statistics (i.e.,  Z-scores) when  the association  analysis was
#' done with missing values unfilled, which is a common practices (i.e., GWAS).
#' Therefore, one does not have to  re-run the assocation analysis with variant
#' imputed  by average.   Also, the  correlation calculation  \link{\code{cst}}
#' always uses \link{\code{imp_avg}}.
#'
#' Mode \link{\code{imp_mod}} and median \link{\code{imp_med}} are as simple as
#' the naive average  in that they only use information  in the target variant,
#' but account  for the  fact that  allele dosage is  discrete and  not normal.
#' However, one  must re-run the  assocation analysis with variants  imputed by
#' mode or median, and use the  new Z-scores for decorrelated statistical tests
#' (see \link{dot_sst}), or otherwise risking inflated type 1 error.
#' 
#' Advanced  approach such  as conditional  expectation (\link{\code{imp_cnd}})
#' explore  the  relationship  between  variants and  borrow  information  from
#' variants other  than the  target to  make a  guess.  The  sample correlation
#' among variants imputed  this way is closer  to the true LD,  and may improve
#' power.  Again, be adviced that one must re-run the association analysis with
#' imputed variants to avoid inflated type I error.
#' 
#' @param g genotype matrix, one row per sample, and one column per variant.
#' @param ... additional parameters.
#'
#' @return imputed genotype matrix.
#'
#' @name imp
NULL

#' @describeIn imp
#'
#' imputation by average.
imp_avg <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- mean(x[-i])
        x
    })
}

#' @describeIn imp
#'
#' imputation by median, re-run the association analysis with imputed genotype!
imp_med <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- median(x[-i])
        x
    })
}

#' @describeIn imp
#' 
#' imputation by mode, re-run the association analysis with imputed genotype!
imp_mod <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        u <- unique(x[-i])
        a <- tabulate(match(x[-i], u))
        x[i] <- rep(u[a == max(a)], l=length(i))
        x
    })
}


#' @describeIn imp
#'
#' imputation by conditional expectation,  re-run the association analysis with
#' imputed genotype!
imp_cnd <- function(g, ...)
{
    N <- nrow(g); M <- ncol(g); Z <- is.na(g); C <- !Z

    .zs <- function(g) # z-scores
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
