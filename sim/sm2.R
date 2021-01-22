#' Simulation main function
#'
#' @param N sample size
#' @param M variant count
#' @param hsq heritability
#' @param frq fraction of casual variants
#' @param psd positive semi-definite threshold
#' @param times repeats of simulation
#' @return a data frame of simulation configuration and p-values
#'
#' If one  draw multivariate normal  from the  LD matrix (mvn=1),  the genotype
#' will have continuouse values instead of {0,  1, 2}, and each SNP is centered
#' at mu=0.
sm2 <- function(N=5e2, M=2, L=5, a=0, b=1, e=1, na=.1, times=5e2, ...)
{
    arg <- get.arg(skp=c("times"))
    set.seed(arg$seed)

    mdl <- Y(M) ~ a@G + b@X
    ## replicate many times
    rpt <- list()
    for(i in seq(times))
    {
        cat(sprintf("iter = %4d, ", i))
        G <- kgp(N, L, ...)
        X <- kgp(N, L, ...)
        flood(sim_rsp(mdl)$dat)
        Y <- Y + matrix(rnorm(N * M), N, M) * e # outcomes

        if(na > 0)
        {
            G[sample(N * L, N * L * na)] <- NA
            Y[sample(N * M, N * M * na)] <- NA
        }
        
        ## gwas
        gwa <- gwa_lm(Y ~ X + {G})
        lhc <- cor.std(gwa$rsp) # left  hand correlation
        rhc <- cor.std(gwa$gmx) # right hand correlation
        zsc <- gwa$zsc

        ## test
        r <- list()
        r[['udt']] <- dot_chisq(zsc, rhc, lhc, ...)
        r[['utq']] <- tsq(zsc, rhc, lhc, ...)
        
        p <- sapply(r, `[[`, 'P')
        l <- sapply(r, `[[`, 'L')
        rpt[[i]] <- data.frame(itr=i, mtd=names(r), pvl=p, egv=l)
        cat("\n")
    }
    set.seed(NULL)
    rpt <- cbind(arg, do.call(rbind, rpt))
    rownames(rpt) <- NULL
    pow(rpt)
}
