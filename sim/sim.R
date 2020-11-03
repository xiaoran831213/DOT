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
sim <- function(N=5e2, M=20, naf=.1, cof=cor.std, hsq=.01, frq=.1, msk=1, wgt=0, times=5e2, ...)
{
    arg <- get.arg()
    arg$times <- NULL
    set.seed(arg$seed)
    stopifnot(msk == 0 || frq < 1)

    ## replicate many times
    rpt <- list()
    for(i in seq(times))
    {
        cat(sprintf("iter = %4d, ", i))
        dat <- try.gen(N, M, NAF=naf, hsq=hsq, frq=frq, msk=msk, wgt=wgt, ...)
        flood(dat) # gmx, tld, wgt, ...
        
        ldm <- cof(gmx, tld)
        r <- list()
        r[['wch']] <- dot_chisq(zsc, ldm, w=wgt, abs=0, ...)
        r[['ach']] <- dot_chisq(zsc, ldm, w=wgt, abs=1, ...)
        r[['uch']] <- dot_chisq(zsc, ldm, w=NL, ...)
        ## r[['ufs']] <- dot_fisher(zsc, ldm, w=NL,  ...)
        ## r[['wfs']] <- dot_fisher(zsc, ldm, w=wgt, ...)
        ## r['dot_tpm'] <- dot_tpm(zsc, ldm)$P
        ## r[['dot_art']] <- dot_art(zsc, ldm, k=M, ...)
        ## r[['dot_arta']] <- dot_arta(zsc, ldm, k=M, ...)
        ## r[['dot_rtp']] <- dot_rtp(zsc, ldm, k=M, ...)
        r[['wtq']] <- tsq(zsc, ldm, w=wgt, ...)
        r[['utq']] <- tsq(zsc, ldm, w=NL, ...)
        ## r['mra'] <- mra(rsp, gmx)$P
        ## r['fdr'] <- fdr(pvl)$P
        ## r['bon'] <- bon(pvl)$P
        ## r['skt'] <- skt(rsp, gmx)$P
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

#' summerize statistical power
#'
#' @param rpt a simulation report in data frame
pow <- function(rpt)
{
    rpt <- subset(rpt, se=-itr)
    grp <- subset(rpt, se=-c(pvl, egv))
    rpt <- by(rpt, grp, function(g)
    {
        cfg <- subset(g, se=-c(pvl, egv))[1, ]
        pow <- with(g, mean(pvl <= 0.05))
        egv <- with(g, mean(egv))
        cbind(cfg, pow=pow, egv=egv, rep=nrow(g))
    })
    rpt <- do.call(rbind, rpt)
    rpt
}
