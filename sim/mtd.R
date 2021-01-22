#' SKAT
#'
#' @param y: vector of phenotype
#' @param g: genotype matrix
#' @param int: use intercept (def=no)?
#' @return the P-value
skt <- function(y, g, ...)
{
    obj <- SKAT_Null_Model(y ~ 1)
    skt <- SKAT(g, obj, 'linear', is_check_genotype = TRUE, is_dosage = FALSE)
    list(P=skt$p.value)
}

#' multi-regression analysis
#'
#' @param y: vector of phenotype
#' @param g: genotype matrix
#' @param int: use intercept (def=no)?
#' @return the degree of freedom, the  sum of squre, the F-test statistics, and
#'     the P-value of F test.
mra <- function(y, g, int=1, ...)
{
    g <- imp.avg(g)
    if(int > 0)
        ret <- anova(lm(y ~ g + 1))[1, c(1, 3, 4, 5)]
    else
        ret <- anova(lm(y ~ g - 1))[1, c(1, 3, 4, 5)]
    names(ret) <- c('df', 'ssq', 'F', 'P')
    ret
}


tsq <- function(Z, C, D=NULL, tol.egv=NULL, ...)
{
    tol.egv <- tol.egv %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1
    dim(Z) <- NULL

    C <- eigen(C, TRUE, TRUE)$values
    D <- eigen(D, TRUE, TRUE)$values
    
    d <- kronecker(D, C)
    dim(d) <- NULL

    ## . <- d > max(d) * tol.egv
    ## if(!all(.))
    ##     d <- d[  .]
    L <- length(d)              # effective number of eigen

    ## P <- imhof(sum(Z^2), d, delta=rep(0, L))$Qq
    P <- davies(sum(Z^2), lambda = d)$Qq
    list(d=d, L=L, P=P)
}

#' False Discovery Rate
#'
#' @param P p-values
fdr <- function(P, ...) list(P=min(p.adjust(P, 'fdr')))

#' bonferroni correction
#'
#' @param P p-values
bon <- function(P, ...) list(P=min(p.adjust(P, 'bon')))

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
