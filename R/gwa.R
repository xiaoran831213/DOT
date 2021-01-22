#' Compute Polynomials
#'
#' @param x matrix of basic terms.
#' @param d degrees to compute.
#' @param t number of terms allowed.
#' @param o use orthoganal transformation (def=FALSE).
#' @return polynomial terms expanded from the given basic terms.
.p <- function(x, d=NULL, t=NULL, o=FALSE)
{
    if(!is.matrix(x))
        dim(x) <- c(length(x), 1L)
    d <- d %||% 1
    t <- t %||% seq.int(max(d))

    ## get polynomial
    x <- poly(x, degree=max(d), raw=!o, simple=TRUE)

    ## select by degree
    x <- x[, attr(x, 'degree') %in% d, drop=FALSE]

    ## select by terms
    . <- rowSums(do.call(rbind, strsplit(colnames(x), "[.]")) != "0")
    x <- x[, . %in% t, drop=FALSE]

    ## clean up
    attr(x, 'class') <- NULL
    attr(x, 'degree') <- NULL

    x
}

pcs <- function(x) x %*% svd(x)$v


#' dosage genotype to factors
#'
#' Recode `aa = 0` to `(1, 0)`, `Aa = 1` to `(1, 1)`, and `AA=2` to `(0, 1)`.
#' @param g matrix of genotype.
#' @return factorized genotype.
fac <- function(g)
{
    g <- cbind(g < 2, g > 0)
    m <- apply(g, 2, any)
    if(!all(m))
        g <- g[, m]
    g <- g + 0
    g <- g[ , apply(g, 2, sd) > 0]
    g
}


gwa_frm <- function(model, ...)
{
    ## separate left and right side
    lh <- frm_lhs(model, 1)
    rh <- frm_rhs(model, 1)
    en <- environment(model)
    
    ## from the r-hand side, separate genotype term(s)
    pt <- "[{][^{}]*[}]"
    gt <- regmatches(rh, gregexpr(pt, rh), invert=0)[[1]]
    gt <- paste(collapse="+", gsub("^[{]|[}]$", "", gt[gt != ""]))

    ## remainig r-hand side, no genotype
    rh <- regmatches(rh, gregexpr(pt, rh), invert=1)[[1]]
    rh <- gsub("^[+]|[+]$", "", rh)
    rh <- paste(collapse="+", rh[rh != ""])
    
    ## from the l-hand side, separate residual regresion term(s)
    if(is.null(lh))
    {
        lh <- ""
        rr <- ""
    }
    else
    {
        pt <- "-[^-+]*"
        rr <- regmatches(lh, gregexpr(pt, lh), invert=0)[[1]]
        lh <- regmatches(lh, gregexpr(pt, lh), invert=1)[[1]]
        rr <- paste(collapse="+", gsub("-", "", rr[rr != ""]))
        lh <- paste(collapse="+", gsub("^[+]|[+]$", "", lh[lh != ""]))
    }
    
    ## return
    lh <- if(lh=="") NULL else as.formula(paste("~", lh), en)
    rh <- if(rh=="") NULL else as.formula(paste("~", rh), en)
    rr <- if(rr=="") NULL else as.formula(paste("~", rr), en)
    gt <- if(gt=="") NULL else as.formula(paste("~", gt), en)

    list(lh=lh, rh=rh, gt=gt, rr=rr)
}

#' Get regression residual
#'
#' @param Y matrix of responses
#' @param ... covariates
get_rsd <- function(Y, ..., int=TRUE)
{
    if(int)
        X <- cbind(rep(1, NROW(Y)), ...)
    else
        X <- cbind(...)
    if(is.null(X))
        Y
    else
        qr.resid(qr(X), Y)
}

#' GWAS Linear Model
#' @examples
#'
#' set.seed(1245)
#' N <- 100                          # sample size
#' D <- mvn(N, 0, cmx(9))            # 9 predictors in total
#' Z <- D[, 5:6]                     # 2 covariates observed
#' X <- hwe(D[, 1:4])                # 4 g-variants observed
#' D[, 1:4] <- X
#' U <- D[, 7:9]                     # 3 predictors unseen
#' D <- std(D)
#' 
#' b <- matrix(rnorm(9 * 3), 9, 3)   # effects
#' Y <- std(D %*% b)                 # 
#' Y <- Y + mvn(N, 0, cmx(3))        # noise
#'
#' Y <- std(Y)
#' X <- std(X)
#' Z <- std(Z)
#'
#' y1 <- gwa_rsp(Y - Z ~ {X})        # residual response
#' y2 <- gwa_rsp(Y ~ {X} + Z)        # original response
#'
#' round(cor(y1, Z), 2)              # expecting zeros
#' round(cor(y2, Z), 2)              # expecting non-zeros
#'
#' round(colMeans(y1), 2)            # expecting zeros
#' round(colMeans(y2), 2)            # expecting zeros
gwa_lm <- function(model, ...)
{
    ## extract the rsponse, genotype, and covariate
    flood(gwa_frm(model, ...))
    
    rsp <- frm_mtx(lh)      # response
    gmx <- frm_mtx(gt, ...) # genotype
    cvr <- frm_mtx(rh, ...) # covariate
    rsr <- frm_mtx(rr)      # residual regression
    rsp <- std(rsp, TRUE, FALSE)
    gmx <- std(gmx, TRUE, FALSE)
    rsp[is.na(rsp)] <- 0
    gmx[is.na(gmx)] <- 0

    ## residuals
    cvr <- cbind(rep(1, nrow(gmx)), cvr)
    qrd <- qr(cvr)
    rsp <- qr.resid(qrd, rsp) # residual response (traits)
    gmx <- qr.resid(qrd, gmx) # residual genotype

    ## model                                                                    
    mdl <- y ~ g - 1

    ## z-scores                                                                 
    nms <- list(colnames(gmx), colnames(rsp))
    zsc <- matrix(0, ncol(gmx), ncol(rsp), dim=nms)
    pvl <- matrix(1, ncol(gmx), ncol(rsp), dim=nms)
    for(i in seq(ncol(gmx)))
    {
        for(j in seq(ncol(rsp)))
        {
            g <- gmx[, i]
            y <- rsp[, j]
            m <- lm(mdl)
            zsc[i, j] <- summary(m)$coef[1, 3]
            pvl[i, j] <- summary(m)$coef[1, 4]
        }
    }
    list(rsp=rsp, rsr=rsr, cvr=cvr, gmx=gmx, zsc=zsc, pvl=pvl)
}


#' GWA Partial LD
gwa_pld <- function(model, ...)
{
    ## rsponse(lh), genotype(gt), residual regressor(rr) and covariate(rh)
    flood(gwa_frm(model, ...))

    gno <- frm_mtx(gt, ...) # genotype
    cvr <- frm_mtx(rh, ...) # covariate
    rsr <- frm_mtx(rr, ...) # residual regresser
    
    ld <- cor(cbind(rsr, cvr, gno))
    if(!is.null(cvr))
        ld <- cov2cor(scp(ld, seq(ncol(cvr))))
    if(!is.null(rsr))
        ld <- scp(ld, seq(ncol(rsr)))

    nm <- unlist(sapply(lapply(all.vars(gt), get, environment(gt)), colnames))
    colnames(ld) <- nm
    rownames(ld) <- nm
    ld
}
