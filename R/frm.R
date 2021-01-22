## left hand side
frm_lhs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)[-1]
    if(length(f) < 2)
        return(NULL)

    f <- gsub("[\n ]*", "", f[1])
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

## right hand side
frm_rhs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)[-1]
    f <- gsub("[\n ]*", "", f[length(f)]) # right
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

## left bar side
frm_lbs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)
    f <- f[length(f)]
    f <- strsplit(f, "[|]")[[1]][1]
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

## right bar side
frm_rbs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)
    f <- f[length(f)]
    f <- strsplit(f, "[|]")[[1]]
    if(length(f) < 2)
        f <- 0
    else
        f <- f[2]
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

## no intercept?
frm_is0 <- function(f, ...)
{
    attr(terms(f), "intercept") == 0 && length(all.vars(f)) == 0
}

## has intercept1?
frm_is1 <- function(f, ...)
{
    attr(terms(f), "intercept") == 1 && length(all.vars(f)) == 0
}

## get matrix
frm_mtx <- function(f, ...)
{
    if(is.null(f))
        return(NULL)

    o <- getOption("na.action")
    options(na.action="na.pass")
    x <- model.matrix(update.formula(f, ~ . - 1), ...)
    options(na.action=o)
    
    attr(x, "assign") <- NULL
    if(is.logical(x))
        x <- x + 0
    x
}
