#' Create data frame
#'
#' A simplified  wrapper of R  data.frame(...), which also ensure  the character
#' variables are not turned into factors.
.d <- function(...) data.frame(..., stringsAsFactors=FALSE)

#' Flood named objects in a container to an environment
#'
#' \code{flood} quickly assign items in a list to the calling environment,
#'
#' So syntax like:
#'   a <- result$a
#'   b <- result$b
#'   ...
#'   z <- result$z
#'
#' becomes: \code{flood(result)}
#'
#' Be careful with the silent overwriting of existing object.
flood <- function(., x=parent.frame())
{
    for(n in names(.)) assign(n, .[[n]], x)
    invisible(NULL)
}

#' get function calling arguments
#'
#' Use \code{get.arg} in a function's body body to capture the calling arguments
#' in a  one-row \code{data.frame}. The  get.arg() is useful in  documenting the
#' simulation configuration which is usually passed in as arguments, such as the
#' desired sample size, the number of variables (i.e., SNPs), the size of noise,
#' etc.
#'
#' \code{get.arg}  is   not  recommended  for  functions   accepting  non-scalar
#' arguments such as matrix of genotype or vector of effects.
#'
#' @return a data.frame of function arguments 
get.arg <- function()
{
    a <- as.list(match.call(sys.function(1), sys.call(1), expand.dots=TRUE))
    a <- lapply(a[-1], function(.)
    {
        switch(class(.), call=deparse(.), name=as.character(.), .)
    })
    do.call(data.frame, c(a, list(stringsAsFactors=FALSE)))
}

