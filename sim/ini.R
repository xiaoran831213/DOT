library(CompQuadForm)
## library(SKAT)
library(MASS)
## library(Matrix)
## library(mvtnorm)

if(!exists("c17"))
{
    c17 <- readRDS('17q12.rds')
    n17 <- nrow(c17)
    m17 <- ncol(c17)
    f17 <- colMeans(c17) / 2
}

for(. in dir("R", "[.]R$", ful=TRUE))
{
    source(.)
}

for(. in dir(".", "[.]R$"))
{
    if(. == "ini.R")
        next
    source(.)
}
