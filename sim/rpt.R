library(ggplot2)

.th <- theme(
    axis.title.x=element_blank(), axis.title.y=element_blank(), 
    strip.text.x = element_text(size=12, face="bold"),
    strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    legend.title=element_blank(), legend.position='bottom')

## cap the values
.cp <- function(dat, grp, val='val', cap=0.01, mtd=c('both', 'upper', 'lower'))
{
    grp <- split(dat, dat[, grp])
    mtd <- match.arg(mtd, c('both', 'upper', 'lower'))
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        if(mtd == 'upper')
            v <- pmin(v, quantile(v, 1-cap, na.rm=TRUE))
        else if(mtd == 'lower')
            v <- pmax(v, quantile(v, 0+cap, na.rm=TRUE))
        else
        {
            v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
            v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        }
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

get.rpt <- function(sim, cache=TRUE)
{
    rds=paste0(sim, '.rds')
    if(file.exists(rds) && cache)
        agg <- readRDS(rds)
    else
    {
        agg <- lapply(dir(sim, 'rds$', full=TRUE), function(f)
        {
            print(f)
            readRDS(f)
        })
        agg <- do.call(rbind, agg)
        saveRDS(agg, rds)
    }
    invisible(agg)
}

get.pow <- function(rpt)
{
    grp <- subset(rpt, se=-c(ngv, pvl, seed))
    rpt <- by(rpt, grp, function(g)
    {
        cfg <- subset(g, se=-c(seed, pvl))[1, ]
        pow <- with(g, mean(pvl <= 0.05))
        cbind(cfg, pow=pow)
    })
    rpt <- do.call(rbind, rpt)
    rpt
}

plt.pow <- function(sim, out=paste0(sim, '.pdf'))
{
    rpt <- get.rpt(sim)
    rpt <- get.pow(rpt)
    g <- ggplot(rpt, aes(x=N, y=pow))
    g <- g + geom_line(aes(color=mtd))
    g <- g + facet_grid(key ~ tag)
    g <- g + .th

    nfy <- length(unique(rpt$key))
    ufy <- 10 / nfy
    nfx <- length(unique(rpt$tag))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    options(bitmapType = 'cairo', device = 'pdf')
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy), scale=0.7)
    invisible(rpt)
}
