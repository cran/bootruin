ruinprob <- function(x, param.list, compmethod = c("dg", "exp"), flmethod = c("nonp", "exp", "lnorm", "custom"), reserve, loading, fl = NA, interval = 0.5, ...){
    stopifnot(is.numeric(x))

    if (!missing(param.list) && is.list(param.list)) {
        try({
                compmethod <- param.list$compmethod
                flmethod   <- param.list$flmethod
                reserve    <- param.list$reserve
                loading    <- param.list$loading
                fl         <- param.list$fl
                interval   <- param.list$interval
            },
            silent = TRUE
        )
    }

    stopifnot(reserve  >= 0.0,
              loading  >= 0.0,
              interval >  0.0)

    if (is.array(x)) {
        apply(X          = x,
              MARGIN     = 2L:length(dim(x)),
              FUN        = ruinprob,
              compmethod = compmethod,
              flmethod   = flmethod,
              reserve    = reserve,
              loading    = loading,
              interval   = interval,
              fl         = fl,
              ...)
    } else {
        x <- as.vector(x)
        x <- x[is.finite(x)]

        compmethod <- match.arg(compmethod)
        flmethod   <- match.arg(flmethod)

        if (flmethod == "custom") {
            stopifnot(is.function(fl))
        } else {
            if (is.function(fl)) {
                warning(paste("The option 'fl' is ignored for this method (flmethod = '", flmethod, "').", sep = ""))
            }
        }

        switch(compmethod,
            #exp = 1 / (1 + loading) * exp(-reserve * loading / (mean(x) * (1 + loading))),
            exp = dexp(reserve * loading / mean(x), 1.0 / (1.0 + loading)),
            dg  = {
                psi.0 <- 1.0 / (1.0 + loading)
                num   <- floor(reserve / interval) + 1L

                if(flmethod == "lnorm"){
                    lx      <- log(x)
                    mymu    <- mean(lx)
                    mysigma <- sd(lx)
                    myseq   <- seq.int(from = 0.0, by = interval, length.out = num)
                } else {
                    myseq <- seq.int(from = 0.0, by = interval, length.out = num + 1L)
                }

                h.l <- switch(
                    flmethod,
                    nonp   = diff(sapply(X   = myseq,
                                         FUN = function(myarg){
                                             sum(pmin(myarg, x))
                                         }) / sum(x)),
                    #exp   = diff(pexp(myseq, rate = 1 / mean(x), lower.tail = TRUE, log.p = FALSE)),
                    exp    = diff(pexp(myseq, 1.0 / mean(x))),
                    lnorm  = sapply(X   = myseq,
                                    FUN = function(x) {
                                        exp(-mymu - mysigma^2.0 / 2.0) *
                                        integrate(f          = stats::plnorm,
                                                  lower      = x,
                                                  upper      = x + interval,
                                                  meanlog    = mymu,
                                                  sdlog      = mysigma,
                                                  lower.tail = FALSE)$value
                                    }),
                    custom = diff(sapply(myseq, fl, ...))
                )

                #h.u <- c(0, h.l)

                if (reserve < interval) {
                    if (reserve == 0.0) {
                        return(psi.0)
                    } else {
                        #lower.limit <- 1 - f.l[1]
                        #upper.limit <- psi.0
                        return(psi.0 * (1.0 + (1.0 - h.l[1L])/(1.0 - psi.0 * h.l[1L])) / 2.0)
                    }
                } else {
                    .C(
                        'rplimits',
                        h.l   = as.double(h.l),
                        #h.u  = as.double(h.u),
                        psi.0 = as.double(psi.0),
                        num   = as.integer(num),
                        rp    = double(1L)
                    )$rp
                }
            }
        )
    }
}
