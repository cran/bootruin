ruinprob <- function(x, param.list, compmethod = c("dg", "exp"), flmethod = c("nonp", "exp", "lnorm", "custom"), reserve, loading, fl = NA, interval = 0.5, ...){
    stopifnot(is.numeric(x))

    if(!missing(param.list) && is.list(param.list)){
        try(
            {
                compmethod <- param.list$compmethod
                flmethod <- param.list$flmethod
                reserve <- param.list$reserve
                loading <- param.list$loading
                fl <- param.list$fl
                interval <- param.list$interval
            },
            silent = TRUE
        )
    }

    stopifnot(reserve >= 0, loading >= 0, interval > 0)


    if(is.array(x)){
        apply(
            X = x,
            MARGIN = 2:length(dim(x)),
            FUN = ruinprob,
            compmethod = compmethod,
            flmethod = flmethod,
            reserve = reserve,
            loading = loading,
            interval = interval,
            fl = fl,
            ...
        )
    } else {
        x <- as.vector(x)
        bad <- is.na(x) | is.nan(x) | is.infinite(x)
        x <- x[!bad]

        compmethod <- match.arg(compmethod)
        flmethod <- match.arg(flmethod)

        if(flmethod == "custom"){
            stopifnot(is.function(fl))
        } else {
            if(is.function(fl)){
                warning(paste("The option 'fl' is ignored for this method (flmethod = '", flmethod, "').", sep = ""))
            }
        }

        switch(compmethod,
            #exp = 1 / (1 + loading) * exp(-reserve * loading / (mean(x) * (1 + loading))),
            exp = .Internal(dexp(reserve * loading / mean(x), 1 + loading, FALSE)),
            dg  = {
                psi.0 <- 1 / (1 + loading)
                num <- floor(reserve / interval) + 1

                if(flmethod == "lnorm"){
                    lx <- log(x)
                    mymu <- mean(lx)
                    mysigma <- sd(lx)
                    myseq <- seq(from = 0, by = interval, length = num)
                } else {
                    myseq <- seq(from = 0, by = interval, length = num + 1)
                }

                h.l <- switch(
                    flmethod,
                    nonp = diff(sapply(
                        X = myseq,
                        FUN = function(myarg){
                            sum(pmin(myarg, x))
                        }
                    ) / sum(x)),
                    #exp = diff(pexp(myseq, rate = 1 / mean(x), lower.tail = TRUE, log.p = FALSE)),
                    exp = diff(.Internal(pexp(myseq, mean(x), TRUE, FALSE))),
                    lnorm = sapply(X = myseq, FUN = function(x){
                        integrate(
                            f = stats::plnorm,
                            lower = x,
                            upper = x + interval,
                            meanlog = mymu,
                            sdlog = mysigma,
                            lower.tail = FALSE
                        )$value * exp(- mymu - mysigma^2 / 2)
                    }),
                    custom = diff(sapply(myseq, fl, ...))
                )

                #h.u <- c(0, h.l)

                if(reserve < interval){
                    if(reserve == 0){
                        return(psi.0)
                    } else {
                        #lower.limit <- 1 - f.l[1]
                        #upper.limit <- psi.0
                        return(psi.0 * (1 + (1 - h.l[1])/(1 - psi.0 * h.l[1])) / 2)
                    }
                } else {
                    .C(
                        'rplimits',
                        h.l = as.double(h.l),
                        #h.u = as.double(h.u),
                        psi.0 = as.double(psi.0),
                        num = as.integer(num),
                        rp = double(1)
                    )$rp
                }
            }
        )
    }
}
