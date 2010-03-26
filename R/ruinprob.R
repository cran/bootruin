ruinprob <- function(x, param.list, compmethod = c("dg", "exp"), flmethod = c("nonp", "exp", "lnorm", "custom"), reserve, loading, fl = NA, interval = 0.5, ...){
    stopifnot(is.numeric(x))

    if(!missing(param.list)){
        stopifnot(is.list(param.list))
        listcont <- names(param.list)
        if(missing(compmethod) && "compmethod" %in% listcont){
            compmethod <- param.list$compmethod
        }
        if(missing(flmethod) && "flmethod" %in% listcont){
            flmethod <- param.list$flmethod
        }
        if(missing(reserve) && "reserve" %in% listcont){
            reserve <- param.list$reserve
        }
        if(missing(loading) && "loading" %in% listcont){
            loading <- param.list$loading
        }
        if(missing(fl) && "fl" %in% listcont){
            fl <- param.list$fl
        }
        if(missing(interval) && "interval" %in% listcont){
            interval <- param.list$interval
        }
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
            exp = 1 / (1 + loading) * exp(-reserve * loading / (mean(x) * (1 + loading))),
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
                    exp = diff(pexp(myseq, rate = 1 / mean(x), lower.tail = TRUE, log.p = FALSE)),
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

                h.u <- c(0, h.l)

                f.l <- numeric(num)
                f.u <- numeric(num)

                f.l[1] <- (1 - psi.0)/(1 - psi.0 * h.l[1])
                f.u[1] <- 1 - psi.0

                if(reserve < interval){
                    if(reserve == 0){
                        lower.limit <- psi.0
                        upper.limit <- psi.0
                    } else {
                        lower.limit <- 1 - f.l[1]
                        upper.limit <- psi.0
                    }
                } else {
                    fac.l <- psi.0 / (1 - psi.0 * h.l[1])

                    for(i in 2:num){
                        f.l[i] <- h.l[2:i] %*% f.l[(i-1):1] * fac.l
                        f.u[i] <- h.u[2:i] %*% f.u[(i-1):1] * psi.0
                    }

                    lower.limit <- max(0, min(1 - sum(head(f.l, -1)), 1))
                    upper.limit <- max(0, min(1 - sum(f.u), 1))
                }

                return(mean(c(lower.limit, upper.limit)))
            }
        )
    }
}
