pvaldistance <- function(x, method = c("ks", "cvm"), dist.to = c("uniform")){
    stopifnot(is.numeric(x), is.character(method), is.character(dist.to))

    method <- match.arg(method)
    dist.to <- match.arg(dist.to)

    if(dist.to == "uniform"){
        x <- sort(as.vector(x))
        num <- length(x)

        switch(method,
            ks = max(abs(sweep(
                x = cbind(c(0, x), c(x, 1)),
                MARGIN = 1,
                STATS = seq(0, 1, 1/num)
            ))),
            cvm = sum(((1:num)/num - 1/(2*num) - x)^2) + 1/(12*num)
        )
    } else {
        NA
    }
}
