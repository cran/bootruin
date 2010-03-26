rpdataboot <- function(x, b, method = c("nonp", "exp", "lnorm")){
    stopifnot(is.numeric(x), is.numeric(b), is.character(method))

    method <- match.arg(method)

    x <- as.matrix(x)
    dx <- dim(x)
    xbad <- is.na(x) | is.nan(x) | is.infinite(x)
    is.na(x) <- xbad

    x.boot <- switch(
        method,
        nonp = apply(
            X = array(x, dim = c(dx, b)),
            MARGIN = c(2, 3),
            FUN = function(arg){
                ok.arg <- !is.na(arg)
                arg[ok.arg] <- sample(x = arg[ok.arg], size = sum(ok.arg), replace = TRUE)
                return(arg)
            }
        ),
        exp = array(
            data = rexp(
                n = prod(dx, b),
                rate = rep(1 / colMeans(x, na.rm = TRUE), times = b, each = dx[1])),
            dim = c(dx, b)
        ),
        lnorm = array(
            data = rlnorm(
                n = prod(dx, b),
                meanlog = rep(colMeans(lx <- log(x), na.rm = TRUE), times = b, each = dx[1]),
                sdlog = rep(sd(lx, na.rm = TRUE), times = b, each = dx[1])
            ),
            dim = c(dx, b)
        )
    )

    is.na(x.boot) <- array(xbad, dim = c(dx, b))
    return(x.boot)
}
