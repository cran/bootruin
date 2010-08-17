rpjack <- function(x, ...){
    stopifnot(is.numeric(x))

    if(is.array(x)){
        return(apply(x, 2:length(dim(x)), rpjack, ...))
    } else {
        x <- as.vector(x)
        bad <- is.na(x) | is.nan(x) | is.infinite(x)
        x <- x[!bad]
        num <- length(x)
        X <- matrix(
            rep(x, num)[-(num*(0:(num - 1)) + 1:num)],
            ncol = num,
            nrow = num - 1,
            byrow = FALSE
        )
        return(sd(ruinprob(x = X, ...)) * (num - 1) / sqrt(num))
    }
}
