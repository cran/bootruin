rppvalue <- function(x, method = c("bootstrap", "normal"), x.boot){
    stopifnot(is.numeric(x), is.character(method))

    method <- match.arg(method)
    if(method != "bootstrap" && !missing(x.boot)){
        warning(paste("The argument x.boot is ignored for method = ", method, ".", sep = ""))
    }
    switch(method,
        normal = pnorm(q = x, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),
        bootstrap = {
            stopifnot(!missing(x.boot))
            stopifnot(is.numeric(x.boot))
            if(!is.matrix(x)){
                x <- matrix(x, nrow = 1)
            }
            if(!is.matrix(x.boot)){
                x.boot <- matrix(x.boot, nrow = 1)
            }
            x <- array(x, dim = c(dim(x), dim(x.boot)[2]))
            return(drop(rowMeans(sweep(x, 2:3, x.boot, ">="), dims = 2)))
        }
    )
}
