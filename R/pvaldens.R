pvaldens <- function(x, bw, rho, method = c("jh", "chen")){
    stopifnot(is.numeric(x))

    x <- as.vector(x)
    x <- x[!is.na(x) & is.finite(x)]

    method <- match.arg(method)
    switch(method,
        jh = {
            if(missing(rho)){
                rho <- ifelse(missing(bw), 0.9, 1 - bw^2)
            } else {
                if(any(rho <= 0, rho >= 1)){
                    stop("'rho' must have a value between 0 and 1.")
                }
            }

            K.jh <- function(u, v){
                if(any(u <= 0, v <= 0, u >= 1, v >= 1)){
                    return(0)
                } else {
                    1 / sqrt(1 - rho^2) *
                    exp(-rho / (2 * (1 - rho^2)) *
                    (rho * qnorm(u)^2 + rho * qnorm(v)^2 -
                    2 * qnorm(u) * qnorm(v)))
                }
            }

            return(
                function(y){
                    rowMeans(outer(y, x, Vectorize(K.jh)))
                }
            )
        },
        chen = {
            if(missing(bw)){
                bw <- 0.01
            }
            if(!missing(rho)){
                warning(paste("The parameter 'rho' is not supported by the method '", method, "'.", sep = ""))
            }

            return(
                function(y){
                    colMeans(outer(x, y, function(x, y){
                        dbeta(y, x / bw + 1, (1 - x) / bw + 1)
                    }), na.rm = TRUE)
                }
            )
        }
    )
}
