rpdataconv <- function(x){
    stopifnot(all(sapply(x, is.numeric)))
    if(is.vector(x) && !is.list(x)){
        return(matrix(x, ncol = 1))
    } else {
        if(is.matrix(x)){
            return(x)
        } else {
            x <- lapply(x, as.vector)
            x.len <- sapply(x, length)
            max.len <- max(x.len)
            if(max.len == 1){
                warning("All entries of x have length 1.")
                x <- unlist(unname(x))
                return(matrix(x, nrow = 1))
            } else {
                x.fill <- lapply(max.len - x.len, rep.int, x = NA)
                x <- mapply(FUN = c, x, x.fill)
                return(unname(as.matrix(as.data.frame(x))))
            }
        }
    }
}
