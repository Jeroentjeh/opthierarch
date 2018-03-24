gower_diff <- function(col){
    N <- length(col)
    if(is.numeric(col) || is.logical(col)){
        x <- sapply(seq_len(length(col) - 1), function(y) {
            return(col[y] - col[(y+1):N])
        })
        x <- unlist(x)
        return(abs(x) / (max(col) - min(col)))
    }else{
        distances <- sapply(col, function(y) y != col)
        x <- sapply(seq_len(length(col) - 1), function(y) {
            return(col[y] != col[(y+1):N])
        })
        return(unlist(x))
    }
}
