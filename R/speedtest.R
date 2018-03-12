# up to 600
nc = 100
nr = 4000
ncluster = nc / 10
cluster_col <- function(i) c(rnorm(nr/2), rnorm(nr/2, 10, 1))

x <- c(sapply(seq_len(ncluster), cluster_col), rnorm(nr*(nc-ncluster)))
x <- matrix(x, nrow = nr, ncol = nc, byrow = F)
x <- as.data.frame(x)
x[,nc-1] <- factor(sample(c("a", "b"), size = nr, replace = T))

plot(x$V1, x$V2)

# r <- find_kmedoids_weights(data = x, k = 2, minimal_memory_mode = T, bounds = c(0,1))
system.time(r <- find_weights(data = x, minimal_memory_mode = T, bounds = c(0,1)))
barplot(r$par)

# Oh wait, perhaps there's an intermediary option for the memory issue. What if instead of either keeping everything in memory or computing everything on the fly, we save the matrices in temporary files?
# I'll look into that tonight / next week. Hopefully it can solve the performance issue, or at least help with it.


# var_list <- list(data = x, bounds = c(0,1), method = "average", "MMM" = T)
# system.time(print(cpcc_derivative(rep(1 / nc, nc-1), var_list)))

library(magrittr)
library(data.table)
b <- as.data.table(x)

dist = cluster::daisy(x, metric = "gower")
cop <- cophenetic(hclust(dist, "average"))
cor(cop, dist)

gower_diff <- function(col){
    if(is.numeric(col) || is.logical(col)){
        c_min = min(col, na.rm = T)
        c_max = max(col, na.rm = T)
        # col <- (col - c_min) / (c_max - c_min)
        #Equal to daisy result if rescaling is not applied
        return(sapply(col, function(y) abs(y-col)) / (c_max - c_min))
    }else{
        return(sapply(col, function(y) y != col))
    }
}

gower <- function(data, weights = rep(1, ncol(data))){
    res <- matrix(0, nrow = nrow(data), ncol = nrow(data))

    for (col in seq_along(data)){
        res <- res + gower_diff(data[,col] * weights[col])
    }
    return(res / sum(weights))
}

N <- 10
system.time(a_or <- as.matrix(cluster::daisy(x = x[,1:N,drop=F], 'gower')))
system.time(a_mod <- gower(x[,1:N]))

a_or[1,1:10]
a_mod[1,1:10]

system.time(a_or <- as.matrix(cluster::daisy(x = x[,1,drop=F], 'gower')))
system.time(a_mod <- gower_diff(x[,1]))

a_or[1:3,1:10]
a_mod[1:3,1:10]

hist(abs(a_or - a_mod))
