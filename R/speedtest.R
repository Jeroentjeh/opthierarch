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
