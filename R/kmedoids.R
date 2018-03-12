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


cost_function = function(w, var_list){
    data = var_list$data
    k = var_list$k
    if(1 - sum(w) < var_list$bounds[1] || 1 - sum(w) > var_list$bounds[2]) return(sum(cluster::daisy(data, "gower")) ^ 2)
    w = c(w, 1 - sum(w))

    dist_matrix = cluster::daisy(data, "gower", weights = w)
    message('determining cost...')
    if(!anyNA(dist_matrix)){
        clustering = cluster::pam(x = dist_matrix, k = k, do.swap = F)
        z = 0
        dist_matrix = as.matrix(dist_matrix)
        for(C in 1:k){
            center = clustering$id.med[C]
            assigned = clustering$clustering == C
            z = z + sum(dist_matrix[center,assigned])
        }
        return(z)
    }else return(sum(cluster::daisy(data, "gower"), na.rm = TRUE) ^ 2)
}

cost_function_derivative = function(w, var_list){
    data = var_list$data
    k = var_list$k
    if(1-sum(w)<var_list$bounds[1] || 1 - sum(w) > var_list$bounds[2]) return(numeric(ncol(data) - 1))
    w = c(w, 1 - sum(w))

    message('calcuting derivative...')
    dist_matrix = cluster::daisy(data, "gower", weights = w)
    if(!anyNA(dist_matrix)){
        cluster = cluster::pam(x = dist_matrix, k = k, do.swap = F)
        dist_matrix = as.matrix(dist_matrix)
        fk = var_list$fk

        process_col <- function(Df) {
            sum(sapply(1:k, function(C)
                sum(Df[cluster$id.med[C], cluster$clustering == C]))
            )}

        if(!is.list(fk)){
            message("using custom gower diff")
            z <- numeric(ncol(data) - 1)
            for(col in seq_len(ncol(data) - 1)){
                fki <- gower_diff(data[,col])
                z[col] <- process_col(fki)
            }
        }else{
            z = sapply(fk, process_col)
        }
        dsum = sum(sapply(1:k,function(C)
            sum(dist_matrix[cluster$medoids[C], cluster$clustering == C]))
        )
        z = z - dsum
        return(z)
    }else return(numeric(ncol(data) - 1))
}

find_kmedoids_weights = function(data, k, start_values = rep(1 / ncol(data), ncol(data) - 1),
                                 n_iterate = 3, minimal_memory_mode = T,
                                 bounds=c(1 / (3 * ncol(data)), 1 - (ncol(data) - 1) / (3 * ncol(data)))){
    #CHECK IF THE INPUT IS CORRECT
    #BASIC INPUT: DATA, CORRECT TYPES, AND CLUSTER PACKAGE
    if (!requireNamespace("cluster", quietly = T))    stop("Package cluster not installed")
    if (!is.data.frame(data) && !is.matrix(data))     stop("Error: data was not supplied as a dataframe or matrix")

    #STARTING VALUES
    if (1 - sum(start_values) < bounds[1] )           stop("Error: starting values cant sum to one. ")
    if (any(start_values < 0))                        stop("Error: negative starting values")

    #QUASI NEWTON ARGUMENTS
    if (n_iterate < 1)                                stop("Error: invalid number of iterations")

    #BOUNDARY CHECKS
    if (bounds[1] < 0)                                stop("Error: lower bound must be positive or 0")
    if (1 - (ncol(data) - 1) * bounds[1] < bounds[1]) stop("Error: lower bound too high. It must be below or equal to 1/ncol(data)")
    if (bounds[2] < bounds[1])                        stop("Error: lower bound must be lower than upper bound")

    #SETUP NECESSITIES
    method = "L-BFGS-B"
    bounds[2] = min(1 - (ncol(data) - 1) * bounds[1], bounds[2])
    control_list = list("maxit" = n_iterate, "fnscale" = 50, factr = .01, pgtol = .001)

    #CALCULATE FK
    if(minimal_memory_mode) fk <- NA
    else {
        fk = list()
        for(i in 1:(ncol(data) - 1)) fk[[i]] = gower_diff(data[,i])
    }

    #INITIATE
    var_list = list(data, k, bounds, fk)
    names(var_list) = c("data","k","bounds","fk")
    message('Starting optimisation')
    RES = optim(par = start_values, fn = cost_function, gr = cost_function_derivative,
                var_list, control = control_list, method = method, lower = bounds[1], upper = bounds[2], hessian = F)
    RES$par = c(RES$par, 1 - sum(RES$par))
    return(RES)
}

opt_kmedoids = function(data, k = k, start_values = rep(1/ncol(data), ncol(data)-1),n_iterate = 10,
                        bounds = c(1 / (3 * ncol(data)), 1 - (ncol(data) - 1)/(3 * ncol(data))) # this equals: [1/3n, (n-1)/3n]
){

    #WRAPPER FOR FIND_KMEDOIDS_WEIGHTS: RETURNS A HCLUST OBJECT AND THE RESULTS OF FIND_WEIGHTS.
    if(!requireNamespace("cluster",quietly = T))    stop("Package cluster not installed")

    RES = find_kmedoids_weights(data, k = k, start_values = start_values,
                                n_iterate = n_iterate, bounds = bounds)

    D = cluster::daisy(data, "gower", weights = RES$par)
    kmedoids = cluster::pam(x = D, k = k)
    res = list(kmedoids, RES)
    names(res) = c("clustering","opt_result")
    return(res)
}

