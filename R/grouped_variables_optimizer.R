grouped_cpcc = function(x,var_list){
    if (!requireNamespace("cluster", quietly = T))    stop("Package cluster not installed")
    data = var_list$data
    method = var_list$method
    weights = numeric(ncol(data))
    combined_indices = var_list$combined_indices
    i = 0
    if(length(combined_indices) > 1){
        for(i in 1:(length(combined_indices) - 1)) {
            weights[combined_indices[[i]]] = x[i]
        }
    }

    weights[combined_indices[[i + 1]]] = (1 - sum(weights))/length(combined_indices[[i + 1]])

    if(any(weights < var_list$bounds[1]) || any(weights > var_list$bounds[2])) return(0)

    dist_matrix = cluster::daisy(data, metric = "gower", weights = weights)
    z = 0
    if(!any(is.na(dist_matrix))) z = - cor(dist_matrix, cophenetic(hclust(dist_matrix, method = method)))
    return(z)
}

grouped_cpcc_derivative = function(x,var_list){
    if (!requireNamespace("cluster", quietly = T))    stop("Package cluster not installed")
    data = var_list$data
    combined_indices = var_list$combined_indices
    weights = numeric(ncol(data))
    i = 0
    if(length(combined_indices) > 1) {
        for(i in 1:(length(combined_indices) - 1)) {
            weights[combined_indices[[i]]] = x[i]
        }
    }

    weights[combined_indices[[i + 1]]] = (1 - sum(weights)) / length(combined_indices[[i + 1]])
    if(any(weights < var_list$bounds[1]) || any(weights > var_list$bounds[2])) {
        return(numeric(length(var_list$combined_indices) - 1))
    }
    dist_matrix = cluster::daisy(data, metric = "gower", weights = weights)

    if(!any(is.na(dist_matrix))){
        method = var_list$method
        minimal_memory_mode = var_list$MMM

        cophenetic_matrix = cophenetic(hclust(dist_matrix, method = method))
        cophenetic_matrix = cophenetic_matrix - mean(cophenetic_matrix, na.rm = T)

        dist_matrix = dist_matrix - mean(dist_matrix, na.rm = T)

        sumsquared = sum(dist_matrix ^ 2, na.rm = T)
        root = sqrt(sumsquared * sum(cophenetic_matrix ^ 2))
        sumsquared = sum(cophenetic_matrix * dist_matrix) / sumsquared
        fin = (dist_matrix * sumsquared - cophenetic_matrix) / root
        rm(root, sumsquared, dist_matrix, cophenetic_matrix)


        if(!minimal_memory_mode) {
            fk = var_list$fk
            return(sapply(fk, function(x)(sum(fin * x, na.rm = T))))
        }else{
            z = numeric(length(combined_indices) - 1)
            for(i in 1:(length(combined_indices) - 1)){
                fk = cluster::daisy(as.data.frame(data[, combined_indices[[i]][1]]), metric  = "gower")
                if(length(combined_indices[[i]]) > 1){
                    for(j in 2:length(combined_indices[[i]])) {
                            fk = fk + cluster::daisy(as.data.frame(data[, combined_indices[[i]][j]]), metric  = "gower")
                        }
                }
                z[i] = sum(fin * fk, na.rm = T)
        }
        return(-z)
        }
    }else return(numeric(length(var_list$combined_indices) - 1))
}

find_grouped_weights = function(data, combined_indices = as.list(1:ncol(data)),
                                start_values = rep(1 / ncol(data), length(combined_indices) - 1),
                                n_iterate = 10, clust_method = "average",
                                bounds = c(1 / (3 * ncol(data)),1 - (ncol(data) - 1) / (3 * ncol(data))),
                                minimal_memory_mode = F){
    #CHECK IF THE INPUT IS CORRECT
    #BASIC INPUT: DATA, CORRECT TYPES, AND CLUSTER PACKAGE
    if (!requireNamespace("cluster", quietly = T))    stop("Package cluster not installed")
    if (!is.data.frame(data) && !is.matrix(data))     stop("Error: data was not supplied as a dataframe or matrix")
    if (!is.logical(minimal_memory_mode))             stop("Error: minimal memory mode must be logical")
    if(length(combined_indices) <= 1) stop("Provide more than 1 group of variables.")

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
    control_list = list(0, 1, c(1,1), 0.001, n_iterate, 0.01, 0.01, NA, 10, NA,5, 0.01, 0.001, NA, NA)
    var_list = list(data, clust_method, combined_indices, bounds, minimal_memory_mode)
    names(var_list) = c("data", "method", "combined_indices", "bounds", "MMM")

    Co = list(0,  1,   c(1,1),  0.001,  n_iterate, 0.01,  0.01,   NA,  10, NA,5,  0.01, 0.001,  NA,   NA)
    if(!minimal_memory_mode){
        print("Computing Fk(i,j)...")
        fk = list(1:(length(combined_indices) - 1))

        #THE ONLY DIFFERENCE IN CALCULATION LIES HERE. FOR EACH GROUP WE SUM THE FK'S, THE REST OF THE COMPUTATION IS THE SAME
        for(i in 1:(length(combined_indices) - 1)){
            fki = 0
                for(j in 1:length(combined_indices[[i]])){
                    fki = fki + cluster::daisy(as.data.frame(data[,combined_indices[[i]][j]]),metric  = "gower")
                }
            fk[[i]] = fki
        }
        print("Finished computing Fk(i,j)...")
        var_list$fk = fk
        rm(fk)
    }
    RES = optim(par = start_values, fn = grouped_cpcc, gr = grouped_cpcc_derivative, var_list, control = control_list,
        method = method, lower = bounds[1], upper = bounds[2], hessian = F)
    i = 0
    weights = numeric(ncol(data))
    if(length(combined_indices) > 1){
        for(i in 1:(length(combined_indices) - 1)) {
            weights[combined_indices[[i]]] = RES$par[i]
        }
    }
    weights[combined_indices[[i + 1]]] = (1 - sum(weights)) / length(combined_indices[[i + 1]])
    RES$par = weights
    RES$value = -RES$value
    RES
}

opt_grouped_hierarchical = function(data, combined_indices = as.list(1:ncol(data)),
                                    start_values = rep(1 / ncol(data), length(combined_indices) - 1),
                                    n_iterate = 10, clust_method = "average",
                                    bounds = c(1 / (3 * ncol(data)),1 - (ncol(data) - 1) / (3 * ncol(data))),
                                    minimal_memory_mode = F){

    #WRAPPER FOR FIND_GROUPED_WEIGHTS: RETURNS A HCLUST OBJECT AND THE RESULTS OF FIND_WEIGHTS.
    if(!requireNamespace("cluster",quietly = T))    stop("Package cluster not installed")

    RES = find_grouped_weights(data, combined_indices = combined_indices, start_values = start_values,
                       n_iterate = n_iterate, clust_method = clust_method,
                       bounds = bounds, minimal_memory_mode = minimal_memory_mode)

    D = cluster::daisy(data, "gower", weights = RES$par)
    hierarch = hclust(D, clust_method)
    res = list(hierarch, RES)
    names(res) = c("clustering","opt_result")
    return(res)
}
