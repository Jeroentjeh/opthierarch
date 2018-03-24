calculate_cpcc = function(x, var_list){
    message('run fn...')

    data = var_list$data
    method = var_list$method
    dist_matrix = cluster::daisy(data, metric = "gower", weights = x)

    #CHECK IF THERE ARE NO PROBLEMS WITH THE DISTANCE MATRIX
    if (!any(is.na(dist_matrix))){
        correlation = cor(dist_matrix, cophenetic(hclust(dist_matrix, method = method)))
    } else correlation = 0

    return(correlation)
}

cpcc_derivative = function(x, var_list){
    data = var_list$data
    message('run derivative...')
    dist_matrix = cluster::daisy(data, metric = "gower", weights = x)


    #CHECK IF THERE ARE NO PROBLEMS WITH THE DISTANCE MATRIX
    if (!any(is.na(dist_matrix))){
        #SETUP NECESSARY VARIABLES BEFOREHAND
        method = var_list$method
        minimal_memory_mode = var_list$MMM

        #CALCULATE NECESSARY VARIABLES BEFOREHAND
        cophenetic_matrix = cophenetic(hclust(dist_matrix, method = method))
        cophenetic_matrix = cophenetic_matrix - mean(cophenetic_matrix, na.rm = T)
        dist_matrix = dist_matrix - mean(dist_matrix, na.rm = T)
        sumsquared = sum(dist_matrix ^ 2, na.rm = T)
        root = sqrt(sumsquared * sum(cophenetic_matrix ^ 2))
        sumsquared = sum(cophenetic_matrix * dist_matrix) / sumsquared
        fin = (dist_matrix * sumsquared - cophenetic_matrix) / root
        rm(root, sumsquared, dist_matrix, cophenetic_matrix)

        #COMPUTE DERIVATIVE
        if (!minimal_memory_mode) {
            fk = var_list$fk
            return(-sapply(fk, function(x)(sum(fin * x, na.rm = T))))
        }else{
            cl <- var_list$cluster
            if(!is.null(cl)){
                return(-parallel::parSapply(cl, data, function(x){
                    sum(fin * gower_diff(x), na.rm = T)
                }))
            } else {
                return(-sapply(data, function(x){
                    sum(fin * gower_diff(x), na.rm = T)
                }))
            }
        }
    }else return(numeric(ncol(data)))
}

find_weights = function(data, start_values = rep(1 / ncol(data), ncol(data)),
                        n_iterate = 10, clust_method = "average",
                        bounds = c(1 / (3 * ncol(data)), 1 - (ncol(data) - 1)/(3 * ncol(data))), # this equals: [1/3n, (n-1)/3n]
                        minimal_memory_mode = FALSE, use_cluster = FALSE){

    #CHECK IF THE INPUT IS CORRECT
    #BASIC INPUT: DATA, CORRECT TYPES, AND CLUSTER PACKAGE
    if (!requireNamespace("cluster", quietly = T))    stop("Package cluster not installed")
    if (!is.data.frame(data) && !is.matrix(data))     stop("Error: data was not supplied as a dataframe or matrix")
    if (!is.logical(minimal_memory_mode))             stop("Error: minimal memory mode must be logical")

    #STARTING VALUES
    # if (1 - sum(start_values) < bounds[1] )           stop("Error: starting values cant sum to one. ")
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
    control_list = list(maxit = n_iterate, fnscale = -1, pgtol = .05, factr = 1e-7)
    var_list = list('data' = data, 'method' = clust_method, 'bounds' = bounds, 'MMM' = minimal_memory_mode)

    #SETUP CLUSTER
    if (class(use_cluster) != "cluster")
        c1 = parallel::makeCluster(parallel::detectCores() - 1, type = "FORK")
    else c1 = use_cluster
    if(use_cluster[1] != F) {
        var_list$cluster <- c1
    }

    #START
    if (!minimal_memory_mode){
        print("Computing Fk(i,j)...")

        if (class(use_cluster) == "cluster" || (class(use_cluster) == "logical" && use_cluster == T)){
            print("Using cluster...")
            if(!requireNamespace("parallel", quietly = T)) stop("Package parallel not installed")


            if (class(use_cluster) == "SOCKcluster"){
                parallel::clusterExport(cl = c1, varlist = "data", envir = environment())
            }

            #COMPUTE FK IN PARALLEL
            fk = parallel::parLapply(c1, data, gower_diff)
        } else {
            #COMPUTE FK NOT PARALLEL
            fk = lapply(data, gower_diff)
        }
        var_list$fk <- fk
        print("Finished computing Fk(i,j)...")
        print("Starting optimization...")
        RES = optim(par = start_values, fn = calculate_cpcc, gr = cpcc_derivative, var_list, control = control_list,
                    method = method, lower = bounds[1], upper = bounds[2], hessian = F)
        rm(fk)
    }else {
        print("Starting optimization...")
        RES = optim(par = start_values, fn = calculate_cpcc, gr = cpcc_derivative, var_list, control = control_list, method = method,
                    lower = bounds[1], upper = bounds[2], hessian = F)
    }

    # Shut down cluster
    if(class(use_cluster) == "cluster" || (class(use_cluster) == "logical" && use_cluster == T)) parallel::stopCluster(c1)
    names(RES$par) = colnames(data)
    RES$par <- RES$par / sum(RES$par)
    # RES$value = -RES$value
    return(RES)
}

opt_hierarchical = function(data, start_values = rep(1/ncol(data), ncol(data)-1),
                            n_iterate = 10, clust_method = "average",
                            bounds = c(1 / (3 * ncol(data)), 1 - (ncol(data) - 1)/(3 * ncol(data))), # this equals: [1/3n, (n-1)/3n]
                            minimal_memory_mode = FALSE, use_cluster = FALSE){

    #WRAPPER FOR FIND_WEIGHTS: RETURNS A HCLUST OBJECT AND THE RESULTS OF FIND_WEIGHTS.
    if(!requireNamespace("cluster",quietly = T))    stop("Package cluster not installed")

    RES = find_weights(data, start_values = start_values,
                       n_iterate = n_iterate, clust_method = clust_method,
                       bounds = bounds, minimal_memory_mode = minimal_memory_mode,
                       use_cluster = use_cluster)

    D = cluster::daisy(data, "gower", weights = RES$par)
    hierarch = hclust(D, clust_method)
    res = list(hierarch, RES)
    names(res) = c("clustering","opt_result")
    return(res)
}
