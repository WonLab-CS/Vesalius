###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################ HORIZONTAL INTEGRATION ###########################
#-----------------------------------------------------------------------------#


#' Aling and integrate spatial assay from the same modality using super pixels
#' @param seed_assay vesalius_assay object - data to be mapped to
#' @param query_assay vesalius_assay objecy - data to map
#' @param neighborhood character - how should the neighborhood be selected?
#' "knn", "radius", "graph"(See details)
#' @param k int ]2, n_points] number of neareset neighbors to be considered for
#' neighborhodd computation.
#' @param radius numeric ]0,1[ proportion of max distance between points 
#' to consider for the neighborhood
#' @param depth int [1, NA] graph depth from cell to consider from neighborhood
#' (See details)
#' @param dimensions Int vector containing latent space dimensions to use
#' @param batch_size number of points per batch in query during assignment
#' problem solving
#' @param signal character (variable_features, all_features, embeddings, custom)
#' - What should  be used as cell signal to generate the cost matrix.
#' Seed details 
#' @param use_norm character - which count data to use
#' @param scale logical - should signal be scaled 
#' @param threshold score threshold below which indicices should be removed.
#' Scores will always be between 0 and 1
#' @param use_cost character string defining how should total cost be computer
#' Available: feature, niche, territory, composition (See details for combinations
#' and custom matrices)
#' @param custom_cost matrix - matrix of size n (query cells) by p (seed cells)
#' containing custom cost matrix. Used instead of vesalius cost matrix
#' @param verbose logical - should I be a noisy boy?
#' @details The goal is to assign the best matching point between a seed set and
#' a query set.
#' 
#' To do so, \code{map_assays} will first extract a
#' biological signal. This can be latent space embeddings per cell, or by using
#' gene counts (or any other modality).
#'
#' If using gene counts, there are a few more options available to
#' you. First, you can select "variable_features" and vesalius will find the
#' intersection between the variable features in your seed_assay and your
#' query_assay. "all_features" will find the intersection of all genes across
#' assays (even if they are not highly variable). Finally, you can also select
#' a custom gene vector, containing only the gene set you are interested in.
#' 
#' The second step is to create a cost matrix. The creation of a cost matrix
#' is achieved by pair-wise sum of various cost matrices. By default, 
#' the map_assays function will use "feature" and "niche" cost matrices. 
#' The feature matrix computes the pearson correlation between the seed and query
#' using which ever signal was defined by the signal argument (variable_features)
#' will compute the correlation between shared variable features in seed 
#' and query).
#' The niche matrix will be computed by using the pearson correlation between
#' niche expression profiles (based on signal). Niche are defined using the
#' neighborhood argument where knn represent the k nearest neighbors algorithm
#' (with k defining the number of nearest neighbors), depth represents the 
#' graph depth of a local neighborhood graph, and radius defining a spatial
#' radius surrunding a center cell. The singal (expression or embedding) is
#' average across all cells in the niche.
#' The territory matrix will compare the average signal of vesalius 
#' territories between seed and query. 
#' The composition matrix will compute a frequency aware jaccard index
#' between cell types present in a niche. Cell types must be assigned 
#' to seed and query vesalius objects  (See add_cells function)
#' Total cost matrix will be computed by computing the pairwise sum 
#' of the complement (1 - p ) of each cost matrix. 
#'
#' This cost matrix is then parsed to a
#' Kuhn–Munkres algorithm that will generate point pairs that minimize
#' the overall cost. 
#' 
#' Since the algorithm complexity is O(n3), it can be time consuming to
#' to run on larger data sets. As such, mapping will be approximated by
#' dividing seed and query into batches defined by batch size. For an
#' exact mapping ensure that batch_size is larger than the number of cells
#' in both query and seed.
#' 
#' Finaly once the matches are found, the coordinates are mapped to its
#' corresponding point and a new object is returned.
#'
#'
#' 
#' @return vesalius_assay
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Create Vesalius object for processing
#' vesalius <- build_vesalius_assay(coordinates, counts)
#' jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)
#' mapped <- map_assays(vesalius, jitter_ves)
#'}
#' @export
map_assays <- function(seed_assay,
    query_assay,
    neighborhood = "knn",
    k = 20,
    radius = 0.05,
    depth = 1,
    dimensions = seq(1, 30),
    batch_size = 1000,
    signal = "variable_features",
    use_norm = "raw",
    scale = FALSE,
    threshold = 0.3,
    use_cost = c("feature","niche"),
    custom_cost = NULL,
    merge = FALSE,
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # making sure we are formatted to accomodate lapplys and mapplys
    #-------------------------------------------------------------------------#
    custom_cost <- check_cost_matrix_validity(custom_cost)
    
    #-------------------------------------------------------------------------#
    # First let's get singal
    # to minimise the memory print - only one seed signal and list for query
    #-------------------------------------------------------------------------#
    signal <- get_signal(seed_assay = seed_assay,
        query_assay = query_assay,
        signal = signal,
        dimensions = dimensions,
        use_norm = use_norm,
        scale = scale,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # Next we map points in the query assay onto the seed assay
    #-------------------------------------------------------------------------#
    mapped <- point_mapping(
        query_signal = signal$query,
        query_assay = query_assay,
        cost = custom_cost,
        seed_signal = signal$seed,
        seed_assay = seed_assay,
        neighborhood = neighborhood,
        k = k,
        radius = radius,
        depth = depth,
        batch_size = batch_size,
        use_cost = use_cost,
        threshold = threshold,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # This is where we actually create a return object that can be used
    # for integrate analysis 
    # this functions needs to be updated
    #-------------------------------------------------------------------------#
    integrated <- integrate_map(
        aligned_graph = mapped,
        query = query_assay,
        seed = seed_assay,
        merge = merge,
        verbose = verbose)
    
    simple_bar(verbose)
    return(integrated)
}
#' get cell signal from vesalius assays
#' @param seed_assay  vesalius_assay object
#' @param query_assay vesalius_assay object
#' @param signal character string where the signal should be taken from
#' @param dimensions int vector - if signal is embeddings which 
#' embeddings should be selected
#' @param use_norm charcater string which counts should be use when
#' extracting signal
#' @param scale logical - should signal be scaled 
#' @param verbose logical - should progress messages be outputed.
#' @return list contain seed signal, query signal and features used
get_signal <- function(seed_assay,
    query_assay,
    signal,
    dimensions = seq(1:30),
    use_norm = "raw",
    scale = FALSE,
    verbose = TRUE) {
    message_switch("signal", verbose = verbose)
    #-------------------------------------------------------------------------#
    # First we check if we are using embedding values 
    #-------------------------------------------------------------------------#
    if (any(grepl(pattern = "embeddings", x = signal))) {
        seed_signal <- t(seed_assay@active)
        query_signal <- t(query_assay@active)
        features <- NA
    } else {
        #---------------------------------------------------------------------#
        # Let's check what feature input we have and then filter 
        #---------------------------------------------------------------------#
        seed_signal <- check_signal(seed_assay, signal, type = use_norm)
        query_signal <- check_signal(query_assay, signal, type = use_norm)
        features <- intersect(seed_signal, query_signal)
        if (length(features) == 0) {
            stop("No common features between seed and query data sets!")
        }
        seed_signal <- get_counts(seed_assay, type = use_norm)[features, ]
        query_signal <- get_counts(query_assay, type = use_norm)[features, ]
        if (scale) {
            seed_signal <- t(scale(t(as.matrix(seed_signal))))
            query_signal <- t(scale(t(as.matrix(query_signal))))
        }
    }
    return(list("seed" = seed_signal,
        "query" = query_signal,
        "features" = features))

}

#' mapping points between data sets
#' @param query_signal processed query signal from query assay
#' @param query vesalius_assay object
#' @param cost matrix - matrix of size n (query cells) by p (seed cells)
#' containing custom cost matrix. 
#' @param seed_signal processed seed signal from seed assay
#' @param seed vesalius_assay object
#' @param k int size of niche (knn)
#' @param radius 0.05 proportion of max distance to use as radius for 
#' neighborhood
#' @param depth graph path depth to condsider for neighborhood. 
#' @param batch_size int number of points in each query batch
#' @param use_cost character string defining how should total cost be computer
#' Available: feature, niche, territory, composition (See details for combinations
#' @param threshold score threshold below which indicices should be removed.
#' Scores will always be between 0 and 1
#' @param verbose logical - out progress to console
#' @return list of matched and aligned coordinates in query
#' @importFrom future.apply future_lapply
point_mapping <- function(query_signal,
    query_assay,
    cost,
    seed_signal,
    seed_assay,
    neighborhood = "knn",
    k = 20,
    radius = 0.05,
    depth = 1,
    batch_size = 1000,
    use_cost = c("feature", "niche"),
    threshold = 0.5,
    verbose = TRUE) {
    assay <- get_assay_names(query_assay)
    check_cost_validity(cost,
        seed_assay,
        seed_signal,
        query_assay,
        query_signal,
        use_cost)
    #--------------------------------------------------------------------------#
    # Correlation between individual cells 
    #--------------------------------------------------------------------------#
    if (any(grepl("feature", use_cost))){
        message_switch("feature_cost", verbose, assay = assay)
        cost <- c(cost, feature_dist(seed_signal,
            query_signal,
            batch_size))
        names(cost)[length(cost)] <-  "feature"

    }
    #--------------------------------------------------------------------------#
    # Correlation between the cellular niche centered around the cell
    #--------------------------------------------------------------------------#
    if (any(grepl("niche", use_cost))) {
        message_switch("get_neigh", verbose, assay = assay)
        seed_signal_niche <- get_neighborhood(seed,
            seed_signal,
            neighborhood,
            k,
            depth,
            radius)
        query_signal_niche <- get_neighborhood(query,
            query_signal,
            neighborhood,
            k,
            depth,
            radius)
        message_switch("neighbor_cost", verbose, assay = assay)
        cost <- c(cost, feature_dist(seed_signal_niche,
            query_signal_niche,
            batch_size))
        names(cost)[length(cost)] <-  "niche"
    }
    #--------------------------------------------------------------------------#
    # Correlation between territories centered around the cell
    # Note: this can be sped up using much smaller matrices and then dispatchin
    # score after
    #--------------------------------------------------------------------------#
    if (any(grepl("territory", use_cost))) {
        message_switch("territory_cost", verbose, assay = assay)
        seed_signal_niche <- get_neighborhood(seed_assay,
            seed_signal,
            "territory")
        query_signal_niche <- get_neighborhood(query_assay,
            query_signal,
            "territory")
        cost <- c(cost, feature_dist(seed_signal_niche,
            query_signal_niche,
            batch_size))
        names(cost)[length(cost)] <-  "territory"
    }
    #--------------------------------------------------------------------------#
    # Computing nich compisition
    #--------------------------------------------------------------------------#
    if (any(grepl("composition", use_cost))) {
        message_switch("composition_cost", verbose, assay = assay)
        seed_niche <- niche_composition(seed,
            seed_assay,
            method = neighborhood,
            k = k,
            depth = depth,
            radius = radius)
        query_niche <- niche_composition(query,
            query_assay,
            method = neighborhood,
            k = k,
            depth = depth,
            radius = radius)

        cost <- c(cost, list(compare_niches(seed_niche, query_niche)))
        names(cost)[length(cost)] <-  "composition"
    }
    #--------------------------------------------------------------------------#
    # filtering and pairwise addition of cost matrices and
    #--------------------------------------------------------------------------#
    cost <- filter_cost(cost, threshold, use_cost)
    cost <- c(cost, concat_cost(cost, use_cost))
    names(cost)[length(cost)] <- "total_cost"
    #--------------------------------------------------------------------------#
    # devide cost matrix
    #--------------------------------------------------------------------------#
    map_cost <- divide_and_conquer(cost$total_cost, batch_size, verbose)
    matched_indices <- match_index(seq(1, length(map_cost)), map_cost)
    scores <- score_matches(matched_indices,
        cost,
        use_cost = use_cost)
    aligned <- align_index(matched_indices,
        seed,
        query,
        verbose)
    return(list("aligned" = aligned, "prob" = scores))
}

#' concat cost - pairwise sum of score complement
#' @param cost list - named list contained score matrices
#' @param use_cost character - which cost matrices to use 
#' @return list with cost matrix 
concat_cost <- function(cost, use_cost) { 
    if (length(use_cost) == 0) {
        stop("Please specify at least one score matrix to use")
    }
    cost <- cost[use_cost]
    if (length(cost) == 0) {
        stop(paste(paste(use_cost, collapse = " "), ": not available in score matrix list"))
    } else if (length(cost) == 1) {
        return(list(1 - cost[[1]]))
    } else {
        buffer <- 1 - cost[[1]]
        for (i in seq(2, length(cost))){
            buffer <- buffer + (1 - cost[[i]])
        }
        return(list(buffer))
    }

}

#' remove scores that are below a certain threshold 
#' @param cost list - score list to be filtered 
#' @param threshold numeric - score threshold [0,1]
#' @param use_cost - character - which score matrices to filter
#' @return filtered cost list 
filter_cost <- function(cost, threshold = 0.3, use_cost) {
    if (length(use_cost) == 0) {
        stop("Please specify at least one score matrix to use")
    }
    
    if (length(cost[use_cost]) == 0) {
        stop(paste(paste(use_cost, collapse = " "), ": not available in score matrix list"))
    } else {
        #---------------------------------------------------------------------#
        # First check which barcodes have at least one barcodes above 
        # threshold - checking seed and query
        # using a for loop for cleaner replacement 
        #---------------------------------------------------------------------#
        for (i in seq_along(use_cost)) {
            tmp <- cost[[use_cost[i]]]
            above <- apply(tmp, 2, function(i, threshold) {
                return(sum(i >= threshold))
            }, threshold = threshold)
            tmp <- tmp[, above > 0]
            # then check query cost
            above <- apply(tmp, 1, function(i, threshold) {
                return(sum(i >= threshold))
            }, threshold = threshold)
            tmp <- tmp[above > 0, ]
            cost[[use_cost[i]]] <- tmp
           
        }
        #---------------------------------------------------------------------#
        # check which barcodes I should keep across all cost that will be
        # used
        #---------------------------------------------------------------------#
        seed_barcodes <- Reduce(intersect, lapply(cost[use_cost], colnames))
        query_barcodes <- Reduce(intersect, lapply(cost[use_cost], rownames))
        if (length(seed_barcodes) == 0 || length(query_barcodes) == 0) {
            stop("Filtering threhsold to High! No indices retained!")
        }
        for (i in seq_along(use_cost)){
            tmp <- cost[[use_cost[i]]]
            tmp <- tmp[rownames(tmp) %in% query_barcodes,
                colnames(tmp)  %in% seed_barcodes]
            cost[[use_cost[i]]] <- tmp
        }
    }
    
    return(cost)

}


#' compute the distance between seed and query signals
#' @param seed seed signal
#' @param query query signal
#' @param batch_size int describing size of each batch parse to the loops
#' @details Chunking cost and signal into smaller chunks to run the 
#' correlation score in paralell. There is room for improvement here.
#' First we could dispatch the longer list to future_lapply
#' but cannot know which one it is and we need to know that so we can 
#' subset the cost. 
#' Also the functions calls feature_cost which is a R wrapper for a 
#' c++ function (cost.cpp).
#' @return matrix with query as rows and seed as colmuns
feature_dist <- function(seed, query, batch_size) {
    #-------------------------------------------------------------------------#
    # First we chunk seed and query 
    # running a parallel lapply internally
    #-------------------------------------------------------------------------#
    seed_batch <- chunk(seq(1, ncol(seed)), batch_size)
    query_batch <- chunk(seq(1, ncol(query)), batch_size)
    total_cost <- vector("list", length(seed_batch))
    #-------------------------------------------------------------------------#
    # Loop over seed batch - idealy we would use the loop over the 
    # smallest batch 
    #-------------------------------------------------------------------------#
    for (i in seq_along(seed_batch)) {
        #---------------------------------------------------------------------#
        # Creating ref to seed matrix
        # why am I not using as.matrix here?
        #---------------------------------------------------------------------#
        local_seed <- matrix(seed[,seed_batch[[i]]],
            ncol = length(seed_batch[[i]]))
        colnames(local_seed) <- colnames(seed)[seed_batch[[i]]]
        #---------------------------------------------------------------------#
        # computing score in batches
        #---------------------------------------------------------------------#
        local_cost <- future.apply::future_lapply(query_batch,
            function(query_batch, seed, query) {
                #-------------------------------------------------------------#
                # Same as above - not sure why I am using this instead of
                # as.matrix 
                #-------------------------------------------------------------#
                local_query <- matrix(query[, query_batch],
                    ncol = length(query_batch))
                colnames(local_query) <- colnames(query)[query_batch]
                cost <- feature_cost(seed, local_query)
                #-------------------------------------------------------------#
                # this can return Nan when SD is 0 - happens when all counts 
                # are 0. Can happen with the overlap between variable features
                # Will replace with 0 instead 
                #-------------------------------------------------------------#
                cost[which(is.na(cost), arr.ind = TRUE)] <- 0
                colnames(cost) <- colnames(seed)
                rownames(cost) <- colnames(query)[query_batch]
                return(cost)
            }, seed = local_seed,
            query = query,
            future.seed = TRUE)
        #---------------------------------------------------------------------#
        # rebuild slice 
        #---------------------------------------------------------------------#
        cost <- do.call("rbind", local_cost)
        total_cost[[i]] <- cost
    }
    #-------------------------------------------------------------------------#
    # Rebuild total matrix
    #-------------------------------------------------------------------------#
    total_cost <- do.call("cbind", total_cost)
    colnames(total_cost) <- colnames(seed)
    rownames(total_cost) <- colnames(query)
    total_cost <- total_cost[order(rownames(total_cost)),
        order(colnames(total_cost))]
    return(list(total_cost))
}


#' Method dispatch function for neighborhood selection
#' @param coord data.frame - coordinates of assay (barcodes, x, y)
#' @param signal matrix - matrix or sparse matrix containing assay 
#' signal for all spatial indices contained in coord
#' @param method character - which method should be use to collect 
#' neighborhood - switch matches
#' @param k int - how many nearest neighbors from KNN algorithm
#' @param depth int - graph path depth to define neighborhood 
#' 0 = self, 1 = direct neigbors, 2 = neighbors of neighbors, etc
#' @param radius - numeric - radius around center cell 
#' @return matrix of average signals for each spatial index and its 
#' neighborhood.
get_neighborhood <- function(coord,
    signal,
    method,
    k = 20,
    depth = 3,
    radius = 20) {
    niches <- switch(method,
        "knn" = knn_neighborhood(coord, k),
        "radius" = radius_neighborhood(coord, radius),
        "graph" = graph_neighborhood(coord, depth),
        "territory" = territory_neighborhood(coord))
    niches <- neighborhood_signal(niches, signal)
    return(niches)
}

#' k nearest neighbors - niche selection
#' @param coord data.frame - coordinates of spatial indices in assay
#' @param k int - number of nearest neighbors to select
#' @return list containing barcodes of nearest neighbors for each 
#' spatial index.
#' @importFrom RANN nn2
knn_neighborhood <- function(coord, k) {
    coord_dist <- RANN::nn2(coord[, c("x", "y")],
        k = k)
    coord_dist <- lapply(seq(1, nrow(coord_dist$nn.idx)),
        function(i,x){return(x[i,])}, x = coord_dist$nn.idx)
    coord_dist <- lapply(coord_dist,
        function(i, coord){return(coord$barcodes[i])}, coord = coord)
    names(coord_dist) <- coord$barcodes
    return(coord_dist)
}

#' Radius based method to select niche
#' @param coord data.frame - coordinates of spatial indices in assay
#' @param radius - numeric - radius around center cell 
#' @return list containing barcodes of nearest neighbors for each 
#' spatial index.
#' @importFrom RANN nn2
radius_neighborhood <- function(coord, radius) {
    coord_dist <- RANN::nn2(coord[, c("x", "y")],
        k = nrow(coord))
    coord_dist <- lapply(seq(1, nrow(coord_dist$nn.idx)),
        function(i,x,r){
            tmp <- x$nn.idx[i, x$nn.dists[i, ] <= r]
            return(tmp)
    }, x = coord_dist, r = radius)
    coord_dist <- lapply(coord_dist,
        function(i, coord){return(coord$barcodes[i])}, coord = coord)
    names(coord_dist) <- coord$barcodes
    return(coord_dist)
}
#' Graph depth method based method to select niche
#' @param coord data.frame - coordinates of spatial indices in assay
#' @param depth int - graph path depth to define neighborhood 
#' 0 = self, 1 = direct neigbors, 2 = neighbors of neighbors, etc
#' @return list containing barcodes of  neighbors for each 
#' spatial index.
graph_neighborhood <- function(coord, depth) {
    coord_dist <- graph_from_voronoi(coord)
    coord_dist <- graph_path_length(coord_dist)
    coord_dist <- lapply(seq(1, nrow(coord)),
        function(i, g, coord, d) {
        tmp <- g[ ,i] <= d
        return(coord$barcodes[tmp])
    }, g = coord_dist, coord = coord, d = depth)
    names(coord_dist) <- coord$barcodes
    return(coord_dist)
}

#' Territory selection for large scale niche
#' @param coord data.frame - coordinates of spatial indices in assay
#' @return list containing barcodes of spatial indices in each territory
territory_neighborhood <- function(coord) {
    # for now only last trial 
    coord <- check_territory_trial(coord, "last")
    coord_dist <- lapply(seq(1, nrow(coord)),function(i, coord) {
        bars <- coord$barcodes[coord$trial == coord$trial[i]]
        return(bars)
    }, coord = coord)
    names(coord_dist) <- coord$barcodes
    return(coord_dist)
}

#' compute average expression of local neighborhood
#' @param neighbors list of local neighbors
#' @param signal count matrix or feature matrix to average
#' @return average feature matrix. The expression of each cell 
#' is replace by the average expression of the k nearest neighbors
neighborhood_signal <- function(neighbors, signal) {
    n_signal <- lapply(neighbors, function(i, signal){
        tmp <- signal[, colnames(signal) %in% i]
        if (is.null(dim(tmp))){
            return(matrix(tmp, ncol = 1))
        } else {
            return(rowMeans(tmp))
        }
    }, signal = signal)
    n_signal <- do.call("cbind", n_signal)
    colnames(n_signal) <- names(neighbors)
    return(n_signal)
}

#' Method dispatch function for neighborhood selection - added
#' flavor specific to composition
#' @param coord data.frame - coordinates of assay (barcodes, x, y)
#' @param signal matrix - matrix or sparse matrix containing assay 
#' signal for all spatial indices contained in coord
#' @param method character - which method should be use to collect 
#' neighborhood - switch matches
#' @param k int - how many nearest neighbors from KNN algorithm
#' @param depth int - graph path depth to define neighborhood 
#' 0 = self, 1 = direct neigbors, 2 = neighbors of neighbors, etc
#' @param radius - numeric - radius around center cell 
#' @return list of all cells and cell types for each niche
niche_composition <- function(coord,
    vesalius_assay,
    method,
    k = 20,
    depth = 3,
    radius = 20) {
    niches <- switch(method,
        "knn" = knn_neighborhood(coord, k),
        "radius" = radius_neighborhood(coord, radius),
        "graph" = graph_neighborhood(coord, depth),
        "territory" = territory_neighborhood(coord))
    niches <- lapply(niches, function(n, assay) {
            composition <- get_territories(assay) %>%
                filter(barcodes %in% n)
            if (!any(grepl("Cells", colnames(composition)))) {
                stop("No cell labels have been added! See add_cells function")
            }
            composition <- composition[, tail(grep("Cells", colnames(composition)),1)]
            return(make.unique(composition))
    }, assay = vesalius_assay)
    return(niches)
}

#' Wrapper for fast niche composition 
#' @param seed list - cell type list in seed
#' @param query list - cell type list in query
#' @return matrix - jaccard index between nich composition in 
#' seed and query
compare_niches <- function(seed, query) {
    # Currently using jaccard but might change towards shannon diversity
    jacques <- compare_niche_fast(seed, query)
    colnames(jacques) <- names(seed)
    rownames(jacques) <- names(query)
    return(jacques)
}

#' splitting cost matrix into batches for batch processing
#' @param cost_mat matrix 
#' @param batch_size integer - size of each batch
#' @param verbose logical out message to console
#' @return list containing cost matrix batches 
divide_and_conquer <- function(cost_mat, batch_size, verbose) {
    #-------------------------------------------------------------------------#
    # checking case
    # case 1 : large data sets both need batching
    #-------------------------------------------------------------------------#
    if (ncol(cost_mat) > batch_size && nrow(cost_mat) > batch_size) {
        message_switch("div_hungarian", verbose)
        query_batch <- chunk(sample(seq(1, nrow(cost_mat)), nrow(cost_mat)),
                batch_size)
        if (ncol(cost_mat) >= nrow(cost_mat)) {
            seed_batch <- chunk(sample(seq(1, ncol(cost_mat)), ncol(cost_mat)),
                batch_size,
                l = length(query_batch))
            coms <- "dispatch"
        } else {
            seed_batch <- lapply(query_batch, function(query, max_idx) {
                    return(sample(seq(1, max_idx), length(query)))
                }, max_idx = ncol(cost_mat))
            coms <- "oversample"
        }
        
    } else if (nrow(cost_mat) > batch_size && nrow(cost_mat) > ncol(cost_mat)) {
        message_switch("div_hungarian", verbose)
        query_batch <- chunk(sample(seq(1, nrow(cost_mat)), nrow(cost_mat)),
                batch_size)
        seed_batch <- rep(list(seq_len(ncol(cost_mat))),
            times = length(query_batch))
        coms <- "subsample"
    } else if (ncol(cost_mat) > batch_size && ncol(cost_mat) > nrow(cost_mat)) {
        message_switch("div_hungarian", verbose)
        seed_batch <- ceiling(ncol(cost_mat) / batch_size)
        seed_batch <- lapply(seq(1, seed_batch),
            function(subs, query_len, max_idx) {
                    return(sample(seq(1, max_idx), query_len))
            },query_len = nrow(cost_mat), max_idx = ncol(cost_mat))
        query_batch <- rep(list(seq_len(nrow(cost_mat))),
            times = length(seed_batch))
        coms <- "reduce"
    } else if (nrow(cost_mat) > ncol(cost_mat)){
        message_switch("hungarian", verbose)
        query_batch <- list(seq_len(nrow(cost_mat)))
        seed_batch <- list(sample(seq_len(ncol(cost_mat)),
            nrow(cost_mat), replace = TRUE))
        coms <- "oversample"
    } else {
        message_switch("hungarian", verbose)
        query_batch <- list(seq_len(nrow(cost_mat)))
        seed_batch <- list(seq_len(ncol(cost_mat)))
        coms <- "exact"
    }
    cost_list <- mapply(function(q, s, cost) {
                return(cost[q, s])
            },
            query_batch,
            seed_batch,
            MoreArgs = list(cost_mat),
            SIMPLIFY = FALSE)   
    comment(cost_list) <- coms
    return(cost_list)
}


#' Finding optimal mapping between query and seed using a mini batch 
#' approach
#' @param idx integer - index of batch 
#' @param cost_matrix list of matrices with cost between query (rows) and 
#' seed (columns)
#' @details Uses a Hungarian solver to minimize overall cost of 
#' pair assignement.
#' @return data.frame with point in query (from) mapped to 
#' which point in seed (to) as well as cost of mapping for that pair
match_index <- function(idx, cost_matrix) {
    coms <- comment(cost_matrix)
    mapping <- lapply(idx, function(idx,
        cost_matrix, n_cost) {
        cost_matrix <- cost_matrix[[idx]]
        map <- RcppHungarian::HungarianSolver(cost_matrix)
        mapping <- as.data.frame(map$pairs)
        colnames(mapping) <- c("from", "to")
        scores <- mapply(function(i, j, cost) {
                return(cost[i, j])
            },
            mapping$from,
            mapping$to,
            MoreArgs = list(cost_matrix))

        mapping$cost <- scores
        mapping$to <- colnames(cost_matrix)[mapping$to]
        mapping$from <- rownames(cost_matrix)[mapping$from]
        return(mapping)
    },
    cost_matrix = cost_matrix)#,
    #future.seed = TRUE)
    mapping <- concat_matches(mapping, coms)
    return(mapping)
}



#' merging batch matches together
#' @param matched_indices list containing matched points in each batch
#' @param coms character - potential changes to be applied to data
#' Attached to list as comment attribute
#' @return data.frame with matched points
concat_matches <- function(matched_indices, coms) {
    matched_indices <- lapply(matched_indices, function(x) {
        return(x[, c("from", "to", "cost")])
    })
    if (coms == "reduce") {
        loc <- sapply(matched_indices, function(x){
                return(sum(x$cost))
            })
        loc <- which(loc == min(loc))
        matched_indices <- matched_indices[[loc]]
    } else {
        matched_indices <- do.call("rbind", matched_indices)
    }
    return(matched_indices)
}

#' Extract score from best matches 
#' @param matched_indices matrix - best matches between seed and query
#' @param cost list - cost list essentially scores 
#' @param use_cost character - which scores should be added to assay
#' @return data frame containing matched indices and associated cost and scores 
score_matches <- function(matched_indices,
    cost,
    use_cost) {
    for (i in seq_along(use_cost)) {
        tmp_cost <- cost[[use_cost[i]]]
        inter_score <- sapply(seq(1, nrow(matched_indices)),
            function(idx,matched, tmp_cost) {
                return(tmp_cost[matched$from[idx],matched$to[idx]])
            },matched = matched_indices, tmp_cost = tmp_cost)
        matched_indices$inter_score <- inter_score
        colnames(matched_indices) <- gsub("inter_score", use_cost[i], colnames(matched_indices))
    }
    return(matched_indices)
}

#' assign coordinates to matched indices
#' @param matched_index data.frame containing matching pairs of 
#' coordinates
#' @param seed data.frame containing seed coordinates 
#' @param query data.frame containing quert cooridates
#' @param verbose logical - should progress message be outputed to the 
#' console
#' @return adjusted query coordinate data.frame where each point
#' receives the coordinates of its best matche in the seed. 
align_index <- function(matched_index,
    seed,
    query,
    threshold = 0.3,
    verbose = TRUE) {
    query <- query[match(matched_index$from, query$barcodes), ]
    seed <- seed[match(matched_index$to, seed$barcodes), ]
    query$x <- seed$x
    query$y <- seed$y
    query$x <- jitter(query$x, factor = 0.1)
    query$y <- jitter(query$y, factor = 0.1)
    return(query)
}



integrate_map <- function(aligned_graph,
    query,
    seed,
    merge = FALSE,
    verbose = TRUE) {
    prob <- list(aligned_graph$prob)
    names(prob) <- "mapping_scores"
    aligned_graph <- aligned_graph$aligned
    local_counts <- get_counts(query, type = "raw")
    local_counts <- local_counts[,
        colnames(local_counts) %in% aligned_graph$barcodes]
    query <- build_vesalius_assay(
        coordinates = aligned_graph[, c("barcodes", "x", "y")],
        counts = local_counts,
        assay = "integrated",
        verbose = FALSE)
    query@meta <- c(query@meta, prob)
    if (!merge) {
        return(query)
    } else {
        return(query)
    }
    
}
#' @importFrom kohonen supersom somgrid map
intergrate_counts <- function(seed,
    query,
    signal,
    dimensions = seq(1:30),
    use_norm = "raw",
    scale = FALSE,
    grid = "auto",
    verbose = TRUE) {
    signal <- get_signal(seed,
        query,
        signal,
        dimensions,
        use_norm,
        scale,
        verbose = FALSE)
    if (grid == "auto") {
        grid <- ceiling(sqrt(ceiling(5 * sqrt(ncol(signal$seed)))))
    } else {
        # should add a type check 
        grid <- ceiling(sqrt(grid))
    }
    grid <- kohonen::somgrid(xdim = grid,
        ydim = grid,
        topo = "rectangular",
        neighbourhood.fct = "gaussian")
    som_seed <- kohonen::supersom(t(as.matrix(signal$seed)),
        grid = grid,
        mode = "pbatch",
        dist.fcts = "euclidean",
        cores = future::nbrOfWorkers())
    som_query <- kohonen::map(som_seed, t(as.matrix(signal$query)))
    som_clusters <- match_som_clusters(som_seed, som_query, signal$seed, signal$query)
    seed_idx <- match(som_clusters$nn, colnames(signal$seed)) - 1
    query_idx <- match(som_clusters$barcodes, colnames(signal$query)) - 1
    query <- rescale_to_seed(seed_idx, query_idx, as.matrix(signal$seed), as.matrix(signal$query))
    browser()
}


match_som_clusters <- function(som_seed, som_query, seed, query) {
    som_seed <- data.frame("barcodes" = colnames(seed),
        "som" = som_seed$unit.classif) 
    som_seed <- split(som_seed, som_seed$som)
    som_query <- data.frame("barcodes" = colnames(query),
        "som" = som_query$unit.classif)
    som_query <- split(som_query, som_query$som)
    for (i in seq_along(som_query)) {
        if (names(som_query)[i] %in% names(som_seed)){
            tmp_seed <- seed[ ,colnames(seed) %in% som_seed[[names(som_query)[i]]]$barcodes]
            tmp_seed <- matrix(tmp_seed, ncol = nrow(som_seed[[names(som_query)[i]]]))
            tmp_query <- query[ ,colnames(query) %in% som_query[[i]]$barcodes]
            tmp_query <- matrix(tmp_query, ncol = nrow(som_query[[i]]))
            nn <- RANN::nn2(data = t(as.matrix(tmp_seed)), query = t(as.matrix(tmp_query)), k = 1)
            som_query[[i]]$nn <- som_seed[[names(som_query)[i]]]$barcodes[nn$nn.idx]
        } else {
            next
        }
        
    }
    som_query <- som_query[!sapply(som_query, ncol) != 3]
    som_query <- do.call("rbind", som_query)
    return(som_query)
}





#-----------------------------------------------------------------------------#
############################### Joint Measure #################################
#-----------------------------------------------------------------------------#

#' Jointly measured spatial omic assays territories
#' @param mod1 vesalius_assay object containing first modality
#' @param mod2 vesalius_assay objecty containing second modality
#' @param dimensions numeric vector describing latent space dimensions 
#' to use during intergration
#' @param method character - integration method. interlace - mean - concat 
#' are available options
#' @param norm_method character - which count values should be use 
#' for integration when using concat method
#' @param dim_reduction characater - which dim reduction methods should be 
#' used for concat integration (PCA,PCA_L,UMAP,LSI,LSI_UMAP,NMF)
#' @param verbose logical - should progress message be outputed to the 
#' console.
#' @return vesalius object containing new image embeddings
#' @export
joint_territories <- function(mod1,
    mod2,
    dimensions = seq(1, 30),
    embedding = "last",
    method = "interlace",
    norm_method = "log_norm",
    dim_reduction = "PCA",
    signal = "variable_features",
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # check place holder 
    #-------------------------------------------------------------------------#

    #-------------------------------------------------------------------------#
    # Get embeddings - need make some changes here
    # for now we assume we get the active 
    #-------------------------------------------------------------------------#
    mod1_embed <- check_embedding_selection(mod1, embedding, dimensions)
    mod2_embed <- check_embedding_selection(mod2, embedding, dimensions)
    #-------------------------------------------------------------------------#
    # Get counts from each vesalius object
    #-------------------------------------------------------------------------#
    mod1_counts <- get_counts(mod1, type = "all")
    names(mod1_counts) <- paste0("mod1_", names(mod1_counts))
    mod2_counts <- get_counts(mod2, type = "all")
    names(mod2_counts) <- paste0("mod2_", names(mod2_counts))
    #-------------------------------------------------------------------------#
    # method switch - which method is best 
    #-------------------------------------------------------------------------#
    integrated_embeds <- switch(EXPR = method,
        "interlace" = interlace_embeds(mod1_embed, mod2_embed, dimensions),
        "mean" = average_embed(mod1_embed, mod2_embed, dimensions),
        "concat" = concat_embed(mod1,
            mod2,
            dimensions,
            norm_method,
            dim_reduction,
            signal = signal))
    integrated_embeds <- list(integrated_embeds)
    names(integrated_embeds) <- method
    #-------------------------------------------------------------------------#
    # get tile overlap - it is possible that there is not a perfect overlap
    # so we will filter
    #-------------------------------------------------------------------------#
    tiles <- mod1@tiles
    tiles <- tiles[tiles$barcodes %in% rownames(integrated_embeds[[1]]), ]
    integrated <- new("vesalius_assay",
        assay = "integrated",
        embeddings = integrated_embeds,
        active = integrated_embeds[[1]],
        tiles = tiles)
    integrated_counts <- c(mod1_counts, mod2_counts)
    integrated <- update_vesalius_assay(vesalius_assay = integrated,
      data = integrated_counts,
      slot = "counts",
      append = FALSE)
    #--------------------------------------------------------------------------#
    # we can update the comment on the count slot list 
    # this comment will indicate which count matrix is set as default
    #--------------------------------------------------------------------------#
    integrated <- add_active_count_tag(integrated, norm = "joint")
    #--------------------------------------------------------------------------#
    # Finally we update the vesalius commit log
    #--------------------------------------------------------------------------#
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(joint_territories))
    integrated <- commit_log(integrated,
      commit,
      "integrated")
    simple_bar(verbose)
    return(integrated)

}

#' interlace image stack between seed and query
#' @param seed matrix - seed embedding image stack
#' @param query matrix - query embedding image stack
#' @param dimensions int vector describing which embeddings
#' should be selected
#' @details Takes selected embedding from seed and query and 
#' creates and interlaced embedding matrix starting with the 
#' first seed embedding
#' @return embedding matrix containing seed embeddings + query 
#' embeddings.
interlace_embeds <- function(seed, query, dimensions) {
    locs <- intersect(rownames(seed), rownames(query))
    seed <- seed[rownames(seed) %in% locs,
        dimensions]
    query <- query[rownames(query) %in% locs,
        dimensions]
    seed <- seed[order(rownames(seed)), ]
    query <- query[order(rownames(query)), ]

    interlaced_embed <- matrix(0,
        ncol = ncol(seed) + ncol(query),
        nrow = nrow(seed))
    rownames(interlaced_embed) <- rownames(seed)
    dimensions <- rep(dimensions, each = 2)
    for (i in seq(1, ncol(interlaced_embed), by = 2)) {
        interlaced_embed[, i] <- seed[, dimensions[i]]
        interlaced_embed[, i + 1] <- query[, dimensions[i + 1]]
    }
    return(interlaced_embed)
}

#' average image stack between seed and query
#' @param seed matrix - seed embedding image stack
#' @param query matrix - query embedding image stack
#' @param dimensions int vector describing which embeddings
#' should be selected
#' @details Takes select embedding from seed and query and 
#' creates and avarage the grey scale pixel values for each spatial
#' location
#' @return embedding matrix containing average pixel value for both seed
#' query
average_embed <- function(seed, query, dimensions) {
    locs <- intersect(rownames(seed), rownames(query))
    seed <- seed[rownames(seed) %in% locs,
        dimensions]
    query <- query[rownames(query) %in% locs,
        dimensions]
    seed <- seed[order(rownames(seed)), ]
    query <- query[order(rownames(query)), ]
    averaged_embed <- matrix(0,
        ncol = length(dimensions),
        nrow = nrow(seed))
    rownames(averaged_embed) <- rownames(seed)
    for (i in seq(1, ncol(averaged_embed))) {
        averaged_embed[, i] <- apply(cbind(seed[, dimensions[i]],
            query[, dimensions[i]]),
            MARGIN = 1,
            mean)
    }
    return(averaged_embed)
}


#' create new embedding from jointly measure spatial omics 
#' @param seed vesalius_assay object of the first modality 
#' @param query vesalius_assay object of the second modality
#' @param dimensions int - number of gray scale images to create
#' @param norm_method string describing which normalisation 
#' method to use. One of the following "log_norm", "SCT", "TFIDF", "raw"
#' @param dim_reduction string describing which dimensionality
#' reduction method should be used. One of the following:
#' "PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP"
#' @details Create new latent space using feature from both modalities
#' Creates a new feature matrix than in nromalized and converted to latent
#' space image stack.
#' @return embedding matrix used as grey scale image stack
concat_embed <- function(seed,
    query,
    dimensions,
    norm_method,
    dim_reduction,
    signal) {
    #-------------------------------------------------------------------------#
    # first we get signal - i.e. these features that should be used
    #-------------------------------------------------------------------------#
    seed_features <- check_signal(seed, signal, type = norm_method)
    query_features <- check_signal(query, signal, type = norm_method)
    #-------------------------------------------------------------------------#
    # next we extract counts and make sure that they are int he right order
    # before rbidning the whole thing 
    #-------------------------------------------------------------------------#
    seed_counts <- get_counts(seed, type = "raw")
    query_counts <- get_counts(query, type = "raw")
    locs <- intersect(colnames(seed_counts), colnames(query_counts))
    seed_counts <- seed_counts[rownames(seed_counts) %in% seed_features,
        colnames(seed_counts) %in% locs]
    query_counts <- query_counts[rownames(query_counts) %in% query_features,
        colnames(query_counts) %in% locs]
    seed_counts <- seed_counts[, order(colnames(seed_counts))]
    query_counts <- query_counts[, order(colnames(query_counts))]
    integrated_counts <- rbind(seed_counts, query_counts)
    #-------------------------------------------------------------------------#
    # Now we can process the counts and create a new embedding
    #-------------------------------------------------------------------------#
    integrated_counts <- process_counts(integrated_counts,
        assay = "integrated",
        method = norm_method,
        use_count = "raw",
        nfeatures = sum(c(length(seed_features), length(query_features))))
    integrated_embeds <- embed_latent_space(integrated_counts$SO,
        assay = "integrated",
        dim_reduction = dim_reduction,
        dimensions = max(dimensions),
        remove_lsi_1 = FALSE,
        verbose = FALSE)
    
    return(integrated_embeds[[1]])
}

