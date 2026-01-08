###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################### MAPPING ASSAYS ################################
#-----------------------------------------------------------------------------#

#' Aling and integrate spatial assay from the same modality using super pixels
#' @param seed_assay vesalius_assay object - data to be mapped to
#' @param query_assay vesalius_assay objecy - data to map
#' @param use_cost character string defining how should total cost be computer
#' Available: feature, niche, territory, composition (See details for combinations
#' and custom matrices)
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

#' @param custom_cost matrix - matrix of size n (query cells) by p (seed cells)
#' containing custom cost matrix. Used instead of vesalius cost matrix
#' @param method character - correlation method to use
#' @param epochs numeric - number of epochs
#' @param allow_duplicates logical - allow duplicate mappings
#' @param filter_cells logical - filter cells
#' @param seed_territory_labels character - territory labels for seed
#' @param query_territory_labels character - territory labels for query
#' @param seed_meta_labels character - meta labels for seed
#' @param query_meta_labels character - meta labels for query
#' @param jitter numeric - jitter amount
#' @param digits numeric - number of digits
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
    signal = "variable_features",
    use_cost = c("feature","niche"),
    method = "pearson",
    neighborhood = "knn",
    k = 20,
    radius = 0.05,
    depth = 1,
    dimensions = seq(1, 30),
    batch_size = 10000,
    epochs = 1,
    allow_duplicates = TRUE,
    threshold = 0.3,
    filter_cells = FALSE,
    use_norm = "raw",
    scale = FALSE,
    custom_cost = NULL,
    seed_territory_labels = "Territory",
    query_territory_labels = "Territory",
    seed_meta_labels = NULL,
    query_meta_labels = NULL,
    jitter = 0,
    digits = 5,
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
    signal_out <- get_signal(seed_assay = seed_assay,
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
        query_signal = signal_out$query,
        query_assay = query_assay,
        cost = custom_cost,
        seed_signal = signal_out$seed,
        seed_assay = seed_assay,
        method = method,
        neighborhood = neighborhood,
        k = k,
        radius = radius,
        depth = depth,
        batch_size = batch_size,
        epochs = epochs,
        use_cost = use_cost,
        threshold = threshold,
        filter_cells = filter_cells,
        seed_territory_labels = seed_territory_labels,
        query_territory_labels = query_territory_labels,
        seed_meta_labels = seed_meta_labels,
        query_meta_labels = query_meta_labels,
        digits = digits,
        verbose = verbose)
    mapped <- filter_maps(mapped,
        allow_duplicates = allow_duplicates,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # Rebuild a base obejct - we will not integrate here. 
    #-------------------------------------------------------------------------#
    vesalius_assay <- build_mapped_assay(mapped,
        seed_assay = seed_assay,
        query_assay = query_assay,
        meta_labels = query_meta_labels,
        jitter = jitter)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(map_assays))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay = get_assay_names(vesalius_assay))
    simple_bar(verbose)
    return(vesalius_assay) 
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
#' @param query_assay vesalius_assay object
#' @param cost matrix - matrix of size n (query cells) by p (seed cells)
#' containing custom cost matrix.
#' @param seed_signal processed seed signal from seed assay
#' @param seed_assay vesalius_assay object
#' @param k int size of niche (knn)
#' @param radius 0.05 proportion of max distance to use as radius for 
#' neighborhood
#' @param depth graph path depth to condsider for neighborhood. 
#' @param batch_size int number of points in each query batch
#' @param use_cost character string defining how should total cost be computer
#' Available: feature, niche, territory, composition (See details for combinations
#' @param method character - correlation method
#' @param neighborhood character - neighborhood method
#' @param epochs numeric - number of epochs
#' @param filter_cells logical - filter cells
#' @param seed_territory_labels character - seed territory labels
#' @param query_territory_labels character - query territory labels
#' @param seed_meta_labels character - seed meta labels
#' @param query_meta_labels character - query meta labels
#' @param digits numeric - number of digits
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
    method = "pearson",
    neighborhood = "knn",
    k = 20,
    radius = 0.05,
    depth = 1,
    batch_size = 10000,
    epochs = 1,
    use_cost = c("feature", "niche"),
    threshold = 0.5,
    filter_cells = FALSE,
    seed_territory_labels = "Territory",
    query_territory_labels = "Territory",
    seed_meta_labels = NULL,
    query_meta_labels = NULL,
    digits = 4,
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
    if ("feature" %in% use_cost){
        
        message_switch("feature_cost", verbose, assay = assay)
        cost <- c(cost, signal_similarity(seed_signal,
            query_signal,
            method = method,
            digits = digits))
        names(cost)[length(cost)] <-  "feature"

    }
     
    #--------------------------------------------------------------------------#
    # Correlation between the cellular niche centered around the cell
    #--------------------------------------------------------------------------#
    if ("niche" %in% use_cost) {
        message_switch("get_neigh", verbose, assay = assay)
        seed_signal_niche <- get_neighborhood_signal(seed,
            seed_signal,
            neighborhood,
            k,
            depth,
            radius)
        query_signal_niche <- get_neighborhood_signal(query,
            query_signal,
            neighborhood,
            k,
            depth,
            radius)
        message_switch("neighbor_cost", verbose, assay = assay)
        cost <- c(cost, signal_similarity(seed_signal_niche,
            query_signal_niche,
            method = method,
            digits = digits))
        names(cost)[length(cost)] <-  "niche"
    }
    #--------------------------------------------------------------------------#
    # Correlation between territories centered around the cell
    # Note: this can be sped up using much smaller matrices and then dispatchin
    # score after
    # This can get a bit weird if there are cells that have been added and no
    # territory is specified. Not sure how I can make this better without 
    # breaking everything else. Essentially, if cells are added after, ves
    # will consider the Cells as territories. Could be find of interesting to 
    # to do though - meta cells? 
    #--------------------------------------------------------------------------#
    if ("territory" %in% use_cost) {
        message_switch("territory_cost", verbose, assay = assay)

        #----------------------------------------------------------------------#
        # Since this is a block of cells - needs to be filter just is case
        # there are drops outs from the custom cost
        #----------------------------------------------------------------------#
        seed_coord  <- check_territory_trial(seed_assay, seed_territory_labels)
        seed_coord <- seed_coord[seed_coord$barcodes %in% seed$barcodes,]
        query_coord  <- check_territory_trial(query_assay, query_territory_labels)
        query_coord <- query_coord[query_coord$barcodes %in% query$barcodes,]
        
        seed_signal_niche <- get_neighborhood_signal(seed_coord,
            seed_signal,
            "territory")
        query_signal_niche <- get_neighborhood_signal(query_coord,
            query_signal,
            "territory")
        cost <- c(cost, signal_similarity(seed_signal_niche,
            query_signal_niche,
            method = method,
            digits = digits))
        names(cost)[length(cost)] <-  "territory"
    }
    #--------------------------------------------------------------------------#
    # Computing nich compisition
    #--------------------------------------------------------------------------#
    if ("composition" %in% use_cost) {
        message_switch("composition_cost", verbose, assay = assay)
        seed_niche <- niche_composition(seed,
            seed_assay,
            method = neighborhood,
            cell_label = seed_meta_labels,
            k = k,
            depth = depth,
            radius = radius)
        query_niche <- niche_composition(query,
            query_assay,
            method = neighborhood,
            cell_label = query_meta_labels,
            k = k,
            depth = depth,
            radius = radius)

        cost <- c(cost, signal_similarity(seed_niche,
            query_niche,
            method = "jaccard",
            digits = digits))
        names(cost)[length(cost)] <-  "composition"
    }
    #--------------------------------------------------------------------------#
    # cell type label comparison => if same label =1 / if differenct label = 0
    #--------------------------------------------------------------------------#
    if ("cell_type" %in% use_cost) {
        message_switch("cell_cost", verbose, assay = assay)
        seed_labels <- check_cell_labels(seed_assay, seed_meta_labels)
        query_labels <- check_cell_labels(query_assay, query_meta_labels)
        cost <- c(cost, cell_type_match(seed_labels, query_labels))
        names(cost)[length(cost)] <- "cell_type"
    }
    #--------------------------------------------------------------------------#
    # filtering and pairwise addition of cost matrices and
    #--------------------------------------------------------------------------#
    cost <- filter_cost(cost, threshold, filter_cells, verbose)
    cost <- c(cost, concat_cost(cost, use_cost))
    names(cost)[length(cost)] <- "total_cost"
    #--------------------------------------------------------------------------#
    # devide cost matrix
    matched_indices <- optimize_matching(cost$total_cost,
        batch_size,
        epochs)
    scores <- score_matches(matched_indices$matched,
        cost,
        use_cost = use_cost)
    return(list("prob" = scores,
        "cost" = cost,
        "cost_by_epoch" = matched_indices$cost_by_epoch))
}

#' concat cost - pairwise sum of score complement
#' @param cost list - named list contained score matrices
#' @param use_cost character - which cost matrices to use
#' @param complement logical - whether to use complement of scores
#' @return list with cost matrix 
concat_cost <- function(cost, use_cost, complement = TRUE) { 
    if (length(use_cost) == 0) {
        stop("Please specify at least one score matrix to use")
    }
    cost <- cost[use_cost]
    if (length(cost) == 0) {
        stop(paste(paste(use_cost, collapse = " "), ": not available in score matrix list"))
    } else if (length(cost) == 1) {
        if (complement) {
            return(list(1 - cost[[1]]))
        } else {
            return(list(cost[[1]]))
        }
    } else {
        if (complement) {
            buffer <- 1 - cost[[1]]
            for (i in seq(2, length(cost))){
                buffer <- buffer + (1 - cost[[i]])
            }
            return(list(buffer))
        } else {
            buffer <- cost[[1]]
            for (i in seq(2, length(cost))){
                buffer <- buffer + (cost[[i]])
            }
        return(list(buffer))
        }
        
    }

}

#' compute the similarity between seed and query signals
#' @param seed seed signal
#' @param query query signal
#' @param method character - correlation method to use
#' @param digits numeric - number of digits for rounding
#' @details Chunking cost and signal into smaller chunks to run the 
#' correlation score in paralell. There is room for improvement here.
#' First we could dispatch the longer list to future_lapply
#' but cannot know which one it is and we need to know that so we can 
#' subset the cost. 
#' Also the functions calls feature_cost which is a R wrapper for a 
#' c++ function (cost.cpp).
#' @return matrix with query as rows and seed as colmuns
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
signal_similarity <- function(seed, query, method = "pearson", digits = 4) {
    batch_size_seed <- ceiling(ncol(seed) / future::nbrOfWorkers())
    batch_size_query <- ceiling(ncol(query) / future::nbrOfWorkers())
    #-------------------------------------------------------------------------#
    # First we chunk seed and query 
    # running a parallel lapply internally
    #-------------------------------------------------------------------------#
    seed_batch <- chunk(seq(1, ncol(seed)), batch_size_seed)
    query_batch <- chunk(seq(1, ncol(query)), batch_size_query)
    seed_bar <- colnames(seed)
    query_bar <- colnames(query)
    seed <- listify(seed, seed_batch)
    query <- listify(query, query_batch)
    total_cost <- vector("list", length(seed_batch))
    #-------------------------------------------------------------------------#
    # Loop over seed batch - idealy we would use the loop over the 
    # smallest batch 
    #-------------------------------------------------------------------------#
    for (i in seq_along(seed)) {
        #---------------------------------------------------------------------#
        # Splitting into sub lists
        #---------------------------------------------------------------------#
        local_seed <- as.matrix(seed[[i]])
        #---------------------------------------------------------------------#
        # computing score in batches
        #---------------------------------------------------------------------#
        local_cost <- future_lapply(query,
            function(local_query, local_seed) {
                local_query <- as.matrix(local_query)
                cost <- switch(EXPR = method,
                    "jaccard" = jaccard_cost(local_seed, local_query),
                    "pearson" = pearson_approx(local_seed, local_query),
                    "pearson_fast" = pearson_fast(local_seed, local_query),
                    "pearson_exact" = pearson_exact(local_seed, local_query))
                cost[which(is.na(cost), arr.ind = TRUE)] <- 0
                colnames(cost) <- colnames(local_seed)
                rownames(cost) <- colnames(local_query)
                cost <- signif(cost, digits = digits)
                return(cost)
            }, local_seed = local_seed, future.seed = TRUE)
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
    colnames(total_cost) <- seed_bar
    rownames(total_cost) <- query_bar
    total_cost <- total_cost[order(rownames(total_cost)),
        order(colnames(total_cost))]
    #-------------------------------------------------------------------------#
    # normalize if distance used
    # we want the max distance the be cosidered as the heighest cost
    #-------------------------------------------------------------------------#
    if (method == "distance") {
        total_cost <- total_cost / max(total_cost) 
        total_cost <- 1 - total_cost
    }
    
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
get_neighborhood_signal <- function(coord,
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
    # avoid int overflow 
    # keeping thsi for now - but not super efficient.
    # Might need to find a better wat of doing this
    if (nrow(coord) > 1000) {
        k <- 1000
    } else {
        k <- nrow(coord)
    }
    coord_dist <- RANN::nn2(coord[, c("x", "y")],
        k = k)
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
#' @param vesalius_assay vesalius_assay object
#' @param method character - which method should be use to collect 
#' neighborhood - switch matches
#' @param k int - how many nearest neighbors from KNN algorithm
#' @param depth int - graph path depth to define neighborhood 
#' 0 = self, 1 = direct neigbors, 2 = neighbors of neighbors, etc
#' @param radius - numeric - radius around center cell
#' @param cell_label character - cell label column name
#' @return list of all cells and cell types for each niche
niche_composition <- function(coord,
    vesalius_assay,
    method,
    cell_label = NULL,
    k = 20,
    depth = 3,
    radius = 20) {
    niches <- switch(method,
        "knn" = knn_neighborhood(coord, k),
        "radius" = radius_neighborhood(coord, radius),
        "graph" = graph_neighborhood(coord, depth),
        "territory" = territory_neighborhood(coord))
    cell_labels <- check_cell_labels(vesalius_assay, cell_label)
    niches <- lapply(niches, function(n, cell_labs) {
            composition <- cell_labs$trial[cell_labs$barcodes %in% n]
            composition <- make.unique(composition, sep = "_")
            return(composition)
    }, cell_labs = cell_labels)
    niches <- make_composition_matrix(niches)
    return(niches)
}

#' Compute score based on cell type labels
#' @param seed_labels cell type labels in seed (reference) data
#' @param query_labels cell type labels in query data 
#' @details Return 1 if same label and 0 if different.
#' @return Score matrix based on cell label similarity
cell_type_match <- function(seed_labels, query_labels) {
    cell_labels <- future_lapply(seq(1, nrow(query_labels)),
        match_cells,
        query_labels = query_labels,
        seed_labels = seed_labels,
        future.seed = TRUE)
    cell_labels <- do.call("rbind", cell_labels)
    colnames(cell_labels) <- seed_labels$barcodes
    rownames(cell_labels) <- query_labels$barcodes
    cell_labels <- cell_labels[order(query_labels$barcodes),
        order(seed_labels$barcodes)]
    return(list(cell_labels))
}

#' Internal function returning cell type label match as numiric
#' @param idx numeric - used to iterate down the labels
#' @param query_labels cell type labels for query
#' @param seed_labels cell type labels for seed
#' @return binary numeric vector 
match_cells <- function(idx, query_labels, seed_labels) {
    return(as.numeric(query_labels[idx, "trial"] == seed_labels$trial))
}


filter_cost <- function(costs,
    threshold = 0,
    filter_cells = TRUE,
    verbose) {
    message_switch("filter_cost", verbose)
    keep <- rep(TRUE, nrow(costs[[1]]))
    for (cost in seq_along(costs)){
        tmp <- costs[[cost]]
        name <- names(costs)[cost]
        if (name == "cell_type" && filter_cells){
            keep <- keep & apply(tmp, 1, sum) > 0
        } else {
            keep <- keep & apply(tmp, 1, max) > threshold
        }
    }
    if (sum(keep) == 0){
        stop("No cells retained at current filter threshold!")
    }

    costs <- lapply(costs, function(mat, keep) {
        return(mat[keep, ])
    },keep = keep)
    return(costs)
}

#' optimize matching scores through batching
#' @param cost_matrix matrix containing mapping cost for each cell
#' @param batch_size int - number of cells to be assigned to each batch
#' @param epochs number of epochs to run the optimization
#' @param verbose logical - output progress messages 
#' @return list with best matching cell pairs (data.frame) and 
#' total cost at each epoch 
optimize_matching <- function(cost_matrix,
    batch_size = 10000,
    epochs = 1,
    verbose = TRUE) {

    matched <- initialize_matches(cost_matrix)
    current_epoch <- 1
    current_cost <- data.frame("cost" = rep(0, epochs),
        "epoch" =  rep(0, epochs))
    while (current_epoch <= epochs) {
        message_switch("mapping", verbose , epoch = current_epoch)
        batch <- dispatch_batch(cost_matrix, batch_size)
        mapped <- future_lapply(batch, map_index, future.seed = TRUE)
        matched <- update_matches(matched, mapped, current_epoch)
        current_cost$epoch[current_epoch] <- current_epoch
        current_cost$cost[current_epoch] <- mean(matched$cost)
        current_epoch <- current_epoch + 1
    }
    matched <- check_for_unmatched(matched)
    return(list("matched" = matched, "cost_by_epoch" = current_cost))

}


#' initialize matched data frame
#' @param cost_matrix matrix containing cost of each cell pair
#' @return data.frame templated which will contain the best matching pairs 
initialize_matches <- function(cost_matrix) {
    if (ncol(cost_matrix) > nrow(cost_matrix)) {
        n_row <- ncol(cost_matrix)
        matches <- data.frame("from" = rep(NA, n_row),
            "to" = sort(colnames(cost_matrix)),
            "epoch" = rep(1, n_row),
            "cost" = rep((max(cost_matrix) + 1), n_row),
            "init" = rep("to", n_row))
    } else {
        n_row <- nrow(cost_matrix)
        matches <- data.frame("from" = sort(rownames(cost_matrix)),
            "to" = rep(NA, n_row),
            "epoch" = rep(1, n_row),
            "cost" = rep((max(cost_matrix) +1), n_row),
            "init" = rep("from", n_row))
    }
    return(matches)
}

#' Dispatch cells into batches
#' @param cost_matrix matrix containing cost for each cell pair
#' @param batch_size int size of batch
#' @details Create cell batches that will dynamically adapt to the size of
#' the data set with respect to batch size. Smalled data sets, cells
#' will be sampled to match the size of the larger data set. This
#' allows for multiple to multiple matching. All cells will be selected
#' at least once.
#' @return Nested list. Each element of the list will contain 
#' a batched cost matrix and the mapping pairs 
dispatch_batch <- function(cost_matrix, batch_size = 5000) {
    seed_barcodes <- seed_current <- colnames(cost_matrix)
    query_barcodes <- query_current <- rownames(cost_matrix)
    #-------------------------------------------------------------------------#
    # effective batch size
    #-------------------------------------------------------------------------#
    batch_seed <- min(c(batch_size, length(seed_barcodes)))
    batch_query <- min(c(batch_size, length(query_barcodes)))
    batch_size <- min(c(max(c(batch_seed, batch_query)), batch_size))
    #-------------------------------------------------------------------------#
    # creating list for batch
    # we can compute the what is the most amount of batch to cover all barcodes
    # at least once
    #-------------------------------------------------------------------------#
    seed_length <- ceiling(length(seed_barcodes) / batch_size)
    query_length <- ceiling(length(query_barcodes) / batch_size)
    batch_length <- max(c(seed_length, query_length))
    batch <- vector("list", batch_length)
    for (i in seq_along(batch)){
        len_seed <- length(seed_current)
        len_query <- length(query_current)
        if (len_seed < batch_size){
            pad <- batch_size - len_seed
            seed <- c(sample(seed_current, size = len_seed),
                sample(seed_barcodes, size = pad, replace = TRUE))
        } else {
            seed <- sample(seed_current, size = batch_size)
            seed_current <- seed_current[!seed_current %in% seed]
        }

        if (len_query < batch_size){
            pad <- batch_size - len_query
            query <- c(sample(query_current, size = len_query),
                sample(query_barcodes, size = pad, replace = TRUE))
        } else {
            query <- sample(query_current, size = batch_size)
            query_current <- query_current[!query_current %in% query]
        }

        batch[[i]] <- list("cost" = cost_matrix[query, seed],
            "match" = data.frame("from" = query, "to" = seed))

    }
    
    return(batch)

}


#' LAPVJ solver 
#' @param batch cost matrix batch to be solved
#' @return data frame with best matches for this batch
map_index <- function(batch) {
    mapped <- TreeDist::LAPJV(batch$cost)$matching
    match <- batch$match
    match$to <- match$to[mapped]
    cost <- mapply(function(i, j, cost) {
               return(cost[i, j])
        },
        match$from,
        match$to,
        MoreArgs = list(batch$cost))
    match$cost <- cost
    return(match)
}

#' Update best matched cell pairs with new mapping costs
#' @param matched data frame containing mapping pairs template
#' @param mapped best mapping pairs for each batch
#' @param epoch int - which epoch was the optimal match found
#' @return updated matched data frame
update_matches <- function(matched, mapped, epoch) {
    if (all(matched$init == "from")) {
        init_col <- "from"
        map_col <- "to"
    } else {
        init_col <- "to"
        map_col <- "from"
    }
    for (i in seq_along(mapped)) {
        loc <- match(mapped[[i]][, init_col], matched[, init_col])
        loc <- loc[!is.na(loc)]
        cost <- mapped[[i]][, "cost"] < matched[loc, "cost"]
        matched[loc[cost], map_col] <- mapped[[i]][cost, map_col]
        matched[loc[cost], "cost"] <- mapped[[i]]$cost[cost]
        matched[loc[cost], "epoch"] <- epoch
    }
    
    return(matched)
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



#' filter mapped cells 
#'@param mapped data frame containing mapped cells
#'@param allow_duplicates logical - if duplicated matches are to be reatained
#'@param verbose logical - output progress messages
#'@return filtered mappped data frame
filter_maps <- function(mapped, allow_duplicates, verbose) {
    message_switch("post_map_filter",verbose)
    map_score <- mapped$prob
    map_score <- map_score[order(map_score$cost), ]
    if (!allow_duplicates) {
        duplicates <- duplicated(map_score$from) | duplicated(map_score$to)
        map_score <- map_score[!duplicates, ]
    }
    mapped$prob <- map_score
    return(mapped)
}



#' assign coordinates to matched indices
#' @param matched_index data.frame containing matching pairs of
#' coordinates
#' @param coord data.frame containing coordinate data
#' @param jitter numeric - amount of jitter to add to duplicated coordinates
#' @return adjusted coordinate data.frame where each point
#' receives the coordinates of its best match. 
align_index <- function(matched_index,
    coord,
    jitter = 0) {
    coord$barcodes <- matched_index$from
    if (jitter){
        locs <- duplicated(paste0(coord$x,"_", coord$y))
        coord$x[locs] <- jitter(coord$x[locs], amount = jitter)
        coord$y[locs] <- jitter(coord$y[locs], amount = jitter)
    }
    coord$barcodes <- make.unique(coord$barcodes, sep = "-")
    return(coord)
}




