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
#' "knn", "radius", "depth" (See details)
#' @param k int ]2, n_points] number of neareset neighbors to be considered for
#' neighborhodd computation.
#' @param radius numeric ]0,1[ proportion of max distance between points 
#' to consider for the neighborhood
#' @param depth int [1, NA] graph depth from cell to consider from neighborhood
#' (See details)
#' @param dimensions Int vector containing latent space dimensions to use
#' @param mapping character string (div - exact) - mapping strategy.
#' Divide and conquer (approximate) or Exact mapping.
#' @param batch_size number of points per batch in query during assignment
#' problem solving
#' @param signal character (variable_features, all_features, embeddings, custom)
#' - What should  be used as cell signal to generate the cost matrix.
#' Seed details 
#' @param use_norm character - which count data to use
#' @param custom_cost matrix - matrix of size n (query cells) by p (seed cells)
#' containing custom cost matrix. Used instead of vesalius cost matrix
#' @param overwrite logical - if custom_cost is not NULL, should this matrix
#' be added to the vesalius matrix or should it overwrite it completely.
#' @param verbose logical - should I be a noisy boy?
#' @details The goal is to assign the best matching point between a seed set and
#' a query set.
#' 
#' To do so, \code{integrate_horizontally} will first extract a
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
#' The second step is to create a cost matrix. The cost matrix computes the
#' distance (using the signal) betwwen each point in the seed and query. It
#' also includes the distance in signal between each points neighborhoods.
#' The neighborhood singal is computed by averaging the signal across k
#' nearest neighbors in space if using Knn. If using radius, the neighborhood
#' will be computed using all cells within a definined radius. 
#' If using depth, coordinates will be used to create a voronoi graph. 
#' The neighborhood will be composed of all cells with "depth" steps i.e.
#' how many nodes in the graph do you need to go through to get to taget cell?
#' A depth of one will contain all cells in direct contact with center cell. 
#'
#'
#' This cost matrix is then parsed to a
#' Kuhn–Munkres algorithm that will generate point pairs that minimize
#' the overall cost. 
#' 
#' Since the algorithm complexity is O(n3), it can be time consuming to
#' to run on larger data sets. We recomned using "exact" mapping when
#' there are less then 1500 points (a few minutes) and use "div" 
#' with larger data sets. The "div" option will split the query data
#' into random batches of pre-defined size (batch_size). The optimization
#' will be run on each batch serately and the results will be concatenated
#' together. Note this will not be an exact match but an approximation.
#' 
#' Finaly once the matches are found, the coordinates are mapped to its
#' corresponding point and a new object is returned.
#'
#'
#' @return vesalius_assay
#' @export


map_assays <- function(seed_assay,
    query_assay,
    neighborhood = "knn",
    k = 20,
    radius = 0.05,
    depth = 1,
    dimensions = seq(1, 30),
    mapping = "div",
    batch_size = 1000,
    signal = "variable_features",
    use_norm = "raw",
    custom_cost = NULL,
    overwrite = FALSE,
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # First let's get singal
    #-------------------------------------------------------------------------#
    if (!is(query_assay, "list")) {
        query_assay <- list(query_assay)
    }
    signal <- get_features(seed_assay = seed_assay,
        query_assay = query_assay,
        signal = signal,
        dimensions = dimensions,
        use_norm = use_norm,
        verbose = verbose)
    mapped <- point_mapping(seed = seed_assay,
        seed_signal = signal$seed,
        query = query_assay,
        query_signal = signal$query,
        k = k,
        mapping = mapping,
        batch_size = batch_size,
        custom_cost = custom_cost,
        overwrite = overwrite,
        verbose = verbose)
    integrated <- integrate_graph(mapped,
        query_assay,
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
#' @param verbose logical - should progress messages be outputed.
#' @return list contain seed signal, query signal and features used
get_features <- function(seed_assay,
    query_assay,
    signal,
    dimensions = seq(1:30),
    use_norm = "raw",
    verbose = TRUE) {
    message_switch("signal", verbose = verbose)
    #-------------------------------------------------------------------------#
    # First we check if we are using embedding values 
    #-------------------------------------------------------------------------#
    if (any(grepl(pattern = "embeddings", x = signal))) {
        seed_signal <- t(seed_assay@active)
        query_signal <- lapply(query_assay, slot, "active") %>%
            lapply(.,t)
        features <- NA
    } else {
        #---------------------------------------------------------------------#
        # Let's check what feature input we have and then filter 
        #---------------------------------------------------------------------#
        seed_signal <- check_signal(seed_assay, signal, type = use_norm)
        query_signal <- lapply(query_assay, check_signal,
            signal = signal,
            type = use_norm)
        features <- intersect(seed_signal,
            unlist(query_signal, recursive = TRUE))
        if (length(features) == 0) {
            stop("No common features between seed and query data sets!")
        }
        seed_signal <- get_counts(seed_assay, type = use_norm)[features, ]
        query_signal <- lapply(query_assay, function(assays, features, type) {
            return(get_counts(assays, type)[features, ])
        }, features = features, type = use_norm)
        seed_signal <- t(scale(t(as.matrix(seed_signal))))
        query_signal <- lapply(query_signal,function(mat){
            return(t(scale(t(as.matrix(mat)))))})
        }
    return(list("seed" = seed_signal,
        "query" = query_signal,
        "features" = features))

}

#' mapping points between data sets
#' @param seed vesalius_assay object
#' @param seed_signal processed seed signal from seed assay
#' @param query vesalius_assay object
#' @param query_signal processed query signal from query assay
#' @param mapping exact mapping or div (divide and conquer)
#' @param k int size of niche (knn)
#' @param batch_size int number of points in each query batch
#' @param custom_cost matrix - matrix of size n (query cells) by p (seed cells)
#' containing custom cost matrix. Used instead of vesalius cost matrix
#' @param overwrite logical - if custom_cost is not NULL, should this matrix
#' be added to the vesalius matrix or should it overwrite it completely.
#' @param verbose logical - out progress to console
#' @return list of matched and aligned coordinates in query
#' @importFrom future.apply future_lapply
point_mapping <- function(seed,
    seed_signal,
    query,
    query_signal,
    mapping = "div",
    k = 20,
    batch_size = 1000,
    custom_cost = NULL,
    overwrite = FALSE,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # first we compute super pixels from embedding values
    # Note this could be updated to parallel
    #-------------------------------------------------------------------------#
    seed <- get_tiles(seed) %>% filter(origin == 1)
    query <- lapply(query, get_tiles) %>% lapply(., filter, origin == 1)
    #-------------------------------------------------------------------------#
    # size based filter??
    # just checking at this point
    #-------------------------------------------------------------------------#
    custom_cost <- check_cost_matrix_validity(custom_cost,
        seed,
        seed_signal,
        query,
        query_signal)
    #-------------------------------------------------------------------------#
    # Next we create a scoring table for each spix
    #-------------------------------------------------------------------------#
    aligned_indices <- vector("list", length(query))
    for (i in seq_along(query)) {
        if (!is.null(custom_cost) && !overwrite) {
            cost_mat <- compute_cell_cost(custom_cost$seed_signal[[i]],
                custom_cost$query_signal[[i]],
                i,
                verbose)
            
            cost_mat <- compute_neighbor_cost(cost_mat,
                custom_cost$seed[[i]],
                custom_cost$query[[i]],
                custom_cost$seed_signal[[i]],
                custom_cost$query_signal[[i]],
                k,
                i,
                verbose)
            cost_mat <- cost_mat + custom_cost$cost[[i]]
        } else if (!is.null(custom_cost) && overwrite) {
            #-----------------------------------------------------------------#
            # Check if we parse a custom matrix or let vesalius compute the 
            #cost
            #-----------------------------------------------------------------#
           
            cost_mat <- custom_cost$cost[[i]]
        } else {
            cost_mat <- compute_cell_cost(seed_signal,
                query_signal[[i]],
                i,
                verbose)
            cost_mat <- compute_neighbor_cost(cost_mat,
                seed,
                query[[i]],
                seed_signal,
                query_signal[[i]],
                k,
                i,
                verbose)
        }
        #---------------------------------------------------------------------#
        # devide cost matrix
        #---------------------------------------------------------------------#
        if (mapping == "div") {
            message_switch("div_hungarian", verbose)
            cost_mat <- divide_and_conquer(cost_mat, batch_size)
            matched_indices <- match_index_batch(seq_along(cost_mat),
                cost_mat)
            matched_indices <- concat_matches(matched_indices)
        } else {
            message_switch("hungarian", verbose)
            matched_indices <- match_index(cost_mat)
        }
        if (!is.null(custom_cost)) {
             aligned_indices[[i]] <- align_index(matched_indices,
                custom_cost$seed[[i]],
                custom_cost$query[[i]],
                verbose)
        } else {
           aligned_indices[[i]] <- align_index(matched_indices,
                seed,
                query[[i]],
                verbose)
        }
        
    }
    return(aligned_indices)
}

#' compute cell matching cost wrapper
#' @param seed_signal count matrix taken from seed data set
#' @param query_signal count matrix taken from query data set
#' @param i numeric - iteration value in query list
#' @param verbose logical - should progress message be outputed
#' @return cost matrix with nrow = # cells in query and 
#' ncol = # cells in seed
compute_cell_cost <- function(seed_signal, query_signal, i, verbose) {
    message_switch("feature_cost", verbose,
            assay = paste("Query Item", i))
    # seed_local <- t(scale(t(as.matrix(seed_signal))))
    # query_local <- t(scale(t(as.matrix(query_signal))))
    cost_mat <- feature_dissim(seed_signal, query_signal)
    return(cost_mat)
}

#' compute cell neighborhood matching cost
#' @param cost_matrix matrix containing assignement cost 
#' from each cell in query to each cell in seed (rows x cols)
#' @param seed data.frame of coordinates for seed
#' @param query data.frame of coordinates for query
#' @param seed_signal matrix - count matrix for seed
#' @param query_signal matrix - count matrix for query
#' @param k integer - how many neighbors should be considered as 
#' part of the local neighborhood
#' @param i numeric - iteration value in query list
#' @param verbose logical - should progress message be outputed
#' @return cost matrix with nrow = # cells in query and 
#' ncol = # cells in seed
#' @importFrom RANN nn2
compute_neighbor_cost <- function(cost_mat,
    seed,
    query,
    seed_signal,
    query_signal,
    k,
    i,
    verbose) {
    message_switch("neighbor_cost", verbose,
            assay = paste("Query Item", i))
    
    seed_local <- neighbor_expression(seed_dist$nn.idx,
        seed_signal)
    query_dist <- RANN::nn2(query[, c("x", "y")],
            k = k)
    rownames(query_dist$nn.idx) <- query$barcodes
    query_local <- neighbor_expression(query_dist$nn.idx,
        query_signal)

    
    cost_mat <- feature_dissim(seed_local, query_local) +
            cost_mat
    return(cost_mat)
}

knn_neighborhood <- function(coord, k) {
    coord_dist <- RANN::nn2(coord[, c("x", "y")],
        k = k)
    rownames(coord_dist$nn.idx) <- coord$barcodes
    return(coord_dist)
}

radius_neighborhood <- function(coord, radius) {

}

depth_neighborhood <- function(coord, detph) {

}

#' compute the distance between seed and query signals
#' @param seed seed signal
#' @param query query signal
#' @return matrix with query as rows and seed as colmuns
#' @importFrom RANN nn2
feature_dissim <- function(seed, query) {
    cores <- nbrOfWorkers()
    seed_index <- chunk(seq_len(ncol(seed)), ceiling(ncol(seed) / cores))
    query_index <- chunk(seq_len(ncol(query)), ceiling(ncol(query) / cores))

    box <- vector("list", length(query_index))
    for (i in seq_along(query_index)){
        local_query <- t(query[, query_index[[i]]])
        seeds <- future_lapply(seed_index, function(idx, seed, query) {
            local_seed <- t(seed[, idx])
            feature_dist_mat <- RANN::nn2(data = local_seed,
                query = query,
                k = nrow(local_seed))
            feature_dist_mat <- arrange_knn_matrix(feature_dist_mat)
            feature_dist_mat <-
                (feature_dist_mat - min(feature_dist_mat)) /
                (max(feature_dist_mat) - min(feature_dist_mat))
            rownames(feature_dist_mat) <- rownames(query)
            colnames(feature_dist_mat) <- rownames(local_seed)
            return(feature_dist_mat)
        }, seed = seed, query = local_query, future.seed = TRUE)
        box[[i]] <- do.call("cbind", seeds)
    }
    box <- do.call("rbind", box)
    return(box)
}

#' compute average expression of local neighborhood
#' @param distance knn neighbor matrix
#' @param features count matrix or feature matrix to average
#' @return average feature matrix. The expression of each cell 
#' is replace by the average expression of the k nearest neighbors
neighbor_expression <- function(distance, features) {
    neigborhood <- matrix(0,
        ncol = ncol(features),
        nrow = nrow(features))
    colnames(neigborhood) <- colnames(features)
    rownames(neigborhood) <- rownames(features)
    for (i in seq_len(nrow(distance))) {
        local_expression <- features[,
            colnames(features) %in% rownames(distance)[distance[i, ]]]
        local_expression <- rowMeans(local_expression)
        neigborhood[, colnames(features) %in% rownames(distance)[i]] <-
            local_expression
    }
    return(neigborhood)
}


#' splitting cost matrix into batches for batch processing
#' @param cost_mat matrix 
#' @param batch_size integer - size of each batch 
#' @return list containing cost matrix batches 
divide_and_conquer <- function(cost_mat, batch_size) {
    query_batch <- chunk(seq(1, nrow(cost_mat)), batch_size)
    query_idx <- sample(seq(1, nrow(cost_mat)), size = nrow(cost_mat))
    seed_batch <- chunk(seq(1, ncol(cost_mat)),
        ceiling(ncol(cost_mat) / length(query_batch)))
    seed_idx <- sample(seq(1, ncol(cost_mat)), size = ncol(cost_mat))
    cost_list <- vector("list", length(query_batch))
    for (i in seq_along(query_batch)) {
        cost_list[[i]] <- cost_mat[query_idx[query_batch[[i]]],
            seed_idx[seed_batch[[i]]]]
    }
    return(cost_list)
}

#' Finding optimal mapping between query and seed
#' @param cost_matrix matrix with cost between query (rows) and 
#' seed (columns)
#' @details Uses a Hungarian solver to minimize overall cost of 
#' pair assignement.
#' @return data.frame with point in query (from) mapped to 
#' which point in seed (to) as well as cost of mapping for that pair
match_index <- function(cost_matrix) {
    if (nrow(cost_matrix) > ncol(cost_matrix)) {
        cost_matrix <- t(cost_matrix)
        with_transpose <- TRUE
    } else {
        with_transpose <- FALSE
    }
    map <- RcppHungarian::HungarianSolver(cost_matrix)
    mapping <- as.data.frame(map$pairs)
    if (with_transpose){
        colnames(mapping) <- c("to", "from")
        scores <- mapply(function(i, j, cost) {
            return(cost[i, j])
        }, mapping$to, mapping$from, MoreArgs = list(cost_matrix))
        mapping$score <- scores
        mapping$to <- rownames(cost_matrix)[mapping$to]
        mapping$from <- colnames(cost_matrix)[mapping$from]
    } else {
        colnames(mapping) <- c("from", "to")
        scores <- mapply(function(i, j, cost) {
            return(cost[i, j])
        }, mapping$from, mapping$to, MoreArgs = list(cost_matrix))
        mapping$score <- scores
        mapping$to <- colnames(cost_matrix)[mapping$to]
        mapping$from <- rownames(cost_matrix)[mapping$from]
    }
    return(mapping)
}


#' Finding optimal mapping between query and seed using a mini batch 
#' approach
#' @param idx integer - index of batch 
#' @param cost_matrix list of matrices with cost between query (rows) and 
#' seed (columns)
#' @param verbose logical - should progress messages be outpute to the console
#' @details Uses a Hungarian solver to minimize overall cost of 
#' pair assignement.
#' @return data.frame with point in query (from) mapped to 
#' which point in seed (to) as well as cost of mapping for that pair
match_index_batch <- function(idx, cost_matrix) {
    mapping <- future_lapply(idx, function(idx,
        cost_matrix) {
        cost_matrix <- cost_matrix[[idx]]
        if (nrow(cost_matrix) > ncol(cost_matrix)) {
            cost_matrix <- t(cost_matrix)
            with_transpose <- TRUE
        } else {
            with_transpose <- FALSE
        }
        map <- RcppHungarian::HungarianSolver(cost_matrix)
        mapping <- as.data.frame(map$pairs)
        if (with_transpose){
            colnames(mapping) <- c("to", "from")
            scores <- mapply(function(i, j, cost) {
                return(cost[i, j])
            }, mapping$to, mapping$from, MoreArgs = list(cost_matrix))
            mapping$score <- scores
            mapping$to <- rownames(cost_matrix)[mapping$to]
            mapping$from <- colnames(cost_matrix)[mapping$from]
        } else {
            colnames(mapping) <- c("from", "to")
            scores <- mapply(function(i, j, cost) {
                return(cost[i, j])
            }, mapping$from, mapping$to, MoreArgs = list(cost_matrix))
            mapping$score <- scores
            mapping$to <- colnames(cost_matrix)[mapping$to]
            mapping$from <- rownames(cost_matrix)[mapping$from]
        }
        return(mapping)
    }, cost_matrix = cost_matrix,
    future.seed = TRUE)
    return(mapping)
}

#' merging batch matches together
#' @param matched_indices list containing matched points in each batch
#' @return data.frame with matched points
concat_matches <- function(matched_indices) {
    matched_indices <- lapply(matched_indices, function(x) {
        return(x[, c("from", "to", "score")])
    })
    matched_indices <- do.call("rbind", matched_indices)
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
    verbose = TRUE) {
    query$norm_with <- 0
    query$x[match(matched_index$from, query$barcodes)] <-
        seed$x[match(matched_index$to, seed$barcodes)]
    query$y[match(matched_index$from, query$barcodes)] <-
        seed$y[match(matched_index$to, seed$barcodes)]
    query$x <- jitter(query$x, factor = 1)
    query$y <- jitter(query$y, factor = 1)
    return(query)
}



integrate_graph <- function(aligned_graph,
    query,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # rebuild vesalius
    #-------------------------------------------------------------------------#
    for (i in seq_along(query)) {
        local_counts <- get_counts(query[[i]], type = "raw")
        local_counts <- local_counts[,
                colnames(local_counts) %in% aligned_graph[[i]]$barcodes]
        query[[i]] <- build_vesalius_assay(
            coordinates = aligned_graph[[i]][, c("barcodes", "x", "y")],
            counts = local_counts,
            assay = "integrated",
            layer = max(query[[i]]@meta$orig_coord$z) + 1,
            verbose = FALSE)
    }
    if (length(query) == 1) {
        return(query[[1]])
    } else {
        return(query)
    }
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

