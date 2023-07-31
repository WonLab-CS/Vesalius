###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################ HORIZONTAL INTEGRATION ###########################
#-----------------------------------------------------------------------------#


#' Aling and integrate spatial assay from the same modality using super pixels
#' @param seed_assay vesalius_assay object - data to be mapped to
#' @param query_assay vesalius_assay objecy - data to map
#' @param k int ]2, n_points] number of neareset neighbors to be considered for
#' neighborhodd computation.
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
#' nearest neighbors in space. This cost matrix is then parsed to a
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
#' @return vesalius_assay
#' @export


integrate_horizontally <- function(seed_assay,
    query_assay,
    k = 20,
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
    spix_mnn <- point_mapping(seed = seed_assay,
        seed_signal = signal$seed,
        query = query_assay,
        query_signal = signal$query,
        k = k,
        mapping = mapping,
        batch_size = batch_size,
        custom_cost = custom_cost,
        overwrite = overwrite,
        verbose = verbose)
    integrated <- integrate_graph(spix_mnn,
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
    # Next we create a scoring table for each spix
    #-------------------------------------------------------------------------#
    aligned_indices <- vector("list", length(query))
    for (i in seq_along(query)) {
        #---------------------------------------------------------------------#
        # Check if we parse a custom matrix or let vesalius compute the cost
        #---------------------------------------------------------------------#
        if (!is.null(custom_cost) && overwrite) {
            cost_mat <- custom_cost
        } else {
            #-----------------------------------------------------------------#
            # compute cost matrix using euclidean distance between cells 
            # and euclidean distance between local neighbohoods 
            #-----------------------------------------------------------------#
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
            #-----------------------------------------------------------------#
            # If a custom matrix is provided and the user wants to add this
            # this cost to cost matrix prodivided by vesalius 
            #-----------------------------------------------------------------#
            if (!is.null(custom_cost) && !overwrite) {
                custom_cost <- check_cost_matrix(custom_cost, cost_mat)
                cost_mat <- custom_cost$custom + custom_cost$cost
            }
        }
        #---------------------------------------------------------------------#
        # devide cost matrix
        #---------------------------------------------------------------------#
        if (mapping == "div") {
            message_switch("div_hungarian", verbose)
            cost_mat <- divide_and_conquer(cost_mat, batch_size)
            matched_indices <- match_index_batch(seq_along(cost_mat),
                cost_mat,
                verbose = verbose)
            matched_indices <- concat_matches(matched_indices)
        } else {
            message_switch("hungarian", verbose)
            matched_indices <- match_index(cost_mat)
        }
        aligned_indices[[i]] <- align_index(matched_indices,
            seed,
            query[[i]],
            verbose)
    }
    return(aligned_indices)
}


compute_cell_cost <- function(seed_signal, query_signal, i, verbose) {
    message_switch("feature_cost", verbose,
            assay = paste("Query Item", i))
    seed_local <- t(scale(t(as.matrix(seed_signal))))
    query_local <- t(scale(t(as.matrix(query_signal))))
    cost_mat <- feature_dissim(seed_local, query_local)
    return(cost_mat)
}

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
    seed_dist <- RANN::nn2(seed[, c("x", "y")],
        k = k)
    rownames(seed_dist$nn.idx) <- seed$barcodes
    seed_local <- neighbor_expression(seed_dist$nn.idx,
        seed_signal)
    query_dist <- RANN::nn2(query[, c("x", "y")],
            k = k)
    rownames(query_dist$nn.idx) <- query$barcodes
    query_local <- neighbor_expression(query_dist$nn.idx,
        query_signal)

    seed_local <- t(scale(t(as.matrix(seed_local))))
    query_local <- t(scale(t(as.matrix(query_local))))
    cost_mat <- feature_dissim(seed_local, query_local) +
            cost_mat
    return(cost_mat)
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

match_index_batch <- function(idx, cost_matrix, verbose) {
    mapping <- future_lapply(idx, function(idx,
        cost_matrix,
        verbose) {
        tot <- length(cost_matrix)
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
        if (verbose) {
            dyn_message_switch("hung", verbose,
            prog = round(idx / tot, digits = 4) * 100)
        }
        return(mapping)
    }, cost_matrix = cost_matrix,
    verbose = verbose)
    return(mapping)
}

concat_matches <- function(matched_indices) {
    matched_indices <- lapply(matched_indices, function(x) {
        return(x[, c("from", "to", "score")])
    })
    matched_indices <- do.call("rbind", matched_indices)
    return(matched_indices)
}


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
    vesalius_assay <- build_vesalius_assay(coordinates = aligned_graph[[1]][ ,1:3],
        counts = get_counts(query[[1]], type = "raw"),
        assay = "integrated",
        layer = 1,
        verbose = FALSE)
    return(vesalius_assay)

}

build_integrated_matrix <- function(...) {
    data_sets <- list(...)
    gene_set <- unique(unlist(lapply(data_sets, rownames)))
    cells <- unlist(lapply(seq(1, length(data_sets)), function(i, d) {
        return(paste0("data", i, "_", colnames(d[[i]])))
    }, d = data_sets))
    integrated_counts <- Matrix(0,
        ncol = length(cells),
        nrow = length(gene_set))
    colnames(integrated_counts) <- cells
    rownames(integrated_counts) <- gene_set
    return(integrated_counts)
}





#-----------------------------------------------------------------------------#
############################ VERTICAL INTEGRATION #############################
#-----------------------------------------------------------------------------#

#' Integrate jointly measured spatial omic assays
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
#' @export
integrate_vertically <- function(mod1,
    mod2,
    dimensions = seq(1, 30),
    embedding = "last",
    method = "interlace",
    norm_method = "log_norm",
    dim_reduction = "PCA",
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
    # method switch - which method is best 
    #-------------------------------------------------------------------------#
    integrated_embeds <- switch(EXPR = method,
        "interlace" = interlace_embeds(mod1_embed, mod2_embed, dimensions),
        "mean" = average_embed(mod1_embed, mod2_embed, dimensions),
        "concat" = concat_embed(mod1,
            mod2,
            dimensions,
            norm_method,
            dim_reduction))
    integrated_embeds <- list(integrated_embeds)
    names(integrated_embeds) <- method
    integrated <- new("vesalius_assay",
        assay = "integrated",
        embeddings = integrated_embeds,
        active = integrated_embeds[[1]],
        tiles = mod1@tiles)
    simple_bar(verbose)
    return(integrated)

}

interlace_embeds <- function(seed, query, dimensions) {
    seed <- seed[, dimensions]
    query <- query[match(rownames(seed), rownames(query)), dimensions]
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

average_embed <- function(seed, query, dimensions) {
    seed <- seed[, dimensions]
    query <- query[match(rownames(seed), rownames(query)), dimensions]
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

concat_embed <- function(seed,
    query,
    dimensions,
    norm_method,
    dim_reduction) {

    seed_features <- seed@meta$variable_features
    query_features <- query@meta$variable_features
    seed_counts <- get_counts(seed, type = "raw")
    seed_counts <- seed_counts[rownames(seed_counts) %in% seed_features, ]
    query_counts <- get_counts(query, type = "raw")
    query_counts <- query_counts[rownames(query_counts) %in% query_features,
        match(colnames(seed_counts), colnames(query_counts))]
    integrated_counts <- rbind(seed_counts, query_counts)
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

