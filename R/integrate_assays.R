###############################################################################
################################   Vesalius      ##############################
###############################################################################

#--------------------------/Integrating & Alinging/---------------------------#


#' Aling and integrate spatial assay from the same modality using super pixels
#' @param seed_assay vesalius_assay object - data to be mapped to
#' @param query_assay vesalius_assay objecy - data to map
#' @param seed_trial name of embedding to use for super pixel generation
#' @param query_trial name of embedding to use for super pixel generation
#' @param scoring_method method used to score the similarity between super 
#' pixels (pearson, spearman, kendall, coherence, index)
#' @param dimensions integer vector - which latent space dimensions should be
#' used for super pixel generation
#' @param scaling numeric ]0,1] describing image scale to consider 
#' during super pixel selection
#' @param compactness numeric ]0,Inf] - importance of the spatial component 
#' in relation to color similarity. See details
#' @param n_centers integer [3, max(spatial_index)] - number of super pixels
#' to generate in each image. See details
#' @param iter integer [1, Inf] - Number of iteration during graph matching
#' phase.
#' @param index_selection character (random, bubble) - how should initial super
#' pixel locations be selected. See detials
#' @param threshold numeric [0,1[ - similarity score threshold. Only super pixel
#' that score above this threshold will be used for graph matching
#' @param n_anchors integer [3, n_centers] - Number of graph anchors used during
#' graph matching
#' @param mut_extent numeric [0,1] - extent of alignment graph that can be 
#' subjected to mutations
#' @param mut_prob numeric [0,1] - probability of alignment graph mutation.
#' @param allow_vertex_merge logical - Should graph vertices be merged or should
#' them be repelled to new location. See details
#' @param signal character (features, counts, embeddings, "custom") - What should 
#' be used as cell signal for super pixel scoring. Seed details 
#' @param verbose logical - should I be a noisy boy?
#' 
#' @export
#' 

integrate_horizontally <- function(seed_assay,
    query_assay,
    seed_trial = "last",
    query_trial = "last",
    scoring_method = "pearson",
    dimensions = seq(1, 30),
    scaling = 0.2,
    compactness = 1,
    n_centers = 2000,
    iter = 10000,
    index_selection = "random",
    threshold = 0.7,
    n_anchors = 20,
    mut_extent = 0.1,
    mut_prob = 0.2,
    allow_vertex_merge = FALSE,
    signal = "features",
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # compute slic for both assays
    #-------------------------------------------------------------------------#
    message_switch("seg", verbose = verbose, method = "slic")
    seed_trial <- slic_segmentation(seed_assay,
        dimensions = dimensions,
        col_resolution = n_centers,
        embedding = seed_trial,
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    query_trial <- slic_segmentation(query_assay,
        dimensions = dimensions,
        col_resolution = n_centers,
        embedding = query_trial,
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    
    #-------------------------------------------------------------------------#
    # get signal either as counts or as embedding values
    #-------------------------------------------------------------------------#
    message_switch("signal", verbose = verbose)
    seed_signal <- check_signal(signal, seed_assay, type = "raw")
    query_signal <- check_signal(signal, query_assay, type = "raw")
    if (grepl(pattern = "embeddings", x = signal)) {
        seed_signal <- seed_trial$active
        query_signal <- query_trial$active
    } else {
        seed_counts <- get_counts(seed_assay, type = "raw")
        query_counts <- get_counts(query_assay, type = "raw")
        seed_genes <- intersect(seed_signal, query_signal)
        if (length(seed_genes) == 0) {
            stop("No common features between seed and query data sets!")
        }
        seed_signal <- seed_counts[seed_genes, ]
        query_signal <- query_counts[seed_genes, ]
    }
    #-------------------------------------------------------------------------#
    # compact signal so we only keep average signal per spix
    #-------------------------------------------------------------------------#
    seed_signal <- compress_signal(seed_signal, seed_trial$segments)
    query_signal <- compress_signal(query_signal, query_trial$segments)
    #-------------------------------------------------------------------------#
    # Get estimated super pixel centers
    # we also generate a graph and score this graph. 
    # score the correlation between each vertex in seed/query graph
    #-------------------------------------------------------------------------#
    mapped <- som_map(seed_trial$segments,
        seed_signal,
        query_trial$segments,
        query_signal,
        anchors = n_anchors)
    
    # message_switch("slic_graph", verbose = verbose, data = "seed")
    # seed_graph <- generate_slic_graph(seed_trial$segments)

    # message_switch("slic_graph", verbose = verbose, data = "query")
    # query_graph <- generate_slic_graph(query_trial$segments)
    # #-------------------------------------------------------------------------#
    # # Now we can compute the same thing but between each graph
    # # For now - we check the number of centeres in case of mismatch
    # # for now I am getting mismatch even with random sampling which 
    # # does not make any sense but to check and fix
    # #-------------------------------------------------------------------------#
    # q_centers <- length(unique(query_graph$from))
    # s_centers <- length(unique(seed_graph$from))
    # integrated_graph <- data.frame(
    #     "from" = rep(unique(seed_graph$from), each = q_centers),
    #     "to" = rep(unique(query_graph$from), times = s_centers))


    # spix_score <- score_graph(integrated_graph,
    #     signal = list(seed_signal, query_signal),
    #     verbose = verbose)
    # # matched_graph <- match_vertex(seed_graph,
    # #         query_graph,
    # #         spix_score,
    # #         depth = depth,
    # #         threshold = threshold,
    # #         verbose = verbose)
    # matched_graph <- match_graph(seed_graph = seed_graph,
    #     seed_trial = seed_trial$segments,
    #     query_graph = query_graph,
    #     query_trial = query_trial$segments,
    #     score = spix_score,
    #     scoring_method = "pearson",
    #     threshold = threshold,
    #     iter = iter,
    #     n_anchors = n_anchors,
    #     mut_extent = mut_extent,
    #     mut_prob = mut_prob,
    #     verbose = verbose)
    aligned_graph <- align_graph(mapped,
        seed_trial$segments,
        seed_graph,
        query_trial$segments,
        query_graph,
        verbose = verbose)
    integrated <- integrate_graph(aligned_graph,
        seed_assay,
        seed_trial$segments,
        query_assay,
        verbose = verbose)
    simple_bar(verbose)
    return(integrated)
}

#' @export
integrate_vertically <- function(seed,
    query,
    dimensions = seq(1, 30),
    embedding = "last",
    method = "interlace",
    norm_method = "raw",
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
    seed_embed <- check_embedding_selection(seed, embedding, dimensions)
    query_embed <- check_embedding_selection(query, embedding, dimensions)
    #-------------------------------------------------------------------------#
    # method switch - which method is best 
    #-------------------------------------------------------------------------#
    integrated_embeds <- switch(EXPR = method,
        "interlace" = interlace_embeds(seed_embed, query_embed, dimensions),
        "mean" = average_embed(seed_embed, query_embed, dimensions),
        "concat" = concat_embed(seed,
            query,
            dimensions,
            norm_method,
            dim_reduction))
    integrated_embeds <- list(integrated_embeds)
    names(integrated_embeds) <- method
    integrated <- new("vesalius_assay",
        assay = "integrated",
        embeddings = integrated_embeds,
        active = integrated_embeds[[1]],
        tiles = seed@tiles)
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



compress_signal <- function(signal, segments) {
    segments <- split(segments, segments$Segment)
    compressed_signal <- vector("list", length(segments))
    names(compressed_signal) <- names(segments)
    for (i in seq_along(segments)) {
        local_signal <- signal[, segments[[i]]$barcodes]
        if (is.null(ncol(local_signal)) || ncol(local_signal) == 1) {
            compressed_signal[[i]] <- as.vector(local_signal)
        } else {
            compressed_signal[[i]] <- apply(local_signal, 1, mean)
        }
    }
    return(compressed_signal)
}

#' importFrom kohonen som map scale somgrid
som_map <- function(seed_trial,
    seed_signal,
    query_trial,
    query_signal,
    anchors){
    seed_spix <- get_super_pixel_centers(seed_trial)
    query_spix <- get_super_pixel_centers(query_trial)
    spixs <- seed_spix$center
    anchors <- sqrt(length(spixs))
    seed_signal <- cbind(seed_spix[,c("x","y")],
        do.call("rbind", seed_signal))
    seed_som <- som(scale(seed_signal),
        grid = somgrid(x = floor(anchors),
            y = ceiling(anchors), "hexagonal"))
    query_signal <- cbind(query_spix[,c("x","y")],
        do.call("rbind", query_signal))
    mapped <- kohonen::map(x = seed_som, newdata = scale(query_signal))
    anchor_map <- cbind(seed_spix, mapped$unit.classif)
    anchor_map$anchor <- 1
    colnames(anchor_map) <- c("x","y","to","from","anchor")
    return(anchor_map)
}


generate_slic_graph <- function(spix,
    k = "auto") {
    #-------------------------------------------------------------------------#
    # first we estimate the pixel centers 
    # we will use this as an estimate for nearest neighbor calculation 
    #-------------------------------------------------------------------------#
    centers <- get_super_pixel_centers(spix) %>% select(c("x", "y"))
    #-------------------------------------------------------------------------#
    # Next we get the nearest neighbors
    # If we use the auto option we compute the nearest neighbors based on 
    # delauney traingulation - not that this will results in an uneven number 
    # of neighbors per point. 
    #-------------------------------------------------------------------------#
    if (k != "auto") {
        k <- min(c(k, nrow(centers)))
        knn <- RANN::nn2(centers, k = k)$nn.idx
        rownames(knn) <- rownames(centers)
        graph <- populate_graph(knn)
    } else {
       graph <- graph_from_voronoi(centers)
    }
    return(graph)
}

graph_from_voronoi <- function(centers) {
    voronoi <- deldir::deldir(x = as.numeric(centers$x),
        y = as.numeric(centers$y))$delsgs
    center <- seq_len(nrow(centers))
    graph <- lapply(center, function(idx, voronoi){
        tri <- voronoi %>% filter(ind1 == idx | ind2 == idx)
        tri <- unique(c(tri$ind1, tri$ind2))
        graph <- data.frame("from" = rep(idx, length(tri)),
            "to" = tri)
        return(graph)
    }, voronoi = voronoi) %>%
    do.call("rbind", .)
    return(graph)
}

get_super_pixel_centers <- function(spix) {
    center_pixels <- sort(unique(spix$Segment))
    centers <- future_lapply(center_pixels, function(center, segments) {
        x <- median(segments$x[segments$Segment == center])
        y <- median(segments$y[segments$Segment == center])
        df <- data.frame("x" = x,
            "y" = y,
            "center" = center)
        rownames(df) <- center
        return(df)
    }, segments = spix) %>% do.call("rbind", .)
    return(centers)
}

#' importFrom future.apply future_lapply
score_graph <- function(graph,
    signal,
    scoring_method = "spearman",
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # assuming that if the input to signal is a list
    # we are comparing 2 data sets
    #-------------------------------------------------------------------------#
    if (is(signal, "list") && length(signal) == 2) {
        seed_signal <- signal[[1]]
        query_signal <- signal[[2]]
    } else {
        seed_signal <- query_signal <- signal
    }
    #-------------------------------------------------------------------------#
    # Create a score graph 
    # and compute correlation between super pixels
    #-------------------------------------------------------------------------#
    graph <- data.frame(graph,
        "score" = rep(0, nrow(graph)))
    for (i in seq_len(nrow(graph))) {
        dyn_message_switch("score_graph", verbose,
            prog = round(i / nrow(graph), 4) * 100)
        c1 <- seed_signal[[graph$from[i]]]
        c2 <- query_signal[[graph$to[i]]]
        graph$score[i] <- cor(c1, c2, method = scoring_method)
    }
    if (verbose){cat("\n")}
    return(graph)
}


match_graph <- function(seed_graph,
    seed_trial,
    query_graph,
    query_trial,
    score,
    scoring_method = "pearson",
    threshold = 0.7,
    iter = 10000,
    n_anchors = 25,
    mut_extent = 0.1,
    mut_prob = 0.3,
    allow_vertex_merge = FALSE,
    verbose = verbose) {
    #-------------------------------------------------------------------------#
    # Initialize optimisation 
    #-------------------------------------------------------------------------#
    
    if (length(unique(score$from)) < n_anchors ||
        length(unique(score$to)) < n_anchors) {
        n_anchors <- min(c(length(unique(score$from)),
            length(unique(score$to))))
        message_switch("anchors_found", verbose, anchors = n_anchors)
    }
    sub_sample <- sample(seq(1,max(score$from)),
        size = n_anchors,
        replace = FALSE)
    seed_anchors <- get_best_vertex(score)$from[sub_sample]
    query_anchors <- get_best_vertex(score)$to[sub_sample]
    centers_1 <- get_super_pixel_centers(seed_trial)
    centers_1$x <- min_max(centers_1$x)
    centers_1$y <- min_max(centers_1$y)
    centers_2 <- get_super_pixel_centers(query_trial)
    centers_2$x <- min_max(centers_2$x)
    centers_2$y <- min_max(centers_2$y)

    indiv_seed <- list("chromosome" = rep(1, n_anchors * 2),
        "indiv" = seed_anchors,
        "score" = sum(rep(1, n_anchors * 2)))
    indiv_query <- list("chromosome" = rep(0, n_anchors * 2),
        "indiv" = query_anchors,
        "score" = sum(rep(1, n_anchors * 2)))
    #-------------------------------------------------------------------------#
    # iterating over random graph and finding best match 
    # at the moment it is super basic
    #-------------------------------------------------------------------------#
    for (i in seq_len(iter)) {
        dyn_message_switch("graph_matching", verbose,
            prog = round(i / iter, 4) * 100)
        cent_1 <- centers_1[match(seed_anchors, centers_1$center), c("x", "y")]
        cent_2 <- centers_2[match(query_anchors, centers_2$center), c("x", "y")]
        dist_1 <- RANN::nn2(data = cent_1,
            k = 2)$nn.dist[, 2]
        dist_2 <- RANN::nn2(data = cent_2,
            k = 2)$nn.dist[, 2]
        loc <- paste0(score$from, "_", score$to) %in%
            paste0(seed_anchors, "_", query_anchors)
        #scores <- 1 - score[loc, "score"]
        scores <- c(abs(dist_1 - dist_2))#, scores)
        if (indiv_seed$score > sum(scores) || indiv_query$score > sum(scores)) {
            indiv_seed$chromosome <- scores
            indiv_query$chromosome <- scores
            indiv_seed$indiv <- seed_anchors
            indiv_query$indiv <- query_anchors
            indiv_seed$score <- sum(scores)
            indiv_query$score <- sum(scores)
            if (runif(1) <= mut_prob) {
                locs <- sample(seq(1, length(seed_anchors)),
                    size = n_anchors * mut_extent)
                new_seed <- sample(centers_1$center[!centers_1$center
                    %in% seed_anchors[-locs]],
                    size = length(locs))
                seed_anchors[locs] <- new_seed
                new_query <- sample(centers_2$center[!centers_2$center
                    %in% query_anchors[-locs]],
                    size = length(locs))
                query_anchors[locs] <- new_query
            } else {
                seed_anchors <- sample(seed_anchors,
                    size = n_anchors,
                    replace = FALSE)
                query_anchors <- sample(query_anchors,
                    size = n_anchors,
                    replace = FALSE)
            }
        }
    }
    if (verbose){cat("\n")}
    anchors <- data.frame("from" = indiv_seed$indiv,
        "to" = indiv_query$indiv,
        "score" = indiv_query$score / n_anchors,
        "anchor" = rep(1, n_anchors))
    return(anchors)
}


match_vertex <- function(seed_graph,
    query_graph,
    scores,
    depth = 1,
    threshold = 0.8,
    verbose = TRUE) {
    message_switch("matching_graphs", verbose)
    seed_paths <- graph_path_length(seed_graph)
    query_paths <- graph_path_length(query_graph)
    max_depth <- min(c(max(seed_paths), max(query_paths), depth))
    best_match <- get_best_vertex(scores, rank = depth)
    best_match <- data.frame(best_match,
        matrix(0, ncol = max_depth, nrow = nrow(best_match)))
    colnames(best_match) <- c(colnames(best_match)[1:3],
        paste0("depth_", seq(1, max_depth)))
    for (d in seq_len(max_depth)) {
        for (i in seq_len(nrow(best_match))) {
            seed_niche <- seed_paths[, as.character(best_match$from[i])]
            seed_niche <- names(seed_niche)[seed_niche <= d &
                seed_niche > 0]
            query_niche <- query_paths[, as.character(best_match$to[i])]
            query_niche <- names(query_niche)[query_niche <= d &
                query_niche > 0]
            niche <- scores[scores$to %in% query_niche, ] %>%
                get_best_vertex()
            overlap <- min(c(sum(niche$from %in% seed_niche),
                min(c(length(query_niche), length(seed_niche))))) /
                min(c(length(query_niche), length(seed_niche)))
            depth_loc <- paste0("depth_", d)
            best_match[i, depth_loc] <- overlap
        }
    }
    anchors <- apply(best_match[, 3:ncol(best_match)], 1, mean) >= threshold
    best_match$anchor <- 0
    best_match$anchor[anchors] <- 1
    message_switch("anchors_found", verbose, anchors = sum(best_match$anchor))
    if (sum(best_match$anchor) == 0) {
        stop("No anchors found! Consider relaxing selection threshold.")
    }
    return(best_match)
}



align_graph <- function(matched_graph,
    seed,
    seed_graph,
    query,
    query_graph,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # First get all anchor trajectories 
    # we used the seed coordinates as center pixel
    #-------------------------------------------------------------------------#
    
    seed_centers <- get_super_pixel_centers(seed)
    query_centers <- get_super_pixel_centers(query)
    anchors <- matched_graph %>%
        filter(anchor == 1) %>%
        select(c("from", "to"))
    anchors$angle <- 0
    anchors$distance <- 0
    message_switch("get_traj", verbose)
    for (i in seq_len(nrow(anchors))) {
        seed_point <- seed_centers[seed_centers$center == anchors$from[i], ]
        query_point <- query_centers[query_centers$center == anchors$to[i], ]
        #angle <- polar_angle(seed_point$x, seed_point$y,
        #    query_point$x, query_point$y)
        angle <- atan2(seed_point$y - query_point$y,
            seed_point$x - query_point$x)
        anchors$angle[i] <- angle
        distance <- matrix(c(seed_point$x, query_point$x,
            seed_point$y, query_point$y), ncol = 2)
        distance <- as.numeric(dist(distance))
        anchors$distance[i] <- distance
    }
    anchors <- anchors[order(anchors$to), ]

    #-------------------------------------------------------------------------#
    # get closest anchor point for all un assigned spix points 
    #-------------------------------------------------------------------------#
    #unassinged <- matched_graph %>% filter(anchor == 0)
    #query_point <- query_centers[query_centers$center %in% unassinged$to, ]
    anchor_point <- query_centers[query_centers$center %in% anchors$to, ]
    #-------------------------------------------------------------------------#
    # Apply compound trajectories to individual points points
    #-------------------------------------------------------------------------#
    message_switch("apply_traj", verbose)
    nn <- RANN::nn2(data = anchor_point[, c("x", "y")],
        query = query_centers[, c("x", "y")],
        k = 1)
    angle <- anchors$angle[nn$nn.idx[, 1]]
    distance <- anchors$distance[nn$nn.idx[, 1]]
    browser()
    
    query_centers$x_new <- query_centers$x + (distance * cos(angle))
    query_centers$y_new <- query_centers$y + (distance * sin(angle))
    #-------------------------------------------------------------------------#
    # get closest seed point for all un assigned points
    #-------------------------------------------------------------------------#
    seed_nn <- RANN::nn2(data = seed_centers[, c("x", "y")],
        query = query[, c("x", "y")],
        k = 1)
    query$norm_with <- seed_centers$center[seed_nn$nn.idx[, 1]]
    
    return(query)
}




get_best_vertex  <- function(score, rank = 1) {
    from <- split(score, score$to)
    best_match <- lapply(from, function(f) {
        return(f[order(f$score, decreasing = TRUE)[rank], ])
    }) %>% do.call("rbind", .)
    return(best_match)
}

#'@importFrom igraph graph_from_data_frame distances
#'@importFrom igraph E
graph_path_length <- function(graph) {
    gr <- igraph::graph_from_data_frame(graph, directed = FALSE)
    #igraph::E(gr)$weight <- graph$score
    path_length <- igraph::distances(gr)
    return(path_length)
}

#' @importFrom dplyr select
integrate_graph <- function(aligned_graph,
    seed,
    seed_spix,
    query,
    norm = "logNorm",
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # first let's get the raw counts from each 
    #-------------------------------------------------------------------------#
    seed_counts <- get_counts(seed)
    query_counts <- get_counts(query)
    genes <- intersect(rownames(seed_counts), rownames(query_counts))
    seed_counts <- seed_counts[rownames(seed_counts) %in% genes, ]
    query_counts <- query_counts[rownames(query_counts) %in% genes, ]
    #-------------------------------------------------------------------------#
    # next split the aligned graph by spix
    # and we can scale the counts in each super pixel 
    #-------------------------------------------------------------------------#
    spix <- split(aligned_graph,
        list(aligned_graph$Segment, aligned_graph$norm_with))
    spix <- spix[sapply(spix, function(x){return(nrow(x)>0)})]
    
    integrated_counts <- vector("list", length(spix))
    query_names <- c()
    for (i in seq_along(spix)) {
        dyn_message_switch("integrate_graph", verbose,
            prog = round(i / length(spix), 4) * 100)
        query_names <- c(query_names,spix[[i]]$barcodes)
        query_local <- query_counts[, spix[[i]]$barcodes]
        seed_local <- seed_counts[, seed_spix$barcodes[
            seed_spix$Segment %in% spix[[i]]$norm_with]]
        max_count <- max(seed_local)
        if (is.null(ncol(seed_local))) {
            sd_count <- sd(query_local)
        } else {
            sd_count <- max(apply(seed_local, 1, sd))
        }
        query_local <- query_local * (max_count / sd_count)

        if (is.null(ncol(query_local))) {
           query_local <- query_local[order(names(query_local))]
        } else {
            query_local <- query_local[order(rownames(query_local)), ]
        }
        integrated_counts[[i]] <- query_local
    }
    if (verbose){cat("\n")}
    #-------------------------------------------------------------------------#
    # Bind everything together - this needs to be cleaned
    #-------------------------------------------------------------------------#
    integrated_counts <- do.call("cbind", integrated_counts)
    # colnames(integrated_counts) <- make.unique(paste0("Query_",
    #     query_names))
    # colnames(seed_counts) <- paste0("Seed_",
    #     colnames(seed_counts))
    # integrated_names <- make.unique(c(colnames(integrated_counts),
    #     colnames(seed_counts)), sep = "_")
    # integrated_counts <- cbind(integrated_counts, seed_counts)
    # colnames(integrated_counts) <- integrated_names
    
    #browser()
    #-------------------------------------------------------------------------#
    # Now re can rebuild the coordinates 
    #-------------------------------------------------------------------------#
    # seed_coordinates <- seed@tiles %>%
    #     filter(origin == 1) %>%
    #     select(c("barcodes", "x", "y"))
    # seed_coordinates$barcodes <- paste0("Seed_", seed_coordinates$barcodes)
    query_coordinates <- aligned_graph[, c("barcodes", "x", "y")]
    # query_coordinates$barcodes <- make.unique(paste0("Query_", query_coordinates$barcodes))
    common_loc <- intersect(query_coordinates$barcodes,
        colnames(integrated_counts))
    query_coordinates <- query_coordinates[query_coordinates$barcodes %in% common_loc, ]
    integrated_counts <- integrated_counts[, common_loc]
    
    query_coordinates$barcodes <- gsub("_et_", "_l_", common_loc)
    colnames(integrated_counts) <- gsub("_et_", "_l_", common_loc)
    # integrated_coordinates <- rbind(query_coordinates, seed_coordinates)
    # integrated_coordinates$barcodes <- make.unique(
    #     integrated_coordinates$barcodes,
    #     sep = "_")
    
    vesalius_assay <- build_vesalius_assay(coordinates = query_coordinates,
        counts = integrated_counts,
        assay = "integrated",
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




#' compute signal similarity between 2 sets of territories
#' @param seed_path list of x and y coordinates of seed territories
#' @param query_path list of x and y coordinates of query territories
#' @param domain similarity method to use (index or dtw)
#' @param partials int - number of FFT partials to use
#' @param threshold numeric time/freq threshold
#' @param weight vector of 3 numeric weights - weight used in index calcultion
#' @returns a matrix containing similarity scores between each territory in each 
#' data set. 
signal_similiarity <- function(seed_path,
    query_path,
    domain = "index",
    partials = 200,
    threshold = 0.7,
    weight = rep(0.33, 3)) {
    sim <- switch(EXPR = domain,
        "index" = index(seed_path,
            query_path,
            partials,
            threshold,
            weight),
        "pearson" = correlation(seed_path, query_path, "pearson"),
        "spearman" = correlation(seed_path, query_path, "spearman"),
        "kendall" = correlation(seed_path, query_path, "kendall"),
        "coherence" = spectral_coherence(seed_path, query_path))
    return(sim)
}

#' circular cross correlation
#' @param seed seed signal as vector 
#' @param query query signal as vector
#' @returns circularised cross correlation
#' @importFrom stats ccf
circular_xcorr <- function(seed, query) {
    corr <- ccf(seed, query, plot = FALSE)$acf
    mid <- length(corr) / 2
    corr <- corr[seq(1, mid)] + corr[seq(mid + 1, length(corr))]
    return(corr)

}



#' calculate weighted index of similarity between 2 signals
#' @param seed_path gene signal of seed territories
#' @param query_path gene signal of query territories
#' @param partials int - number of FFT partials to use
#' @param threshold numeric time/freq threshold
#' @param weight vector of 3 numeric weights - weight used in index calcultion
#' @returns matrix of similarity indices between each territory
#' @importFrom stats fft
index <- function(seed_path,
    query_path,
    partials = 500,
    threshold = 0.7,
    weight = rep(0.3, 3)) {
    #-------------------------------------------------------------------------#
    # Remove Isolated territories 
    #-------------------------------------------------------------------------#
    seed_path <- seed_path[names(seed_path) != "isolated"]
    query_path <- query_path[names(query_path) != "isolated"]
    indexed <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    for (seed in seq_along(seed_path)){
        for (query in seq_along(query_path)){
            #-----------------------------------------------------------------#
            # Get partial max and subset paths
            #-----------------------------------------------------------------#
            seed_sub <- seed_path[[seed]]
            query_sub <- query_path[[query]]
            #-----------------------------------------------------------------#
            # compute signal similarity
            #-----------------------------------------------------------------#
            partial_min <- min(c(length(seed_sub), partials))
            thresh <- threshold * partial_min
            time_stat <- max(seed_sub * query_sub)
            time_shifted <- circular_xcorr(seed_sub, query_sub)
            time_shifted <- sum(time_shifted[time_shifted > thresh])
            freq_static <- max(Re(fft(seed_sub))[seq(1, partial_min)] *
                Re(fft(query_sub))[seq(1, partial_min)])

            freq_shifted <- Re(fft(seed_sub * query_sub))
            freq_shifted <- sum(freq_shifted[freq_shifted > thresh])
            indexed[query, seed] <- weight[1] * time_stat +
                weight[1] * time_shifted +
                weight[2] * freq_static +
                weight[3] * freq_shifted
        }
    }
    colnames(indexed) <- names(seed_path)
    rownames(indexed) <- names(query_path)
    return(indexed)
}

#' importFrom gsignal mscohere
spectral_coherence <- function(seed_path, query_path) {
    seed_path <- seed_path[names(seed_path) != "isolated"]
    query_path <- query_path[names(query_path) != "isolated"]
    coh <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    for (seed in seq_along(seed_path)){
        for (query in seq_along(query_path)){
            #-----------------------------------------------------------------#
            # Get partial max and subset paths
            #-----------------------------------------------------------------#
            seed_sub <- seed_path[[seed]]
            query_sub <- query_path[[query]]
            signal <- cbind(seed_sub, query_sub)
            #-----------------------------------------------------------------#
            # coherence 
            #-----------------------------------------------------------------#
            coh[query, seed] <- max(gsignal::mscohere(signal)$coh)
        }
    }
    colnames(coh) <- names(seed_path)
    rownames(coh) <- names(query_path)
    return(coh)
}

correlation <- function(seed_path, query_path, method) {
    seed_path <- seed_path[names(seed_path) != "isolated"]
    query_path <- query_path[names(query_path) != "isolated"]
    co <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    for (seed in seq_along(seed_path)){
        for (query in seq_along(query_path)){
            #-----------------------------------------------------------------#
            # Get partial max and subset paths
            #-----------------------------------------------------------------#
            seed_sub <- seed_path[[seed]]
            query_sub <- query_path[[query]]
            #-----------------------------------------------------------------#
            # correlation  
            #-----------------------------------------------------------------#
            co[query, seed] <- cor(seed_sub, query_sub, method = method)
        }
    }
    colnames(co) <- names(seed_path)
    rownames(co) <- names(query_path)
    return(co)
}
