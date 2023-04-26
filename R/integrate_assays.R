################################################################################
################################   Vesalius      ###############################
###############################################################################

#------------------------------/Fourier Transform/----------------------------#



#' @export
#' 

integrate_assays <- function(seed_assay,
    query_assay,
    seed_trial = "last",
    query_trial = "last",
    scoring_method = "pearson",
    dimensions = seq(1, 30),
    scaling = 0.2,
    compactness = 1,
    n_centers = 2000,
    max_iter = 1000,
    index_selection = "random",
    threshold = 0.9,
    k = "auto",
    signal = "features",
    expand_genes = TRUE,
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # compute slic for both assays
    #-------------------------------------------------------------------------#
    message_switch("seg", verbose = verbose, method = "slic")
    seed_trial <- slic_segmentation(seed_assay,
        dimensions = dimensions,
        col_resolution = n_centers,
        embedding = "last",
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    query_trial <- slic_segmentation(query_assay,
        dimensions = dimensions,
        col_resolution = n_centers,
        embedding = "last",
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
    message_switch("slic_graph", verbose = verbose, data = "seed")
    seed_graph <- generate_slic_graph(seed_trial$segments, k = k)

    message_switch("slic_graph", verbose = verbose, data = "query")
    query_graph <- generate_slic_graph(query_trial$segments, k = k)
    #-------------------------------------------------------------------------#
    # Now we can compute the same thing but between each graph
    # For now - we check the number of centeres in case of mismatch
    # for now I am getting mismatch even with random sampling which 
    # does not make any sense but to check and fix
    #-------------------------------------------------------------------------#
    q_centers <- length(unique(query_graph$from))
    s_centers <- length(unique(seed_graph$from))
    integrated_graph <- data.frame(
        "from" = rep(unique(seed_graph$from), each = q_centers),
        "to" = rep(unique(query_graph$from), times = s_centers))
    message_switch("score_graph", verbose = verbose)
    spix_score <- score_graph(integrated_graph,
        signal = list(seed_signal, query_signal))
    message_switch("matching_graphs", verbose = verbose)
    matched_graph <- match_vertex(seed_graph, query_graph, spix_score)
    aligned_graph <- align_graph(matched_graph,
        seed_trial$segments,
        seed_graph,
        query_trial$segments,
        query_graph)
    integrated <- integrate_graph(aligned_graph,
        seed_assay,
        seed_trial$segments,
        query_assay,
        expand_genes = expand_genes)
    simple_bar(verbose)
    return(integrated)
}

compress_signal <- function(signal, segments) {
    segments <- split(segments, segments$Segment)
    compressed_signal <- vector("list", length(segments))
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
    scoring_method = "spearman") {
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
    cat("\n")
    for (i in seq_len(nrow(graph))) {
        cat(paste0(i, "\r"))
        c1 <- seed_signal[[graph$from[i]]]
        c2 <- query_signal[[graph$to[i]]]
        graph$score[i] <- cor(c1, c2, method = scoring_method)
    }
    return(graph)
}


match_vertex <- function(seed_graph,
    query_graph,
    scores,
    depth = 1,
    threshold = 0.8) {
    seed_paths <- graph_path_length(seed_graph)
    query_paths <- graph_path_length(query_graph)
    max_depth <- min(c(max(seed_paths), max(query_paths)))
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
    best_match$anchor <- 0
    anchors <- best_match$score >= threshold &
        apply(best_match[, grep("depth", colnames(best_match))], 1, mean) >=
        threshold
    best_match$anchor[anchors] <- 1
    return(best_match)
}



align_graph <- function(matched_graph,
    seed,
    seed_graph,
    query,
    query_graph) {
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
    for (i in seq_len(nrow(anchors))) {
        seed_point <- seed_centers[seed_centers$center == anchors$from[i], ]
        query_point <- query_centers[query_centers$center == anchors$to[i], ]
        angle <- polar_angle(seed_point$x, seed_point$y,
            query_point$x, query_point$y)
        anchors$angle[i] <- angle
        distance <- matrix(c(seed_point$x, query_point$x,
            seed_point$y, query_point$y), ncol = 2)
        distance <- as.numeric(dist(distance))
        anchors$distance[i] <- distance
    }

    #-------------------------------------------------------------------------#
    # get closest anchor point for all un assigned spix points 
    #-------------------------------------------------------------------------#
    unassinged <- matched_graph %>% filter(anchor == 0)
    query_point <- query_centers[query_centers$center %in% unassinged$to, ]
    anchor_point <- query_centers[query_centers$center %in% anchors$to, ]
    #-------------------------------------------------------------------------#
    # Apply compound trajectories to individual points points
    #-------------------------------------------------------------------------#
    nn <- RANN::nn2(data = anchor_point[, c("x", "y")],
        query = query[, c("x", "y")],
        k = 2)
    angle <- anchors$angle[nn$nn.idx[, 1]] - anchors$angle[nn$nn.idx[, 2]]
    distance <- sqrt(((anchors$distance[nn$nn.idx[, 1]])^2 + 
        (anchors$distance[nn$nn.idx[, 2]])^2))
    query$x <- query$x + distance * cos(angle * pi / 180)
    query$y <- query$y + distance * sin(angle * pi / 180)
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


integrate_graph <- function(aligned_graph,
    seed,
    seed_spix,
    query,
    norm = "logNorm",
    expand_genes = TRUE) {
    #-------------------------------------------------------------------------#
    # first let's get the raw counts from each 
    #-------------------------------------------------------------------------#
    seed_counts <- get_counts(seed)
    query_counts <- get_counts(query)
    #-------------------------------------------------------------------------#
    # next split the aligned graph by spix
    # and we can scale the counts in each super pixel 
    #-------------------------------------------------------------------------#
    spix <- split(aligned_graph, aligned_graph$Segment)
    integrated_counts <- vector("list", length(spix))
    for (i in seq_along(spix)) {
        query_local <- query_counts[, spix[[i]]$barcodes]
        if (is.null(ncol(query_local))){
            query_local <- matrix(query_local, ncol = 1)
            rownames(query_local) <- names(query_local)
        }
        colnames(query_local) <- paste0("data2_", colnames(query_local))
        seed_local <- seed_counts[, seed_spix$barcodes[
            seed_spix$Segment %in% spix[[i]]$norm_with]]
        if (is.null(ncol(seed_local))){
            seed_local <- matrix(seed_local, ncol = 1)
            rownames(seed_local) <- names(seed_local)
        }
        
        colnames(seed_local) <- paste0("data1_", colnames(seed_local))
        max_count <- max(c(max(seed_local), max(query_local)))
        sd_count <- max(c(apply(seed_local, 2, sd), apply(query_local, 2, sd)))
        seed_local <- seed_local * (max_count / sd_count)
        query_local <- query_local * (max_count / sd_count)
        genes <- intersect(rownames(query_local), rownames(seed_local))
        print(i)
        counts <- cbind(seed_local[genes, ], query_local[genes, ])
        if (expand_genes) {
            diff_d1 <- setdiff(rownames(query_local), rownames(seed_local))
            diff_d2 <- setdiff(rownames(seed_local), rownames(query_local))
            expanded_counts <- matrix(0, ncol = ncol(counts),
                nrow = length(c(diff_d1, diff_d2)))
            colnames(expanded_counts) <- colnames(counts)
            rownames(expanded_counts) <- c(diff_d1, diff_d2)
            expanded_counts[diff_d1, colnames(seed_local)] <- 
                seed_local[diff_d1, ]
            expanded_counts[diff_d2, colnames(query_local)] <- 
                query_local[diff_d2, ]
            counts <- rbind(counts, expanded_counts)
        }
        counts <- counts[order(rownames(counts)), ]
        #counts <- log(x+1, base = 10)
        integrated_counts[[i]] <- counts
    }
    integrated_counts <- do.call("cbind", integrated_counts)
    colnames(integrated_counts) <- make.unique(colnames(integrated_counts),
        sep = "_")
    seed_coordinates <- seed@tiles %>%
        filter(origin == 1) %>%
        select(c("barcodes", "x", "y"))
    seed_coordinates$barcodes <- paste0("data1_", seed_coordinates$barcodes)
    query_coordinates <- aligned_graph[, c("barcodes", "x", "y")]
    query_coordinates$barcodes <- paste0("data1_",query_coordinates$barcodes)
    integrated_coordinates <- rbind(seed_coordinates, query_coordinates)
    integrated_coordinates$barcodes <- make.unique(
        integrated_coordinates$barcodes,
        sep = "_")
    vesalius_assay <- build_vesalius_assay(coordinates = integrated_coordinates,
        counts = integrated_counts,
        assay = "integrated",
        verbose = FALSE)
    return(vesalius_assay)

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
