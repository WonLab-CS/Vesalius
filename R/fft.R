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
    use_slic = TRUE,
    dimensions = seq(1, 30),
    scaling = 0.2,
    compactness = 1,
    n_centers = 2000,
    threshold = 0.9,
    k = "auto",
    use_counts = FALSE,
    use_norm = "raw",
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # compute slic for both assays 
    #-------------------------------------------------------------------------#
    if (use_slic) {
        seed_trial <- segment_image(seed_assay,
            method = "slic",
            dimensions = dimensions,
            col_resolution = n_centers,
            compactness = compactness,
            scaling = scaling,
            verbose = FALSE)
        query_trial <- segment_image(query_assay,
            method = "slic",
            dimensions = dimensions,
            col_resolution = n_centers,
            compactness = compactness,
            scaling = scaling,
            verbose = FALSE)
    }
    #-------------------------------------------------------------------------#
    # get signal either as counts or as embedding values
    #-------------------------------------------------------------------------#
    if (use_counts) {
        seed_signal <- get_counts(seed_assay, type = use_norm)
        query_signal <- get_counts(query_assay, type = use_norm)
        seed_genes <- intersect(rownames(seed_counts), rownames(query_counts))
        seed_signal <- seed_counts[seed_genes, ]
        query_signal <- query_counts[seed_genes, ]
    } else {
        seed_signal <- t(get_embeddings(seed_assay))
        query_signal <- t(get_embeddings(query_assay))
    }
    #-------------------------------------------------------------------------#
    # Get estimated super pixel centers
    # we also generate a graph and score this graph. 
    # score the correlation between each vertex in seed/query graph
    #-------------------------------------------------------------------------#
    seed_centers <- check_segment_trial(seed_trial)
    seed_graph <- generate_slic_graph(seed_centers,
        signals = seed_signal,
        k = k,
        scoring_method = scoring_method)
    colnames(seed_graph) <- c("seed_1", "seed_2", "seed_score")
    query_centers <- check_segment_trial(query_trial)
    query_graph <- generate_slic_graph(query_centers,
        signal = query_signal,
        k = k,
        scoring_method = scoring_method)
    colnames(query_graph) <- c("query_1", "query_2", "query_score")
    #-------------------------------------------------------------------------#
    # Now we can compute the same thing but between each graph
    #-------------------------------------------------------------------------#
    
    spix_score <- score_graph(
        g1 = rep(unique(seed_graph$seed_1), times = n_centers),
        g2 = rep(unique(query_graph$query_1), each = n_centers),
        signal = list(seed_signal, query_signal),
        centers = list(seed_centers, query_centers))
    colnames(spix_score) <- c("seed_spix", "query_spix", "spix_score")
    
    simple_bar(verbose)
    return(list("seed_score" = seed_graph,
        "query_score" = query_graph,
        "spix_score" = spix_score))
}

territory_signal <- function(counts, territories) {
    territories <- split(territories, territories$trial)
    for (i in seq_along(territories)) {
        territories[[i]] <- as.vector(
            apply(counts[, territories[[i]]$barcodes],
            MARGIN = 1,
            mean))
    }
    return(territories)
}



# integrate_by_territory <- function(seed_assay,
#     query_assay,
#     seed_trial = "last",
#     query_trial = "last",
#     method = "coherence",
#     k = 5,
#     use_counts = TRUE,
#     use_norm = "raw",
#     verbose = TRUE) {
#     simple_bar(verbose)
#     #-------------------------------------------------------------------------#
#     # Get counts 
#     #-------------------------------------------------------------------------#
#     if (use_counts) {
#         seed_counts <- get_counts(seed_assay, type = use_norm)
#         query_counts <- get_counts(query_assay, type = use_norm)
#         seed_genes <- intersect(rownames(seed_counts), rownames(query_counts))
#         seed_counts <- seed_counts[seed_genes, ]
#         query_counts <- query_counts[seed_genes, ]
#     } else {
#         seed_counts <- t(get_embeddings(seed_assay))
#         query_counts <- t(get_embeddings(query_assay))
#     }
    
#     #-------------------------------------------------------------------------#
#     # get territory information
#     #-------------------------------------------------------------------------#
#     seed_assay <- vesalius:::check_territory_trial(seed_assay, seed_trial)
#     query_assay <- vesalius:::check_territory_trial(query_assay, query_trial)
#     #-------------------------------------------------------------------------#
#     # Get variable features - place holder for now 
#     #-------------------------------------------------------------------------#
#     seed_signals <- territory_signal(seed_counts, seed_assay)
#     query_signals <- territory_signal(query_counts, query_assay)
#     #-------------------------------------------------------------------------#
#     # Comapring signals 
#     #-------------------------------------------------------------------------#
#     sim <- signal_similiarity(seed_signals,
#         query_signals,
#         domain = method)
    
#     seed_rank <- generate_territory_graph(seed_assay[seed_assay$trial !=
#         "isolated", ], k)
#     query_rank <- generate_territory_graph(query_assay[query_assay$trial !=
#         "isolated", ], k)
#     neighborhood_sim <- score_neighbor_graph(seed_rank, query_rank, sim)
#     simple_bar
#     return(list("sim" = sim, "n_sim" = neighborhood_sim))
# }

# territory_signal <- function(counts, territories) {
#     territories <- split(territories, territories$trial)
#     for (i in seq_along(territories)) {
#         territories[[i]] <- as.vector(
#             apply(counts[, territories[[i]]$barcodes],
#             MARGIN = 1,
#             mean))
#     }
#     return(territories)
# }




# #' retrieve the points contained in the edge of each territory
# #' @param trial name of territory trial that should be selected
# #' @param tiles vesalius tiles 
# #' @param start string - which point should be used as starting point
# #' @param verbose logical - should progress message be print to console
# #' @details Here we are using convex as start point. Essentially, we 
# #' order the coordinates based on their polar coordinates using the 
# #' median coordinate as the center point. 
# #' @returns a data frame containing edge of each territory.
# unpack_territory_path <- function(trial,
#     tiles,
#     start = "convex",
#     verbose = TRUE) {
#     #-------------------------------------------------------------------------#
#     # First we convert to pixset and detect territory edge 
#     #-------------------------------------------------------------------------#
#     trial_split <- vector("list", length(unique(trial$trial)))
#     names(trial_split) <- unique(trial$trial)
#     for (i in seq_along(unique(trial$trial))) {
#         territory <- unique(trial$trial)[i]
#         ter <- right_join(trial, tiles, by = "barcodes") %>%
#             filter(trial %in% territory) %>%
#             mutate(value = 1) %>%
#             select(c("barcodes", "x.y", "y.y", "value", "origin", "trial"))
#         colnames(ter) <- c("barcodes", "x", "y", "value", "origin", "trial")
#         edge <- extend_boundary(ter, 1) %>%
#             detect_edges() %>%
#             grow(1) %>%
#             as.data.frame()
#         edge <- inner_join(edge, ter, by = c("x", "y")) %>%
#             select("barcodes") %>% unique()
#         edge <- tiles %>% filter(barcodes %in% edge$barcodes & origin == 1)
#         trial_split[[i]] <- edge

#     }

#     #-------------------------------------------------------------------------#
#     # next we remove NULLs - this happens when no edge can be deteced 
#     #-------------------------------------------------------------------------#
#     nulls <- sapply(trial_split, nrow) == 0
#     trial_split <- trial_split[!nulls]
#     if (sum(nulls) > 0) {
#         message_switch("edge_detect", verbose,
#             nulls = paste(names(trial_split)[nulls]))
#     } else if (length(trial_split) == 0) {
#         stop("No edge can be detect in territories! Granularity too high.
#          Consider increasing smoothing and/or decreasing segmentationd depth")
#     }
#     #-------------------------------------------------------------------------#
#     # select starting point for path 
#     # ATM we convert edge to ordered shape using polar coordinates 
#     #-------------------------------------------------------------------------#
#     trial <- switch(EXPR = start,
#         "convex" = lapply(trial_split, function(trial) {
#             ord <- convexify(trial$x,
#                 trial$y,
#                 median(trial$x),
#                 median(trial$y),
#                 order = TRUE)
#             trial <- trial[ord$x, ]
#             return(trial)
#         }),
#         "connected" = lapply(trial_split, connected_points))
#     return(trial)
# }

# #' create path from neighboring points 
# #' @param trail data frame containign x and y path
# #' @return order x y coordinates
# #' @importFrom RANN nn2
# #' 
# connected_points <- function(trial) {
#     knn <- RANN::nn2(trial[, c("x","y")], k = nrow(trial))$nn.idx
#     ord <- rep(0, nrow(knn))
#     ord[1] <- 1
#     for (i in seq(2, nrow(knn))) {
#         if (i == 2) {
#             ord[2] <- knn[(i - 1), 2]
#         } else {
#             ord[i] <- knn[ord[(i - 1)], min(which(!knn[ord[(i - 1)], ] %in% ord))]
#         }
#     }
#     return(trial[ord, ])
# }




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


# #' territory graph 
# #' @param territories territory coordinates 
# #' @param k number of neigbors to create graph 
# #' @return nearest neighbors matrix 
# #' 
# #' 

# generate_territory_graph <- function(territories, k) {
#     #-------------------------------------------------------------------------#
#     # Computing distance matrix between all points and initialising 
#     # rank matrix for nearest neighbors 
#     #-------------------------------------------------------------------------#
#     coordinates <- as.matrix(dist(territories[, c("x", "y")]))
#     territory_list <- unique(territories$trial)
#     if (k > length(territory_list)){
#         warning("Value for k nearest neighbors is to high.
#         Returning all possible neighbors")
#         k <- length(territory_list)
#     }
#     rank_matrix <- matrix(0, nrow = length(territory_list), ncol = k + 1)
#     rownames(rank_matrix) <- territory_list
#     colnames(rank_matrix) <- seq(0, k)
#     rank_matrix[, 1] <- territory_list
#     #-------------------------------------------------------------------------#
#     # For each territory, check which coordinates from a different 
#     # territory are the closest. Assign nearest neighbor of rank 1.
#     # repeat the process until k nearest neighbor are found 
#     #-------------------------------------------------------------------------#
#     for (i in seq_along(territory_list)) {
#         #---------------------------------------------------------------------#
#         # intialise buffer storing barcodes that have already be
#         # assgined to nearest neighbor territory 
#         #---------------------------------------------------------------------#
#         init_territory <- which(colnames(coordinates) %in%
#                 territories$barcodes[territories$trial == territory_list[i]])
#         buffer <- colnames(coordinates)[init_territory]
#         for (j in seq(2, k + 1)) {
#             #-----------------------------------------------------------------#
#             # For clarity,  create variable for subsets
#             #-----------------------------------------------------------------#
#             not_in_buffer <- !rownames(coordinates) %in% buffer
#             tmp <- coordinates[not_in_buffer, init_territory]
#             tmp <- rownames(which(tmp == min(tmp), arr.ind = TRUE))[1L]
#             rank_matrix[i, j] <- territories$trial[territories$barcodes ==
#                 tmp]
#             buffer <- c(buffer, territories$barcodes[territories$trial ==
#                 rank_matrix[i, j]])
#         }
#     }
#     return(rank_matrix)
# }

# score_neighbor_graph <- function(seed_rank, query_rank, score_matrix) {
#     best_rank <- matrix(0, ncol = ncol(score_matrix),
#         nrow = nrow(score_matrix))
#     colnames(best_rank) <- colnames(score_matrix)
#     rownames(best_rank) <- rownames(score_matrix)
#     ranked_score <- apply(score_matrix, 2, order, decreasing = TRUE)
#     rownames(ranked_score) <- rownames(score_matrix)
#     for (i in seq_len(ncol(score_matrix))) {
#         for (j in seq_len(nrow(score_matrix))) {
#             seed <- seed_rank[rownames(seed_rank) ==
#                 colnames(score_matrix)[i], ]
#             query <- query_rank[rownames(query_rank) ==
#                 rownames(score_matrix)[j], ]
#             score <- best_neighbor_match(ranked_score[query, seed])
#             best_rank[j, i] <- score
#         }
#     }
#     return(best_rank)
# }

# best_neighbor_match <- function(score) {
#     vec <- as.vector(score)
#     ord <- ((order(vec, decreasing = TRUE) - 1) %/% ncol(score)) + 1
#     index <- c()
#     selected <- c()
#     for (i in seq_along(ord)){
#         if (ord[i] %in% unique(ord) && !ord[i] %in% selected) {
#             index <- c(index, i)
#             selected <- c(selected, ord[i])
#         }
#     }
#     index <- (sort(vec, decreasing = TRUE))[index]
#     return(mean(index))
# }

