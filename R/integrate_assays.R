###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################ HORIZONTAL INTEGRATION ###########################
#-----------------------------------------------------------------------------#


#' Aling and integrate spatial assay from the same modality using super pixels
#' @param seed_assay vesalius_assay object - data to be mapped to
#' @param query_assay vesalius_assay objecy - data to map
#' @param seed_trial name of embedding to use for super pixel generation
#' @param query_trial name of embedding to use for super pixel generation
#' @param scoring_method method used to score the similarity between super
#' pixels (pearson, spearman, kendall, coherence, index)
#' @param index_selection character (random, bubble, hex)
#' how should initial super pixel locations be selected. See detials
#' @param mapping character - method used for mapping query data onto see data
#' (mesh, spix, morph)
#' @param dimensions integer vector - which latent space dimensions should be
#' used for super pixel generation
#' @param scaling numeric ]0,1] describing image scale to consider
#' during super pixel selection
#' @param compactness numeric ]0,Inf] - importance of the spatial component
#' in relation to color similarity. See details
#' @param grid integer [2, Inf[]
#' How many grid elements should be use to create inital hex grid (see details)
#' @param iter integer [1, Inf] - Number of iteration during graph matching
#' phase.
#' @param n_landmarks [1, N_spatial_index] - number of super pixels to use as
#' landmarks for mesh adjustement.

#' @param threshold numeric [0,1[ - similarity score threshold. Only super pixel
#' that score above this threshold will be used for graph matching
#' @param signal character (features, counts, embeddings, "custom") - What should 
#' be used as cell signal for super pixel scoring. Seed details 
#' @param verbose logical - should I be a noisy boy?
#' @details Hex grid is expanded to include all points and then reduced to exclude
#' exmpty triangles.
#' @export
#' @importFrom Morpho computeTransform applyTransform

integrate_horizontally <- function(seed_assay,
    query_assay,
    scoring_method = "pearson",
    index_selection = "bubble",
    mapping = "mesh",
    trial = "last",
    dimensions = seq(1, 30),
    scaling = 0.2,
    compactness = 5,
    grid = 100,
    n_landmarks = 50,
    threshold = 0.7,
    depth = 1,
    signal = "variable_features",
    use_norm = "raw",
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # First let's get singal
    #-------------------------------------------------------------------------#
    if (!is(query_assay,"list")) {
        query_assay <- list(query_assay)
    }
    signal <- get_features(seed_assay = seed_assay,
        query_assay = query_assay,
        signal = signal,
        use_norm = use_norm,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # Compute landmarks
    #-------------------------------------------------------------------------#
    # landmarks <- get_landmarks(seed = seed_assay,
    #     seed_signal = signal$seed,
    #     query = query_assay,
    #     query_signal = signal$query,
    #     trial = trial,
    #     scoring_method = scoring_method,
    #     index_selection = index_selection,
    #     dimensions = dimensions,
    #     scaling = scaling,
    #     compactness = compactness,
    #     n_landmarks = n_landmarks,
    #     depth = depth,
    #     threshold = threshold,
    #     verbose)
    spix_mnn <- get_spix_MNN(seed = seed_assay,
        seed_signal = signal$seed,
        query = query_assay,
        query_signal = signal$query,
        trial = trial,
        scoring_method = scoring_method,
        index_selection = index_selection,
        dimensions = dimensions,
        scaling = scaling,
        compactness = compactness,
        n_landmarks = n_landmarks,
        depth = depth,
        threshold = threshold,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # Mapping points from landmarks
    #-------------------------------------------------------------------------#
    # aligned <- switch(EXPR = mapping,
    #     "mesh" = mesh_mapping(landmarks,
    #         signal,
    #         grid,
    #         depth,
    #         threshold,
    #         verbose),
    #     "morph" = morph_mapping(landmarks),
    #     "spix" = spix_mapping(landmarks),
    #     "spix_nn" = spix_nn_mapping(landmarks,
    #         signal))
    # message_switch("mesh", verbose)
    # seed_mesh <- generate_mesh(seed_assay,
    #     n_centers = grid)
    # query_mesh <- generate_mesh(query_assay,
    #     n_centers = grid)
    # #browser()
    # #-------------------------------------------------------------------------#
    # # Adjust grid to match landmark spix 
    # #-------------------------------------------------------------------------#
    # seed_mesh <- adjust_mesh(seed_mesh, seed_landmarks)
    # query_mesh <- adjust_mesh(query_mesh, query_landmarks)
    
    
   
    #-------------------------------------------------------------------------#
    # Get estimated super pixel centers
    # we also generate a graph and score this graph. 
    # score the correlation between each vertex in seed/query graph
    #-------------------------------------------------------------------------#
    
    # message_switch("slic_graph", verbose = verbose, data = "seed")
    # seed_graph <- generate_slic_graph(seed_mesh$segments)

    # message_switch("slic_graph", verbose = verbose, data = "query")
    # query_graph <- generate_slic_graph(query_mesh$segments)
    # #-------------------------------------------------------------------------#
    # # Now we can compute the same thing but between each graph
    # # For now - we check the number of centeres in case of mismatch
    # # for now I am getting mismatch even with random sampling which 
    # # does not make any sense but to check and fix
    # #-------------------------------------------------------------------------#
    # q_centers <- length(unique(query_graph$from))
    # s_centers <- length(unique(seed_graph$from))
    # integrated_graph <- data.frame(
    #     "to" = rep(unique(seed_graph$from), each = q_centers),
    #     "from" = rep(unique(query_graph$from), times = s_centers))


    # spix_score <- score_graph(integrated_graph,
    #     signal = list(seed_signal, query_signal),
    #     verbose = verbose)
    # best_vertex <-get_best_vertex(spix_score)
    # land <-  get_best_vertex(spix_score) %>% filter(score > 0.5)
    # land_1 <- get_super_pixel_centers(seed_landmarks$segments)[land$to, ]
    # land_2 <- get_super_pixel_centers(query_landmarks$segments)[land$from, ]
    # trasfo <- computeTransform(as.matrix(land_1[,c("x","y")]),
    #     as.matrix(land_2[,c("x","y")]), reflection = TRUE)
    # transformed <- applyTransform(as.matrix(query_landmarks$segments[,c("x","y")]),trasfo)
    # aligned_graph <- query_landmarks$segments
    # aligned_graph[,c("x","y")] <- transformed
    # aligned_graph$norm_with <- aligned_graph$Segment
    #browser()
    # matched_graph <- match_graph(seed_graph = seed_graph,
    #     seed_trial = seed_mesh$segments,
    #     query_graph = query_graph,
    #     query_trial = query_mesh$segments,
    #     score = spix_score,
    #     depth = depth,
    #     verbose = verbose)
    # aligned_graph <- align_graph(matched_graph,
    #     seed_mesh$segments,
    #     seed_mesh$mesh,
    #     query_mesh$segments,
    #     query_mesh$mesh,
    #     verbose = verbose)
    integrated <- integrate_graph(spix_mnn,
        query_assay,
        verbose = verbose)
    # simple_bar(verbose)
    return(integrated)
}

get_features <- function(seed_assay,
    query_assay,
    signal,
    use_norm = "raw",
    verbose = TRUE) {
    message_switch("signal", verbose = verbose)
    #-------------------------------------------------------------------------#
    # First we check if we are using embedding values 
    #-------------------------------------------------------------------------#
    if (any(grepl(pattern = "embeddings", x = signal))) {
        seed_signal <- seed@active
        query_signal <- lapply(query, slot, "active")
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
        seed_counts <- get_counts(seed_assay, type = use_norm)[features, ]
        query_counts <- lapply(query_assay, function(assays, features, type) {
            return(get_counts(assays, type)[features, ])
        }, features = features, type = use_norm)
    }
    return(list("seed" = seed_counts,
        "query" = query_counts,
        "features" = features))

}

get_landmarks <- function(seed,
    seed_signal,
    query,
    query_signal,
    trial = "last",
    scoring_method = "pearson",
    index_selection = "random",
    dimensions = seq(1,30),
    scaling = 0.2,
    compactness = 5,
    n_landmarks = 50,
    depth = 1,
    threshold = 0.5,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # first we compute super pixels from embedding values
    # Note this could be updated to parallel 
    #-------------------------------------------------------------------------#
    message_switch("landmarks", verbose, assay = get_assay_names(seed))
    seed_spix <- slic_segmentation(seed,
        dimensions = dimensions,
        col_resolution = n_landmarks,
        embedding = trial,
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    message_switch("landmarks", verbose, assay = "Query List")
    query_spix <- lapply(query, slic_segmentation,
        dimensions = dimensions,
        col_resolution = n_landmarks,
        embedding = trial,
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    #-------------------------------------------------------------------------#
    # Next we create a scoring table for each spix
    #-------------------------------------------------------------------------#
    message_switch("optimal_land", verbose)
    matched_graph <- vector("list", length(query_spix))
    for (i in seq_along(query_spix)) {
        seed_nodes <- sort(unique(seed_spix$segments$Segment))
        query_nodes <- sort(unique(query_spix[[i]]$segments$Segment))
        score_table <- data.frame(
            "from" = rep(query_nodes, each = length(seed_nodes)),
            "to" = rep(seed_nodes, times = length(query_nodes)))
        seed_local <- compress_signal(seed_signal,
            seed_spix$segments)
        query_local <- compress_signal(query_signal[[i]],
            query_spix[[i]]$segments)
        spix_score <- score_graph(score_table,
            signal = list(seed_local, query_local),
            assay = paste("Query List", i),
            verbose = verbose)
        seed_graph <- generate_slic_graph(seed_spix$segments)
        query_graph <- generate_slic_graph(query_spix[[i]]$segments)
        matched_graph[[i]] <- match_graph(seed_graph = seed_graph,
            seed_trial = seed_spix$segments,
            query_graph = query_graph,
            query_trial = query_spix[[i]]$segments,
            score = spix_score,
            depth = depth,
            threshold = threshold,
            verbose = verbose)
    }
    return(list("seed_spix" = seed_spix,
        "query_spix" = query_spix,
        "matched_graph" = matched_graph))

}

get_spix_MNN <- function(seed,
    seed_signal,
    query,
    query_signal,
    trial = "last",
    scoring_method = "pearson",
    index_selection = "random",
    dimensions = seq(1,30),
    scaling = 0.2,
    compactness = 5,
    n_landmarks = 50,
    depth = 1,
    threshold = 0.5,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # first we compute super pixels from embedding values
    # Note this could be updated to parallel 
    #-------------------------------------------------------------------------#
    message_switch("landmarks", verbose, assay = get_assay_names(seed))
    seed_spix <- slic_segmentation(seed,
        dimensions = dimensions,
        col_resolution = n_landmarks,
        embedding = trial,
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    seed_spix <- merge_spix(seed_spix)
    message_switch("landmarks", verbose, assay = "Query List")
    query_spix <- lapply(query, slic_segmentation,
        dimensions = dimensions,
        col_resolution = n_landmarks,
        embedding = trial,
        index_selection = index_selection,
        compactness = compactness,
        scaling = scaling,
        verbose = FALSE)
    query_spix <- lapply(query_spix, merge_spix)
    #-------------------------------------------------------------------------#
    # Next we create a scoring table for each spix
    #-------------------------------------------------------------------------#
    message_switch("spix_MNN", verbose)
    aligned_indices <- vector("list", length(query_spix))
    for (i in seq_along(query_spix)) {
        #---------------------------------------------------------------------#
        # Preparing matrices to generate a cost matrix for hugarian solver
        #---------------------------------------------------------------------#
        seed_nodes <- sort(unique(seed_spix$segments$Segment))
        query_nodes <- sort(unique(query_spix[[i]]$segments$Segment))
        score_table <- data.frame(
            "from" = rep(query_nodes, each = length(seed_nodes)),
            "to" = rep(seed_nodes, times = length(query_nodes)))
        #---------------------------------------------------------------------#
        # First we create a transcriptional  matrix
        #---------------------------------------------------------------------#
        seed_local <- t(scale(t(as.matrix(seed_signal))))
        query_local <- t(scale(t(as.matrix(query_signal[[i]]))))
        feature_dist_mat <- RANN::nn2(t(seed_local),
            t(query_local),
            k = min(ncol(seed_local), ncol(query_local)))
        ordered_mat <- feature_dist_mat$nn.idx
        rownames(ordered_mat) <- colnames(query_local)
        colnames(ordered_mat) <- colnames(seed_local)
        feature_dist_mat <- 
            (feature_dist_mat$nn.dists - min(feature_dist_mat$nn.dists)) /
            (max(feature_dist_mat$nn.dists) - min(feature_dist_mat$nn.dists))
        rownames(feature_dist_mat) <- colnames(query_local)
        colnames(feature_dist_mat) <- colnames(seed_local)
        #---------------------------------------------------------------------#
        # Next we create a relative spatial distance matrix
        #---------------------------------------------------------------------#
        seed_dist <- RANN::nn2(seed_spix$segments[, c("x","y")],
            k = nrow(seed_spix$segments))
        query_dist <- RANN::nn2(query_spix[[i]]$segments[, c("x","y")],
            k = nrow(query_spix[[i]]$segments))
        
        #---------------------------------------------------------------------#
        # Now we create a spix neighbohood matrix - generates a score 
        # for each super pixel 
        #---------------------------------------------------------------------#
        spix_score <- mutual_spix_score(ordered_mat,
            graph = score_table,
            seed = seed_spix$segments,
            query = query_spix[[i]]$segments,
            assay = paste("Query list" , i),
            verbose = TRUE)
        #---------------------------------------------------------------------#
        # Create a cost matrix
        #---------------------------------------------------------------------#
        cost_mat <- feature_dist_mat + spix_score
        matched_indices <- match_index(cost_mat)
        #browser()
        aligned_indices[[i]] <- align_index(matched_indices,
            seed_spix$segments,
            query_spix[[i]]$segments,
            verbose)
       
    }
    return(aligned_indices)
}

spix_mapping <- function(landmarks) {
    seed <- landmarks$seed_spix
    query <- landmarks$query_spix
    matched_graph <- landmarks$matched_graph
    for (i in seq_along(query)){
        seed_centers <- get_super_pixel_centers(seed$segments)
        query_centers <- get_super_pixel_centers(query[[i]]$segments)
        anchors <- matched_graph[[i]] %>%
            filter(anchors == 1)
        query_local <- query[[i]]$segments %>%
            filter(Segment %in% anchors$from)
        query_local$norm_with <- NA
        for (j in seq_len(nrow(anchors))) {
            seed_point <- seed_centers[seed_centers$center ==
                anchors$to[j], ]
            query_point <- query_centers[query_centers$center ==
                anchors$from[j], ]
            angle <- atan2(seed_point$y - query_point$y,
                seed_point$x - query_point$x)
            distance <- matrix(c(seed_point$x, query_point$x,
                seed_point$y, query_point$y), ncol = 2)
            distance <- as.numeric(dist(distance))
            query_local$x[query_local$Segment == anchors$from[j]] <-
                query_local$x[query_local$Segment == anchors$from[j]] +
                (distance * cos(angle))
            query_local$y[query_local$Segment == anchors$from[j]] <-
                query_local$y[query_local$Segment == anchors$from[j]] +
                (distance * sin(angle))
            query_local$norm_with[query_local$Segment == anchors$from[j]] <-
                anchors$to[j]
        }
        query[[i]] <- query_local
    }
    return(query)
}

mesh_mapping <- function(landmarks,
    signal,
    grid,
    depth = 1,
    threshold = 0.7,
    verbose = TRUE) {
    seed <- landmarks$seed_spix
    query <- landmarks$query_spix
    seed_signal <- signal$seed
    query_signal <- signal$query
    seed_mesh <- generate_mesh(seed$segments, grid) %>%
            adjust_mesh()
    aligned_mesh <- vector("list", length(query))
    for (i in seq_along(query)) {
        query_mesh <- generate_mesh(query[[i]]$segments, grid) %>%
            adjust_mesh()
        seed_nodes <- sort(unique(seed_mesh$coord$Segment))
        query_nodes <- sort(unique(query_mesh$coord$Segment))
        score_table <- data.frame(
            "from" = rep(query_nodes, each = length(seed_nodes)),
            "to" = rep(seed_nodes, times = length(query_nodes)))
        seed_local <- compress_signal(seed_signal,
            seed_mesh$coord)
        query_local <- compress_signal(query_signal[[i]],
            query_mesh$coord)
        spix_score <- score_graph(score_table,
            signal = list(seed_local, query_local),
            assay = paste("Query Mesh List", i),
            verbose = verbose)
        seed_graph <- generate_slic_graph(seed_mesh$coord)
        query_graph <- generate_slic_graph(query_mesh$coord)
        matched_graph <- match_graph(seed_graph = seed_graph,
            seed_trial = seed_mesh$coord,
            query_graph = query_graph,
            query_trial = query_mesh$coord,
            score = spix_score,
            depth = depth,
            threshold = threshold,
            verbose = verbose)
        aligned_mesh[[i]] <- align_graph(matched_graph,
            seed = seed_mesh$coord,
            seed_graph = seed_mesh$mesh,
            query = query_mesh$coord,
            query_graph = query_mesh$mesh,
            verbose)
        
    }
    return(aligned_mesh)
}
#' @importFrom Morpho computeTransform applyTransform
morph_mapping <- function(landmarks) {
    seed <- landmarks$seed_spix
    query <- landmarks$query_spix
    landmarks_points <- landmarks$matched_graph
    for (i in seq_along(query)) {
        land <- landmarks_points[[i]]
        seed_spix <- get_super_pixel_centers(seed$segments)[land$to, ]
        query_spix <- get_super_pixel_centers(query[[i]]$segments)[land$from, ]
        transform <- Morpho::computeTransform(
            as.matrix(seed_spix[, c("x", "y")]),
            as.matrix(query_spix[, c("x", "y")]),
            reflection = TRUE,
            type = "rigid",
            weights = land$score)
        transform <- Morpho::applyTransform(
            as.matrix(query[[i]]$segments[, c("x", "y")]),
            transform)
        aligned_graph <- query[[i]]$segments
        aligned_graph[, c("x","y")] <- transform
        # To UPDATE - this is not how it should be done 
        aligned_graph$norm_with <- aligned_graph$Segment
        query[[i]] <- aligned_graph
    }
    return(query)
}


mutual_spix_score <- function(ordered_mat,
    graph,
    seed,
    query,
    assay,
    verbose = TRUE) {
    graph <- data.frame(graph,
        "score" = rep(0, nrow(graph)))
    for (i in seq_len(nrow(graph))) {
        dyn_message_switch("score_graph", verbose,
            prog = round(i / nrow(graph), 4) * 100,
            assay = assay)
        seed_bar <- seed$barcodes[seed$Segment == graph$to[i]]
        query_bar <- query$barcodes[query$Segment == graph$from[i]]
        if (length(seed_bar) == 0 || length(query_bar) == 0) {
            #browser()
            graph$score <- NA
        } else {
            #browser()
            local_ord_mat <- ordered_mat[rownames(ordered_mat) %in% seed_bar,
            colnames(ordered_mat) %in% query_bar]
            dim(local_ord_mat) <- c(length(seed_bar), length(query_bar))
            ordered_mat[rownames(ordered_mat) %in% seed_bar,
                colnames(ordered_mat) %in% query_bar]   <- sum(local_ord_mat) /
                    (nrow(local_ord_mat) * ncol(local_ord_mat))
        }
    }
    if (verbose){cat("\n")}
    ordered_mat <- (ordered_mat - min(ordered_mat)) /
        (max(ordered_mat) - min(ordered_mat))
    return(ordered_mat)
}

#' importFrom future.apply future_lapply
score_graph <- function(graph,
    signal,
    scoring_method = "spearman",
    assay,
    verbose = TRUE){
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
            prog = round(i / nrow(graph), 4) * 100,
            assay = assay)
        c1 <- seed_signal[[graph$to[i]]]
        c2 <- query_signal[[graph$from[i]]]
        if (length(c1) == 0 || length(c2) == 0){
            graph$score <- 0
        } else {
            graph$score[i] <- switch(EXPR = scoring_method,
                "pearson" = cor(c1, c2, method = scoring_method),
                "spearman" = cor(c1, c2, method = scoring_method),
                "kendall" = cor(c1, c2, method = scoring_method))
        }
    }
    if (verbose){cat("\n")}
    return(graph)
}


merge_spix <- function(spix) {
    spixs <- spix$segments
    n_bar <- table(spixs$Segment)
    low_n <- names(n_bar)[n_bar < 3]
    if (length(low_n) > 0) {
        knn <- RANN::nn2(spixs[!spixs$Segment %in% low_n, c("x", "y")],
            spixs[spixs$Segment %in% low_n, c("x", "y")],
            k = 1)
        spixs$Segment[spixs$Segment %in% low_n] <- 
            spixs$Segment[which(!spixs$Segment %in% low_n)[knn$nn.idx[,1]]]
        #spixs$Segment <- spixs$Segment[]
        spix$segments <- spixs
    }
    return(spix)
}

match_index <- function(cost_matrix) {
    mapping <- RcppHungarian::HungarianSolver(cost_matrix)$pairs
    mapping <- as.data.frame(mapping)
    colnames(mapping) <- c("to", "from")
    scores <- mapply(function(i, j, cost) {
        return(cost[i, j])
    }, mapping$to, mapping$from, MoreArgs = list(cost_matrix))
    mapping$score <- scores
    mapping$to <- rownames(cost_matrix)[mapping$to]
    mapping$from <- colnames(cost_matrix)[mapping$from]
    return(mapping)
}


#' @importFrom RcppHungarian HungarianSolver
match_graph <- function(seed_graph,
    seed_trial,
    query_graph,
    query_trial,
    score,
    depth = 1,
    threshold = 0.5,
    verbose = verbose) {
    #-------------------------------------------------------------------------#
    # Initialize optimisation 
    #-------------------------------------------------------------------------#
    query_paths <- graph_path_length(query_graph)
    seed_path <- graph_path_length(seed_graph)
    seed_trial <- get_super_pixel_centers(seed_trial)
    seed_trial$x <- min_max(seed_trial$x)
    seed_trial$y <- min_max(seed_trial$y)
    #-------------------------------------------------------------------------#
    # Create a cost matrix fro optimizer
    #-------------------------------------------------------------------------#
    cost_matrix <- lapply(seq(1, nrow(score)), build_cost_matrix,
        score = score,
        query_paths = query_paths,
        seed_path = seed_path,
        seed_trial = seed_trial,
        depth = depth,
        verbose = verbose)
    if (verbose){cat("\n")}
    cost <- sapply(cost_matrix, "[[", 1)
    #-------------------------------------------------------------------------#
    # Optimize using a hungarian solver
    #-------------------------------------------------------------------------#
    message_switch("hungarian", verbose)
    cost_matrix <- matrix(cost,
        ncol = length(unique(score$from)),
        nrow = length(unique(score$to)),
        byrow = TRUE)
    mapping <- RcppHungarian::HungarianSolver(cost_matrix)$pairs
    mapping <- as.data.frame(mapping)
    colnames(mapping) <- c("to", "from")
    #-------------------------------------------------------------------------#
    # define anchors
    #-------------------------------------------------------------------------#
    scores <- mapply(function(i, j, cost){
        return(cost[i, j])
    }, mapping$to, mapping$from, MoreArgs=list(cost_matrix))
    scores <- 1 - ((max(cost) - scores) / (max(cost)))
    mapping$score <- scores
    mapping$anchors <- 0
    mapping$anchors[scores > threshold] <- 1
    message_switch("landmarks_found", verbose, anchors = sum(mapping$anchors))
    return(mapping)
}




build_cost_matrix <- function(i,
    score,
    query_paths,
    seed_path,
    seed_trial,
    depth,
    verbose) {
    dyn_message_switch("cost_mat", verbose,
            prog = round(i / nrow(score), 4) * 100)
    #---------------------------------------------------------------------#
    # First we get the score of for the center spix
    # since we want to minimize the cost we do 1 - score
    #---------------------------------------------------------------------#
    feature_score <- 1 - score$score[i]
    #---------------------------------------------------------------------#
    # next we want to minimize the score difference betwee neighborhoods
    # what is the best score for the query neighbor hood if mapped to 
    # the seed 
    #---------------------------------------------------------------------#
    query_niche <- query_paths[, as.character(score$from[i])]
    query_niche <- names(query_niche)[query_niche <= depth &
            query_niche > 0]
    niche <- score[score$from %in% query_niche, ] %>%
                get_best_vertex()
    niche_score <- 1 - niche$score
    #---------------------------------------------------------------------#
    # For the best matches we want to minimize the distance between 
    # query neighborhood and best matches
    # we already know that we have the best matches for the query niche
    #---------------------------------------------------------------------#
    mapped_to <- seed_trial[seed_trial$center == score$to[i], ]
    seed_niche <- seed_trial[seed_trial$center %in% niche$to, ]
        
    dist <- RANN::nn2(mapped_to[, c("x", "y")],
        query = seed_niche[, c("x", "y")])$nn.dist
    out <- list("cost" = sum(c(feature_score, niche_score, mean(dist))),
        "score" = feature_score,
        "niche" = niche_score,
        "distance" = dist)
    return(out)
}


generate_mesh <- function(coord,
    grid) {
    #-------------------------------------------------------------------------#
    # get coordinates
    #-------------------------------------------------------------------------#
    coord$mesh <- 0
    mesh <- hex_grid(coord, grid, return_index = FALSE)
    mesh <- create_triangle_mesh(mesh[, 1], mesh[, 2])
    names(mesh) <- seq_along(mesh)
    for (i in seq_along(mesh)){
        in_triangle <- sp::point.in.polygon(coord$x, coord$y,
            mesh[[i]]$x, mesh[[i]]$y) != 0
        if (sum(in_triangle) == 0) {
            mesh[[i]] <- NA
            next()
        }
        coord$mesh[in_triangle] <- i
    }
    mesh <- mesh[!is.na(mesh)]
    coord$mesh <- seq_along(mesh)[
        match(as.character(coord$Segment), names(mesh))]
    names(mesh) <- seq_along(mesh)
    return(list("coord" = coord, "mesh" = mesh))
}

adjust_mesh <- function(mesh) {
    #-------------------------------------------------------------------------#
    # first get points to from mesh to shift
    # and shift to spix center
    #-------------------------------------------------------------------------#
    coord <- mesh$coord
    mesh <- mesh$mesh
    landmark_points <- get_super_pixel_centers(coord)
    buffer_mesh <- do.call("rbind", mesh)
    closest <- RANN::nn2(data = landmark_points[, c("x", "y")],
        query = buffer_mesh[, c("x","y")], k = 6)
    xlim <- sort(unique(buffer_mesh$x))
    xlim <- c(head(xlim,2), tail(xlim, 2))
    ylim <- range(buffer_mesh$y)
    non_edge <- TRUE
    for (i in seq_along(unique(landmark_points$center))) {
        nn <- 1
        node <- landmark_points$center[i]
        while (non_edge) {
            loc <- which(closest$nn.idx[, nn] == node)
            if (length(loc) == 0) {
                break()
            }
            min_loc <- loc[which(closest$nn.dist[loc, nn] ==
                min(closest$nn.dist[loc, nn]))]
            local_coord <- buffer_mesh[min_loc, ]
            if (any(local_coord$x %in% xlim) || any(local_coord$y %in% ylim)) {
                nn <- nn + 1
            } else {
                non_edge <- FALSE
            }
            if (nn > 6) {
                non_edge <- FALSE
                min_loc <- NULL
            }
        }
        non_edge <- TRUE
        buffer_mesh[min_loc, ] <- landmark_points[node, c("x", "y")]
    }
    #-------------------------------------------------------------------------#
    # Rebuild mesh list 
    #-------------------------------------------------------------------------#
    updated_mesh <- lapply(seq(1, nrow(buffer_mesh), by = 3),
        function(i, mesh) {
            return(mesh[i:(i + 2), ])
        }, buffer_mesh)
    names(updated_mesh) <- seq_along(updated_mesh)
    #-------------------------------------------------------------------------#
    # Add new mesh annotations to segment for scoring
    #-------------------------------------------------------------------------#
    new_seg <- 1
    for (i in seq_along(updated_mesh)){
        in_triangle <- sp::point.in.polygon(coord$x, coord$y,
            updated_mesh[[i]]$x, updated_mesh[[i]]$y) != 0
        if (sum(in_triangle) == 0) {
            updated_mesh[[i]] <- NA
            next()
        }
        coord$Segment[in_triangle] <- new_seg
        new_seg <- new_seg + 1
    }
    updated_mesh <- updated_mesh[!is.na(updated_mesh)]
    coord$Segment <- seq_along(updated_mesh)[
        match(as.character(coord$Segment), names(updated_mesh))]
    names(updated_mesh) <- seq_along(updated_mesh)
    return(list("coord" = coord, "mesh" = updated_mesh))

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
        tri <- voronoi %>% filter(ind2 == idx)
        tri <- c(tri$ind1, idx)
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



#' 
#' @importFrom geometry cart2bary bary2cart
align_graph <- function(matched_graph,
    seed,
    seed_graph,
    query,
    query_graph,
    verbose = TRUE) {
    query$norm_with <- 0
    for (i in seq_len(nrow(matched_graph))) {
        pts <- query[query$Segment == matched_graph$from[i], c("x", "y")]
        tile <- query_graph[[as.character(matched_graph$from[i])]]
        bary <- geometry::cart2bary(X = as.matrix(tile),
            P = as.matrix(pts))
        new_tile <- seed_graph[[as.character(matched_graph$to[i])]]
        new_coord <- geometry::bary2cart(as.matrix(new_tile), bary)
        query$x[query$Segment == matched_graph$from[i]] <- new_coord[, 1]
        query$y[query$Segment == matched_graph$from[i]] <- new_coord[, 2]
        query$norm_with[query$Segment == matched_graph$from[i]] <- 
            as.character(matched_graph$to[i])

    }
    
    #-------------------------------------------------------------------------#
    # get closest seed point for all un assigned points
    #-------------------------------------------------------------------------#
    return(query)
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
    return(query)
}


# create_triangle_mesh <- function(x, y) {
#     voronoi <- deldir::deldir(x = as.numeric(x),
#         y = as.numeric(y))
#     anchors <- sort(unique(c(voronoi$delsgs$ind1, voronoi$delsgs$ind2)))
#     graph <- lapply(anchors, function(idx, voronoi){
#         triangle <- voronoi$delsgs %>% 
#             filter(ind1 == idx | ind2 == idx)
#         loc <- which(triangle[, c("ind1", "ind2")] != idx, arr.ind = TRUE)
#         triangle <- apply(loc, 1, function(loc, tri){
#             if (loc[2] == 1){
#                 tmp <- tri[loc[1], c("x1", "y1")]
#             } else {
#                 tmp <- tri[loc[1], c("x2", "y2")]
#             }
#             colnames(tmp) <- c("x", "y")
#             return(tmp)
#         }, triangle) %>%
#         do.call("rbind", .)
#         center <- voronoi$summary[idx, c("x", "y")]
#         triangle <- convexify(triangle$x,triangle$y, center$x, center$y)
#         all_tiles <- vector("list", nrow(triangle))
#         for (i in seq_len(nrow(triangle))){
#             if (i < nrow(triangle)){
#                 all_tiles[[i]] <- data.frame(rbind(center,
#                     triangle[i, ],
#                     triangle[i + 1, ]),
#                     "anchor" = rep(idx, 3))
#             } else {
#                 all_tiles[[i]] <- data.frame(rbind(center,
#                     triangle[i, ],
#                     triangle[1, ]),
#                     "anchor" = rep(idx, 3))
#             }
#         }
#         return(all_tiles)
#     }, voronoi = voronoi)
#     graph <- unlist(graph, recursive = FALSE)
#     return(graph)
# }

create_triangle_mesh <- function(x, y) {
    voronoi <- deldir::deldir(x = as.numeric(x),
        y = as.numeric(y))
    triangles <- seq_len(nrow(voronoi$delsgs))
    graph <- lapply(triangles, function(idx, voronoi){
        p1 <- voronoi$delsgs$ind1[idx]
        p2 <- voronoi$delsgs$ind2[idx]
        sub_1 <- voronoi$delsgs %>% filter(ind2 == p1)
        sub_2 <- voronoi$delsgs %>% filter(ind2 == p2)
        p3s <- intersect(sub_1$ind1, sub_2$ind1)
        triangles <- lapply(p3s, function(p3, p1, p2, coord){
            tmp <- rbind(coord[p1, c("x", "y")],
                coord[p2, c("x", "y")],
                coord[p3, c("x", "y")])
        },p1, p2, voronoi$summary)
        return(triangles)
    }, voronoi = voronoi)
    graph <- graph[sapply(graph, length) > 0]
    return(unlist(graph, recursive = FALSE))
}

get_triangular_idx <- function(nn) {
    coord_loc <- nn$nn.idx[, 1]
    i <- 1
    while (length(unique(coord_loc)) != 3) {
        coord_loc[duplicated(coord_loc)] <- nn$nn.idx[duplicated(coord_loc), i]
        i <- i + 1
    }
    return(coord_loc)
}





get_best_vertex  <- function(score, rank = 1, decreasing = TRUE) {
    from <- split(score, score$from)
    best_match <- lapply(from, function(f) {
        return(f[order(f$score, decreasing = decreasing)[rank], ])
    }) %>% do.call("rbind", .)
    return(best_match)
}

#'@importFrom igraph graph_from_data_frame distances
#'@importFrom igraph E
graph_path_length <- function(graph) {
    gr <- igraph::graph_from_data_frame(graph, directed = FALSE)
    path_length <- igraph::distances(gr)
    return(path_length)
}

#' @importFrom dplyr select
integrate_graph <- function(aligned_graph,
    query,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # rebuild vesalius
    #-------------------------------------------------------------------------#
    browser()
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

# #' importFrom gsignal mscohere
# spectral_coherence <- function(seed_path, query_path) {
#     seed_path <- seed_path[names(seed_path) != "isolated"]
#     query_path <- query_path[names(query_path) != "isolated"]
#     coh <- matrix(0,
#         ncol = length(seed_path),
#         nrow = length(query_path))
#     for (seed in seq_along(seed_path)){
#         for (query in seq_along(query_path)){
#             #-----------------------------------------------------------------#
#             # Get partial max and subset paths
#             #-----------------------------------------------------------------#
#             seed_sub <- seed_path[[seed]]
#             query_sub <- query_path[[query]]
#             signal <- cbind(seed_sub, query_sub)
#             #-----------------------------------------------------------------#
#             # coherence 
#             #-----------------------------------------------------------------#
#             coh[query, seed] <- max(gsignal::mscohere(signal)$coh)
#         }
#     }
#     colnames(coh) <- names(seed_path)
#     rownames(coh) <- names(query_path)
#     return(coh)
# }

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

