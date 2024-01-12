###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################## COUNT INTEGRATION ##############################
#-----------------------------------------------------------------------------#
integrate_map <- function(aligned_graph,
    query,
    seed,
    return_cost = FALSE,
    merge = FALSE,
    verbose = TRUE) {
    prob <- list(aligned_graph$prob)
    names(prob) <- "mapping_scores"
    cost <- list(aligned_graph$cost)
    names(cost) <- "score_matrices"
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
    query@meta <- c(query@meta, cost)
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
    return(query)
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



