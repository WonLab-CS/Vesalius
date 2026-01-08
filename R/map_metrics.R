###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################ Mapping Matrics ##################################
#-----------------------------------------------------------------------------#

#' Cluster query cells based on which reference cells they tend to mapped to
#' @param vesalius_assay vesalius_assay object after mapping a query onto a 
#' reference.
#' @param use_cost character vector describing which cost matrices should be
#' used to compare cells
#' @param cluster_method character string - which method should be used for
#' clustering (hclust, louvain, leiden)
#' @param trial character string defining which trial should be used for
#' clustering if any. If NULL, will search for "Cells".
#' @param group_identity character vector - which specific substes of trial 
#' should be used for clustering By default will use all labels present.
#' @param ref_cells character vector with reference cell barcodes
#' (by default will use all barcodes)
#' @param query_cells character with query cell barcodes
#' (by default will use all barcodes)
#' @param top_nn int - how many cells should be used to define clustering
#' similarity (see details)
#' @param h numeric - normalized height to use as hclust cutoff [0,1]
#' @param k int - number of cluster to obtain from hclust
#' @param nn int - number of nearest neighbors to use when creating graph
#' for community clustering algorithms
#' @param resolution numeric - clustering resolution to be parsed to 
#' community clustering algorithms
#' @param verbose logical - print output message
#' @param ... additional arguments
#' @details Once we have mapped cells between sample, we can identify which cells
#' tend to map to the same group of cells. To achieve this, we first create a 
#' cost matrix that will serve as a basis to find similar-mapping instances.
#' The cost matrix can be constructed from any cost matrix that was used during the 
#' mapping phase. 
#' Next, for each query cell we extract the top_nn cells in the reference with lowest
#' cost. Using the ordered index as a character label, we compute a jaccard index 
#' between overlapping labels. Query cells with a high jaccard index tend to map 
#' to the same reference cells. We then use the reciprocal to define a distance between
#' cells and cluster cells based on this distance.
#' The same approach is used for every clustering method provided. 
#' This function will add a new column with the metric clustering results.
#' @return vesalius_assay with clustering results 
#' @export
get_metric_clusters <- function(vesalius_assay,
    use_cost = "feature",
    cluster_method = "hclust",
    trial = NULL,
    group_identity = NULL,
    ref_cells = NULL,
    query_cells = NULL,
    top_nn = 30,
    h = 0.75,
    k = NULL,
    nn = 30,
    resolution = 1,
    verbose = TRUE,
    ...) {
    simple_bar(verbose)
    # add this correctly
    args <- list(...)
    #-------------------------------------------------------------------------#
    # set up cost matrics 
    #-------------------------------------------------------------------------#
    cost <- get_cost(vesalius_assay, use_cost)
    score <- concat_cost(cost, use_cost, complement = FALSE)[[1L]]
    #-------------------------------------------------------------------------#
    # get barcodes representing cells / sub group of cells 
    #-------------------------------------------------------------------------#
    cells <- dispatch_cost_groups(vesalius_assay,
        cost = score,
        trial = trial,
        group_identity = group_identity,
        ref_cells = ref_cells,
        query_cells = query_cells)
    query_cells <- match(cells$query, rownames(score))
    ref_cells <- match(cells$ref, colnames(score))
    score <- score[query_cells[!is.na(query_cells)],ref_cells[!is.na(ref_cells)]]
    score <- overlap_distance_matrix(score, top_nn, verbose)
    #-------------------------------------------------------------------------#
    # clustering 
    #-------------------------------------------------------------------------#
    clusters <- switch(EXPR = cluster_method,
        "hclust" = hclust_scores(score,  h, k,  verbose),
        "louvain" = louvain_scores(score, resolution, nn, verbose),
        "leiden" = leiden_scores(score, resolution, nn, verbose))
    #-------------------------------------------------------------------------#
    # rebuilding and adding to map slot 
    # for this we need to create a tag for both 
    #-------------------------------------------------------------------------#
    trial <- vesalius_assay@territories
    trial$trial <- "Not Selected"
    locs <- match(make.unique(names(clusters),sep = "-"), trial$barcodes)
    trial$trial[locs[!is.na(locs)]] <- unname(clusters)
    new_trial <- create_trial_tag(colnames(trial),
        "Map_cluster") %>%
        tail(1)
    colnames(trial) <- gsub("trial", new_trial, colnames(trial))
    #-------------------------------------------------------------------------#
    # update vesalius assay with 
    #-------------------------------------------------------------------------#
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
        data = trial,
        slot = "territories",
        append = TRUE)
    commit <- create_commit_log(arg_match = as.list(match.call()),
        default = formals(get_metric_clusters))
    vesalius_assay <- commit_log(vesalius_assay,
        commit,
        get_assay_names(vesalius_assay))
  simple_bar(verbose)
  return(vesalius_assay)
}

#' Hiearchal clustering of mapping scores for co mapping events
#' @param score matrix containing mapping scores
#' @param h numeric representing the normalized height that will be 
#' used for clustering (normalized to ensure that it will be between 0 and 1)
#' @param k int - number of clusters
#' @param verbose logical - print progress messages
#' @importFrom stats hclust cutree as.dist sd var IQR
hclust_scores <- function(score, h, k, verbose) {
    message_switch("hclust_scores", verbose)
    clusters <- hclust(as.dist(score))
    h <- h * (max(clusters$height) - min(clusters$height) / max(clusters$height))
    clusters <- cutree(clusters, h = h, k = k)
    return(clusters)
}



#' Using Louvain community clustering to for co-mapping events
#' @param score matrix containing mapping scores
#' @param resolution numeric - clustering resolution
#' @param nn int - number of nearest neighbors to use for graph construction
#' @param verbose logical - print output messages
#' @importFrom igraph graph_from_data_frame cluster_louvain
louvain_scores <- function(score, resolution, nn, verbose) {
    message_switch("louvain_scores", verbose)
    nn <- min(c(nn, nrow(score)))
    graph <- lapply(seq_len(ncol(score)),function(idx, score, nn) {
        return(order(score[, idx], decreasing = TRUE)[seq(1, nn)])
    })
    graph <- do.call("rbind",graph)
    graph <- populate_graph(graph)
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
    clusters <- igraph::cluster_louvain(graph, resolution = resolution)
    cluster <- clusters$membership
    names(cluster) <- clusters$names
    cluster <- cluster[!names(cluster) %in% colnames(score)]
    return(cluster)
}

#' Using Leiden community clustering to for co-mapping events
#' @param score matrix containing mapping scores
#' @param resolution numeric - clustering resolution
#' @param nn int - number of nearest neighbors to use for graph construction
#' @param verbose logical - print output messages
#' @importFrom igraph graph_from_data_frame cluster_leiden
leiden_scores <- function(score, resolution, nn, verbose) {
    message_switch("louvain_scores", verbose)
    graph <- lapply(seq_len(ncol(score)),function(idx, score, nn) {
        return(order(score[, idx], decreasing = TRUE)[seq(1, nn)])
    })
    graph <- do.call("rbind",graph)
    graph <- populate_graph(graph)
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
    clusters <- igraph::cluster_leiden(graph, resolution_parameter = resolution)
    cluster <- clusters$membership
    names(cluster) <- clusters$names
    cluster <- cluster[!names(cluster) %in% colnames(score)]
    return(cluster)
}

#' create a distance matrix based on mapping scores overlaps
#' @param score matrix - cost matrix 
#' @param top_nn integer nearest neighbors
#' @param verbose logical - print progress messages
#' @importFrom future.apply future_lapply
overlap_distance_matrix <- function(score, top_nn = 10, verbose = TRUE) {
  message_switch("overlap_scores", verbose)

  top_nn <- min(top_nn, ncol(score))
  rows <- rownames(score)
  n_cells <- nrow(score)

  # Create a character matrix with top NN column names for each row
  score_mat <- matrix(NA_character_, nrow = top_nn, ncol = n_cells)
  colnames(score_mat) <- rows
  rownames(score_mat) <- paste0("rank_", seq_len(top_nn))

  for (i in seq_len(n_cells)) {
    top_indices <- order(score[i, ], decreasing = TRUE)[seq_len(top_nn)]
    score_mat[, i] <- colnames(score)[top_indices]
  }

  # Compute the pairwise jaccard cost using the C++ function
  tmp <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    local_seed <- score_mat[, i, drop = FALSE]
    buffer <- future_lapply(seq_len(n_cells), function(idx) {
      local_query <- score_mat[, idx, drop = FALSE]
      local_score <- jaccard_cost(local_seed, local_query)
      colnames(local_score) <- colnames(local_seed)
      rownames(local_score) <- colnames(local_query)
      return(local_score)
    })
    tmp[[i]] <- do.call(rbind, buffer)
  }

  tmp <- do.call(cbind, tmp)
  colnames(tmp) <- rows
  rownames(tmp) <- rows

  return(1 - tmp)
}


#-----------------------------------------------------------------------------#
############################ Cost Contribution ################################
#-----------------------------------------------------------------------------#
#' Compute mapping score contribution to mapping
#' @param vesalius_assay vesalius_assay object after mapping a query 
#' onto a reference assay.
#' @param method character - which method to use to compute contribution
#' currently only dispersion availble
#' @param verbose logical - print progress messages
#' @export
get_cost_contribution <- function(vesalius_assay,
    method = "dispersion",
    verbose = TRUE){
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # check if cost matrix has been returned
    #-------------------------------------------------------------------------#
    cost <- check_cost_contribution(vesalius_assay)
    #-------------------------------------------------------------------------#
    # for now i will use dispersion scores
    #-------------------------------------------------------------------------#
    contribution <- switch(EXPR = method,
        "dispersion" = lapply(cost, dispersion),
        "CV" = lapply(cost, coef_of_var),
        "gini" = lapply(cost, gini),
        "IQR" = lapply(cost, IQR),
        "POC" = poc(cost))

    contribution <- data.frame("metric" = names(contribution),
        "score" = unlist(contribution))
    colnames(contribution) <- sub("score", method, colnames(contribution))
    contribution <- list("contribution_score" = contribution)
    #-------------------------------------------------------------------------#
    # rebuild and add
    #-------------------------------------------------------------------------#
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
        data = contribution,
        slot = "cost",
        append = TRUE)
    commit <- create_commit_log(arg_match = as.list(match.call()),
        default = formals(get_cost_contribution))
    vesalius_assay <- commit_log(vesalius_assay,
        commit,
        get_assay_names(vesalius_assay))
  simple_bar(verbose)
  return(vesalius_assay)
}


dispersion <- function(x) {
    return(var(x) / mean(x))
}

coef_of_var <- function(x) {
    return(sd(x) / mean(x))
}

gini <- function(x) {
    x <- sort(x)
    num <- sum(x * seq_along(x))
    dem <- seq_along(x) * sum(x) - num
    g <- num / dem
    return(g)
}

get_range <- function(x) {
    return(max(x) - min(x))
}

poc <- function(cost) {
    contrib <- sapply(cost,sum)
    contrib <- contrib / cost
    names(contrib) <- names(cost)
    return(contrib)
}
