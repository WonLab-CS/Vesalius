###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################ Mapping Matrics####### ###########################
#-----------------------------------------------------------------------------#

#' @export
compare_assays <- function(seed_assay,
    query_assay,
    compare = "niche",
    method = "wilcox",
    group_1 = NULL,
    group_2 = NULL,
    aggregate = FALSE,
    neighborhood = "knn",
    signal = "variable_features",
    use_norm = "log_norm",
    seed_cells = "Cells",
    query_cells = "Cells",
    k = 10,
    depth = 1,
    radius = 1,
    log_fc = 0.25,
    pval = 0.05,
    min_pct = 0.05,
    min_spatial_index = 10,
    verbose = TRUE,
    ...) {
    simple_bar(verbose)
    args <- list(...)
    #-------------------------------------------------------------------------#
    # add some check here to make sure the assays are mapped assays and 
    # what not
    #-------------------------------------------------------------------------#
    method <- check_deg_method(method)
    if (compare == "feature" && is.null(group_1) && is.null(group_2)) {
        stop("If using 'feature' to compare assays \n 
        please specifiy barcodes to compare in group_1 and group_2")
    }
    #-------------------------------------------------------------------------#
    # dispatch groups
    #-------------------------------------------------------------------------#
    if (grepl("niche|composition", compare)) {
        if (is.null(group_1)){
            group_1 <- query_assay@map$mapping_score$to
        }
        group_1 <- dispatch_niches(assay = seed_assay,
            group = group_1,
            aggregate = aggregate,
            neighborhood = neighborhood,
            k = k,
            depth = depth,
            radius = radius,
            group_id = "1")
        if (is.null(group_2)){
            group_2 <- query_assay@map$mapping_score$from
        }
        group_2 <- dispatch_niches(assay = query_assay,
            group = group_2,
            aggregate = aggregate,
            neighborhood = neighborhood,
            k = k,
            depth = depth,
            radius = radius,
            group_id = "2")
    }
    if (grepl("territory", compare)) {
        group_1 <- dispatch_territories(seed_assay,
            group = group_1,
            aggregate = aggregate,
            group_id = "1")
        group_2 <- dispatch_territories(query_assay,
            group = group_2,
            aggregate = aggregate,
            group_id = "2")
    }
    
    if(length(group_1) != length(group_2)){
        stop("Groups must have the same length for pairwise comparison!
        Use aggregate = TRUE to aggregate groups of different lengths into 1")
    }
    #-------------------------------------------------------------------------#
    # get data to compare
    #-------------------------------------------------------------------------#
    if (grepl("niche|territory|feature", compare)) {
       signal <- get_signal(seed_assay,
            query_assay,
            signal,
            verbose = verbose)
    }
    if(grepl("composition", compare)) {
        signal <- get_composition(seed_assay,
            query_assay,
            seed_cells = seed_cells,
            query_cells = query_cells,
            verbose = verbose)
    }
    #-------------------------------------------------------------------------#
    # get diff expression
    #-------------------------------------------------------------------------#
    deg <- vector("list", length(group_1))
    for (i in seq_along(group_1)) {
        seed <- signal$seed[,group_1[[i]]]
        seed_id <- names(group_1)[i]
        query <- signal$query[,group_2[[i]]]
        query_id <- names(group_2)[i]

        deg[[i]] <- vesalius_deg(as.matrix(seed),
            as.matrix(query),
            seed_id,
            query_id,
            method,
            log_fc,
            pval,
            min_pct,
            min_spatial_index,
            verbose,
            args)
    }
    deg <- do.call("rbind", deg)
    simple_bar(verbose)
    return(deg)

}


#' @export 
get_neighborhood <- function(vesalius_assay,
    neighborhood = "knn",
    k = 10,
    depth = 1,
    radius = 1) {
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    niches <- switch(neighborhood,
        "knn" = knn_neighborhood(coord, k),
        "radius" = radius_neighborhood(coord, radius),
        "graph" = graph_neighborhood(coord, depth))
    return(niches)
}



dispatch_niches <- function(assay,
            group = NULL,
            aggregate = FALSE,
            neighborhood = "knn",
            k = 10,
            depth = 1,
            radius = 1,
            group_id = "1") {
    niches <- get_neighborhood(assay,
        neighborhood = neighborhood,
        k = k,
        depth = depth,
        radius = radius)
    loc <- match(group, names(niches))
    if (aggregate) {
        niches <- list(unique(unlist(niches[loc[!is.na(loc)]])))
        names(niches) <- paste0("group_", group_id)
    } else {
        niches <- niches[loc[!is.na(loc)]]
    }
    return(niches)
}


dispatch_territories <- function(assay,
    group,
    aggregate = FALSE,
    group_id = "1") {
    trial <- check_territory_trial(assay, "last")
    ter <- check_group_value(trial, group)
    if (is.null(ter)) {
        ter <- sort(unique(trial$trial))
    }
    trial <- lapply(ter, function(ter,trial){
        return(trial$barcodes[trial$trial == ter])
    }, trial = trial)
    names(trial) <- ter
    if (aggregate) {
        trial <- list(unlist(trial))
        names(trial) <- paste0("group_", group_id)
    }
    return(trial)
}

get_composition <- function(seed_assay,
    query_assay,
    seed_cells = "Cells",
    query_cells = "Cells",
    verbose = TRUE) {
    seed_cells <- check_territory_trial(seed_assay, trial = seed_cells)
    query_cells <- check_territory_trial(query_assay, trial = query_cells)
    cell_union <- sort(union(unique(seed_cells$trial),unique(query_cells$trial)))

    seed_composition <- matrix(0,
        nrow = length(cell_union),
        ncol = nrow(seed_cells))
    rownames(seed_composition) <- cell_union
    colnames(seed_composition) <- seed_cells$barcodes
    query_composition <- matrix(0,
        nrow = length(cell_union),
        ncol = nrow(query_cells))
    rownames(query_composition) <- cell_union
    colnames(query_composition) <- query_cells$barcodes

    seed_composition[match(seed_cells$trial, cell_union) + seq(0,l= nrow(seed_cells), by = nrow(seed_composition))] <- 1
    query_composition[match(query_cells$trial, cell_union) +seq(0,l= nrow(query_cells), by = nrow(query_composition))] <- 1
    return(list("seed" = seed_composition,"query" = query_composition,"feature" = cell_union))
    
}

#' @export
get_metric_clusters <- function(vesalius_assay,
    use_cost = "feature",
    cluster_method = "hclust",
    distance = "euclidean",
    trial = NULL,
    group_identity = NULL,
    by_similarity = TRUE,
    ref_cells = NULL,
    query_cells = NULL,
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
    score <- concat_cost(cost, use_cost, complement = !by_similarity)[[1L]]
    #-------------------------------------------------------------------------#
    # get barcodes representing cells / sub group of cells 
    #-------------------------------------------------------------------------#
    cells <- dispatch_cost_groups(vesalius_assay,
        cost = score,
        trial = trial,
        group_identity = group_identity,
        ref_cells = ref_cells,
        query_cells = query_cells)
    score <- score[match(cells$query, rownames(score)), cells$ref]
    #-------------------------------------------------------------------------#
    # clustering 
    #-------------------------------------------------------------------------#
    clusters <- switch(EXPR = cluster_method,
        "hclust" = hclust_scores(score, h, k, distance, verbose),
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

#' hclust of mapping scores
#' @param score matrix containing mapping scores
#' @param h cutree height will take precedent over k
#' @param k number of cluster to get from hclust
#' @param distance which distance metric should be used to compute hclust
#' @param verbose logical for progress message output
#' @return named vector - vector = clusters and names = barcodes
#' @importFrom stats hclust cutree
hclust_scores <- function(score, h, k, distance, verbose) {
    message_switch("hclust_scores", verbose)
    clusters <- hclust(dist(score, method = distance))
    h <- h * (max(clusters$height) - min(clusters$height) / max(clusters$height))
    clusters <- cutree(clusters, h = h, k = k)
    return(clusters)
}


#' louvain clustering applied to mapping scores
#' @param scores matrix containing mapping scores
#' @param resolution numeric - clustering resolution
#' @param nn integer - number of nearest neighbors based on scores
#' @param verbose logical - if progress message should be printed
#' @details Unweighted louvain custering could be weighted by score
#' this could help discriminate of scores are low but still above 
#' filter criteria
#' @return named vector - vector = clusters and names = barcodes
#' @importFrom igraph graph_from_data_frame cluster_louvain
louvain_scores <- function(score, resolution, nn, verbose) {
    message_switch("louvain_scores", verbose)
    locs <- lapply(seq_len(nrow(score)), function(idx, scores, nn){
        locs <- order(scores[idx, ],decreasing = TRUE)[seq(1,nn)]
        local_graph <- data.frame("from" = rownames(scores)[idx],
            "to" = colnames(score)[locs])
        return(local_graph)
    }, scores = score, nn = nn)
    graph <- do.call("rbind", locs)
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
    clusters <- igraph::cluster_louvain(graph, resolution = resolution)
    cluster <- clusters$membership
    names(cluster) <- clusters$names
    cluster <- cluster[!names(cluster) %in% colnames(score)]
    return(cluster)
}

#' leiden clustering applied to mapping scores
#' @param scores matrix containing mapping scores
#' @param resolution numeric - clustering resolution
#' @param nn integer - number of nearest neighbors based on scores
#' @param verbose logical - if progress message should be printed
#' @details Unweighted leiden custering could be weighted by score
#' this could help discriminate of scores are low but still above 
#' filter criteria
#' @return named vector - vector = clusters and names = barcodes
#' @importFrom igraph graph_from_data_frame cluster_leiden
leiden_scores <- function(score, resolution, nn, verbose) {
    message_switch("louvain_scores", verbose)
    locs <- lapply(seq_len(nrow(score)), function(idx, scores, nn){
        locs <- order(scores[idx, ],decreasing = TRUE)[seq(1,nn)]
        local_graph <- data.frame("from" = rownames(scores)[idx],
            "to" = colnames(score)[locs])
        return(local_graph)
    }, scores = score, nn = nn)
    graph <- do.call("rbind", locs)
    graph <- graphigraph::graph_from_data_frame(graph, directed = FALSE)
    clusters <- igraph::cluster_leiden(graph, resolution_parameter = resolution)
    cluster <- clusters$membership
    names(cluster) <- clusters$names
    cluster <- cluster[!names(cluster) %in% colnames(score)]
    return(cluster)
}



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
        "dispersion" = lapply(cost, dispersion))
    contribution <- data.frame("metric" = names(contribution),
        "score" = unlist(contribution))
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
    return(var(x)/mean(x))
}