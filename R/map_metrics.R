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
            group_1 <- query_assay@meta$mapping_score$to
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
            group_2 <- query_assay@meta$mapping_score$from
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




# #' remove scores that are below a certain threshold 
# #' @param cost list - score list to be filtered 
# #' @param threshold numeric - score threshold [0,1]
# #' @param use_cost - character - which score matrices to filter
# #' @return filtered cost list 
# #' @export
# filter_maps <- function(vesalius_assay, by = NULL, threshold = 0.3) {
#     if (is.null(by)) {
#         warning("No mapping score requested - returning assay as is")
#         return(vesalius_assay)
#     }
#     maps <- check_map_selection(vesalius_assay, by = by)
#     map_score <- vesalius_assay@meta$mapping_scores
#     for (i in seq_along(maps)) {
#         map_score <- map_score[map_score[, maps[i]] >= threshold, ]
#     }
#     # vesalius_assay <- filter_assay(vesalius_assay,
#     #     cells = map_score$from)
#     # vesalius_assay@meta$mapping_scores <- map_score
#     return(map_score)
# }

