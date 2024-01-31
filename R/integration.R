###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################## COUNT INTEGRATION ##############################
#-----------------------------------------------------------------------------#
#' @export 
integrate_map <- function(mapped,
    query,
    seed,
    signal = "variable_features",
    return_cost = FALSE,
    threshold = 0.3,
    allow_duplicates = TRUE,
    seed_cell_labels = NULL,
    query_cell_labels = NULL,
    merge = TRUE,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # filter based on threshold and duplicates
    #-------------------------------------------------------------------------#
    mapped <- filter_maps(mapped, threshold, allow_duplicates, verbose)
    #-------------------------------------------------------------------------#
    # prepare counts for integration 
    # Vesalius won't like having duplicates so got to make sure that 
    # all barcodes have unique names
    # we also filter any cells that need to be removed
    #-------------------------------------------------------------------------#
    seed_counts <- get_counts(seed, type = "raw")
    #seed_counts <- seed_counts[, colnames(seed_counts) %in% mapped$prob$to]
    colnames(seed_counts) <- make.unique(colnames(seed_counts), sep = "_")
    query_counts <- get_counts(query, type = "raw")
    query_counts <- query_counts[, colnames(query_counts) %in% mapped$prob$from]
    colnames(query_counts) <- make.unique(colnames(query_counts), sep = "_")
    intergrated <- intergrate_counts(seed_counts, query_counts, signal, verbose)
    comment(intergrated) <- "intergrated"
    #-------------------------------------------------------------------------#
    # Prepare coordinates
    #-------------------------------------------------------------------------#
    seed_coord <- get_coordinates(seed,
        original = TRUE)[, c("barcodes","x_orig","y_orig","z")]
    colnames(seed_coord) <- c("barcodes", "x","y","z")
    query_coord <- align_index(mapped$matched, seed_coord)
    query_coord$z <- max(seed_coord$z) + 1
    query_coord <- query_coord[query_coord$barcodes %in% mapped$prob$from,]
    #-------------------------------------------------------------------------#
    # check for images
    #-------------------------------------------------------------------------#
    seed_image <- check_image(seed, image_name = NULL, as_is = TRUE)
    query_image <- check_image(query, image_name = NULL, as_is = TRUE)

    #-------------------------------------------------------------------------#
    # check scales
    #-------------------------------------------------------------------------#
    seed_scale <- seed@meta[c("scale","unit")]
    query_scale <- query@meta[c("scale","unit")]
    #-------------------------------------------------------------------------#
    # check for cell type annotations
    #-------------------------------------------------------------------------#
    seed_cells <- check_cell_labels(seed, seed_cell_labels)
    seed_cells <- cbind(
        seed_coord[match(seed_cells$barcodes, seed_coord$barcodes), ],
        seed_cells$trial)
    colnames(seed_cells) <- c("barcodes","x", "y", "Cells")
    query_cells <- check_cell_labels(seed, seed_cell_labels)
    query_cells <- cbind(
        query_coord[match(query_cells$barcodes, query_coord$barcodes),],
        query_cells$trial)
    colnames(query_cells) <- c("barcodes","x", "y", "Cells")
    #-------------------------------------------------------------------------#
    # rebuild object
    #-------------------------------------------------------------------------#
    if (merge) {
        vesalius_assay <- new("vesalius_assay")
        vesalius_assay@assay <- "mapped"
        vesalius_assay@tiles <- check_coordinates(
            rbind(seed_coord, query_coord),
            "mapped",
            verbose = FALSE)
        vesalius_assay@counts <- intergrated
        vesalius_assay@territories <- rbind(seed_cells, query_cells)
        vesalius_assay@image <- list("image" = seed_image)
        vesalius_assay@meta <- c(seed_scale,
            "orig_coord" = list(rbind(seed_coord, query_coord)),
            "mapping_scores" = list(mapped$prob))

    } else {
        vesalius_assay <- new("vesalius_assay")
        vesalius_assay@assay <- "mapped"
        vesalius_assay@tiles <- check_coordinates(
            query_coord,
            "mapped",
            verbose = FALSE)
        intergrated <- intergrated[c("raw_query","log_norm_query","intergrated")]
        names(intergrated) <- c("raw","log_norm","intergrated")
        intergreated$intergrated <- intergreated$intergrated [,query_coord$barcodes]
        vesalius_assay@counts <- intergreated
        vesalius_assay@territories <- query_cells
        vesalius_assay@image <- list("image" = query_image)
        vesalius_assay@meta <- c(query_scale,
            "orig_coord" = list(query_coord),
            "mapping_scores" = list(mapped$prob))
    }
    #-------------------------------------------------------------------------#
    # adding cost matrices
    #-------------------------------------------------------------------------#
    if (return_cost){
        vesalius_assay@meta$cost <- mapped$cost
    }

    return(vesalius_assay) 
}


filter_maps <- function(mapped, threshold, allow_duplicates, verbose) {
    message_switch("post_map_filter",verbose)
    map_score <- mapped$prob
    #-------------------------------------------------------------------------#
    # First we remove points that have a score below threshold
    #-------------------------------------------------------------------------#
    locs <- apply(
        X = map_score[,seq((grep("cost", colnames(map_score)) + 1),
            ncol(map_score))],
        MARGIN = 1,
        FUN = function(r, t) {sum(r < t)}, t = threshold)
    map_score <- map_score[locs == 0, ]
    mapped$prob <- map_score
    #-------------------------------------------------------------------------#
    # Next we check for dupliactes and gert the best ones
    #-------------------------------------------------------------------------#
    if (!allow_duplicates) {
        map_score <- map_score[order(map_score$cost), ]
        duplicates <- duplicated(map_score$from) | duplicated(map_score$to)
        map_score <- map_score[!duplicates, ]
        mapped$prob <- map_score
    }
    
    #-------------------------------------------------------------------------#
    # Remove thos points from cost matrices
    #-------------------------------------------------------------------------#
    cost <- mapped$cost
    cost <- lapply(cost,
        function(cost, row, col) {
            return(cost[row, col])
        },row = map_score$from,
        col = map_score$to)
    mapped$cost <- cost
    return(mapped)
}

intergrate_counts <- function(seed,
    query,
    signal,
    verbose) {
    message_switch("integrate",verbose)
    seed_cells <- colnames(seed)
    seed_genes <- rownames(seed)
    seed <- Seurat::CreateSeuratObject(seed)
    query_cells <- colnames(query)
    query_genes <- rownames(query)
    query <- Seurat::CreateSeuratObject(query)
    intergrated <- merge(seed, query)
    intergrated <- intergrated %>% 
        Seurat::NormalizeData(verbose = FALSE) %>%
        Seurat::FindVariableFeatures(verbose = FALSE, n_features = 2000) %>%
        Seurat::ScaleData(verbose = FALSE) %>%
        Seurat::RunPCA(verbose = FALSE)
    features <- check_feature_integration(signal, intergrated)
    intergrated <- Seurat::IntegrateLayers(intergrated,
        method = "CCAIntegration",
        new.reduction = "integrated",
        features = features,
        verbose = FALSE)
    counts <- intergrated@assays$RNA@layers
    names(counts) <- c("raw_seed","raw_query","log_norm_seed","log_norm_query","intergrated")
    counts <- rename_counts(counts,
        seed_cells,
        seed_genes,
        query_cells,
        query_genes,
        features)
    
    return(counts)
}


rename_counts <- function(counts,
    seed_cells,
    seed_genes,
    query_cells,
    query_genes,
    features) {
    seed <- lapply(counts[c(1,3)],function(counts, names){
            colnames(counts) <- names
            return(counts)
        }, names = seed_cells)
    seed <- lapply(seed,function(counts, names){
            rownames(counts) <- names
            return(counts)
        }, names = seed_genes)
    query <- lapply(counts[c(2,4)],function(counts, names){
            colnames(counts) <- names
            return(counts)
        }, names = query_cells)
    query <- lapply(query,function(counts, names){
            rownames(counts) <- names
            return(counts)
        }, names = query_genes)
    intergrated <- counts[[5]]
    colnames(intergrated) <- c(seed_cells, query_cells)
    rownames(intergrated) <- features
    counts <- c(seed,query, list(intergrated))
    names(counts) <- c("raw_seed","raw_query","log_norm_seed","log_norm_query","intergrated")
    return(counts)

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
    seed) {
    seed <- seed[match(matched_index$to, seed$barcodes), ]
    matched_index$x <- seed$x
    matched_index$y <- seed$y
    matched_index$x <- jitter(matched_index$x, factor = 0.1)
    matched_index$y <- jitter(matched_index$y, factor = 0.1)
    matched_index <- matched_index[, c("from", "x", "y")]
    colnames(matched_index) <-  c("barcodes", "x", "y")
    return(matched_index)
}
