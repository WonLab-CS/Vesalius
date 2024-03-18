###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################## COUNT INTEGRATION ##############################
#-----------------------------------------------------------------------------#

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
    seed_counts <- seed_counts[, match(mapped$prob$to, colnames(seed_counts))] 
    colnames(seed_counts) <- make.unique(mapped$prob$to, sep = "_")
    query_counts <- get_counts(query, type = "raw")
    query_counts <- query_counts[, match(mapped$prob$from, colnames(query_counts))]
    colnames(query_counts) <- make.unique(mapped$prob$from, sep = "_")
    integrated <- integrate_counts(seed_counts, query_counts, signal, verbose)
    counts <- integrated$counts
    embeds <- integrated$embeds
    comment(counts) <- "integrated"
    comment(embeds) <- "CCA"
    #-------------------------------------------------------------------------#
    # Prepare coordinates
    #-------------------------------------------------------------------------#
    seed_coord <- get_coordinates(seed,
        original = TRUE)[, c("barcodes","x_orig","y_orig","z")]
    colnames(seed_coord) <- c("barcodes", "x","y","z")
    seed_coord <- seed_coord[match(mapped$prob$to, seed_coord$barcodes), ]
    seed_coord <- adjust_duplicated_coord(seed_coord)
    query_coord <- align_index(mapped$prob, seed_coord)
    query_coord$z <- max(seed_coord$z) + 1
    #query_coord <- query_coord[match(mapped$prob$from, query_coord$barcodes), ]
    
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
    locs <- match(mapped$prob$to, seed_cells$barcodes)
    seed_cells <- cbind(
        seed_coord,
        seed_cells$trial[locs[!is.na(locs)]])
    colnames(seed_cells) <- c("barcodes", "x", "y", "z", "Cells")
    seed_cells$barcodes <- make.unique(mapped$prob$to, sep = "_")
    
    query_cells <- check_cell_labels(query, query_cell_labels)
    locs <- match(mapped$prob$from, query_cells$barcodes)
    query_cells <- cbind(
        query_coord,
        query_cells$trial[locs[!is.na(locs)]])
    
    colnames(query_cells) <- c("barcodes", "x", "y", "z", "Cells")
    query_cells$barcodes <- make.unique(mapped$prob$from, sep = "_")
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
        vesalius_assay@active <- embeds
        vesalius_assay@embeddings <- setNames(list(embeds),"CCA")
        vesalius_assay@counts <- counts
        vesalius_assay@territories <- rbind(seed_cells, query_cells)
        vesalius_assay@image <- list("image" = seed_image)
        vesalius_assay@meta <- c(seed_scale,
            "orig_coord" = list(rbind(seed_coord, query_coord)),
            "mapping_scores" = list(mapped$prob))
        if (return_cost){
            vesalius_assay@meta$cost <- mapped$cost
        }

    } else {
        vesalius_assay_2 <- new("vesalius_assay")
        vesalius_assay_2@assay <- "query"
        vesalius_assay_2@tiles <- check_coordinates(
            query_coord,
            "query",
            verbose = FALSE)
        counts_2 <- counts[c("raw_2","log_norm_2","integrated")]
        names(counts_2) <- c("raw","log_norm","integrated")
        comment(counts_2) <- "integrated"
        counts_2$integrated <- counts_2$counts[,query_coord$barcodes]
        vesalius_assay_2@counts <- counts_2
        vesalius_assay_2@territories <- query_cells
        vesalius_assay_2@image <- list("image" = query_image)
        vesalius_assay_2@meta <- c(query_scale,
            "orig_coord" = list(query_coord),
            "mapping_scores" = list(mapped$prob))
        vesalius_assay_1 <- new("vesalius_assay")
        vesalius_assay_1@assay <- "reference"
        vesalius_assay_1@tiles <- check_coordinates(
            seed_coord,
            "reference",
            verbose = FALSE)
        counts_1 <- counts[c("raw_1","log_norm_1","integrated")]
        names(counts_1) <- c("raw","log_norm","integrated")
        comment(counts_1) <- "integrated"
        counts_1$integrated <- counts_1$counts[,seed_coord$barcodes]
        vesalius_assay_1@counts <- counts_1
        vesalius_assay_1@territories <- seed_cells
        vesalius_assay_1@image <- list("image" = seed_image)
        vesalius_assay_1@meta <- c(seed_scale,
            "orig_coord" = list(seed_coord),
            "mapping_scores" = list(mapped$prob))
        if (return_cost){
            vesalius_assay_1@meta$cost <- mapped$cost
            vesalius_assay_2@meta$cost <- mapped$cost
        }
        vesalius_assay <- list("seed" = vesalius_assay_1,
            "query" = vesalius_assay_2)
    }
    return(vesalius_assay)
}


filter_maps <- function(mapped, threshold, allow_duplicates, verbose) {
    message_switch("post_map_filter",verbose)
    map_score <- mapped$prob
    #-------------------------------------------------------------------------#
    # First we remove points that have a score below threshold
    #-------------------------------------------------------------------------#
    cols <- seq((grep("init", colnames(map_score)) + 1), ncol(map_score))
    if ( length(cols) == 1) {
        tmp <- matrix(map_score[,cols], ncol = length(cols))
    } else {
        tmp <- map_score[,cols]
    }
    
    locs <- apply(
        X = tmp,
        MARGIN = 1,
        FUN = function(r, t) {sum(r < t)}, t = threshold)
    map_score <- map_score[locs == 0, ]
    if (nrow(map_score) == 0){
        stop("No cells retained under current filter threshold")
    }
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

#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Seurat ScaleData RunPCA IntegrateLayers
integrate_counts <- function(seed,
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
    integrated <- merge(seed, query)
    integrated <- integrated %>% 
        Seurat::NormalizeData(verbose = FALSE) %>%
        Seurat::FindVariableFeatures(verbose = FALSE, nfeatures = 2000) %>%
        Seurat::ScaleData(verbose = FALSE) %>%
        Seurat::RunPCA(verbose = FALSE)
    features <- check_feature_integration(signal, integrated)
    integrated <- Seurat::IntegrateLayers(integrated,
        method = "CCAIntegration",
        new.reduction = "integrated",
        features = features,
        verbose = FALSE)
    counts <- integrated@assays$RNA@layers
    counts <- rename_counts(counts,
        seed_cells,
        seed_genes,
        query_cells,
        query_genes,
        features)
    counts <- back_infer(counts, integrated@reductions$integrated)
    integrated <- integrated@reductions$integrated@cell.embeddings
    return(list("counts" = counts, "embeds" = integrated))
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
    names(counts) <- c("raw_1", "log_norm_1", "raw_2", "log_norm_2", "scaled")
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
    matched_index$x <- jitter(seed$x, amount = 1)
    matched_index$y <- jitter(seed$y, amount = 1)
    matched_index <- matched_index[, c("from", "x", "y")]
    colnames(matched_index) <-  c("barcodes", "x", "y")
    matched_index$barcodes <- make.unique(matched_index$barcodes, sep = "_")
    return(matched_index)
}

adjust_duplicated_coord <- function(coord) {
    dups <- duplicated(coord$barcodes)
    coord$x[dups] <- jitter(coord$x[dups], amount = 1)
    coord$y[dups] <- jitter(coord$y[dups], amount = 1)
    coord$barcodes <- make.unique(coord$barcodes, sep = "_")
    return(coord)
}

#' @importFrom Matrix Matrix
back_infer <- function(counts, embeds) {
    integrate <- t(as.matrix(cbind(counts[["log_norm_1"]], counts[["log_norm_2"]])))
    integrate <- integrate[match(rownames(integrate), rownames(embeds@cell.embeddings)),
        match(rownames(embeds@feature.loadings),colnames(integrate))]
    mu <- colMeans(integrate)
    back <- embeds@cell.embeddings %*% t(embeds@feature.loadings)
    back <- t(scale(back, center = -mu, scale = FALSE))
    back[which(back <= 0, arr.ind = TRUE)] <- 0
    counts <- c(counts, setNames(list(Matrix(back)), "integrated"))
    return(counts)
}
