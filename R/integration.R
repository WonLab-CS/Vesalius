###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################## COUNT INTEGRATION ##############################
#-----------------------------------------------------------------------------#


integrate_assays <- function(seed,
    query,
    method = "harmony",
    infer = TRUE,
    use_counts = "raw",
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # First check if mapping has been done
    # If not we can just do the count integration
    #-------------------------------------------------------------------------#
    query_map <- check_maps(query)
    if (is.null(query_map)) {
        message_switch("count_only",
            verbose = verbose)
        seed_counts <- get_counts(seed, type = use_counts)
        query_counts <- get_counts(query, type = use_counts)
    } else {
        seed_counts <- get_counts(seed, type = use_counts)
        seed_counts <- seed_counts[, match(query_map$mapping, colnames(seed_counts))] 
        colnames(seed_counts) <- make.unique(mapped$prob$to, sep = "_")
        query_counts <- get_counts(query, type = use_counts)
        query_counts <- query_counts[, match(mapped$prob$from, colnames(query_counts))]
        colnames(query_counts) <- make.unique(mapped$prob$from, sep = "-")
    }

    counts <- integrate_counts(seed = seed,
            query = query,
            method = method,
            infer = infer,
            verbose = verbose)
}


#' @export
integrate_assays <- function(mapped,
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
    colnames(query_counts) <- make.unique(mapped$prob$from, sep = "-")
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


#' @importFrom Matrix Matrix
back_infer <- function(counts, embeds) {
    genes <- intersect(rownames(counts[["log_norm_1"]]), rownames(counts[["log_norm_2"]]))
    log_1 <- counts[["log_norm_1"]][genes, ]
    log_2 <- counts[["log_norm_2"]][genes, ]
    integrate <- t(as.matrix(cbind(log_1, log_2)))
    integrate <- integrate[match(rownames(integrate), rownames(embeds@cell.embeddings)),
        match(rownames(embeds@feature.loadings),colnames(integrate))]
    mu <- colMeans(integrate)
    back <- embeds@cell.embeddings %*% t(embeds@feature.loadings)
    back <- t(scale(back, center = -mu, scale = FALSE))
    back[which(back <= 0, arr.ind = TRUE)] <- 0
    counts <- c(counts, setNames(list(Matrix(back)), "integrated"))
    return(counts)
}
