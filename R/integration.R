###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################## COUNT INTEGRATION ##############################
#-----------------------------------------------------------------------------#

#' integrate counts from 2 vesalius assays
#' @param mapped vesalius_assay object - matched veslius assay (map_assays)
#' @param reference vesalius_assay object - reference vesalius_assay (seed)
#' @param method character - count integration method (methods provided by 
#' Seurat v5)
#' @param infer logical - back infer original counts by reversing reduced 
#' dimensional space roations. 
#' @param use_counts character - which count matrix to use during integration
#' @param verbose logical - should progressed message be printed
#' @export 
integrate_assays <- function(mapped,
    reference,
    method = "CCAIntegration",
    nfeatures = 2000,
    signal = "variable_features",
    dimensions = 30,
    infer = TRUE,
    use_counts = "raw",
    cell_label_mapped = "Cells",
    cell_label_reference = "Cells",
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # check map and get counts - use only barcodes that were actually mapped
    #-------------------------------------------------------------------------#
    query_map <- check_maps(mapped)

    query_counts <- get_counts(mapped, type = use_counts)
    #from <- remove_suffix(query_map$from)
    from <- query_map$from
    query_counts <- query_counts[, match(from, colnames(query_counts))]

    reference_counts <- get_counts(reference, type = use_counts)
    #to <- remove_suffix(query_map$to)
    to <- query_map$to
    reference_counts <- reference_counts[, match(to, colnames(reference_counts))] 

    #-------------------------------------------------------------------------#
    # Count integration using Seurat methods 
    # Using harmony as default since it seems that it is the best method
    #-------------------------------------------------------------------------#
    integrated <- integrate_counts(matched = query_counts,
        reference = reference_counts,
        method = method,
        nfeatures = nfeatures,
        infer = infer,
        signal = signal,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # Construct vesalius object 
    #-------------------------------------------------------------------------#
    coordinates <- merge_coordinates(mapped,
        reference,
        rownames(integrated$integrated))
    territories <- merge_territories(mapped,
        reference,
        coordinates,
        cell_label_mapped,
        cell_label_reference)
    cells <- merge_cells(mapped,
        reference,
        rownames(integrated$integrated),
        cell_label_mapped,
        cell_label_reference)
    vesalius_assay <- build_vesalius_assay(coordinates = coordinates, verbose = FALSE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay,
        data = integrated$counts,
        slot = "counts",
        append = FALSE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay,
        data = integrated["integrated"],
        slot = "embeddings",
        append = TRUE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay,
        data = integrated[["integrated"]],
        slot = "active",
        append = FALSE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay,
        data = territories,
        slot = "territories",
        append = FALSE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay,
        data = mapped@map,
        slot = "map",
        append = FALSE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay,
        data = mapped@cost,
        slot = "cost",
        append = FALSE)
    vesalius_assay <- add_cells(vesalius_assay, cells = cells, verbose = FALSE)
    message_switch("integrated_embed", verbose, tag = "integrated")
    vesalius_assay <- add_active_embedding_tag(vesalius_assay, "integrated")
    message_switch("integrated_counts", verbose,
        tag = ifelse(infer, "inferred","scaled"))
    vesalius_assay <- add_active_count_tag(vesalius_assay,
        norm = ifelse(infer, "inferred","scaled"))
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(integrate_assays))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      "Integrated")
    simple_bar(verbose)
    return(vesalius_assay)
}

#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Seurat ScaleData RunPCA IntegrateLayers
integrate_counts <- function(matched,
    reference,
    features,
    method = "HarmonyIntegration",
    nfeatures = 2000,
    dimensions = 30,
    infer = FALSE,
    signal = "variable_features",
    verbose) {
    
    message_switch("integrate",verbose)
    matched_cells <- colnames(matched)
    matched_genes <- rownames(matched)
    matched <- Seurat::CreateSeuratObject(matched)
    reference_cells <- colnames(reference)
    reference_genes <- rownames(reference)
    reference <- Seurat::CreateSeuratObject(reference)
    integrated <- merge(reference, matched)
    integrated <- integrated %>% 
        Seurat::NormalizeData(verbose = FALSE) %>%
        Seurat::FindVariableFeatures(verbose = FALSE, nfeatures = nfeatures) %>%
        Seurat::ScaleData(verbose = FALSE) %>%
        Seurat::RunPCA(verbose = FALSE, npcs = dimensions)
    features <- check_feature_integration(signal, integrated)
    integrated <- Seurat::IntegrateLayers(integrated,
        method = method,
        new.reduction = "integrated",
        features = features,
        verbose = FALSE)
    counts <- integrated@assays$RNA@layers
    counts <- rename_counts(counts,
        reference_cells,
        reference_genes,
        matched_cells,
        matched_genes,
        features)
    if (infer) {
        back <- back_infer(counts, integrated@reductions$integrated)
        tags <- c(names(counts), "inferred")
        counts <- c(counts, back)
        names(counts) <- tags
    }
    integrated <- integrated@reductions$integrated@cell.embeddings
    return(list("counts" = counts, "integrated" = integrated))
}


rename_counts <- function(counts,
    seed_cells,
    seed_genes,
    query_cells,
    query_genes,
    features) {
    seed <- lapply(counts[c("counts.1","data.1")],function(counts, names){
            colnames(counts) <- names
            return(counts)
        }, names = seed_cells)
    seed <- lapply(seed,function(counts, names){
            rownames(counts) <- names
            return(counts)
        }, names = seed_genes)
    query <- lapply(counts[c("counts.2","data.2")],function(counts, names){
            colnames(counts) <- names
            return(counts)
        }, names = query_cells)
    query <- lapply(query,function(counts, names){
            rownames(counts) <- names
            return(counts)
        }, names = query_genes)
    intergrated <- counts[["scale.data"]]
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
    return(list(back))
}


merge_coordinates <- function(matched, reference, barcodes) {
    matched <- get_coordinates(matched, original = TRUE)
    colnames(matched) <- c("barcodes", "x","y")
    reference <- get_coordinates(reference, original = TRUE)
    colnames(reference) <- c("barcodes", "x","y")
    merged <- rbind(matched, reference)
    locs <- paste0(merged$x,"_",merged$y)
    dups <- duplicated(locs)
    while (sum(dups) > 0) {
        warning("No duplicated coordinates allowed - Adding noise!")    
        merged$x[dups] <- jitter(merged$x[dups], amount = 1)
        merged$y[dups] <- jitter(merged$y[dups], amount = 1)
        locs <- paste0(merged$x,"_",merged$y)
        dups <- duplicated(locs)
    }
    colnames(merged) <- c("barcodes", "x","y")
    merged <- merged[merged$barcodes %in% barcodes, c("barcodes", "x","y")]
    return(merged)
}
#' @importFrom dplyr right_join
merge_territories <- function(matched,
    reference,
    coordinates,
    cell_label_mapped,
    cell_label_reference) {
    matched <- get_territories(matched)
    matched_barcodes <- matched$barcodes
    reference <- get_territories(reference)
    reference_barcodes <- reference$barcodes
    if (!is.null(matched)){
        matched <- right_join(matched, coordinates, by = "barcodes")
        matched <- matched[,grep(x = colnames(matched),pattern = cell_label_mapped, invert = TRUE)]
        matched <- data.frame(matched[ ,!colnames(matched) %in% c("barcodes","x.x","y.x","x.y","y.y")])
        if (ncol(matched) > 0){
            colnames(matched) <- paste0("mapped_", colnames(matched))
        }
        
    }
    if (!is.null(reference)) {
        reference <- right_join(reference, coordinates, by = "barcodes")
        reference <- reference[,grep(x = colnames(reference),pattern = cell_label_reference, invert = TRUE)]
        reference <- data.frame(reference[ ,!colnames(reference) %in% c("barcodes","x.x","y.x","x.y","y.y")])
        if (ncol(reference) > 0){
            colnames(reference) <- paste0("ref_", colnames(reference))
        }
    }
    coordinates$sample <- NA
    coordinates$sample[coordinates$barcodes %in% matched_barcodes] <- "matched"
    coordinates$sample[coordinates$barcodes %in% reference_barcodes] <- "reference"
    merged <- cbind(coordinates, matched,reference)
    rownames(merged) <- NULL
    return(merged)
}

merge_cells <- function(matched,
    reference,
    barcodes,
    cell_labels_matched,
    cell_labels_reference) {
    matched <- get_territories(matched)
    matched_cells <- which(colnames(matched) %in% cell_labels_matched)
    if (length(matched_cells) > 0){
        matched_cells <- matched[, matched_cells]
        names(matched_cells) <- matched$barcodes
    } else {
        matched_cells <- NULL
    }
    reference <- get_territories(reference)
    reference_cells <- which(colnames(reference) %in% cell_labels_reference)
    if (length(reference_cells) > 0){
        reference_cells <- reference[, reference_cells]
        names(reference_cells) <- reference$barcodes
    } else {
        reference_cells <- NULL
    }
    merged_cells <- c(matched_cells,reference_cells)
    merged_cells <- merged_cells[names(merged_cells) %in% barcodes]
    return(merged_cells)
}
