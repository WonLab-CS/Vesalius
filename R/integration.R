###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################## COUNT INTEGRATION ##############################
#-----------------------------------------------------------------------------#

#' integrate 2 vesalius assays
#' @param mapped vesalius_assay object - mapped veslius assay (map_assays)
#' @param reference vesalius_assay object - reference vesalius_assay (seed)
#' @param method character - count integration method (methods provided by 
#' Seurat v5)
#' @param nfeatures integer - number of variable features to sue during 
#' integration.
#' @param signal character - defining which signal should be returned:
#' variable_features, all_features or custom gene list.
#' @param dimensions interger - number of dimensions integrated latent space
#' dimensions.
#' @param infer logical - back infer original counts by reversing reduced 
#' dimensional space roations. 
#' @param use_counts character - which count matrix to use during integration
#' @param labels_mapped  character - which columns in the mapped data assay
#' should be merged with the reference data (see details)
#' @param labels_reference character - which columns in the reference data assay
#' should be merged with the mapped data (see details)
#' @param verbose logical - should progressed message be printed
#' @details After mapping coordinates from a query onto a reference, vesalius
#' provides a way to then integrate the assays together. This function will:
#'
#' * Integrate Counts using Seurat 
#' * Merge coordinates (adding a jitter to avoid overlapping coordinayes)
#' * Merge territories (pair-wise merging using labels_mapped and labels_reference 
#' - everything else will have a separate column).
#'
#' We also infer log nomalized counts from CCA latent space. 
#'
#' The final output of this function is a vesalius_assay object containing 
#' coordinates from the mapped and reference, integrated latent space (e.g.CCA)
#' integrated counts, merged territories, an additional meta data.
#'
#' It should be noted that this object does not contain tiles. To use this 
#' vesalius_assay as any other, add tiles by using \code{\link{generate_tiles}}
#' @return a vesalius_assay object
#' @export 
integrate_assays <- function(mapped,
    reference,
    method = "CCAIntegration",
    nfeatures = 2000,
    signal = "variable_features",
    dimensions = 30,
    infer = TRUE,
    use_counts = "raw",
    labels_mapped = NULL,
    labels_reference = NULL,
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # check map and get counts - use only barcodes that were actually mapped
    # NOTE update this - need to be able to integrate assays without needing 
    # mapped
    #-------------------------------------------------------------------------#
    query_map <- check_maps(mapped)
    query_counts <- get_counts(mapped, type = use_counts)
    from <- query_map$from
    query_counts <- query_counts[, match(from, colnames(query_counts))]
    reference_counts <- get_counts(reference, type = use_counts)
    to <- query_map$to
    reference_counts <- reference_counts[, match(to, colnames(reference_counts))] 
    #-------------------------------------------------------------------------#
    # Count integration using Seurat methods 
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
        labels_mapped,
        labels_reference)
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

#' Integrate counts using Seurat
#' @param matched matrix - matrix containing counts from matched/mapped assay
#' @param reference matrix - matrix containing counts from reference assay
#' @param method character - Seurat integration method to use
#' @param nfeatures integer - number of features to use during integration
#' @param dimensions interger - number of dimensions integrated latent space
#' dimensions.
#' @param infer logical - back infer original counts by reversing reduced 
#' dimensional space roations. 
#' @param signal character - defining which signal should be returned:
#' variable_features, all_features or custom gene list.
#' @param verbose logical - print output messages
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#' @importFrom Seurat ScaleData RunPCA IntegrateLayers
integrate_counts <- function(matched,
    reference,
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
    
    counts <- rename_counts(integrated,
        reference_cells,
        reference_genes,
        matched_cells,
        matched_genes)
    if (infer) {
        back <- back_infer(counts, integrated@reductions$integrated)
        tags <- c(names(counts), "inferred")
        counts <- c(counts, back)
        names(counts) <- tags
    }
    integrated <- integrated@reductions$integrated@cell.embeddings
    integrated <- apply(integrated, 2, min_max)
    return(list("counts" = counts, "integrated" = integrated))
}

#' renaming counts to remain consistent with vesalius nomenclature
#' @param integrated seurat object containing integrated count data
#' @param seed_cells character vector contain seed barcodes 
#' @param seed_genes character vector contain seed genes 
#' @param query_cells character vector contain query barcodes 
#' @param query_genes character vector contain query genes
#' @return list of count matrices with cell and gene names added and 
#' each matrix is renamed to follow the vesalius nomencalture  
rename_counts <- function(integrated,
    seed_cells,
    seed_genes,
    query_cells,
    query_genes) {
    features <- Seurat::VariableFeatures(integrated)
    counts <- integrated@assays$RNA@layers
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
    integrated <- counts[["scale.data"]]
    colnames(integrated) <- c(seed_cells, query_cells)
    rownames(integrated) <- features
    counts <- c(seed,query, list(integrated))
    names(counts) <- c("raw_1", "log_norm_1", "raw_2", "log_norm_2", "scaled")
    return(counts)

}

#' Infer normalized counts from shared embedding values 
#' @param counts list of count matrices containing normalized gene expression
#' @param embeds shared embedding matrix produced by Seurat.
#' @return list containing inferred count values from embeddings  
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

#' Merge coordinates from 2 vesalius assays
#' @param matched matched vesalius_assay object (query object)
#' @param reference reference vesalius_assay object
#' @param barcodes barcodes to use for matching and merging
#' @details Simple merge followed a duplication check. Duplicated
#' coordinates produce errors when running tesselation. Just
#' adding some noise
#' @return merged coordinate data frame (barcodes, x, y)
merge_coordinates <- function(matched, reference, barcodes) {
    matched <- get_coordinates(matched, original = FALSE)
    colnames(matched) <- c("barcodes", "x","y")
    matched  <- matched[, c("barcodes", "x","y")]
    reference <- get_coordinates(reference, original = FALSE)
    colnames(reference) <- c("barcodes", "x","y")
    reference <- reference[, c("barcodes", "x","y")]
    merged <- rbind(matched, reference)
    locs <- paste0(merged$x,"_",merged$y)
    dups <- duplicated(locs)
    while (sum(dups) > 0) {
        warning("No duplicated coordinates allowed - Adding noise!")
        merged$x[dups] <- merged$x[dups] + runif(length(merged$x[dups]), min = 1, max = 2)
        merged$y[dups] <- merged$y[dups] + runif(length(merged$y[dups]), min = 1, max = 2)
        locs <- paste0(merged$x,"_",merged$y)
        dups <- duplicated(locs)
    }
    colnames(merged) <- c("barcodes", "x","y")
    merged <- merged[merged$barcodes %in% barcodes, c("barcodes", "x","y")]
    return(merged)
}

#' Merge territories using specified labels
#' @param matched matched (query) vesalius_assay object
#' @param reference reference vesalius_assay
#' @param coordinates coord data frame after merging (see merge_coordinates)
#' @param labels_mapped Territory columns to merge in mapped object
#' @param labels_reference Territory columns to merge in reference object
#' @details Territory data frame generally contain a lot of information 
#' depending on how much analysis the user has run. 
#' (Segments, Cells, territories, morphology). 
#' Using pariwise matching of labels parse in labels_mapped and labels_reference
#' the data will be merged into a single column. 
#' A new column will be added to distinguish samples. All other columns that are not
#' present in the labels argument will retain their name and NA will be added 
#' to cell present in the other data set (similar to left and right joins).
#' @return territory data frame with merged columns 
#' @importFrom dplyr right_join
merge_territories <- function(matched,
    reference,
    coordinates,
    labels_mapped,
    labels_reference) {
    matched <- get_territories(matched)
    matched_barcodes <- matched$barcodes
    reference <- get_territories(reference)
    reference_barcodes <- reference$barcodes
    coordinates$sample <- NA
    coordinates$sample[coordinates$barcodes %in% matched_barcodes] <- "matched"
    coordinates$sample[coordinates$barcodes %in% reference_barcodes] <- "reference"
    #-------------------------------------------------------------------------#
    # Add collumns that should not ne merged
    #-------------------------------------------------------------------------#
    if (!is.null(matched)){
        matched_cols <- colnames(matched)[!colnames(matched) %in%
            c("barcodes","x","y", labels_mapped)]
        for (i in seq_along(matched_cols)) {
            local_name <- paste0("mapped_",matched_cols[i])
            tmp <- match(coordinates$barcodes, matched$barcodes)
            coordinates[local_name] <- tmp
            coordinates[!is.na(tmp), local_name] <- matched[!is.na(tmp),matched_cols[i]]
        }
    }
    if (!is.null(reference)){
        reference_cols <- colnames(reference)[!colnames(reference) %in%
            c("barcodes","x","y", labels_reference)]
        for (i in seq_along(reference_cols)) {
            local_name <- paste0("reference_",reference_cols[i])
            tmp <- match(coordinates$barcodes, reference[reference_cols[i]])
            coordinates[local_name] <- tmp
            coordinates[!is.na(tmp), local_name] <- reference[!is.na(tmp),reference_cols[i]]
        }
    }
    #-------------------------------------------------------------------------#
    # Merge columns 
    # Could add this as sanity check function
    #-------------------------------------------------------------------------#
    if (!is.null(labels_mapped) && !is.null(labels_reference)) {
        if (length(labels_mapped) != length(labels_reference)) {
            stop("label vectors are not same length - Only pairwise merging possible!")
        }
        for (i in seq_along(labels_mapped)) {
            mapped_local <- match(matched$barcodes, coordinates$barcodes)
            ref_local <- match(reference$barcodes, coordinates$barcodes)
            name_local <- paste0(labels_mapped[i],"_",labels_reference[i])
            coordinates[name_local] <- NA
            coordinates[mapped_local[!is.na(mapped_local)], name_local] <-
                matched[!is.na(mapped_local),labels_mapped[i]]
            coordinates[ref_local[!is.na(ref_local)], name_local] <-
                reference[!is.na(ref_local),labels_reference[i]]
        }
    }
    rownames(coordinates) <- NULL
    return(coordinates)
}


