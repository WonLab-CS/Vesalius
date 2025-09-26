
#------------------------/ Preprocessing counts /------------------------------#
# NOTE: this might change if we decide do move away from Seurat as dependancy.

#' process counts 
#' 
#' pre-process count matrices
#' @param counts count matrix in the form of a sparse matrix 
#' @param assay character string describing the assay that is being 
#' pre-processed in the vesaliusObject or vesalius_assay
#' @param method character string describing which normalisation method to use.
#' One of the following "log_norm", "SCT", "TFIDF", "none".
#' @param use_count string describing which counts should be used for the 
#' generating emebddings. Default = "raw".
#' @param nfeatures numeric describing the number of variable features to use.
#' @param min_cutoff only used when dimensionality reduction method is
#' LSI or LSI_UMAP cutoff for feature to be included in the 
#' VariableFeatures for the object.
#' @param verbose logical - progress messages outputed or not
#' @details The `use_count` argument specifies which count matrix should be used
#' for normalization. This argument is only necessary if you use a custom
#' normalised count matrix. In this case, set this argument to the name
#' you gave your count matrix (see \code{\link{add_counts}}) and
#' `generate_embeddings` will skip the normalization and use your custom
#' count matrix to generate image embeddings. 
#' @importFrom Seurat CreateSeuratObject CreateAssayObject
process_counts <- function(counts,
  assay,
  method = "log_norm",
  use_count = "raw",
  nfeatures = 2000,
  min_cutoff = "q5",
  verbose = TRUE) {
     message_switch("in_assay",
      verbose,
      comp_type = "Pre-processing counts",
      assay = assay)
    #--------------------------------------------------------------------------#
    # We are still going to use Seurat for now
    # rememmber that if we do decide to change things
    # we have to change things in the embbeddings as well
    #--------------------------------------------------------------------------#
    #counts <- Seurat::CreateAssayObject(counts = counts)
    counts <- suppressWarnings(Seurat::CreateSeuratObject(counts))
    counts <- switch(method,
                    "log_norm" = log_norm(counts, nfeatures),
                    "SCTransform" = int_sctransform(counts, nfeatures),
                    "TFIDF" = tfidf_norm(counts, min_cutoff = min_cutoff),
                    "none" = no_norm(counts, use_count))
    return(counts)
}

#' no norm
#' 
#' no normalisation applied simply return raw counts
#' @param counts seurat object containing counts
#' @param use_count string describing name that needs to be added to
#' list element. This list will be appended to the count slot in
#' the vesalius_assay. 
#' @details Here, either the user doesn't want to normalise the data or
#' they provide their custom count matrix. In this case, we parse it 
#' as "none" to avoid writing another function and add the custom name. 
#' @return list with seurat object used later and raw counts to be stored in
#' the vesalius objects 
#' @importFrom Seurat GetAssayData
#' @importFrom SeuratObject Features 
no_norm <- function(counts, use_count = "raw") {
    #--------------------------------------------------------------------------#
    # Essentially we want people to be able to parse their matrix
    # If they want to use a different type of norm method that is not present
    # or play around with parameters not provided by vesalius
    # they can do that and just always call norm
    # We are using this just for formating at the moment
    # We have to be a bit hacky with the Seurat object
    # This will become an issue when change are made to Seurat objects
    # Really need to move aways from seurat for internal computing
    #--------------------------------------------------------------------------#
    counts <- Seurat::NormalizeData(object = counts, verbose = FALSE)
    counts <- Seurat::ScaleData(counts, verbose = FALSE)
    counts <- suppressWarnings(Seurat::FindVariableFeatures(counts,
    nfeatures = length(Features(counts)),
    verbose = FALSE))
    counts@assays$RNA@layers$data <- as.matrix(Seurat::GetAssayData(counts,
        layer = "counts"))
    counts@assays$RNA@layers$scale.data <- as.matrix(Seurat::GetAssayData(counts,
        layer = "counts"))
    
    return(list("SO" = counts, "norm" = NULL))
}
#' log norm
#' 
#' log normalisation, scaling and top variable features
#' @param counts seurat object containing counts
#' @param nfeatures number of top variable features to select
#' @return list with seurat object used later and normalised counts to be stored
#' in a vesalius object
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat GetAssayData
log_norm <- function(counts, nfeatures) {
  counts <- Seurat::NormalizeData(object = counts, verbose = FALSE)
  counts <- Seurat::ScaleData(counts, verbose = FALSE)
  counts <- suppressWarnings(Seurat::FindVariableFeatures(counts,
    nfeatures = nfeatures,
    verbose = FALSE))
  norm_counts <- list(Seurat::GetAssayData(counts, layer = "data"))
  names(norm_counts) <- "log_norm"
  return(list("SO" = counts, "norm" = norm_counts))
}

#' SCTransform
#' 
#' SCTransform pre-processing from Seurat
#' @param counts seurat object containing counts
#' @param nfeatures number of top variable features to select
#' @return list with seurat object used later and normalised counts to be stored
#' in a vesalius object
#' @importFrom Seurat SCTransform
#' @importFrom Seurat GetAssayData
int_sctransform <- function(counts, nfeatures) {
    counts <- suppressWarnings(Seurat::SCTransform(counts,
      variable.features.n = nfeatures, verbose = FALSE))
    norm_counts <- list(GetAssayData(counts, layer = "data"))
    names(norm_counts) <- "SCTransform"
    return(list("SO" = counts, "norm" = norm_counts))
}

#' tf idf normalisation
#' 
#' nornalising count using TF IDF
#' @param counts Seurat object containing counts
#' @param min_cutoff min cutoff of for top features
#' list with seurat object used later and normalised counts to be stored
#' in a vesalius object
#' @importFrom Signac RunTFIDF
#' @importFrom Seurat ScaleData
#' @importFrom Signac FindTopFeatures
#' @importFrom Seurat GetAssayData
tfidf_norm <- function(counts, min_cutoff) {
  counts <- Signac::RunTFIDF(counts, verbose = FALSE)
  counts <- Seurat::ScaleData(counts, verbose = FALSE)
  counts <- Signac::FindTopFeatures(counts, min.cutoff = min_cutoff)
  norm_counts <- list(Seurat::GetAssayData(counts, layer = "data"))
  names(norm_counts) <- "TFIDF"
  return(list("SO" = counts, "norm" = norm_counts))
}


