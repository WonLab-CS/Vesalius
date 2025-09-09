###############################################################################
################################   Vesalius      ##############################
###############################################################################

#--------------------------------/ From Vesalius /----------------------------#


as.Seurat.vesalius_assay <- function(obj) {
    coord <- get_coordinates(obj)
    counts <- get_counts(obj, type = "all")
    raw <- counts[names(counts) == "raw"]
    counts <- counts[names(counts) != "raw"]
    embeddings <- get_embeddings(obj, active = FALSE)
    img <- obj@img
}


as.SingleCellExperiment.vesalius_assay <- function(obj, assay = NULL) {

}

as.SpatialExperiment.vesalius_assay <- function(obj, assay = NULL) {

}

as.anndata.vesalius_assay <- function(obj, assay = NULL) {

}

#--------------------------------/ To Vesalius /------------------------------#