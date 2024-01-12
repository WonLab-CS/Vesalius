###############################################################################
################################   Vesalius      ##############################
###############################################################################

#-----------------------------------------------------------------------------#
############################ Mapping Matrics####### ###########################
#-----------------------------------------------------------------------------#
#' @export 
get_neighborhood_matrix <- function(vesalius_assay,
    neighborhood = "knn",
    signal = "variable_features",
    use_norm = "raw",
    k = 10,
    depth = 1,
    radius = 1,
    dimensions = seq(1, 30),
    verbose = TRUE) {
    features <- check_signal(vesalius_assay, signal, type = use_norm)
    signal <- get_counts(vesalius_assay, type = "use_norm")[features, ]
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    niche <- get_neighborhood(coord,
        signal,
        neighborhood,
        k,
        depth,
        radius)
    return(niche)
}


get_composition_matrix <- function(vesalius_assay,
    as_graph = FALSE,
    verbose = TRUE){

}