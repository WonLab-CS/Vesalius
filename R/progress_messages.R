################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Progress Message/--------------------------------#
# Section for progress messages. Essentially this will print out
# whatever message we need. Just add element to switch to add new message.
# Prettify is only for graphical stuff. Could add more stuff.

#-----------------------------/Prettify/-----------------------------------#
# Simple bar output used to start any progress message output
simple_bar <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    width_display <- round(options()$width)
    bar <- paste0("#",
        paste0(rep("-", width_display), collapse = ""),
        "#")
    cat(bar, "\n")

}

#-----------------------------/Messages/-----------------------------------#

#' switch message output based in input
#' @param type name of message output to produce
#' @param verbose logical if message should be outputed
#' @param ... any other parameter
#' @details Essentially, all message types are listed in this giant switch.
#' The type defines what message you want to output. We select named
#' arguments from `...` and parse them when required.
#' Adding message is straight forward by simply using the template in 
#' other switch types. This also means that you can add and remove messages
#' without causing errors. While this might lead to dead code, better dead
#' code than buggy code with undefined functions...
message_switch <- function(type, verbose = TRUE, ...) {
    args <- list(...)
    t <- format(Sys.time(), format = "%Y-%M-%D %H:%M:%S")
    if (verbose) {
        switch(EXPR = type,
        "in_assay" = cat(paste(t,
            "===>", args$comp_type, "in", args$assay, "<===\n")),
        "pca_tensor" = cat(paste(t,
            " Running Principal Component Analysis \n")),
        "nmf_tensor" = cat(paste(t,
           " Running Non-Negative Matrix Factorization \n")),
        "pca_rgb_tensor" = cat(paste(t,
            " Converting PCA Embedding Values to gray scale \n")),
        "pcal_rgb_tensor" = cat(paste(t,
            " Converting Loading Values to gray scale in PC", args$pc, "\n")),
        "svd_tensor" = cat(paste(t,
            " Running Single Value Decomposition \n")),
        "svd_rgb_tensor" = cat(paste(t,
            " Converting Embeddings to gray scale \n")),
        "umap_rgb_tensor" = cat(paste(t,
            " Converting UMAP Embeddings to gray scale \n")),
        "distance_beads" = cat(paste(t,
            " Filtering outlier beads \n")),
        "tess" = cat(paste(t,
            " Generating Voronoi Tesselation \n")),
        "raster" = cat(paste(t,
            " Rasterising Tiles \n")),
        "tensor_res" = cat(paste(t,
            " Reducing tensor resolution\n")),
        "filter_tiles" = cat(paste(t,
            " Filtering Tiles\n")),
        "adj_counts" = cat(paste(t,
            " Adjusting count matrix\n")),
        "reg" = cat(paste(t,
            " Regularising Image \n")),
        "seg" = cat(paste(t,
            " Segmenting Image using", args$method, "\n")),
        "smooth" = cat(paste(t,
            " Smoothing Image Arrays \n")),
        "ter_pool" = cat(paste(t,
            " Pooling Segment ", args$ter, "\n")),
        "eq" = cat(paste(t,
            " Equalizing Histogram \n")),
        "morph" = cat(paste(t,
            " Converting to pixset and morphing territory \n")),
        "layer" = cat(paste(t,
            " Generating layers\n")),
        "norm_check" = cat(paste(t,
            " Ignoring norm method for", args$method, "- Using Raw counts\n")),
        "deg_dispatch_all_null" = cat(paste(t,
            " No territory Specified - Comparing all territories\n")),
        "deg_prog_each" = cat(paste(t,
            "===>", args$seed, "VS", args$query, "<===", "\n")),
        "check_coord" = cat(paste(t,
            " Checking Coordinates in", args$assay, "\n")),
        "check_counts" = cat(paste(t,
            " Checking Counts in", args$assay, "\n")),
        "scale" = cat(paste(t,
            " Calculating Assay scale from coordinates\n")),
        "vtc" = cat(paste(t,
            " Converting Vesalius to Image\n")),
        "ctv" = cat(paste(t,
            " Converting Images to Vesalius\n")),
        "force_count" = cat(paste(t,
            " Force setting normalized counts as raw counts\n")),
        "add_counts" = cat(paste(t,
            " Adding counts to", args$assay, "\n")),
        "new_cells" = cat(paste(t,
            " Adding cell labels to",args$new_trial,"territory column in", args$assay,"\n")),
        "raw_count" = cat(paste(t,
            " Raw count matrix already present. Adding", args$count_type, "tag\n")),
        "add_embeds" = cat(paste(t,
            " Adding embeddings to", args$assay, "\n")),
        "signal" = cat(paste(t,
            " Extracting assay signal\n")),
        "custom_cost" = cat(paste(t,
            " Checking custom cost matrix",args$cost,"\n")),
        "mapping" = cat(paste(t,
            " Mapping query to seed - Epoch = ",args$epoch,"\n")),
        "feature_cost" = cat(paste(t,
            " Computing feature cost in", args$assay, "\n")),
        "neighbor_cost" = cat(paste(t,
            " Computing neighborhood cost in", args$assay, "\n")),
        "composition_cost" = cat(paste(t,
            " Computing niche composition index in", args$assay, "\n")),
        "territory_cost" = cat(paste(t,
            " Computing Territory cost in", args$assay, "\n")),
        "get_neigh" = cat(paste(t,
            " Getting Neighborhoods in", args$assay, "\n")),
        "cell_cost" = cat(paste(t,
            " Computing Cell Type cost in", args$assay, "\n")),
        "integrate" = cat(paste(t,
            " Intergrating Counts \n")),
        "integrated_embed" = cat(paste(t,
            " Setting",args$tag, "as active embedding \n")),
        "integrated_counts" = cat(paste(t,
            " Setting",args$tag, "as active count matrix \n")),
        "post_map_filter" = cat(paste(t,
            " Filtering Mapped Indices \n")),
        "count_only" = cat(paste(t,
            " No barcode maps found - Integrating counts only \n")),
        "overlap_scores" = cat(paste(t,
            " Computing mapping overlaps between query cells \n")),
        "hclust_scores" = cat(paste(t,
            " Applying hierarchical clustering to mapping scores \n")),
        "louvain_scores" = cat(paste(t,
            " Applying Louvain clustering to mapping scores \n")),
        "leiden_scores" = cat(paste(t,
            " Applying Leiden clustering to mapping scores \n"))
        )
    } else {
        return(NULL)
    }
}
