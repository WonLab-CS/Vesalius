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
    t <- Sys.time()
    if (verbose) {
        switch(EXPR = type,
        "in_assay" = cat(paste(t,
            "===>", args$comp_type, "in", args$assay, "<===\n")),
        "pca_tensor" = cat(paste(t,
            " Running Principal Component Analysis \n")),
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
            " Ignoring norm method for", args$method, "Using Raw counts\n")),
        "deg_dispatch_all_null" = cat(paste(t,
            " No territory Spacified - Comparing all territories\n")),
        "deg_prog_each" = cat(paste(t,
            "===>", args$seed, "VS", args$query, "<===", "\n")),
        "check_coord" = cat(paste(t,
            " Checking Coordinates in", args$assay, "\n")),
        "check_counts" = cat(paste(t,
            " Checking Counts in", args$assay, "\n")),
        "vtc" = cat(paste(t,
            " Converting Vesalius to Image\n")),
        "ctv" = cat(paste(t,
            " Converting Images to Vesalius\n")),
        "force_count" = cat(paste(t,
            " Force setting normalized counts as raw counts\n")),
        "add_counts" = cat(paste(t,
            " Adding counts to", args$assay, "\n")),
        "raw_count" = cat(paste(t,
            " Raw count matrix already present. Adding", args$count_type, "tag\n")),
        "add_embeds" = cat(paste(t,
            " Adding embeddings to", args$assay, "\n"))
        )
    } else {
        return(NULL)
    }
}
