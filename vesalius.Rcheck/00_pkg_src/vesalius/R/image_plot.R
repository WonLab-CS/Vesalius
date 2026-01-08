
#' image_plot - plotting vesalius embeddings
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions which dimensions to use to generate image (see details)
#' @param embedding character string descrining which embedding should be 
#' used for image generation.
#' @param cex numeric - font and point size resizing factor.
#' @details Once you have generated your embeddings, 
#' you can visualise these embeddings using the \code{image_plot} function.
#' This function will generate a ggplot object representing the embedding 
#' image containing all pixels. You can select any dimension and in any 
#' combination you desire, however you can only select 1 or 3 dimensions 
#' to visualise at a time. This will either generate grey scale image 
#' or RGB images. 
#' 
#' This function is always applied to the active embedding. By default,
#' this is the last you generated. This means that you can also use
#' this function to visualise the effect if smoothing, equalization,
#' regularisation or segmentation. 
#' 
#' Please note that if you are using "louvain" or "leiden" for 
#' image segmentation, this will always generate grey scale images
#' even if you select multiple dimensions. "Kmeans" on the other hand 
#' will still produce RGB color segments.
#' 
#' The usage of this function remains the same in after any processing steps.
#' @return ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run 
#' ves <- build_vesalius_embeddings(ves)
#' # Plot 1st 3 PCs
#' p <- image_plot(ves)
#'
#'}
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @importFrom dplyr right_join
#' @importFrom stats na.exclude
#' @importFrom grDevices rgb gray

image_plot <- function(vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  cex = 10) {
    #--------------------------------------------------------------------------#
    # Checking dims - only gray scale or 3 coloured images for now
    # What approach could we take here to use different number of dims?
    #--------------------------------------------------------------------------#
    if (length(dimensions) > 3) {
       stop("To many dims selected")
    } else if (length(dimensions) != 3 && length(dimensions) != 1) {
       stop("Non RGB /gray scale images not supported")
    }
    #--------------------------------------------------------------------------#
    # First get image data from vesalius object
    #--------------------------------------------------------------------------#
    coordinates <- check_tiles(vesalius_assay)
    #--------------------------------------------------------------------------#
    # get last embedding
    #--------------------------------------------------------------------------#
    tile_colour <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    if (embedding == "last") {
      embed_name <- get_active_embedding_tag(vesalius_assay)
    } else {
      embed_name <- embedding
    }
    #--------------------------------------------------------------------------#
    # reformat to a data frame for easy use with ggplot
    # rebalence colors - min-max norm colors jjust in case 
    # and reformats to clean data frame for hex color generation
    #--------------------------------------------------------------------------#
    tile_colour <- as.data.frame(tile_colour)
    tile_colour$barcodes <- rownames(tile_colour)
    coordinates <- right_join(coordinates, tile_colour, by = "barcodes") %>%
      na.exclude()
    coordinates <- rebalence_colors(coordinates,
      length(dimensions),
      method = "minmax")
    #--------------------------------------------------------------------------#
    # Generate colors and plots
    #--------------------------------------------------------------------------#
    if (length(dimensions) == 3) {
      cols <- rgb(coordinates$R, coordinates$G, coordinates$B)
      title <- paste0(embed_name, " - Dims = ",
        paste0(dimensions, collapse = "-"))
    } else {
      # to do
      cols <- gray(coordinates$Gray)
      title <- paste0(embed_name, " - Dim = ",
        paste0(dimensions, collapse = "-"))
    }
    image <- ggplot(coordinates, aes(x, y)) +
        geom_raster(aes(fill = cols)) +
        scale_fill_identity() +
        theme_classic() +
        theme(axis.text = element_text(size = cex * 1.2),
          axis.title = element_text(size = cex * 1.2),
          plot.title = element_text(size = cex * 1.5)) +
        labs(title = title, x = "X coordinates", y = "Y coordinates")
    return(image)
}


#' rebalence colors
#' @param coordinates data frame containing barcodes, x/y coord, origin,
#' and color value. 
#' @param dimensions number dimensions select for plotting
#' @param method character string: min max or truncate
#' @details This function is use to re-bound values between 0 and 1. 
#' Some image processing steps may lead to negative values being introduced.
#' Imager handles this without an issue but for plotting we want to make sure 
#' that all values are indeed bound between 0 and 1. We will try to different 
#' methods either we truncate the value or we min max norm. 
#' @return data frame with [0,1] bound values
#' @importFrom dplyr select
rebalence_colors <- function(coordinates, dimensions, method = "minmax") {
    if (dimensions == 3) {
        template <- select(coordinates, c("barcodes", "x", "y", "origin"))
        colors <- coordinates[, c(5, 6, 7)]
        if (method == "minmax") {
            colors <- apply(colors, 2, min_max)
        } else {
            colors[which(colors < 0, arr.ind = TRUE)] <- 0
            colors[which(colors > 1, arr.ind = TRUE)] <- 1
        }
        template <- data.frame(template, colors)
        colnames(template) <- c("barcodes", "x", "y", "origin", "R", "G", "B")
    } else {
        template <- select(coordinates, c("barcodes", "x", "y", "origin"))
        colors <- coordinates[, c(5)]
        if (method == "minmax") {
            colors <- min_max(colors)
        } else {
            colors[which(colors < 0)] <- 0
            colors[which(colors > 1)] <- 1
        }
        template <- data.frame(template, colors)
        colnames(template) <- c("barcodes", "x", "y", "origin", "Gray")
    }
    return(template)
}