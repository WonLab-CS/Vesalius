################################################################################
############################   ST package        ###############################
################################################################################

#' imagePlot - Vesalius image plotting
#' @param image a Vesalius data frame containing at least
#' barcodes, x, y, cc, value
#' @param as.cimg logical - If TRUE, will plot directly using cimg plotting.
#' If FALSE, will produce a ggplot object.
#' @param cex numeric - font and point size resizing factor.
#' @details Vesalius provides basic plotting functions. This function
#' only shows Vesalius images (i.e. using true embedded colours or segmented
#' colours).
#'
#' Currently, we advise to use as.cimg as false for custom plots.
#'
#' As Vesalius nearly always returns data frames, custom plots are left to
#' the discretion of the user.
#' @return cimg plot or ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' imagePlot(image)
#' # as ggplot
#' g <- imagePlot(image,as.cimg = F)
#' }
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
    coordinates <- check_tiles(vesalius)
    # need to check this -
    tile_colour <- check_embedding(vesalius,
      embedding,
      dimensions)[, dimensions]
    #### format needs to be reworked here
    tile_colour <- as.data.frame(tile_colour)
    tile_colour$barcodes <- rownames(tile_colour)
    coordinates <- right_join(coordinates, tile_colour, by = "barcodes") %>%
      na.exclude()
    #--------------------------------------------------------------------------#
    # Generate plots
    # Ceck to see how you can avoid the hard coding of the dimensions
    #--------------------------------------------------------------------------#
    if (length(dimensions) == 3) {
      cols <- rgb(coordinates[, 5], coordinates[, 6], coordinates[, 7])
      title <- paste0("Vesalius - Dims = ", paste0(dimensions, collapse = "-"))
    } else {
      # to do
      cols <- gray(coordinates[, 5])
      title <- paste0("Vesalius - Dim = ", paste0(dimensions, collapse = "-"))
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




#' territoryPlot plotting Vesalius territories
#' @param territories a Vesalius data frame conatining territory information
#' (i.e. containing the terriotry column)
#' @param split logical - If TRUE, territories will be plotted in
#' separate panels
#' @param randomise logical - If TRUE, colour palette will be randomised.
#' @param cex numeric describing font size multiplier.
#' @param cex.pt numeric describing point size multiplier.
#' @details Territory plots show all territories in false colour after they
#' have been isolated from a Vesalius image.
#' @return a ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' g <- territoryPlot(image,cex = 12, cex.pt = 1)
#' }
#' @export

territory_plot <- function(vesalius_assay,
  trial = "last",
  split = FALSE,
  randomise = TRUE,
  cex = 10,
  cex_pt = 1) {
    #--------------------------------------------------------------------------#
    # Dirty ggplot - this is just a quick and dirty plot to show what it look
    # like
    # At the moment - I thinking the user should make their own
    # Not a prority for custom plotting functions
    # SANITY check here and format
    #--------------------------------------------------------------------------#
    territories <- check_territories(vesalius_assay, trial)
    legend <- sapply(strsplit(trial, "_"), "[[", 1)
    #--------------------------------------------------------------------------#
    # Changing label order because factor can suck ass sometimes
    #--------------------------------------------------------------------------#
    sorted_labels <- order(levels(as.factor(territories$trial)))
    if (any(grepl("isolated", territories$trial))) {
      sorted_labels[length(sorted_labels)] <- "isolated"
    }

    territories$trial <- factor(territories$trial,
      levels = sorted_labels)
    #--------------------------------------------------------------------------#
    # My pure hatred for the standard ggplot rainbow colours has forced me
    # to use this palette instead - Sorry Hadely
    #--------------------------------------------------------------------------#
    ter_col <- length(levels(territories$trial))
    ter_pal <- colorRampPalette(c("#999999",
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7"))

    if (randomise) {
        ter_col <- sample(ter_pal(ter_col), ter_col)
    } else {
        ter_col <- ter_pal(ter_col)
    }
    if (split) {
      ter_plot <- ggplot(territories, aes(x, y, col = trial)) +
          geom_point(size = cex_pt, alpha = 0.65) +
          facet_wrap(~trial) +
          theme_classic() +
          scale_color_manual(values = ter_col) +
          theme(legend.text = element_text(size = cex * 1.2),
            axis.text = element_text(size = cex * 1.2),
            axis.title = element_text(size = cex * 1.2),
            plot.title = element_text(size = cex * 1.5),
            legend.title = element_text(size = cex * 1.2)) +
          guides(colour = guide_legend(
            override.aes = list(size = cex * 0.3))) +
          labs(colour = legend, title = paste("Vesalius", trial),
            x = "X coordinates", y = "Y coordinates")
    } else {
      ter_plot <- ggplot(territories, aes(x, y, col = trial)) +
          geom_point(size = cex_pt, alpha = 0.65) +
          theme_classic() +
          scale_color_manual(values = ter_col) +
          theme(legend.text = element_text(size = cex * 1.2),
            axis.text = element_text(size = cex * 1.2),
            axis.title = element_text(size = cex * 1.2),
            plot.title = element_text(size = cex * 1.5),
            legend.title = element_text(size = cex * 1.2)) +
          guides(colour = guide_legend(
            override.aes = list(size = cex * 0.3))) +
          labs(colour = legend, title = paste("Vesalius", trial),
            x = "X coordinates", y = "Y coordinates")
    }
    return(ter_plot)
}




#' viewGeneExpression - plot gene expression.
#' @param image a Vesalius data frame containing barcodes, x, y, cc, value,
#' cluster, and territory.
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param ter integer - Territory ID in which gene expression will be viewed. If
#' NULL, all territories will be selected.
#' @param genes character - gene of interest (only on gene at a time)
#' @param normalise logical - If TRUE, gene expression values will be min/max
#' normalised.
#' @param cex numeric - font size modulator
#' @details Visualization of gene expression in all or in selected territories.
#' Gene expression is shown "as is". This means that if no transformation
#' is applied to the data then normalized raw count will be shown.
#'
#' If normalization and scaling have been applied, normalized counts will be
#' shown.
#'
#' This also applies to any data transformation applied on the count matrix
#' or the Seurat object.
#'
#' We recommend using Seurat Objects.
#'
#' NOTE : You might be promoted to use geom_tiles instead of geom_raster.
#' However - this will increase file size to unreasonable heights.
#' @return a ggplot object (geom_raster)
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' # In all points
#' g <- viewGeneExpression(image,vesalius, genes = "Cst3")
#' # In a specific territory
#' g1 <- viewGeneExpression(image, vesalius, ter = 1, genes = "Cst3")
#' }
#' @export


view_gene_expression <- function(vesalius,
  genes = NULL,
  norm_method = "last",
  trial = "last",
  territory_1 = NULL,
  territory_2 = NULL,
  cells = NULL,
  norm = TRUE,
  as_layer = FALSE,
  cex = 10) {
    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    counts <- check_norm(vesalius, norm_method)
    territories <- check_territories(vesalius, trial) %>%
      territory_dispatch(territory_1, territory_2, cells)
    #--------------------------------------------------------------------------#
    # Getting genes
    # We will facet wrap if there is more than one.
    #--------------------------------------------------------------------------#
    genes <- genes %||%
      stop("Please specifiy which gene(s) you would like to visualize")
    #----------------------------------------------------------------------#
    # First lets make a list to store our counts and then plots
    #----------------------------------------------------------------------#
    gene_list <- vector("list", length(genes))
    names(gene_list) <- genes
    #----------------------------------------------------------------------#
    # Next we loop,extract, combine
    #----------------------------------------------------------------------#
    for (i in seq_along(genes)) {
      gene <- rownames(counts) == genes[i]
      if (sum(gene) == 0) {
        warning(paste(genes[i], "is not present in count matrix.
          Returning NULL"), immediate. = TRUE)
        gene_list[[i]] <- NULL
        next()
      }
      gene <- counts[gene, ]
      gene <- data.frame(names(gene), gene)
      colnames(gene) <- c("barcodes", "gene")
      gene <- left_join(gene, territories, by = c("barcodes")) %>%
        na.exclude()
      if (norm) {
        gene$gene <- min_max(gene$gene)
        type <- "Norm. Expression"
      } else {
        type <- "Expression"
      }
      gene_list[[i]] <- gene
    }
    #--------------------------------------------------------------------------#
    # Now we can loop of gene list
    #--------------------------------------------------------------------------#
    for (i in seq_along(gene_list)) {
      #------------------------------------------------------------------------#
      # Extract data and convert count to mean count if as_layer=T
      #------------------------------------------------------------------------#
      other <- filter(gene_list[[i]], trial == "other") %>% droplevels()
      territory <- filter(gene_list[[i]], trial != "other") %>% droplevels()

      if (as_layer) {
          territory <- territory %>%
            group_by(trial) %>%
            mutate(gene = mean(gene))
      }
      #------------------------------------------------------------------------#
      # Create background
      #------------------------------------------------------------------------#
      if (nrow(other) > 0) {
        ge <- ggplot() +
          geom_point(data = other, aes(x, y),
            alpha = 0.75,
            fill = "#2e2e2e",
            size = cex * 0.01,
            show.legend = TRUE)
      } else {
        ge <- ggplot()
      }
      #------------------------------------------------------------------------#
      # Create ggplot with foreground
      #------------------------------------------------------------------------#
      if (is.null(cells)) {
        ge <- ge + geom_point(data = territory,
          aes(x = x, y = y, col = gene), size = cex * 0.1)
      } else {
        ge <- ge + geom_point(data = territory,
          aes(x = x, y = y, shape = trial, col = gene), size = cex * 0.1)
      }

      ge <- ge +
        scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
        theme_classic() +
        theme(#panel.background = element_rect(fill = "#474747"),
          legend.title = element_text(size = cex),
          legend.text = element_text(size = cex),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
        labs(col = type, title = names(gene_list)[i]) +
        guides(shape = guide_legend(override.aes = list(size = cex * 0.5)))

      gene_list[[i]] <- ge
    }
    #--------------------------------------------------------------------------#
    # returning ggplot or ggplot list with ggarrange
    #--------------------------------------------------------------------------#
    if (length(gene_list) == 1) {
      gene_list <- gene_list[[1L]]
    } else {
      gene_list <- ggarrange(plotlist = gene_list,
        ncol = floor(sqrt(length(gene_list))))
    }
    return(gene_list)
}







## assuming a list after chunking images
view_transform <- function(ves) {
    for (i in seq_along(ves)) {
        for (j in seq_along(ves[[i]])) {
            plot(ves[[i]][[j]]$img, main = "Image")
            plot(ves[[i]][[j]]$fft, main = "Imaginary")
            sqrt(ves[[i]][[j]]$real^2 + ves[[i]][[j]]$fft^2) %>%
              plot(main = "Power spectrum")
        }
    }
}
