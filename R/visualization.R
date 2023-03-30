################################################################################
############################   ST package        ###############################
################################################################################

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

#' territory_plot - plotting Vesalius territories
#' @param vesalius_assay a vesalius_Assay object
#' @param trial character string describing which segmentation trial
#' to use. Default is "last" which is the last segmentation trial used.
#' @param split logical - If TRUE, territories will be plotted in
#' separate panels
#' @param randomise logical - If TRUE, colour palette will be randomised.
#' @param highlight numeric vector describing which territories should be 
#' highlighted.
#' @param cex numeric describing font size multiplier.
#' @param cex_pt numeric describing point size multiplier.
#' @details Territory plots show all territories in false colour after they
#' have been isolated from a Vesalius image.
#' 
#' Note that this function can be applied to image segments, territories,
#' and layers.
#' @return a ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple smoothing
#' ves <- smooth_image(ves, dimensions = seq(1, 30))
#' 
#' # quick segmentation
#' ves <- segment_image(ves, dimensions = seq(1, 30))
#' 
#' # isolate territories
#' ves <- isolate_territories(ves)
#' 
#' # Plot Territories
#' p <- territory_plot(ves)
#'}
#' @export
#' @importFrom ggplot2 ggplot geom_point aes facet_wrap theme_classic
#' @importFrom ggplot2 scale_color_manual theme element_text
#' @importFrom ggplot2 guides guide_legend labs

territory_plot <- function(vesalius_assay,
  trial = "last",
  split = FALSE,
  highlight = NULL,
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
    territories <- check_territory_trial(vesalius_assay, trial)
    if (any(territories$trial == 0)) {
      territories$trial <- territories$trial + 1
    }
    if (!is.null(highlight)){
        highlight <- check_group_value(territories, highlight)
    }
    
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
    ter_col <- create_palette(territories, randomise)
    ter_alpha <- create_alpha(territories, highlight)
    if (split) {
      ter_plot <- ggplot(territories, aes(x, y, col = trial)) +
          geom_point(size = cex_pt, alpha = ter_alpha) +
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
          geom_point(size = cex_pt, alpha = ter_alpha) +
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



#' create color palette from predefine scheme
#' @param territories vesalius territories taken from a vesalius_assay
#' @param randomise logical describing if colour palette should be 
#' randomised.
#' @details We use a predefined palette that use colour blind friendly
#' base colours. We generate a color palette based on the number of 
#' territories present. If required the colours will be randomly assinged
#' to each territory. Note that as the territory plot 
#' return a ggplot object, you can easily override the color scheme. 
#' @return color vector
#' @importFrom grDevices colorRampPalette
create_palette <- function(territories, randomise) {
  ter_col <- length(levels(territories$trial))
  base_colours <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999")
  if (ter_col < length(base_colours)) {
      ter_pal <- colorRampPalette(base_colours[seq(1, ter_col)])
  } else {
      ter_pal <- colorRampPalette(base_colours)
  }
  if (randomise) {
        ter_col <- sample(ter_pal(ter_col), ter_col)
  } else {
        ter_col <- ter_pal(ter_col)
  }
  return(ter_col)
}

#' create alpha value if territories need to be highlighted
#' @param territories vesalius territories taken from a vesalius_assay
#' @param highlight numeric vector describing which territories should 
#' be highlighted
#' @details If highlight is null, will return the same alpha values 
#' for all territories
#' @return vector of alpha values
create_alpha <- function(territories, highlight) {
  if (!is.null(highlight)){
    ter_col <- rep(0.25, length(levels(territories$trial)))
    loc <- as.character(levels(territories$trial)) %in% highlight
    ter_col[loc] <- 1
  } else {
    ter_col <- rep(0.65, length(levels(territories$trial)))
  }
  return(ter_col[as.integer(territories$trial)])
}

#' view_gene_expression
#' 
#' View gene expression in spatial omics data, in specific territories or
#' the expression of genes in a subset of cells.
#' @param vesalius_assay a vesalius_assay object
#' @param genes character vector containing the genes you wish to 
#' visualise.
#' @param norm_method character string - which count matrix should be used.
#' @param trial character string describing which segmentation trial
#' to use. Default is "last" which is the last trial used.
#' @param territory_1 integer or vector of integers descrbing territories in
#' group 1 (see details)
#' @param territory_2 integer or vector of integers descrbing territories in
#' group 2 (see details)
#' @param cells charactor vector containing barcodes/spatial_indices
#' associated with cell types of interest (see details)
#' @param norm logical indicating if expression should be min/max normalised
#' @param as_layer logical indicating if expression should represented as
#' a territory/ layer.
#' @param cex numeric - font size modulator
#' @param cex_pt numeric point size
#' @param alpha point transparency
#' @details Vesalius offers a plotting function that allows you to 
#' visualise gene expression. 
#' 
#' This function offers multiple options depending on what you provide.
#' 
#' 1. Overall expression
#' You can visualise the overall expression pattern of a set of genes
#' by providing a vesalius_assay object containing counts.
#' If \code{as_layer} is set to FALSE this will show the expression at
#' each sptial index indivdually. If set to TRUE, this will show the
#' average expression on a gene in all territories present.
#' 
#' 2. Expression in a territory
#' You can visualise the expression of a gene in an isolated territory.
#' 
#' 3. Expression of cells in one or more territory
#' If you want to visualise the expression of specific cells, you can
#' parse a character vector containing your cells of interest. This
#' function will automatically subset the relevant territory data and
#' show the expression only in the spatial indeces hat are associated
#' with your cell type of interest. You can use this option to contrast
#' the expression of cells between territories by also providing which
#' territories you wish to contrast (`territory_1` and `territory_2`).
#' If only a single territory is provided, vesalius will only shows cell
#' in that territory.
#' 
#' If you provide more than one gene, `view_gene_expression` will return 
#' a ggarrange list containing all your genes as individual plots. 
#' @return a ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple smoothing
#' ves <- smooth_image(ves, dimensions = seq(1, 30))
#' 
#' # quick segmentation
#' ves <- segment_image(ves, dimensions = seq(1, 30))
#' 
#' # isolate territories
#' ves <- isolate_territories(ves)
#' 
#' # view over all expression
#' p <- view_gene_expression(ves, genes = "Malat1")
#' p1 <- view_gene_expression(ves, genes = "Malat1", as_layer = TRUE)
#' 
#' # view expression in isolated territory 
#' p2 <- view_gene_expression(ves, genes = "Malat1", territory_1 = 5)
#' 
#' # view expression of cells
#' cells <- sample(colnames(get_counts(ves)),300)
#' p3 <- view_gene_expression(ves,
#'  genes = "Malat",
#'  cells = cells,
#'  territory_1 = 5,
#'  terriotry_2 = 8)
#'}
#' @export
#' @importFrom infix %||%
#' @importFrom dplyr left_join filter group_by mutate
#' @importFrom stats na.exclude
#' @importFrom ggplot2 ggplot geom_point aes facet_wrap theme_classic
#' @importFrom ggplot2 scale_color_gradientn theme element_text
#' @importFrom ggplot2 guides guide_legend labs margin
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
view_gene_expression <- function(vesalius_assay,
  genes = NULL,
  norm_method = "last",
  trial = "last",
  territory_1 = NULL,
  territory_2 = NULL,
  cells = NULL,
  norm = TRUE,
  as_layer = FALSE,
  with_background = FALSE,
  cex = 10,
  cex_pt = 1,
  alpha = 0.75) {
    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    counts <- check_norm(vesalius_assay, norm_method)
    territories <- check_territory_trial(vesalius_assay, trial) %>%
      dispatch_territory(territory_1, territory_2, cells)
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
      gene_list[[i]] <- gene
    }
    #--------------------------------------------------------------------------#
    # Now we can loop of gene list
    #--------------------------------------------------------------------------#
    for (i in seq_along(gene_list)) {
      #------------------------------------------------------------------------#
      # Extract data and convert count to mean count if as_layer=T
      #------------------------------------------------------------------------#
      if (is.null(gene_list[[i]])) next()

      other <- filter(gene_list[[i]], trial == "other") %>% droplevels()
      territory <- filter(gene_list[[i]], trial != "other") %>% droplevels()
      if (as_layer) {
          territory <- territory %>%
            group_by(trial) %>%
            mutate(gene = mean(gene)) %>%
            ungroup()
      }
      #------------------------------------------------------------------------#
      # checking if layer and normalising post layering if layering required 
      #------------------------------------------------------------------------#
      if (norm) {
        territory$gene <- min_max(territory$gene)
        type <- "Norm. Expression"
      } else {
        type <- "Expression"
      }
      #------------------------------------------------------------------------#
      # Create background
      #------------------------------------------------------------------------#
      if (nrow(other) > 0) {
        ge <- ggplot() +
          geom_point(data = other, aes(x, y),
            alpha = alpha,
            fill = "#2e2e2e",
            size = cex_pt,
            show.legend = TRUE)
      } else {
        ge <- ggplot()
      }
      #------------------------------------------------------------------------#
      # Create ggplot with foreground
      #------------------------------------------------------------------------#
      if (is.null(cells)) {
        ge <- ge + geom_point(data = territory,
          aes(x = x, y = y, col = gene), size = cex_pt * 2, alpha = alpha)
      } else {
        ge <- ge + geom_point(data = territory,
          aes(x = x, y = y, shape = trial, col = gene),
          size = cex_pt * 2,
          alpha = alpha)
      }

      ge <- ge +
        scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
        theme_classic() +
        theme(legend.title = element_text(size = cex),
          legend.text = element_text(size = cex),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
        labs(col = type, title = names(gene_list)[i]) +
        guides(shape = guide_legend(override.aes = list(size = cex * 0.5)))
      
      if (with_background){
        ge <- ge +
          theme(panel.background = element_rect(
            fill = rev(brewer.pal(11, "Spectral"))[1]))
      }
      gene_list[[i]] <- ge
    }
    #--------------------------------------------------------------------------#
    # returning ggplot or ggplot list with ggarrange
    #--------------------------------------------------------------------------#
    if (is.null(gene_list)) {
      return(NULL)
    } else {
      gene_list <- gene_list[!sapply(gene_list, is.null)]
    }

    if (length(gene_list) == 1) {
      gene_list <- gene_list[[1L]]
    } else {
      gene_list <- ggarrange(plotlist = gene_list,
        ncol = floor(sqrt(length(gene_list))),
        nrow = floor(sqrt(length(gene_list))))
    }
    return(gene_list)
}

