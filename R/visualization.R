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
#' @param contour if territory contours should be added. Availble:
#' "None", "convex", "concave"
#' @param cex numeric describing font size multiplier.
#' @param cex_pt numeric describing point size multiplier.
#' @param alpha opacity factor ]0,1[
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
#' @importFrom ggplot2 guides guide_legend labs geom_polygon
#' @importFrom ggplot2 scale_fill_manual scale_x_continuous scale_y_continuous
#' @importFrom imager imrotate
#' @importFrom ggnewscale new_scale

territory_plot <- function(vesalius_assay,
  trial = "last",
  split = FALSE,
  highlight = NULL,
  contour = "None",
  randomise = TRUE,
  cex = 10,
  cex_pt = 1,
  alpha = 0.65,
  use_image = FALSE) {
    #--------------------------------------------------------------------------#
    # Dirty ggplot - this is just a quick and dirty plot to show what it look
    # like
    # At the moment - I thinking the user should make their own
    # Not a prority for custom plotting functions
    # SANITY check here and format
    #--------------------------------------------------------------------------#
    territories <- check_territory_trial(vesalius_assay, trial)
    tiles <- get_tiles(vesalius_assay)
    # if (any(territories$trial == 0)) {
    #   territories$trial <- territories$trial + 1
    # }
    if (!is.null(highlight)){
        highlight <- check_group_value(territories, highlight)
    }
    
    if (use_image) {
        img <- vesalius_assay@image
        if (length(img) == 0) {
          warning(paste0("No Image found in ",
            get_assay_names(vesalius_assay)))
          img <- NULL
        } else {
          #-------------------------------------------------------------------#
          # NOTE: not ideal - i made it a list to handle multiple images 
          # will need to think about hwo to handle this data 
          # especially when it comes to interoperable object
          #-------------------------------------------------------------------#
          img <- imager::imrotate(img[[1]], 90)
          img <- as.data.frame(img, wide = "c") %>%
            mutate(rgb_val = rgb(c.1, c.2, c.3))
          img$x <- rev(img$x)
          territories <- adjust_cooridnates(territories, vesalius_assay)
          #tiles <- adjust_cooridnates(tiles, vesalius_assay)
        }
    } else {
       img <- NULL
    }
    if (as_contour <- !grepl("none|None|NONE", contour)) {
      territories <- unpack_territory_path(territories,
        tiles,
        method = contour)
    }
    
    legend <- sapply(strsplit(trial, "_"), "[[", 1)
    #--------------------------------------------------------------------------#
    # Changing label order because factor can suck ass sometimes
    #--------------------------------------------------------------------------#
    sorted_labels <- order_labels(territories)

    territories$trial <- factor(territories$trial,
      levels = sorted_labels)
    
    #--------------------------------------------------------------------------#
    # My pure hatred for the standard ggplot rainbow colours has forced me
    # to use this palette instead - Sorry Hadely
    #--------------------------------------------------------------------------#
    ter_col <- create_palette(territories, randomise)
    ter_alpha <- create_alpha(territories, highlight, alpha)

    ter_plot <- ggplot()
    if (!is.null(img)) {
      ter_plot <- ter_plot + geom_raster(data = img,
        aes(x = x, y = y, fill = rgb_val)) +
        scale_fill_identity() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        new_scale("fill")
    }
    if (as_contour) {
      ter_plot <- ter_plot +
        geom_polygon(data = territories,
          aes(x, y, fill = trial),
          colour = ter_col[territories$trial],
          alpha = ter_alpha,
          size = cex_pt) +
        scale_fill_manual(values = ter_col)
    } else {
      ter_plot <- ter_plot +
        geom_point(data = territories,
          aes(x, y, col = trial),
          size = cex_pt,
          alpha = ter_alpha)+
        scale_color_manual(values = ter_col)
    }
    if (split) {
      ter_plot <- ter_plot + facet_wrap(~trial)
    }
    ter_plot <- ter_plot +
      theme_classic() +
      theme(legend.text = element_text(size = cex * 1.2),
        axis.text = element_text(size = cex * 1.2),
        axis.title = element_text(size = cex * 1.2),
        plot.title = element_text(size = cex * 1.5),
        legend.title = element_text(size = cex * 1.2)) +
      guides(colour = guide_legend(
        override.aes = list(size = cex * 0.3))) +
      labs(colour = legend, title = paste("Vesalius", trial),
        x = "X coordinates", y = "Y coordinates")
    return(ter_plot)
}

order_labels <- function(territories) {
      ter <- unique(territories$trial)
      iso <- ter[ter == "isolated"]
      tmp <- ter[ter != "isolated"]
      if (all(is.na(suppressWarnings(as.numeric(tmp))))) {
         labels <- c(sort(tmp), iso)
      } else {
         labels <- c(sort(as.numeric(tmp)), iso)
      }
      return(labels)
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
#' @param alpha tranaparent factor
#' @details If highlight is null, will return the same alpha values 
#' for all territories
#' @return vector of alpha values
create_alpha <- function(territories, highlight, alpha) {
  if (!is.null(highlight)) {
    ter_col <- rep(alpha * 0.25, length(levels(territories$trial)))
    loc <- as.character(levels(territories$trial)) %in% highlight
    ter_col[loc] <- alpha
  } else {
    ter_col <- rep(alpha, length(levels(territories$trial)))
  }
  return(ter_col[as.integer(territories$trial)])
}

#' retrieve the points contained in the edge of each territory
#' @param trial name of territory trial that should be selected
#' @param tiles vesalius tiles 
#' @param method string - method for how territories edges should be
#' selected
#' @details Here we are using convex as start point. Essentially, we 
#' order the coordinates based on their polar coordinates using the 
#' median coordinate as the center point. 
#' @returns a data frame containing edge of each territory.
unpack_territory_path <- function(trial,
    tiles,
    method = "none") {
    #-------------------------------------------------------------------------#
    # First we convert to pixset and detect territory edge 
    #-------------------------------------------------------------------------#
    trial_split <- vector("list", length(unique(trial$trial)))
    names(trial_split) <- unique(trial$trial)
    for (i in seq_along(unique(trial$trial))) {
        territory <- unique(trial$trial)[i]
        path <- switch(method,
          "none" = trial,
          "edge" = territory_edge(trial, tiles, territory),
          "concave" = territory_concave(trial, territory),
          "convex" = territory_convex(trial, territory))
        trial_split[[i]] <- path
    }

    #-------------------------------------------------------------------------#
    # next we remove NULLs - this happens when no edge can be deteced 
    #-------------------------------------------------------------------------#
    nulls <- sapply(trial_split, nrow) == 0
    trial_split <- trial_split[!nulls]
    if (sum(nulls) > 0) {
      warning(paste(names(trial_split)[nulls], "could not be convert to path!
       Skiiping territories in plot\n"))
    } else if (length(trial_split) == 0) {
        stop("No edge can be detect in territories! Granularity too high.
         Consider increasing smoothing and/or decreasing segmentationd depth")
    }
    #-------------------------------------------------------------------------#
    # select starting point for path 
    # ATM we convert edge to ordered shape using polar coordinates 
    #-------------------------------------------------------------------------#
    trial <- lapply(trial_split, function(trial) {
            ord <- convexify(trial$x,
                trial$y,
                median(trial$x),
                median(trial$y),
                order = TRUE)
            trial <- trial[ord$x, ]
            return(trial)
        })
    trial <- do.call("rbind", trial)
    trial <- trial %>% filter(trial != "isolated")
    return(trial)
}

territory_edge <- function(trial, tiles, territory) {
  ter <- right_join(trial, tiles, by = "barcodes") %>%
    filter(trial %in% territory) %>%
    mutate(value = 1) %>%
    select(c("barcodes", "x.y", "y.y", "value", "origin", "trial"))
  colnames(ter) <- c("barcodes", "x", "y", "value", "origin", "trial")
  edge <- extend_boundary(ter, c(10)) %>%
    detect_edges() %>%
    grow(1) %>%
    as.data.frame()
  edge <- inner_join(edge, ter, by = c("x", "y")) %>%
    select("barcodes") %>% 
    unique()
  edge <- tiles %>% filter(barcodes %in% edge$barcodes & origin == 1)
  edge$trial <- territory
  edge <- rbind(edge, edge[1, ])
  return(edge)
}

#' @importFrom grDevices chull
territory_convex <- function(trial, territory) {
  trial <- trial %>%
    filter(trial %in% territory)
  hull <- chull(trial$x, trial$y)
  trial <- trial[hull, ]
  trial <- rbind(trial, trial[1, ])
  return(trial)
}


# #' @importFrom concaveman concaveman
# territory_concave <- function(trial, territory) {
#   trial <- trial %>%
#     filter(trial %in% territory)
#   bars <- trial$barcodes
#   trial_mat <- as.matrix(trial[, c("x", "y")])
#   rownames(trial_mat) <- bars
#   hull <- as.data.frame(concaveman(trial_mat, concavity = 3))
#   colnames(hull) <- c("x", "y")
#   knn <- RANN::nn2(data = trial_mat,
#     query = hull,
#     k = 1)
#   trial <- trial[knn$nn.idx[, 1], ]
#   return(trial)
# }

adjust_cooridnates <- function(trial, vesalius_assay) {
    orig_coord <- vesalius_assay@meta$orig_coord
    #-------------------------------------------------------------------------#
    # First let's split barcodes between adjusted and not 
    #-------------------------------------------------------------------------#
    adj_barcodes <- grep("_et_", trial$barcodes, value = TRUE)
    non_adj_barcodes <- trial$barcodes[!trial$barcodes %in% adj_barcodes]
    #-------------------------------------------------------------------------#
    # next get original coordinates and match them in trial for non adjusted
    #-------------------------------------------------------------------------#
    in_trial <- match(orig_coord$barcodes, non_adj_barcodes) %>%
      na.exclude()
    in_orig <- match(non_adj_barcodes, orig_coord$barcodes) %>%
      na.exclude()
    
    trial$x[in_trial] <- orig_coord$x_orig[in_orig]
    trial$y[in_trial] <- orig_coord$y_orig[in_orig]
    #-------------------------------------------------------------------------#
    # unpack adjusted 
    #-------------------------------------------------------------------------#
    adj_barcodes_sp <- split(adj_barcodes, "_et_")
    for (i in seq_along(adj_barcodes_sp)) {
        x <- orig_coord[orig_coord$barcodes %in% adj_barcodes_sp[[i]],
          "x_orig"]
        y <- orig_coord[orig_coord$barcodes %in% adj_barcodes_sp[[i]],
          "y_orig"]
        trial$x[trial$barcodes == adj_barcodes[i]] <- median(x)
        trial$y[trial$barcodes == adj_barcodes[i]] <- median(y)
    }
    #-------------------------------------------------------------------------#
    # adjuste using scale - note!!! This is dodgy as fuck!
    # mainly because my original use of scale was intended to be used this way
    # Using this like this since it is faster for ad hoc analysis 
    #-------------------------------------------------------------------------#
    scale <- vesalius_assay@meta$scale$scale
    trial$x <- trial$x * scale
    trial$y <- trial$y * scale
    return(trial)

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
#' @param return_as_list logical - should plot be returned as simple list
#' or as a ggplot object (single gene)/ patchwork object (multiple genes)
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
#' @importFrom ggplot2 guides guide_legend labs margin element_rect
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
#' @importFrom patchwork wrap_plots
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
  alpha = 0.75,
  max_size = 5,
  return_as_list = FALSE) {
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

    #----------------------------------------------------------------------#
    # if nothing to plot
    #----------------------------------------------------------------------#
    nulls <- sapply(gene_list, is.null)
    if (all(nulls)) {
        return(NULL)
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
    } else if (!return_as_list) {
      gene_list <- ggarrange(plotlist = gene_list,
        ncol = min(c(floor(sqrt(length(gene_list))), max_size)),
        nrow = min(c(floor(sqrt(length(gene_list))), max_size)))
    }
    return(gene_list)
}


#' @export 
view_mapping_score <- function(vesalius_assay,
    score = "total_score",
    cex_pt = 1,
    cex = 15) {
    coord <- get_tiles(vesalius_assay) %>%
        filter(origin == 1)
    scores <- vesalius_assay@meta$mapping_probability
    scores <- scores[match(coord$barcodes, scores$from), score]
    coord$score <- scores
    g <- ggplot(coord, aes(x,y, col = score)) +
        geom_point(size = cex_pt) +
        scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")),
            limits = c(0,1)) +
        theme_classic() +
        theme(legend.title = element_text(size = cex),
          legend.text = element_text(size = cex),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
        labs(col = "Mapping score", title = score) +
        guides(shape = guide_legend(override.aes = list(size = cex * 0.5)))
    return(g)
}

