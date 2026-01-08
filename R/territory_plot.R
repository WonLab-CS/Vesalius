################################################################################
############################   ST package        ###############################
################################################################################

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
#' @param use_image logical - whether to use image background
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


#need to rework this one
order_labels <- function(territories) {
      ter <- unique(territories$trial)
      iso <- ter[grepl(pattern = "isolated|Not Selected",x = ter)]
      tmp <- ter[!grepl(pattern = "isolated|Not Selected",x = ter)]
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
    if (length(alpha) > 1){
        ter_col <- rep(alpha[2L], length(levels(territories$trial)))
        loc <- as.character(levels(territories$trial)) %in% highlight
        ter_col[loc] <- alpha[1L]
    } else {
        ter_col <- rep(alpha[1L] * 0.15, length(levels(territories$trial)))
        loc <- as.character(levels(territories$trial)) %in% highlight
        ter_col[loc] <- alpha[1L]
    }
    
  } else {
    ter_col <- rep(alpha[1L], length(levels(territories$trial)))
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



