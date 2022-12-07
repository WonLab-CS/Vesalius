################################################################################
############################   ST package        ###############################
################################################################################
# Territory morphology


#' territory_morphing applies morphological operators to a set of territoriees
#' @param vesalius_assay vesalius_assay object
#' @param territory integer or vector of integers desrcining territories 
#' to morph.
#' @param trial character string - which territory trial that 
#' should be used to select
#' territorires. Default is last one computed
#' @param morphology_factor integer or vector of integers describing growth
#' and/or shrink extent.
#' @param verbose logical - progress message output.
#' @details Territory morphing can manipulate territories by growing, shrinking,
#' filling, and cleaning territories.
#' Growing = Positive integers - Territory will be dilated by x number of pixels
#' Shrinking = Negative integers - Territory will be contracted by x number of
#' pixels
#' Filling = grow followed by shrink.
#' Cleaning = shrink followed by grow.
#' @return a vesalius_assay
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
#' # morph territory
#'
#' ves <- territory_morphing(ves, 8, morphology_factor = 30)
#' ves <- terriotry_morphing(ves, 1, morpholgy_factor = c(-15, 15))
#'
#' # view territory morphing
#' territory_plot(ves)
#'}
#' @export
#' @importFrom infix %||%
#' @importFrom dplyr filter left_join mutate select inner_join
#' @importFrom imager grow shrink

territory_morphing <- function(vesalius_assay,
  territory = NULL,
  trial = "last",
  morphology_factor = 0,
  verbose = TRUE) {
    simple_bar(verbose)
    #------------------------------------------------------------------------#
    # lets make a few checks
    # for now we will assume that either segments or territories can be used
    # TOADD methods to access parameters associated with each trial run
    #------------------------------------------------------------------------#
    territory <- territory %||%
      stop("No specified territory for territory morphing!")
    adjusted_territory_values <- min(territory)

    ter <- check_territories(vesalius_assay, trial)
    ter <- filter(ter, trial %in% territory) %>%
               left_join(get_tiles(vesalius_assay), by = "barcodes") %>%
               mutate(value = 1) %>%
               select(c("barcodes", "x.y", "y.y", "trial", "origin", "value"))
    colnames(ter) <- c("barcodes", "x", "y", "trial", "origin", "value")
    message_switch("morph", verbose)
    #------------------------------------------------------------------------#
    # getting last name if any to create new column
    # Later it would be good to write these sections as utility functions
    # they often do the same thing and it make things cleaner
    # not high priority - for now we need to make it work
    #------------------------------------------------------------------------#
    new_trial <- create_trial_tag(colnames(get_territories(vesalius_assay)),
      "Morphology")
    #------------------------------------------------------------------------#
    # First we define territory limits and add a little on each
    # side - this ensures that we won't be clipping any parts of the
    # territory
    #------------------------------------------------------------------------#
    ter <- extend_boundary(ter, morphology_factor)
    #------------------------------------------------------------------------#
    # Now we can do the morphing - grow, erode, clean and fill
    # For some reason I need to create mf seperately
    # maybe due to the fact that the grow/shrink function
    # is parsing weird stuff to imager/C++ ?
    #------------------------------------------------------------------------#
    for (i in seq_along(morphology_factor)) {
      mf <- abs(morphology_factor[i])
      if (morphology_factor[i] >= 0) {
        ter <- imager::grow(ter, mf)
      } else {
        ter <- imager::shrink(ter, mf)
      }
    }
    #------------------------------------------------------------------------#
    # Next we rebuild the image data frame with dilated
    #------------------------------------------------------------------------#
    ter <- ter %>% as.data.frame()
    ter <- inner_join(ter, get_tiles(vesalius_assay), by = c("x", "y")) %>%
      filter(origin == 1)
    buffer$trial[buffer$barcodes %in% ter$barcodes] <- adjusted_territory_values
    colnames(buffer) <-  new_trial

    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = ter,
    slot = "territories",
    append = TRUE)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(territory_morphing))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
    get_assay_names(vesalius_assay))
    simple_bar(verbose)
    return(vesalius_assay)
}


#' extend image boundary around territory 
#' @param territories data frame containing x/y and color value of
#' territories to extend
#' @param morphology_factor integer or vector of integers describing growth
#' and/or shrink extent.
#' @details we want to avoid clipping territory if they sit at the edge of the
#' image. To avoid this we simply extend the image boundary.
#' @importFrom dplyr %>% select
#' @importFrom imager as.cimg
extend_boundary <- function(territories, morphology_factor) {
  ymin <- ifelse((min(territories$y) - max(abs(morphology_factor)) * 2) <= 0, 1,
         min(territories$y) - morphology_factor * 2)
  xmin <- ifelse((min(territories$x) - max(abs(morphology_factor)) * 2) <= 0, 1,
         min(territories$x) - max(abs(morphology_factor)) * 2)
  ymax <- max(territories$y) + max(abs(morphology_factor)) * 2
  xmax <- max(territories$x) + max(abs(morphology_factor)) * 2
  territories <- territories %>% 
    select(c("x", "y", "value")) %>%
    rbind(., c(xmin, ymin, 1), c(xmax, ymax, 1)) %>%
    suppressWarnings(as.cimg())
  return(territories)
}

#' layer_territory generates layer from the outside to the inside of 
#' a territory
#' @param vesalius_assay vesalius_assay object
#' @param territory integer or vector of integers desrcining territories 
#' to morph.
#' @param trial character string - which territory trial that 
#' should be used to select
#' territorires. Default is last one computed
#' @param layer_depth integer describing the number of final layers.
#' @param morphology_factor integer or vector of integers describing growth
#' and/or shrink extent.
#' @param verbose logical - progress message output.
#' @details Each territory can be subdivided into a series of layers. 
#' Each layer will be considered a seperate territory and can be treated as
#' such for functions such as \code{\link{identify_markers}} and 
#' \code{\link{territory_plot}}.
#' 
#' However, all other territories present will be labled as "out". 
#' This means that for the time being you can only work with a single
#' territory at a time. 
#' @return a vesalius_assay
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
#' # morph territory
#'
#' ves <- layer_territory(ves)
#' 
#' # view territory morphing
#' territory_plot(ves)
#'}
#' @export
#' @importFrom dplyr right_join filter mutate select %>% inner_join
#' @importFrom imager grow as.cimg

layer_territory <- function(vesalius_assay,
  territory = NULL,
  trial = "last",
  layer_depth = NULL,
  morphology_factor = 0,
  verbose = TRUE) {
  simple_bar(verbose)
  #--------------------------------------------------------------------------#
  # This might be a bit messier but i'm just going to make a call to
  # territoryMorphing function
  # this has some un necssary transformations but this is cleaner
  # We check first if we need to run this if no morphological operator
  # If we run this then we always take the last trial run
  # this will update the morphology only if it is run
  #--------------------------------------------------------------------------#
  if (any(morphology_factor != 0)) {
    message_switch("morph", verbose)
    vesalius_assay <- territory_morphing(vesalius_assay,
      territory,
      trial,
      morphology_factor,
      verbose = FALSE)
      # adding last trial as last as here we will have added a morph trial
      trial <- "last"
      territory <- min(territory)
  }

  ter <- check_territories(vesalius_assay, trial)
  buffer <- ter
  ter <- right_join(ter, get_tiles(vesalius), by = "barcodes") %>%
    filter(trial %in% territory) %>%
    mutate(value = 1) %>%
    select(c("barcodes", "x.y", "y.y", "value", "origin", "trial"))
  colnames(ter) <- c("barcodes", "x", "y", "value", "origin", "trial")
  ter_for_loop <- ter
  #------------------------------------------------------------------------#
  # getting last name if any to create new column
  # This will yield two new columns
  #------------------------------------------------------------------------#
  new_trial <- create_trial_tag(colnames(get_territories(vesalius_assay)),
    "Layer")
  #------------------------------------------------------------------------#
  # First we define territory limits and add a little on each
  # side - this ensures that we won't be clipping any parts of the
  # territory
  #------------------------------------------------------------------------#
  ter <- extend_boundary(ter, morphology_factor)
  #--------------------------------------------------------------------------#
  # Now we can get edges of shape and compare this to tiles
  # and pool this "edge" into layers
  #--------------------------------------------------------------------------#
  message_switch("layer", verbose)
  counter <- 1
  layer <- list()
  while (nrow(ter_for_loop) > 0) {
    grad <- ter  %>%
      detect_edges() %>%
      grow(1) %>%
      suppressWarnings(as.cimg()) %>%
      as.data.frame() %>%
      filter(value > 0)
    #------------------------------------------------------------------------#
    # getting barcodes from territory
    #------------------------------------------------------------------------#
    edge <- inner_join(grad, ter_for_loop, by = c("x", "y")) %>%
      select(c("barcodes"))

    #------------------------------------------------------------------------#
    # Resizing ter - removing barcodes that are part of the edge
    #------------------------------------------------------------------------#
    ter_for_loop <- filter(ter_for_loop, !barcodes %in% unique(edge$barcodes))

    #------------------------------------------------------------------------#
    # Rebuilding an image but adding a little extra space
    #------------------------------------------------------------------------#
    if (nrow(ter_for_loop) > 0) {
        ter <- extend_boundary(ter_for_loop, morphology_factor)
    }
    #------------------------------------------------------------------------#
    # Adding edge to layer list and counting up
    #------------------------------------------------------------------------#
    layer[[counter]] <- unique(edge$barcodes)
    counter <- counter + 1
  }
  #--------------------------------------------------------------------------#
  # Now we can add the layers to the original territory
  # A rename layer if required
  #--------------------------------------------------------------------------#
  buffer$trial <- "out"
  for (lay in seq_along(layer)) {
    buffer$trial[buffer$barcodes %in% layer[[lay]]] <- lay
  }

  #--------------------------------------------------------------------------#
  # Finally we can split the different layers if we want to combine
  #--------------------------------------------------------------------------#
  layers <- unique(buffer$trial)
  if (!is.null(layer_depth)) {
    if (length(layers) < layer_depth) {
      warning("Layer depth exceeds layers in Territory -
        Using layers in territories", immediate. = TRUE)
    } else {
      idx <- floor(seq(1, length(layers), length.out = layer_depth + 1))
      for (i in seq(1, length.out = layer_depth)) {
        buffer$trial[buffer$trial %in% seq(idx[i], idx[i + 1])] <- i
      }
    }
  }
  #--------------------------------------------------------------------------#
  # rename new column
  colnames(buffer) <-  new_trial

  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = ter,
    slot = "territories",
    append = TRUE)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(layer_territory))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
    get_assay_names(vesalius_assay))
    simple_bar(verbose)
    return(vesalius_assay)
}
