################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#

#' Smooth Image
#'
#' Apply iterative smoothing to Vesalius images
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions numeric vector of latent space dimensions to use.
#' @param embedding character string describing which embedding should
#' be used.
#' @param method character describing smoothing method to use "median" ,
#' "iso"  or "box" or a combination of them.
#' @param iter numeric - number of smoothing iteration
#' @param sigma numeric - standard deviation associated with isoblur (Gaussian)
#' @param box numeric describing box size (centered around center pixel)
#' for smoothing
#' @param threshold numeric - discard pixels that are too low in value (cutoff
#' threshold only applied in box/median blurs).
#' @param neuman logical describing If Neumann boundary conditions should be
#' used, Dirichlet otherwise (default true, Neumann)
#' @param gaussian logical - use gaussian filter
#' @param na.rm logical describing if NA values should be removed
#' @param across_levels character - method used to account for multiple
#' smoothing levels (see details). Select from: "min","mean", "max"
#' @param verbose logical - progress message output.
#' @details The smooth_image function provides a series of options to smooth
#' your grey scale images contained within the vesalius_assay object.
#'
#' You can select any number of dimensions to be smoothed. As default, we take
#' the first 3 dimensions.
#'
#' Vesalius provides 3 smoothing methods:
#'  * median: computes median color in box (box size defined by box)
#'  * box: computes max color value in box (box size defined by box)
#'  * iso: compute gaussian kernel using sigma (Gussian SD defined by sigma)
#'
#' Vesalius can apply the same smoothing paramters over multiple iterations
#' as defined by the iter argument.
#' It is also possible to provide multiple values to box and sigma as numeric
#' vectors. Vesalius will run the smoothing over all provided values and return
#' either the maximum, minimum or mean color value as described by the
#' across_levels argument.
#'
#' Please note that unless specified the smoothing always applied to
#' the active embedding (default is last embedding computed). If you
#' want to start over you can simply set \code{embedding} to which ever
#' embedding you want to work with.
#'
#' This will take the raw embedding values before any image processing
#' has been applied. No need to re-run the embedding if you are not
#' satisfied with the smoothing.
#'
#' For more information, we suggest reading through the imager vignette.
#'
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
#' ves <- smooth_image(ves, embedding = "PCA")
#' # multiple rounds
#' ves <- smooth_image(ves, iter = 3, embedding = "PCA")
#' # accross level
#' ves <- smooth_image(ves, box = seq(3,11),
#'  accross_level = "mean",
#'  embedding = "PCA")
#'}
#' @export
#' @importFrom future.apply future_lapply
smooth_image <- function(vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "iso",
  iter = 1,
  sigma = 1,
  box = 20,
  threshold = 0,
  neuman = TRUE,
  gaussian = TRUE,
  na.rm = FALSE,
  across_levels = "min",
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # shifting format
    # ves_to_c => format.R
    #--------------------------------------------------------------------------#
    simple_bar(verbose)
    if (is(vesalius_assay, "vesalius_assay")) {
      assay <- get_assay_names(vesalius_assay)
      images <- format_ves_to_c(vesalius_assay = vesalius_assay,
        embed = embedding,
        dims = dimensions,
        verbose = verbose)
      method <- check_smoothing_kernel(method)
      across_levels <- check_smoothing_level(across_levels)
    } else {
      stop("Unsupported format to smooth_image function")
    }
    #--------------------------------------------------------------------------#
    # We can smooth over different arrays in parallel - I hope...
    # this might lead to some issue with low level C
    # Leave as is for now - but could be better to set up some check for
    # core allocation. Don't want to waste time assigning threads if not
    # required.
    #--------------------------------------------------------------------------#
    message_switch("smooth", verbose)
    images <- future_lapply(images, internal_smooth,
      method = method,
      iter = iter,
      sigma = sigma,
      box = box,
      threshold = threshold,
      neuman = neuman,
      gaussian = gaussian,
      na.rm = na.rm,
      across_levels = across_levels,
      future.seed = TRUE)
    #--------------------------------------------------------------------------#
    # shifting format
    # c_to_ves => format.R
    #--------------------------------------------------------------------------#
    embeds <- format_c_to_ves(images,
      vesalius_assay,
      dimensions,
      embed = embedding,
      verbose = verbose)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeds,
      slot = "active",
      append = FALSE)
    vesalius_assay <- add_active_embedding_tag(vesalius_assay, embedding)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(smooth_image))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
    simple_bar(verbose)
    return(vesalius_assay)

}


#' internal smoothing function
#' Function that does the smoothing on individual image array
#' @param image cimg array
#' @inheritParams smooth_image
#' @importFrom imager imrotate
#' @importFrom imager medianblur
#' @importFrom imager isoblur
#' @importFrom imager boxblur
#' @importFrom imager map_il
#' @importFrom imager parmin
#' @importFrom imager parmax
#' @importFrom imager average
internal_smooth <- function(image,
  method = c("median", "iso", "box"),
  iter = 1,
  sigma = 1,
  box = 20,
  threshold = 0,
  neuman = TRUE,
  gaussian = TRUE,
  na.rm = FALSE,
  across_levels = "min") {

  for (i in seq_len(iter)) {
        #----------------------------------------------------------------------#
        # This might seem strange but when doing a lot of smoothing there is
        # directionality bias. This introduces a shift in pixels and colour
        # To avoid this problem we can rotate the image.
        #----------------------------------------------------------------------#
    if (i %% 2 == 0) {
      image <- imager::imrotate(image, 180)
    }
    for (j in method) {
      image <- switch(j,
        "median" = imager::map_il(box, ~imager::medianblur(image,
          .,
          threshold)),
        "iso" = imager::map_il(sigma, ~imager::isoblur(image,
          .,
          neuman,
          gaussian,
          na.rm)),
        "box" = imager::map_il(box, ~imager::boxblur(image,
          .,
          neuman)))
      image <- switch(across_levels,
        "min" = imager::parmin(image),
        "max" = imager::parmax(image),
        "mean" = imager::average(image))
    }
    if (i %% 2 == 0) {
      image <- imager::imrotate(image, 180)
    }
  }
  return(image)
}


#' equalise image histogram
#'
#' equalizeHistogram image enhancement via colour histogram equalization.
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions numeric vector of latent space dimensions to use.
#' @param embedding character string describing which embedding should
#' be used.
#' @param type character - histogram EQ type. Select from: BalanceSimplest,
#' EqualizePiecewise, SPE, EqualizeDP, EqualizeADP, ECDF (see details)
#' @param N numeric describing how each colour channel will be mapped back to
#' the image (Higher N = Higher greyscale contrast).
#' Used with EqualizePiecewise
#' @param smax numeric - upper limit if contrast stretching.
#' Used with EqualizePiecewise
#' @param sleft numeric - Range 0 - 100. Percentage of pixel to be saturated on
#' the left side of the histogram. Used with BalanceSimplest
#' @param sright numeric - Range 0 - 100. Percentage of pixel to be saturated on
#' the right side of the histogram. Used with BalanceSimplest
#' @param lambda numeric - strength of background correction.
#' Used with SPE (Screened Poisson Equation).
#' @param up numeric - color value threshold in the upper limit.
#' Used with EqualizeDP.
#' @param down numeric color value threshold in the lower limit.
#' Used with EqualizeDP.
#' @param verbose logical - progress message output.
#' @details Histogram equalization ensures that image details are amplified.
#' In turn, territories may be extract with greater precision. We recommend
#' balancing the histogram prior to smoothing.
#'
#' For further details on each method described here, please refer to
#' \href{imagerExtra Vignette}{https://cran.r-project.org/web/packages/imagerExtra/vignettes/gettingstarted.html}
#' @return a vesalius_assay object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple EQ
#' ves <- equalisz_image(ves, embedding = "PCA")
#'}
#' @export
#' @importFrom imagerExtra EqualizePiecewise
#' @importFrom imagerExtra BalanceSimplest
#' @importFrom imagerExtra SPE
#' @importFrom imagerExtra EqualizeDP
#' @importFrom imagerExtra EqualizeADP
#' @importFrom future.apply future_lapply

equalize_image <- function(vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "BalanceSimplest",
  N = 1,
  smax = 1,
  sleft = 1,
  sright = 1,
  lambda = 0.1,
  up = 100,
  down = 10,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # shifting format
    #--------------------------------------------------------------------------#
    if (is(vesalius_assay, "vesalius_assay")) {
      assay <- get_assay_names(vesalius_assay)
      images <- format_ves_to_c(vesalius_assay = vesalius_assay,
        embed = embedding,
        dims = dimensions,
        verbose = verbose)
      method <- check_eq_method(method)
    } else {
      stop("Unsupported format to equalize_image function")
    }
    message_switch("eq", verbose)
    #--------------------------------------------------------------------------#
    # Equalizing histogram - depending on the image different methods will
    # work differently. It worth keeping in mind that these are made and
    # optimized for real images.
    #--------------------------------------------------------------------------#

    images <- switch(method,
      "EqualizePiecewise" = future_lapply(images,
        imagerExtra::EqualizePiecewise, N,
        future.seed = TRUE),
      "BalanceSimplest" = future_lapply(images,
        imagerExtra::BalanceSimplest, sleft, sright, range = c(0, 1),
        future.seed = TRUE),
      "SPE" = future_lapply(images, imagerExtra::SPE, lambda,
        future.seed = TRUE),
      "EqualizeDP"  = future_lapply(images, imagerExtra::EqualizeDP, down, up,
        future.seed = TRUE),
      "EqualizeADP" = future_lapply(images, imagerExtra::EqualizeADP,
        future.seed = TRUE),
      "ECDF" = future_lapply(images, ecdf_eq,
        future.seed = TRUE))
    #--------------------------------------------------------------------------#
    # shifting format
    # c_to_ves => format.R
    #--------------------------------------------------------------------------#
    embeds <- format_c_to_ves(images,
      vesalius_assay,
      dimensions,
      embed = embedding,
      verbose = verbose)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeds,
      slot = "active",
      append = FALSE)
    vesalius_assay <- add_active_embedding_tag(vesalius_assay, embedding)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(equalize_image))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
    simple_bar(verbose)
    return(vesalius_assay)
}

#' internal ecdf eq 
#' used for formatting
#' @param im image matrix
#' @importFrom imager as.cimg
#' @importFrom stats ecdf
ecdf_eq <- function(im) {
   return(suppressWarnings(imager::as.cimg(stats::ecdf(im)(im), dim = dim(im))))
}





#' regularise image
#' 
#' regularise_image denoise Vesalius images via variance regularization
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions numeric vector of latent space dimensions to use.
#' @param embedding character string describing which embedding should
#' be used.
#' @param lambda numeric - positive real numbers describing regularization
#' parameter (see details)
#' @param niter numeric - number of variance regularization iterations
#' (Default = 100)
#' @param verbose logical - progress message output.
#' @details Image regularization can be seen as a form of image denoising.
#' Details on each method can be found in the tvR package under the denoise2
#' function \href{tvR}{https://cran.r-project.org/web/packages/tvR/tvR.pdf}.
#'
#'
#' @return a vesalius_assay
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple regularisation
#' ves <- regularise_image(ves, embedding = "PCA")
#'}
#' @export
#' @importFrom future.apply future_lapply


regularise_image <- function(vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  lambda = 1,
  niter = 100,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # shifting format
    # ves_to_c => format.R
    #--------------------------------------------------------------------------#
    if (is(vesalius_assay, "vesalius_assay")) {
      assay <- get_assay_names(vesalius_assay)
      images <- format_ves_to_c(vesalius_assay = vesalius_assay,
        embed = embedding,
        dims = dimensions,
        verbose = verbose)
    } else {
      stop("Unsupported format to regularise_image function")
    }
    #--------------------------------------------------------------------------#
    # now we can do some var reg!
    # Will add the imager denoising function as well
    # regularisation is essentially denoising anyway
    #--------------------------------------------------------------------------#
    message_switch("reg", verbose)
    images <- future_lapply(images, regularise,
      lambda,
      niter,
      future.seed = TRUE)

    embeds <- format_c_to_ves(images,
      vesalius_assay,
      dimensions,
      embed = embedding,
      verbose = verbose)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeds,
      slot = "active",
      append = FALSE)
    vesalius_assay <- add_active_embedding_tag(vesalius_assay, embedding)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(regularise_image))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
    simple_bar(verbose)
    return(vesalius_assay)
}

#' Internal regularise function.
#' performs total variance regularisation 
#' @param img image array - generally a matrix
#' @param lambda lambda value for regularisation
#' @param niter numeric number of rounds of regularisation 
#' @importFrom tvR denoise2
#' @importFrom imager as.cimg
regularise <- function(img,
  lambda = 1,
  niter = 100) {
    #--------------------------------------------------------------------------#
    # Also not parsing any of the other tvR methods for denoise 
    # They don't work that well - we will see if it worth adding or not
    # This could do with some in depth optimization 
    # Beyond the scope of current dev cycle 
    # Force norm in all cases. 
    #--------------------------------------------------------------------------#
    img <- tvR::denoise2(as.matrix(img), lambda = lambda, niter = niter,
           method = "TVL2.FiniteDifference", normalize = FALSE)
    img <- (img - min(img)) / (max(img) - min(img))
    return(suppressWarnings(imager::as.cimg(img)))
}




#' segment image
#'
#' segment vesalius images to find initial territories
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions numeric vector of latent space dimensions to use.
#' @param embedding character string describing which embedding should
#' be used.
#' @param method character string for which method should be used for
#' segmentation. Select from "kmeans", "louvain", or "leiden".
#' @param col_resolution numeric colour resolution used for segmentation. 
#' (see details)
#' @param verbose logical - progress message output.
#' @details Applying image segmentation ensures a reduction in colour
#' complexity.
#'
#' Vesalius provides 3 different methods for clustering colours and
#' reducing color complexity: **Kmeans**, **Louvain**, and **Leiden**.
#'
#' In the case of kmeans clustering the \code{col_resolution} argument
#' shows the number of colours that the images should be reduced to.
#' In this case, \code{col_resolution} should be an integer and
#' we suggest first looking at values between 3 and 20.
#' 
#' In the case of **leiden** and **louvain** clustering, the
#' \code{col_resolution} is the granularity of the clustering.
#' In this case, we suggest using values between 0.01 and 1 to start with.
#' We recommned uisng **louvain** clustering over **leiden** in 
#' this context. 
#' 
#' The optimal \code{col_resolution} will depend on your interest and 
#' biological question at hand. You might be interested in more or less
#' granular territories. Along with smoothing, the number of segments is
#' one way to control this granularity. 
#
#' @return a vesalius_assay object
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
#'}
#' @export

segment_image <- function(vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "kmeans",
  col_resolution = 10,
  verbose = TRUE) {
  simple_bar(verbose)
  message_switch("seg", verbose, method = method)
  #----------------------------------------------------------------------------#
  # Parsing vesalius object so we can recontruct it internally and not need to
  # rebuild intermediates and shift between formats...
  #----------------------------------------------------------------------------#
  method <- check_segmentation_method(method)
  segments <- switch(method,
    "kmeans" = kmeans_segmentation(vesalius_assay,
      col_resolution = col_resolution,
      dimensions = dimensions,
      embedding = embedding,
      verbose = verbose),
    "leiden" = leiden_segmentation(vesalius_assay,
      dimensions = dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      verbose = verbose),
    "louvain" = louvain_segmentation(vesalius_assay,
      dimensions = dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      verbose = verbose))

  #vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
  #  data = segments$segments,
  #  slot = "active",
  #  append = FALSE)
  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = segments,
    slot = "territories",
    append = TRUE)
  vesalius_assay <- add_active_embedding_tag(vesalius_assay, embedding)
  commit <- create_commit_log(arg_match = as.list(match.call()),
    default = formals(segment_image))
  vesalius_assay <- commit_log(vesalius_assay,
    commit,
    get_assay_names(vesalius_assay))
  simple_bar(verbose)
  return(vesalius_assay)
}

#' kmeans segmentation function
#' @param vesalius_assay a vesalius_assay
#' @param col_resolution integer or vector of positive integers.
#' Colour depth used for segmentation.
#' @inheritParams segment_image
#' @details Run an interaive kmeans segmentation with the possibility 
#' to run multiple rounds of smoothing. 
#' @return list containing segmented image as an active embedding and
#' territory cluster for all barcodes.
#' @importFrom dplyr right_join
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom stats kmeans
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom imager as.cimg
#' @importFrom utils tail
kmeans_segmentation <- function(vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 10,
  embedding = "last",
  verbose = TRUE) {
  coord <- get_tiles(vesalius_assay) %>%
    filter(origin == 1)
  embeddings <- check_embedding_selection(vesalius_assay,
    embedding,
    dimensions)

  #--------------------------------------------------------------------------#
  # Now lets cluster colours
  # Remember that here colours are your new active embedding values
  #--------------------------------------------------------------------------#
  km <- suppressWarnings(kmeans(embeddings,
    col_resolution,
    iter.max = 100,
    nstart = 10))
  clusters <- km$cluster
  #kcenters <- km$centers
  match_loc <- match(coord$barcodes, names(clusters))
  clusters <- data.frame(coord, "Segment" = clusters[match_loc])
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(clusters)
}



#' leiden segmentation
#'
#' using leiden clustering to cluster colors
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions embedding dimensions used for clustering
#' @param col_resolution clustering resolution used for leiden
#' @param embedding embedding type used for clustering
#' @param verbose logical if progress message should outputed
#' @returns list with updated segmented embedding values 
#' and segment territories.
#' @importFrom igraph cluster_leiden
#' @importFrom dplyr %>%
leiden_segmentation <- function(vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 0.01,
  embedding = "last",
  verbose = TRUE) {
  coord <- get_tiles(vesalius_assay) %>%
    filter(origin == 1)
  embeddings <- check_embedding_selection(vesalius_assay,
    embedding,
    dimensions)
  graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
  clusters <- igraph::cluster_leiden(graph,
    resolution_parameter = col_resolution)
  cluster <- data.frame("cluster" = clusters$membership,
    "barcodes" = clusters$names)

  
  match_loc <- match(coord$barcodes, cluster$barcodes)
  clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(clusters)
}

#' louvain segmentation
#'
#' using leiden clustering to cluster colors
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions embedding dimensions used for clustering
#' @param col_resolution clustering resolution used for leiden
#' @param embedding embedding type used for clustering
#' @param verbose logical if progress message should outputed
#' @returns list with updated segmented embedding values 
#' and segment territories.
#' @importFrom  igraph cluster_louvain
#' @importFrom dplyr %>%
louvain_segmentation <- function(vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 0.01,
  embedding = "last",
  verbose = TRUE) {
  coord <- get_tiles(vesalius_assay) %>%
    filter(origin == 1)
  embeddings <- check_embedding_selection(vesalius_assay,
    embedding,
    dimensions)
  graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
  clusters <- igraph::cluster_louvain(graph, resolution = col_resolution)
  cluster <- data.frame("cluster" = clusters$membership,
    "barcodes" = clusters$names)


  match_loc <- match(coord$barcodes, cluster$barcodes)
  clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(clusters)
}

#' compute and greate nearest neighbor graph
#' @param embeddings embedding matrix
#' @param k numeric describing number of nearest neighbors
#' @return igraph object
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_data_frame
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
compute_nearest_neighbor_graph <- function(embeddings, k = 20) {
    knn <- RANN::nn2(embeddings, k = k + 1)$nn.idx
    rownames(knn) <- rownames(embeddings)
    graph <- populate_graph(knn)
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
}



#' populate graph with network connections
#' @param chunk data frame with nearest neighbors
#' @returns data frame with network connections
populate_graph <- function(chunk) {
  barcodes <- rownames(chunk)
  nxk <- ncol(chunk)
  template <- data.frame("e1" = rep(barcodes, each = nxk),
    "e2" = barcodes[c(t(chunk))])
  return(template)
}

#' isolating territories from vesalius image segments
#' @param vesalius_assay vesalius_Assay object
#' @param method character describing barcode pooling method.
#' Currently, only "distance" availble
#' @param trial character string describing which segmentation trial
#' to use. Default is "last" which is the last segmentation trial used.
#' @param capture_radius numeric - proportion of maximum distance between
#' barcodes that will be used to pool barcodes together (range 0 - 1).
#' @param global logical - If TRUE, territories will be numbered across all
#' colour segments. If FALSE, territories will be numbered within each colour
#' segment.
#' @param min_spatial_index integer - minimum number of barcodes/spots/beads
#' required in each territory
#' @param verbose logical - progress message output.
#' @details Image segments can be further subdivided into 2D
#' seperated territorires. This is accomplished by pooling barcodes
#' that are associated with a colour cluster into territories based on the
#' distance between each barcode.
#'
#' First, \code{isolate_territories} considers the maximum distance
#' between all beads. The \code{capture_radius} will define which 
#' proportion of this distance should be considered.
#'
#' Second, a seed barcode will be selected and all barcodes that are within the
#' capture distance of the seed barcode with be pooled together. This process
#' is then applied on barcodes that have been selected in this manner. The
#' process is repeated until all barcodes have been pooled into a territory.
#' If there are still barcodes remaining, a new seed barcode is selected and the
#' whole process is repeated. NOTE : Territory isolation is applied to each
#' colour segment independently.
#'
#' If a territory does not contain enough barcodes, it will be pooled into the
#' isolated territory. This territory contains all isolated territories
#' regardless of colour cluster of origin.
#'
#' @return a vesalius_assay object
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
#'}
#' @export
isolate_territories <- function(vesalius_assay,
  method = "distance",
  trial = "last",
  capture_radius = 0.05,
  global = TRUE,
  min_spatial_index = 10,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # Get stuff out as usual
    # Not super happy with this check
    # it's a bit messy - we might need to consider to do a whole sanity check
    # of inout data and see if that makes sense - this will include checking log
    #--------------------------------------------------------------------------#
    ter <- check_segment_trial(vesalius_assay, trial) %>%
      na.exclude()
    #--------------------------------------------------------------------------#
    # Compute real capture Radius
    # Only one method for now so this is not necessary 
    # keep it for later
    #--------------------------------------------------------------------------#
    method <- check_isolation_method(method)
    if (method[1L] == "distance") {
      capture_radius <- sqrt(((max(ter$x) - min(ter$x))^2 +
        (max(ter$y) - min(ter$y))^2)) * capture_radius
    }
    #--------------------------------------------------------------------------#
    # Creating new trial column name and adding it to input data
    # The input data here is a subset of the full territory df
    # we at least make sure that we are using the right input
    #--------------------------------------------------------------------------#
    new_trial <- create_trial_tag(colnames(get_territories(vesalius_assay)),
      "Territory") %>%
      tail(1)
    ter$trial <- 0
    #--------------------------------------------------------------------------#
    # Now we can dispatch the necessary information
    # and run the pooling algorithm
    #--------------------------------------------------------------------------#
    segment <- unique(ter$segment)
    for (i in seq_along(segment)) {
      message_switch("ter_pool", verbose, ter = i)
      #------------------------------------------------------------------------#
      # We don't need all data in this case only clusters and locations
      # We can just rebuild everything afterwards
      # At least we don't don't need to compute anything unnecessarily
      ## Note other method not in use for now
      # Might not be worthwile to implement them
      # Argument could be removed
      #------------------------------------------------------------------------#

      tmp <- ter[ter$segment == segment[i], ]
      tmp <- switch(method[1L],
        "distance" = distance_pooling(tmp, capture_radius,
          min_spatial_index))
      #------------------------------------------------------------------------#
      # Skipping if does not meet min cell requirements
      # filtering can be done later
      #------------------------------------------------------------------------#
      if (is.null(tmp)) next()
      #------------------------------------------------------------------------#
      # adding territories
      #------------------------------------------------------------------------#
      ter$trial[ter$segment == segment[i]] <- tmp

    }

    #--------------------------------------------------------------------------#
    # Globalise territories crate a numbering system for all colour clusters
    # If set to false territories only describe territories for any given colour
    # cluster
    #--------------------------------------------------------------------------#
    if (global) {
        ter <- globalise_territories(ter)
    }
    #--------------------------------------------------------------------------#
    #Now we can add it to the full territory df and update vesalius
    #--------------------------------------------------------------------------#

  colnames(ter) <- gsub("trial", new_trial, colnames(ter))
  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = ter,
    slot = "territories",
    append = TRUE)
  commit <- create_commit_log(arg_match = as.list(match.call()),
    default = formals(isolate_territories))
  vesalius_assay <- commit_log(vesalius_assay,
    commit,
    get_assay_names(vesalius_assay))
  simple_bar(verbose)
  return(vesalius_assay)
}



#' distance pooling beads of colour segment into seperate territories
#' @param img data frame contain all barcodes of a single sgement
#' @param capture_radius numeric proportion of max distance between beads
#' to use as distance threshold between beads
#' @param min_spatial_index numeric minimum number of beads that should
#' be contained in a territory. 
#' @details Beads that are too far away or bead cluster that are below
#' the minimum number of spatial indices will all be pooled under the
#' isolated label. Note that this label is used across color segments.
#' @importFrom dplyr %>% distinct
distance_pooling <- function(img, capture_radius, min_spatial_index) {
    #--------------------------------------------------------------------------#
    # Select center point of each tile for only one channel
    # Dont need to run it for all channels
    #--------------------------------------------------------------------------#
    img_copy <- img  %>%
      distinct(barcodes, .keep_all = TRUE) 
    if (nrow(img_copy) < 1) { return(NULL) }
    #--------------------------------------------------------------------------#
    # Compute distances
    #--------------------------------------------------------------------------#
    idx <- seq_len(nrow(img_copy))
    distance_matrix <- lapply(idx, function(idx, mat) {
      xo <- mat$x[idx]
      yo <- mat$y[idx]
      xp <- mat$x
      yp <- mat$y
      distance <- sqrt(((abs(xp - xo))^2 + (abs(yp - yo))^2))
      return(distance)
    }, img_copy)


    #--------------------------------------------------------------------------#
    # Buildan actual matrix
    #--------------------------------------------------------------------------#

    distance_matrix <- do.call("rbind", distance_matrix)
    colnames(distance_matrix) <- img_copy$barcodes
    rownames(distance_matrix) <- img_copy$barcodes

    #--------------------------------------------------------------------------#
    # If there is only one point
    # In this case you only have one territory as well
    # Don't need to any pooling
    #--------------------------------------------------------------------------#
    if (sum(dim(distance_matrix)) == 2) {
        return(1)
    }


    #--------------------------------------------------------------------------#
    # Pooling points together
    #--------------------------------------------------------------------------#
    barcodes <- img_copy$barcodes
    territories <- list()
    count <- 1

    while (length(barcodes) > 0) {
         #---------------------------------------------------------------------#
         # First lets select a random barcode in the colour segment
         # And create a pool of barcodes to select based on capture radius
         # This first while loop checks if there are still barcodes
         # left in the colour segment
         #---------------------------------------------------------------------#
          tmp <- distance_matrix[, sample(barcodes, 1)]
          pool <- names(tmp)[tmp <= capture_radius]
          inter <- pool
          converge <- FALSE

          while (!converge) {
            #------------------------------------------------------------------#
            # This while loop checks if all possible barcodes have been pooled
            # into the current territory
            #------------------------------------------------------------------#
            if (length(inter) == 1) {
              #------------------------------------------------------------#
              # when there is only one barcodes
              # remove barcode from pool and move on
              #------------------------------------------------------------#
              territories[[count]] <- pool
              barcodes <- barcodes[!barcodes %in% pool]
              count <- count + 1
              converge <- TRUE
            } else {
              #------------------------------------------------------------#
              # Get a new pool from the distance matrix
              # and check which ones are within capture radius
              #------------------------------------------------------------#
              new_pool <- distance_matrix[, inter]
              new_pool <- unique(unlist(lapply(seq_len(ncol(new_pool)),
                function(idx, np, capture_radius) {
                  res <- rownames(np)[np[, idx] <= capture_radius]
                  return(res)
              }, new_pool, capture_radius)))
              #------------------------------------------------------------#
              # check which barcodes in the new pool overlap with
              # the ones in the full pool
              # If there is a perfect overlap then there are no new bacodes
              # to pool into a territory
              #------------------------------------------------------------#
              overlap <- new_pool %in% pool
              if (sum(overlap) != length(new_pool)) {
                #------------------------------------------------------#
                # There are still some new barcodes to pool
                # lets do some more looping then
                #------------------------------------------------------#
                pool <- unique(c(pool, new_pool[!overlap]))
                inter <- unique(new_pool[!overlap])
                converge <- FALSE
              } else {
                #------------------------------------------------------#
                # it is done ! no more barcodes for this territory
                #------------------------------------------------------#
                territories[[count]] <- pool
                count <- count + 1
                barcodes <- barcodes[!barcodes %in% pool]
                converge <- TRUE
              }
            }
          }
      }
    #--------------------------------------------------------------------------#
    # Clean up drop outs
    #--------------------------------------------------------------------------#

    all_ters <- img$trial

    for (ter in seq_along(territories)) {
        loc <- img$barcodes %in% territories[[ter]]
        if (length(territories[[ter]]) <= min_spatial_index) {
            all_ters[loc] <- "isolated"
        } else {
            all_ters[loc] <- ter
        }
    }
      return(all_ters)
}


#' create a registration template image
#' @param vesalius_assay a vesalius_assay object
#' @param from character string defining what should be used to create
#' the image template
#' @param dimensions numeric vector defining which dimensions should be used
#' if calling en embedding
#' @details creating a template gray image that we will use for image registration
#' The point here is to make sure that we have something can use to compared the
#' the image to a microscopy image
#' @return a cimg object
generate_image_template <- function(vesalius_assay,
  from = "Territory",
  dimensions = seq(1:3)) {
  #---------------------------------------------------------------------------#
  # check if source for image template is availbale
  #---------------------------------------------------------------------------#
  tiles <- get_tiles(vesalius_assay)
  image <- check_template_source(vesalius_assay, from, dimensions)
  type <- comment(image)

  if (type != "Embedding") {
    #-------------------------------------------------------------------------#
    # Creating evenly space color values from segments
    #-------------------------------------------------------------------------#
    value <- suppressWarnings(as.numeric(image[, 4]))
    value[is.na(value)] <- max(value[!is.na(value)]) + 1
    value <- (value / max(value)) * 255
    image[, 4] <- value
    colnames(image) <- c("barcodes", "x", "y", "value")
    image <- right_join(tiles, image, by = "barcodes") %>%
    select(c("barcodes", "x.x", "y.x", "value"))
    colnames(image) <- c("barcodes", "x", "y", "value")
    image <- na.exclude(image)
  } else {
    #-------------------------------------------------------------------------#
    # Color value based on sum across latent space dimensions
    #-------------------------------------------------------------------------#
    value <- apply(image, 1, sum)
    value <- min_max(value)
    image <- data.frame("barcodes" = names(value), "value" = value)
    image <- right_join(tiles, image, by = "barcodes") %>%
      select(c("barcodes", "x", "y", "value"))
  }
  #---------------------------------------------------------------------------#
  # rebuild image with median colored background
  #---------------------------------------------------------------------------#
  img <- as.cimg(array(median(image$value), c(max(image$x), max(image$y))))
  ind <- imager:::index.coord(img, image[, c("x", "y"), drop = FALSE])
  img[ind] <- image[["value"]]
  return(img)
}