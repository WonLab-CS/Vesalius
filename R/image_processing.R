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
    # c_to_ves => format_and_dispatch.R
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
#' A higher value for lambda will results in a smoother image. It should be noted
#' that in the context of spatial omics the more sparse the points in the data
#' (the more space between coordinates), the more you will need to increase
#' the value of lambda to obtain better denoising. 
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
#' segmentation. Select from "kmeans", "louvain", "leiden", "slic", 
#' "leiden_slic","louvain_slic","som"
#' @param col_resolution numeric colour resolution used for segmentation. 
#' (see details)
#' @param compactness numeric - factor defining super pixel compaction.
#' @param scaling numeric - scaling image ration during super pixel 
#' segmentation.
#' @param k numeric - number of closest super pixel neighbors to consider
#' when generating segments from super pixels
#' @param threshold numeric [0,1] - correlation threshold between 
#' nearest neighbors when generating segments from super pixels.
#' @param verbose logical - progress message output.
#' @details Applying image segmentation ensures a reduction in colour
#' complexity.
#'
#' Vesalius provides 7 different methods for clustering colours and
#' reducing color complexity: **Kmeans**, **Louvain**, **Leiden**,
#' **slic**, **leiden_slic**, **louvain_slic**, and **som**
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
#' In the case of slic, the col_resolution define the number of starting
#' points used to generate super pixels. Depending on the number of
#' points there are in the assay, we suggested using 10% of the total 
#' number of points as starting point. 
#' For example, if you have  1000 spatial indices, you can set 
#' col_resolution to 100.
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
  compactness = 1,
  scaling = 0.5,
  threshold = 0.9,
  index_selection = "bubble",
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
      verbose = verbose),
    "slic" = slic_segmentation(vesalius_assay,
      dimensions = dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      index_selection = index_selection,
      compactness = compactness,
      scaling = scaling,
      threshold = threshold,
      verbose = verbose),
    "louvain_slic" = louvain_slic_segmentation(vesalius_assay,
        dimensions = dimensions,
        col_resolution = col_resolution,
        embedding = embedding,
        compactness = compactness,
        scaling = scaling,
        verbose = verbose),
    "leiden_slic" = leiden_slic_segmentation(vesalius_assay,
        dimensions = dimensions,
        col_resolution = col_resolution,
        embedding = embedding,
        compactness = compactness,
        scaling = scaling,
        verbose = verbose),
    "som" = som_segmentation(vesalius_assay,
      dimensions= dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      compactness = compactness,
      scaling = scaling,
      verbose = verbose))

  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = segments$active,
    slot = "active",
    append = FALSE)
  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = segments$segments,
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
    iter.max = 10000,
    nstart = 10))
  clusters <- km$cluster
  kcenters <- km$centers
  match_loc <- match(coord$barcodes, names(clusters))
  clusters <- data.frame(coord, "Segment" = clusters[match_loc])
  active <- assign_centers(vesalius_assay, clusters, kcenters, dimensions)
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(list("segments" = clusters, "active" = active))
}

#' assign centroid values to active embedding 
#' @param vesalius_assay a vesalius assy object 
#' @param clusters data.frame containing cluster values
#' @param kcenters matrix containing centroid values for each dimension
#' @param dimensions vector (nummeric / int) describin which latent space
#' @param ratio if used in the context of super pixel - spatial ration
#' dimensiuons shouls be used. 
#' @return matrix for the active embedding usiong color segementation 

assign_centers <- function(vesalius_assay,
  clusters,
  kcenters,
  dimensions,
  ratio = NULL) {
  #sp <- map(1:spectrum(images),~ km$centers[km$cluster,2+.]) %>% do.call(c,.) %>% as.cimg(dim=dim(images))
  active <- vesalius_assay@active
  for (d in seq_along(dimensions)) {
      tmp <- kcenters[clusters$Segment, dimensions[d]]
      if (!is.null(ratio)) {
        tmp <- tmp / ratio 
      }
      active[, dimensions[d]] <- tmp
  }
  return(active)
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
    dimensions)[,dimensions]
  
  loc <- match(coord$barcodes, rownames(embeddings))
  embeddings <- cbind(coord[, c("x", "y")], embeddings[loc, ])
  graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
  clusters <- igraph::cluster_leiden(graph,
    resolution_parameter = col_resolution)
  cluster <- data.frame("cluster" = clusters$membership,
    "barcodes" = clusters$names)
  match_loc <- match(coord$barcodes, cluster$barcodes)
  clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
  active <- create_pseudo_centroids(vesalius_assay,
    clusters,
    dimensions)
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(list("segments" = cluster, "active" = active))
}


#' create centroid vlues for louvain and leiden 
#' @param vesalius_assay  vesalius assay object
#' @param clusters data frame containing clusters/segments 
#' @param dimensions vector (nummeric / int) describin which latent space
#' dimensiuons shouls be used.
#' @return matrix for the active embedding usiong color segementation 
create_pseudo_centroids <- function(vesalius_assay, clusters, dimensions) {
  active <- vesalius_assay@active
  for (d in dimensions) {
    for (clust in unique(clusters$Segment)) {
        loc <- match(clusters$barcodes[clusters$Segment == clust],
          rownames(active))
        active[loc, d] <- mean(active[loc, d])
    }
  }
  return(active)
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
    dimensions)[, dimensions]
  loc <- match(coord$barcodes, rownames(embeddings))
  embeddings <- cbind(coord[, c("x", "y")], embeddings[loc, ])
  graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
  clusters <- igraph::cluster_louvain(graph, resolution = col_resolution)
  cluster <- data.frame("cluster" = clusters$membership,
    "barcodes" = clusters$names)


  match_loc <- match(coord$barcodes, cluster$barcodes)
  clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
  active <- create_pseudo_centroids(vesalius_assay,
    clusters,
    dimensions)
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(list("segments" = cluster, "active" = active))
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



#' @importFrom imager imappend imsplit spectrum
louvain_slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    compactness = 1,
    scaling = 0.3,
    verbose) {
    #-------------------------------------------------------------------------#
    # first we get tiles and get images 
    # then we compute the a scaling metric - this is normally based on 
    # empirical data (0.28 is one way) but it essentially a measure of the
    # scaling value between pixel distance and color distance 
    # we measure the max spread of colors as well
    # all this serves as a way to scale the "spread" of the color values 
    # to match the "spread" of the spatial coordinates
    #-------------------------------------------------------------------------#
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    loc <- match(coord$barcodes, rownames(embeddings))
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    embeddings$barcodes <- barcodes <- rownames(embeddings)
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    rownames(embeddings) <- barcodes
    #-------------------------------------------------------------------------#
    # Run louvain
    #-------------------------------------------------------------------------#
    graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
    clusters <- igraph::cluster_louvain(graph, resolution = col_resolution)
    cluster <- data.frame("cluster" = clusters$membership,
      "barcodes" = clusters$names)


    match_loc <- match(coord$barcodes, cluster$barcodes)
    clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
    active <- create_pseudo_centroids(vesalius_assay,
      clusters,
      dimensions)
    active <- active / ratio
    
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
      "Segment") %>%
      tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = active))
}


#' @importFrom imager imappend imsplit spectrum
leiden_slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    compactness = 1,
    scaling = 0.3,
    verbose) {
    #-------------------------------------------------------------------------#
    # first we get tiles and get images 
    # then we compute the a scaling metric - this is normally based on 
    # empirical data (0.28 is one way) but it essentially a measure of the
    # scaling value between pixel distance and color distance 
    # we measure the max spread of colors as well
    # all this serves as a way to scale the "spread" of the color values 
    # to match the "spread" of the spatial coordinates
    #-------------------------------------------------------------------------#
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    loc <- match(coord$barcodes, rownames(embeddings))
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    embeddings$barcodes <- barcodes <- rownames(embeddings)
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    rownames(embeddings) <- barcodes
    #-------------------------------------------------------------------------#
    # Run louvain
    #-------------------------------------------------------------------------#
    graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
    clusters <- igraph::cluster_leiden(graph,
      resolution_parameter = col_resolution)
    cluster <- data.frame("cluster" = clusters$membership,
      "barcodes" = clusters$names)


    match_loc <- match(coord$barcodes, cluster$barcodes)
    clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
    active <- create_pseudo_centroids(vesalius_assay,
      clusters,
      dimensions)
    active <- active / ratio
    
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
      "Segment") %>%
      tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = active))
}



#' @importFrom imager imappend imsplit spectrum
#' @importFrom purrr map_dbl map
slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    index_selection = "bubble",
    compactness = 1,
    scaling = 0.3,
    threshold = 0.9,
    max_iter = 1000,
    verbose) {
    #-------------------------------------------------------------------------#
    # first we get tiles and get images 
    # then we compute the a scaling metric - this is normally based on 
    # empirical data (0.28 is one way) but it essentially a measure of the
    # scaling value between pixel distance and color distance 
    # we measure the max spread of colors as well
    # all this serves as a way to scale the "spread" of the color values 
    # to match the "spread" of the spatial coordinates
    #-------------------------------------------------------------------------#
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    
    #tiles <- get_tiles(vesalius_assay)
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    embeddings$barcodes <- rownames(embeddings)
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    #-------------------------------------------------------------------------#
    # Generate initial centers from a grid
    # This will need to be updated to make sure
    # At the moment this will take barcodes within the data set to use as 
    # inital centers.
    # It would be worth looking into using pixels instead. Mainly for low
    # res data sets - we might want to have pixels between spatial indices 
    #-------------------------------------------------------------------------#
    index <- select_initial_indices(coord,
        #embeddings,
        type = index_selection,
        n_centers = col_resolution,
        max_iter = max_iter)
    
    #-------------------------------------------------------------------------#
    # Run k-means - you should have read the mnual pages you silly goose
    #-------------------------------------------------------------------------#
    km <- suppressWarnings(kmeans(embeddings,
        embeddings[index, ],
        iter.max = max_iter))
    centroids <- purrr::map(seq(1, l = ncol(embeddings) - 2),
        ~ km$centers[km$cluster, .]) %>%
        do.call("cbind", .)
    centroids <- centroids / ratio
    clusters <- cbind(embeddings[, c("x", "y")], km$cluster) %>%
        as.data.frame() %>%
        right_join(coord, by = c("x", "y"))
    colnames(clusters) <- c("x", "y", "value", "barcodes", "origin")
    match_loc <- match(coord$barcodes, clusters$barcodes)
    clusters <- data.frame(coord, "Segment" = clusters[match_loc, "value"])
    
    #-------------------------------------------------------------------------#
    # rescaling back to original values
    # rebuilding everything 
    # ON HOLD: combining super pixels into large scale segments
    #-------------------------------------------------------------------------#
    
    embeddings[, seq(1, ncol(embeddings) - 2)] <- centroids

    embeddings <- lapply(seq_len(ncol(centroids)), function(i, embed) {
        ret <- as.data.frame(embed[, c("x", "y", paste0("dim_", i))])
        colnames(ret) <- c("x", "y", "value")
        return(ret)
    }, embed = embeddings)
    #clusters <- select_similar(embeddings, coordinates = coord)
    embeddings <- format_c_to_ves(embeddings,
      vesalius_assay,
      dimensions,
      embed = embedding,
      verbose = FALSE)
    
    #clusters <- connected_pixels(clusters, embeddings, k, threshold, verbose)
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
        "Segment") %>%
    tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = embeddings))
}


#' @importFrom kohonen som somgrid
som_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    compactness = 1,
    scaling = 0.3,
    verbose) {
    #-------------------------------------------------------------------------#
    # first we get tiles and get images 
    # then we compute the a scaling metric - this is normally based on 
    # empirical data (0.28 is one way) but it essentially a measure of the
    # scaling value between pixel distance and color distance 
    # we measure the max spread of colors as well
    # all this serves as a way to scale the "spread" of the color values 
    # to match the "spread" of the spatial coordinates
    #-------------------------------------------------------------------------#
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    
    #tiles <- get_tiles(vesalius_assay)
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------#
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    bar <- rownames(embeddings)
    embeddings$barcodes <- bar
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    #-------------------------------------------------------------------------#
    # create SOM map of data 
    #-------------------------------------------------------------------------#
    som <- som(as.matrix(embeddings),
      grid = somgrid(xdim = floor(sqrt(col_resolution)),
        ydim = ceiling(sqrt(col_resolution))))
    
    cluster <- data.frame("cluster" = som$unit.classif,
      "barcodes" = bar)
    match_loc <- match(coord$barcodes, cluster$barcodes)
    clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
    active <- create_pseudo_centroids(vesalius_assay,
      clusters,
      dimensions)
    active <- active / ratio
    
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
      "Segment") %>%
      tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = active))

}

#' select inital indices
#' @param coordinates data frame containing spatial coordinates of beads
#' @param embeddings matrix containing embedding values - full pixel image
#' @param n_centers numeric number of beads to select as super pixel centers
#' @param max_iter numeric number of iteration before returning result if 
#' no coveregnce.
#' @return barcodes of starting coordinates
#' @importFrom RANN nn2
select_initial_indices <- function(coordinates,
    embeddings,
    type = "bubble",
    n_centers = 500,
    max_iter = 500) {
    indices <- switch(type,
        # "linear" = uniform_sampling(coordinates,
        #     embeddings,
        #     n_centers),
        "bubble" = bubble_stack(coordinates,
            #embeddings,
            n_centers,
            max_iter),
        # "pixel" = pixel_sampling(coordinates,
        #     embeddings,
        #     n_centers),
        "random" = random_sampling(coordinates,
          n_centers = n_centers),
        "hex" = hex_grid(coordinates,
          n_centers) )
    return(indices)
}

# linear_sampling <- function(coordinates,
#     embeddings,
#     n_centers = 500) {
#     coordinates <- coordinates[order(coordinates$x, decreasing = FALSE),]
#     coordinates <- split(coordinates, coordinates$x) %>%
#         lapply(., function(df) {return(df[order(df$y, decreasing = FALSE),])})
#     coordinates <- do.call("rbind", coordinates)
#     coordinates <- coordinates[seq(1, nrow(coordinates), l = n_centers), ]

#     in_image <- paste0(embeddings[, "x"], "_", embeddings[, "y"])
#     in_background <- paste0(coordinates$x, "_", coordinates$y)
#     in_image <- which(in_image %in% in_background)
#     return(in_image)
# }

random_sampling <- function(coordinates, n_centers) {
    #idx <- which(coordinates$origin == 1)
    return(sample(seq_len(nrow(coordinates)),
      size = n_centers, replace = FALSE))
}

#' importFrom imager nPix
# pixel_sampling <- function(coordinates, embeddings, n_centers) {
#   images <- coordinates[, c("x", "y")]
#   images$value <- coordinates$origin
#   images <- suppressWarnings(as.cimg(images))
#   ind <- round(seq(1, nPix(images) / spectrum(images), l = n_centers))
#   images <- as.data.frame(images)[ind, ]
#   in_background <- paste0(images[, "x"], "_", images[, "y"])
#   in_image <- paste0(embeddings[,"x"], "_", embeddings[,"y"])
#   in_image <- which(in_image %in% in_background)
#   return(in_image)
# } 

bubble_stack <- function(coordinates,
    n_centers = 500,
    max_iter = 500) {
    coordinates <- coordinates %>% filter(origin == 1)
    #-------------------------------------------------------------------------#
    # intialise background grid, active grod and get nearest neighbors
    # we assume that you will have at least 4 inital points
    # This needs to be changed later - I dont think using knn is the best 
    # option here
    #-------------------------------------------------------------------------#
    knn <- RANN::nn2(coordinates[, c("x", "y")],
      k = floor(nrow(coordinates) * 0.2))
    radius <- max(knn$nn.dist)
    #-------------------------------------------------------------------------#
    # iterated a reduced space until you get the desired numbers of circles  
    #-------------------------------------------------------------------------#
    convergence <- FALSE
    iter <- 1
    while (!convergence) {
        #---------------------------------------------------------------------#
        # intialize everything we need 
        #---------------------------------------------------------------------#
        background_grid <- c()
        active <- seq_len(nrow(coordinates))
        nn_index <- knn$nn.idx
        nn_dist <- knn$nn.dist
        #---------------------------------------------------------------------#
        # as long as their points present we repeat the process of selecting 
        # point that are outisde the selection radius
        #---------------------------------------------------------------------#
        while (length(active) > 0) {
            random_start <- active[sample(x =
                seq(1, l = length(active)),
                size = 1)]
            
            background_grid <- c(background_grid, random_start)
            removing <- nn_index[random_start,
                nn_dist[random_start, ] <= radius]
            active <- active[!active %in% unique(removing)]
        }
        if (length(background_grid) == n_centers) {
            convergence <- TRUE
        }
        if (iter >= max_iter) {
            convergence <- TRUE
            warning("Max Iteration with no convergence!
            Returning Approximation")
        }
        #---------------------------------------------------------------------#
        # we change the radius based on how many points there 
        # we assume that if there are 2 few barcodes then the radius is to 
        # large 
        # if there are too many barcodes then the radius is 2 small. 
        #---------------------------------------------------------------------#
        if (length(background_grid) < n_centers) {
            radius <- radius - (radius / 2)
        } else {
            radius <- radius + (radius / 2)
        }
        iter <- iter + 1
    }
    #-------------------------------------------------------------------------#
    # Next we need to find what are the indices of these barcodes 
    # within the full pixel image 
    # could use right_join and the likes but no 
    #-------------------------------------------------------------------------#
    
    return(background_grid)
}


hex_grid <- function(coordinates, n_centers, return_index = TRUE) {
    dims <- ceiling(sqrt(n_centers))
    range_y <- max(coordinates$y)
    range_x <- max(coordinates$x)
    y <- seq(0, range_y + 1, l = dims)
    x_short <- seq(0, range_x, l = dims)
    shift <- (x_short[2] - x_short[1]) / 2
    x_long <- seq(-shift, max(x_short) + shift, l = dims + 1)
    c_x <- c()
    c_y <- c()
    for (i in seq_along(y)){
        if (i %% 2 == 0) {
            c_x <- c(c_x, x_long)
            c_y <- c(c_y, rep(y[i], times = length(x_long)))
        } else {
            c_x <- c(c_x, x_short)
            c_y <- c(c_y, rep(y[i], times = length(x_short)))
        }
    }
    coord <- data.frame("x" = c_x, "y" = c_y)
    if (return_index) {
      background_grid <- RANN::nn2(data = coordinates[, c("x", "y")],
      query = coord,
      k = 1)
      return(background_grid$nn.idx[, 1])
    } else {
      return(coord)
    }
    
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
    segment <- unique(ter$Segment)
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

      tmp <- ter[ter$Segment == segment[i], ]
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
      ter$trial[ter$Segment == segment[i]] <- tmp

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


#' @importFrom imager imsplit  threshold split_connected where
#' @importFrom imagerExtra ThresholdML
#' @importFrom dplyr inner_join
select_similar <- function(img,
  coordinates,
  threshold = 1) {
  img <- img %>%
    imsplit("cc")

  pos <- lapply(img, function(x, threshold){
      ret <- ThresholdML(x, threshold)
      ret <- split_connected(ret)
      return(ret)
  }, threshold = threshold)
  coordinates$Segment <- 1
  count <- 2
  for (i in seq_along(pos)){
    for (j in seq_along(pos[[i]])) {
      tmp <- where(pos[[i]][[j]]) %>%
        inner_join(coordinates, by = c("x", "y"))
         coordinates$Segment[coordinates$barcodes %in% tmp$barcodes &
          coordinates$Segment == 1] <- count
        count <- count + 1
    }
  }
  all_ter <- unique(coordinates$Segment)
  coordinates$Segment <- seq_along(all_ter)[match(coordinates$Segment, all_ter)]
  return(coordinates)
}

#' @importFrom future.apply future_lapply
connected_pixels <- function(clusters,
  embeddings,
  k = 6,
  threshold = 0.90,
  verbose = TRUE) {
    message_switch("connect_pixel", verbose)
    #-------------------------------------------------------------------------#
    # First we need to get super pixel centers 
    #-------------------------------------------------------------------------#
    center_pixels <- sort(unique(clusters$Segment))
    centers <- future_lapply(center_pixels, function(center, segments){
        x <- median(segments$x[segments$Segment == center])
        y <- median(segments$y[segments$Segment == center])
        df <- data.frame("x" = x, "y" = y)
        rownames(df) <- center
        return(df)
    }, segments = clusters) %>% do.call("rbind", .)
    #-------------------------------------------------------------------------#
    # Next we get the nearest neighbors
    #-------------------------------------------------------------------------#
    knn <- RANN::nn2(centers, k = k + 1)$nn.idx
    rownames(knn) <- rownames(centers)
    #-------------------------------------------------------------------------#
    # Next we intialise a graph and then compute correlation
    #-------------------------------------------------------------------------#
    graph <- populate_graph(knn)
    graph$cor <- 0
    for (i in seq_len(nrow(graph))) {
        c1 <- embeddings[clusters$barcodes[clusters$Segment == graph$e1[i]], ]
        if (!is.null(nrow(c1))) {
            c1 <- apply(c1, 2, mean)
        }
        c2 <- embeddings[clusters$barcodes[clusters$Segment == graph$e2[i]], ]
        if (!is.null(nrow(c2))) {
            c2 <- apply(c2, 2, mean)
        }
        graph$cor[i] <- cor(c1, c2, method = "pearson")
    }
    #-------------------------------------------------------------------------#
    # Then we interatively pool pixels together under the a transitive 
    # correlation assumption i.e if A cor B and B cor with C then A cor C
    #-------------------------------------------------------------------------#
    initial_pixels <- unique(graph$e1)
    total_pool <- c()
    segments <-  list()
    count <- 1
    while (length(initial_pixels > 0)) {
        start_pixel <- sample(initial_pixels, size = 1)
        pool <- graph$e2[graph$e1 == start_pixel & graph$cor >= threshold]
        inter <- pool
        total_pool <- c(total_pool, pool)
        converge <- FALSE
        #browser()
        while (!converge) {
            if (length(inter) == 1) {
                segments[[count]] <- pool
                initial_pixels <- initial_pixels[!initial_pixels %in% pool]
                count <- count + 1
                converge <- TRUE
            } else {
                new_pool <- unique(graph$e2[graph$e1 %in% inter &
                  graph$cor >= threshold])
                overlap <- new_pool %in% pool & !new_pool %in% total_pool
                if (sum(overlap) != length(new_pool)) {
                  pool <- unique(c(pool, new_pool[!overlap]))
                  
                  inter <- unique(new_pool[!overlap])
                  converge <- FALSE
                } else {
                  segments[[count]] <- pool
                  total_pool <- c(total_pool, pool)
                  count <- count + 1
                  initial_pixels <- initial_pixels[!initial_pixels %in% pool]
                  converge <- TRUE
                }
            }
        }
    }
    
    for (seg in seq_along(segments)){
        loc <- clusters$Segment %in% segments[[seg]]
        clusters$Segment[loc] <- seg
    }
    return(clusters)
}