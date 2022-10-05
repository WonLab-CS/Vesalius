################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#


#' Apply iterative smoothing to Vesalius images
#' @param image a Vesalius data frame
#' (containing at least barcodes, x, y, cc, value) or a cimg object.
#' @param method character describing smoothing method to use "median" ,
#' "iso"  or "box"
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
#' @param acrossLevels character - method used to account for multiple smoothing
#' levels (see details). Select from: "min","mean", "max"
#' @param invert logical - If TRUE, colours will be inverted i.e. 1 - colorValue
#' (background set to 1 instead of 0).
#' @param verbose logical - progress message output.
#' @details Images produced by Vesalius should be smoothed in order to retrieve
#' territories. The \code{smoothArray} function provides multiple ways to smooth
#' an image array.
#'
#' First, consider the smoothing method : median, iso ,or box. Median and box
#' both require the \code{box} argument to be set according to your needs.
#' Median converts all pixels surronding a center pixel to the median value of
#' the selected box. Box converts all pixels to the value of the center pixel.
#' Iso applies a gaussian blur and is modulated by the \code{sigma} argument.
#'
#' Both the \code{box} and \code{sigma} arguments can take more than one value.
#' Both may take a vector of positive numeric values. During every smoothing
#' iteration, Vesalius will map this vector over multiple copies of the image
#' and summarize the final smoothed image by taking values across different
#' levels. The levels are: min, max, and mean. Min will take the lowest pixel
#' value (in all three colour channels) across all levels. Max will take the
#' the highest pixel value (in all three colour channels) across levels. Mean
#' will take the average pixel value (in all three colour channels) across
#' levels.
#'
#' The optimal smoothing values will depend on the question at hand and which
#' territories you wish to extract.
#'
#' NOTE: this function is applied internally in the
#' \code{iterativeSegmentation.array} function.
#'
#' For more information, we suggest reading through the imager vignette.
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile",...
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
#' image <- smoothArray(image)
#' # across multiple levels
#' image <- smoothArray(image,method = c("iso"),
#' acrossLevels = "mean",sigma = seq(0.5,1.5,l = 10)))
#' imagePlot(image)
#' }
smooth_array <- function(vesalius,
  dimensions = seq(1, 3),
  embedding = "last",
  method = c("median", "iso", "box"),
  iter = 1,
  sigma = 1,
  box = 20,
  threshold = 0,
  neuman = TRUE,
  gaussian = TRUE,
  na.rm = FALSE,
  across_levels = "min",
  cores = 1,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # shifting format
    # ves_to_c => format.R
    #--------------------------------------------------------------------------#
    simple_bar(verbose)
    if (is(vesalius)[1L] == "vesaliusObject") {
        images <- ves_to_c(object = vesalius,
          embed = embedding,
          dims = dimensions,
          verbose = verbose)
    } else {
      stop("Unsupported format to smooth_array function")
    }

    #--------------------------------------------------------------------------#
    # We can smooth over different arrays in parallel - I hope...
    # this might lead to some issue with low level C
    # Leave as is for now - but could be better to set up some check for
    # core allocation. Don't want to waste time assigning threads if not
    # required.
    #--------------------------------------------------------------------------#
    message_switch("smooth", verbose)
    images <- parallel::mclapply(images, internal_smooth,
      method = method,
      iter = iter,
      sigma = sigma,
      box = box,
      threshold = threshold,
      neuman = neuman,
      gaussian = gaussian,
      na.rm = na.rm,
      across_levels = across_levels,
      mc.cores = cores)
    #--------------------------------------------------------------------------#
    # shifting format
    # c_to_ves => format.R
    #--------------------------------------------------------------------------#
    vesalius <- c_to_ves(images,
      vesalius,
      dimensions,
      embed = embedding,
      verbose = verbose)
    vesalius <- update_vesalius(vesalius = vesalius,
      data = vesalius@activeEmbeddings,
      slot = "activeEmbeddings",
      commit = as.list(match.call()),
      defaults = as.list(args(smooth_array)),
      append = FALSE)

    simple_bar(verbose)
    return(vesalius)

}

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
            image <- imrotate(image, 180)
        }
        for (j in method) {
          image <- switch(j,
                     "median" = map_il(box, ~medianblur(image,
                        .,
                        threshold)),
                     "iso" = map_il(sigma, ~isoblur(image,
                        .,
                        neuman,
                        gaussian,
                        na.rm)),
                     "box" = map_il(box, ~boxblur(image,
                        .,
                        neuman)))
          image <- switch(across_levels,
                     "min" = parmin(image),
                     "max" = parmax(image),
                     "mean" = average(image))
        }
        if (i %% 2 == 0) {
            image <- imrotate(image,180)
        }

    }
    return(image)
}



#' equalizeHistogram image enhancement via colour histogram equalization.
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,...)
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
#' @param invert logical - If TRUE, colours will be inverted i.e. 1 - colorValue
#' (background set to 1 instead of 0).
#' @param verbose logical - progress message output.
#' @details Histogram equalization ensures that image details are amplified.
#' In turn, territories may be extract with greater precision. We recommend
#' balancing the histogram prior to smoothing.
#'
#' For further details on each method described here, please refer to
#' \href{imagerExtra Vignette}{https://cran.r-project.org/web/packages/imagerExtra/vignettes/gettingstarted.html}
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile",...
#' @examples
#' \dontrun{
#' data(vesalius)
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' image <- equalizeHistogram(image)
#' imagePlot(image)
#' }

equalize_histogram <- function(vesalius,
  dimensions = seq(1, 3),
  embedding = "last",
  type = "BalanceSimplest",
  N = 1,
  smax = 1,
  sleft = 1,
  sright = 1,
  lambda = 0.1,
  up = 100,
  down = 10,
  cores = 1,
  verbose = TRUE) {
    .simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # shifting format
    # ves_to_c => format.R
    #--------------------------------------------------------------------------#
    if (is(vesalius)[1L] == "vesaliusObject") {
        images <- ves_to_c(object = vesalius,
          embed = embedding,
          dims = dimensions,
          verbose = verbose)
    } else {
      stop("Unsupported format to equalize_histogram function")
    }
    message_switch("eq", verbose)
    #--------------------------------------------------------------------------#
    # Equalizing histogram - depending on the image different methods will
    # work differently. It worth keeping in mind that these are made and
    # optimized for real images.
    #--------------------------------------------------------------------------#

    images <- switch(type,
           "EqualizePiecewise" = parallel::mclapply(images, EqualizePiecewise,
             N, mc.cores = cores),
           "BalanceSimplest" = parallel::mclapply(images, BalanceSimplest,
             sleft, sright, range = c(0, 1), mc.cores = cores),
           "SPE" = parallel::mclapply(image, SPE,
             lambda, mc.cores = cores),
           "EqualizeDP"  = parallel::mclapply(images, EqualizeDP,
             down, up, mc.cores = cores),
           "EqualizeADP" = parallel::mclapply(images, EqualizeADP,
             mc.cores = cores),
           "ECDF" = parallel::mclapply(images, ecdf_eq,
             mc.cores = cores))
    #--------------------------------------------------------------------------#
    # shifting format
    # c_to_ves => format.R
    #--------------------------------------------------------------------------#
    vesalius <- c_to_ves(images,
      vesalius,
      dimensions,
      embed = embedding,
      verbose = verbose)
    vesalius <- update_vesalius(vesalius = vesalius,
      data = vesalius@activeEmbeddings,
      slot = "activeEmbeddings",
      commit = as.list(match.call()),
      defaults = as.list(args(equalize_histogram)),
      append = FALSE)
    simple_bar(verbose)
    return(vesalius)
}

# Internal function related to ECDF historgram eq.
# It is just for formatting really
ecdf_eq <- function(im){
   return(as.cimg(ecdf(im)(im), dim = dim(im)))
}






#' regulariseImage denoise Vesalius images via variance regularization
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,tile,...)
#' @param lambda numeric - positive real numbers describing regularization
#' parameter (see details)
#' @param niter numeric - number of variance regularization iterations
#' (Default = 100)
#' @param method character - regularization method.
#' Select from: "TVL2.FiniteDifference","TVL1.PrimalDual",and "TVL2.PrimalDual"
#' (see details)
#' @param normalise logical - If TRUE, regularized colour values will be
#' min max normalised.
#' @param na.rm logical - If TRUE, NAs are removed
#' @param invert logical - If TRUE, colours will be inverted i.e. 1 - colorValue
#' (background set to 1 instead of 0).
#' @param verbose logical - progress message output.
#' @details Image regularization can be seen as a form of image denoising.
#' Details on each method can be found in the tvR package under the denoise2
#' function \href{tvR}{https://cran.r-project.org/web/packages/tvR/tvR.pdf}.
#'
#' We recommend using the TVL2.FiniteDifference method.
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile",etc
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
#' image <- regulariseImage(image)
#' }


regularise_image <- function(vesalius,
  dimensions = seq(1, 3),
  embedding = "last",
  lambda = 1,
  niter = 100,
  method = "TVL2.FiniteDifference",
  normalise = TRUE,
  na.rm = TRUE,
  cores = 1,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # shifting format
    # ves_to_c => format.R
    #--------------------------------------------------------------------------#
    if (is(vesalius)[1L] == "vesaliusObject") {
        images <- ves_to_c(object = vesalius,
          embed = embedding,
          dims = dimensions,
          verbose = verbose)
    } else {
      stop("Unsupported format to regularise_image function")
    }
    message_switch("reg", verbose)
    #--------------------------------------------------------------------------#
    # now we can do some var reg!
    # Will add the imager denoising function as well
    # regularisation is essentially denoising anyway
    #--------------------------------------------------------------------------#
    images <- parallel::mclapply(images, regularise,
      lambda,
      niter,
      method,
      normalise,
      mc.cores = cores)

    vesalius <- c_to_ves(images,
      vesalius,
      dimensions,
      embed = embedding,
      verbose = verbose)
    vesalius <- update_vesalius(vesalius = vesalius,
      data = vesalius@activeEmbeddings,
      slot = "activeEmbeddings",
      commit = as.list(match.call()),
      defaults = as.list(args(regularise_image)),
      append = FALSE)
    simple_bar(verbose)
    return(vesalius)
}

# Internal regularise function.

regularise <- function(img,
  lambda = 1,
  niter = 100,
  method = "TVL2.FiniteDifference",
  normalise = TRUE) {
    #--------------------------------------------------------------------------#
    # might need to change the normalise part
    # it doesn't seem it works well
    # or at all for that matter
    #--------------------------------------------------------------------------#
    img <- denoise2(as.matrix(img), lambda = lambda, niter = niter,
           method = method, normalize = FALSE)
    if (normalise) {
        img <- (img - min(img)) / (max(img) - min(img))
    }
    return(as.cimg(img))
}




#' iterativeSegmentation.array - image segmentation via iterative application
#' of kmeans clustering.
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,...)
#' @param colDepth integer or vector of positive integers.
#' Colour depth used for segmentation. (see details)
#' @param smoothIter integer - number of smoothing iterations
#' @param method character describing smoothing method to use "median" ,
#' "iso"  or "box"
#' @param acrossLevels character - method used to account for multiple smoothing
#' levels (see details). Select from: "min","mean", "max"
#' @param sigma numeric - standard deviation associated with isoblur (Gaussian)
#' @param box numeric describing box size (centered around a center pixel)
#' for smoothing
#' @param threshold numeric - discard pixels that are too low in value (cutoff
#' threshold only applied in box/median blurs).
#' @param neuman logical describing If Neumann boundary conditions should be
#' used, Dirichlet otherwise (default true, Neumann)
#' @param gaussian logical - use gaussian filter
#' @param useCenter logical - If TRUE, only the center pixel value will be used
#' during segmentation. If FALSE, all pixels will be used (see details)
#' @param na.rm logical describing if NA values should be removed
#' @param invert logical - If TRUE, colours will be inverted
#' i.e. 1 - colourValue (background set to 1 instead of 0).
#' @param verbose logical - progress message output.
#' @details Once images have been produced by Vesalius, applying image
#' segmentation ensure a reduction in colour complexity. The colDepth argument
#' describes the number of colour segments to select in the image. It should
#' be noted that this function can take a vector of positive integers. This
#' becomes handy if you wish to gradually decrease the colour complexity. It is
#' generally recommend to use a vector of decreasing colDepth values.
#'
#' The segmentation process is always proceeded by smoothing and multiple rounds
#' may be applied. The segmentation process uses kmeans clustering with k number
#' of cluster being represented by colDepth values.
#'
#' The segmentation process can be applied on all pixel or only center pixels.
#' Center pixels are described by the value in the "tile" column in a Vesalius
#' data frame. Every pixel with a value of 1 in the "tile" column corresponds
#' to the original barcode location in the Spatial Transcriptomic Assay
#' (i.e. Original coordinates) before tesselation and rasterisation
#' (see \code{buildImageArray})
#
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile","cluster".
#' Cluster represents the colour segement the pixel belongs to.
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
#' # Multiple segmentation rounds & multiple smoothing rounds
#' image <- iterativeSegmentation.array(image,smoothIter = 3,
#'       colDepth = seq(12,8, by = -2))
#' # smoothing across levels
#' image <- iterativeSegmentation.array(image, smoothIter = 5,
#' method = c("iso"), acrossLevels = "mean",sigma = seq(0.5,1.5,l = 10))
#' }

image_segmentation <- function(vesalius,
  method = c("kmeans", "SISKmeans", "SIS"),
  embedding = "last",
  col_depth = 10,
  dimensions = seq(1, 3),
  smooth_iter = 1,
  smooth_type = c("median", "iso", "box"),
  across_levels = "min",
  sigma = 1,
  box = 20,
  threshold = 0,
  neuman = TRUE,
  gaussian = TRUE,
  na.rm = TRUE,
  use_center = TRUE,
  cores = 1,
  verbose = TRUE) {
  simple_bar(verbose)
  #----------------------------------------------------------------------------#
  # Parsing vesalius object so we can recontruct it internally and not need to
  # rebuild intermediates and shift between formats...
  #----------------------------------------------------------------------------#
  image <- switch(method[1L],
    "kmeans" = vesalius_kmeans(vesalius,
      col_depth = col_depth,
      dims = dimensions,
      embedding = embedding,
      smooth_iter = smooth_iter,
      method = smooth_type,
      across_levels = across_levels,
      sigma = sigma,
      box = box,
      threshold = threshold,
      neuman = neuman,
      gaussian = gaussian,
      na.rm = na.rm,
      use_center = use_center,
      verbose = verbose),
    "SISKmeans" = .sis_kmeans(vesalius, dimensions, cores, verbose),
    "SIS" = .sis(vesalius, dimensions, cores, verbose))

  vesalius <- update_vesalius(vesalius = vesalius,
    data = image$vesalius@activeEmbeddings,
    slot = "activeEmbeddings",
    commit = as.list(match.call()),
    defaults = as.list(args(image_segmentation)),
    append = FALSE)
  vesalius <- updateVesalius(vesalius = vesalius,
    data = image$clusters,
    slot = "territories",
    commit = as.list(match.call()),
    defaults = as.list(args(image_segmentation)),
    append = TRUE)
  simple_bar(verbose)
  return(vesalius)
}

vesalius_kmeans <- function(vesalius,
  dims = seq(1, 3),
  col_depth = 10,
  embedding = "last",
  smooth_iter = 1,
  method = c("median", "iso", "box"),
  across_levels = "min",
  sigma = 1,
  box = 20,
  threshold = 0,
  neuman = TRUE,
  gaussian = TRUE,
  na.rm=FALSE,
  use_center = TRUE,
  cores = 1,
  verbose = TRUE) {
  

  if (is(vesalius)[1L] == "vesaliusObject") {
      images <- ves_to_c(object = vesalius,
        embed = embedding,
        dims = dims,
        verbose = verbose)
  } else {
    stop("Unsupported format to image_segmentation function")
  }
  #----------------------------------------------------------------------------#
  # Segmenting image by iteratively decreasing colour depth and smoothing
  # Well that's true only the user parse an array of decreasing values
  # add some sanity checks maybe?
  #----------------------------------------------------------------------------#
  for (j in seq_along(col_depth)) {
    #--------------------------------------------------------------------------#
    # Lets smooooooth this image
    # insert https://www.youtube.com/watch?v=4TYv2PhG89A&ab_channel=SadeVEVO
    #--------------------------------------------------------------------------#
    images <- parallel::mclapply(images, internal_smooth,
      method = method,
      iter = smooth_iter,
      sigma = sigma,
      box = box,
      threshold = threshold,
      neuman = neuman,
      gaussian = gaussian,
      na.rm = na.rm,
      across_levels = across_levels,
      mc.cores = cores)
    #--------------------------------------------------------------------------#
    # Let's get the value vector out of each gray scale image
    #--------------------------------------------------------------------------#
    colours <- lapply(images,function(img) {
        return(as.data.frame(img)$value)
    })
    colours <- as.matrix(do.call("cbind", colours))
    if (use_center) {
        colours <- cbind(as.data.frame(images[[1L]])[, c("x", "y")],
                         as.data.frame(colours)) %>%
                   right_join(vesalius@tiles,by = c("x", "y")) %>%
                   filter(origin == 1)
        coord <- colours[, c("x", "y")]
        colours <- colours[, !colnames(colours) %in%
          c("x", "y", "barcodes", "origin")] %>% as.matrix()
    }
    message_switch("seg", verbose, seg = j)
    cat("\n")
    #--------------------------------------------------------------------------#
    # Now lets cluster colours
    # Remember that here colours are your new active embedding values
    #--------------------------------------------------------------------------#
    km <- kmeans(colours, col_depth[j], iter.max = 200, nstart = 50)

    cluster <- km$cluster
    kcenters <- km$centers

    for (i in seq_len(ncol(colours))) {
          colours[, i] <- kcenters[cluster, i]
    }

    #------------------------------------------------------------------------#
    # This section is so we don't end up having barcode associated with
    # multiple clusters/segments.
    # There are a lot of pixel associated to each barcode. When smoothing
    # it is possible that some pixel will be assigned to a different segment
    # than the center pixel. This is a pain to deal with in later steps
    # instead: take mean value per barcode and take the cluster that is the
    # most represented as the cluster value that we will use
    # this is mostly relevant when using useCentre =F and multiple colDepth
    # values.
    #------------------------------------------------------------------------#
    if (use_center) {
      clusters <- data.frame(coord, colours, cluster) %>%
        right_join(vesalius@tiles, by = c("x", "y"))
      embeds <- colnames(clusters)[!colnames(clusters) %in%
        c("x", "y", "cluster", "barcodes", "origin")]

    } else {
      clusters <- cbind(as.data.frame(
        images[[1L]])[, c("x", "y")], colours, cluster)

      embeds <- colnames(clusters)[!colnames(clusters) %in%
        c("x", "y", "cluster")]

      clusters <- right_join(clusters, vesalius@tiles, by = c("x", "y")) %>%
             group_by(barcodes) %>%
             mutate(across(all_of(embeds),mean),
              cluster = top_cluster(cluster)) %>%
             ungroup
    }


     images <- lapply(embeds, function(idx,cols) {
          tmp <- cols[, c("x", "y", as.character(idx))]
          colnames(tmp) <- c("x", "y", "value")
          return(as.cimg(tmp))
     }, cols = clusters)
    }
    #--------------------------------------------------------------------------#
    # Let's rebuild everything
    #--------------------------------------------------------------------------#

    vesalius <- c_to_ves(images,
      vesalius,
      dims,
      mbed = embedding,
      verbose = verbose)
    clusters <- clusters %>% 
      filter(origin == 1) %>%
      select(c("barcodes", "x", "y", "cluster")) %>%
      as.data.frame()

    if (!is.null(vesalius@territories)) {
      last <- get_last_entry(vesalius@territories)
      colnames(clusters) <- c(colnames(clusters)[seq_len(ncol(clusters) - 1)],
        paste0("Segment_Trial_", last + 1))
    } else {
      colnames(clusters) <- c(colnames(clusters)[seq_len(ncol(clusters) - 1)],
        "Segment_Trial_1")
    }
  
  return(list("vesalius" = vesalius, "clusters" = clusters))
}

top_cluster <- function(cluster) {
    top <- table(cluster)
    top <- names(top)[order(top, decreasing = TRUE)]
    return(as.numeric(top[1]))
}

## NOTE this might become redundant if I implement the other segmentation method
sis_kmeans <- function(image, col_depth, box) {
  #----------------------------------------------------------------------------#
  # First we need to convert cimg array to simple array
  #----------------------------------------------------------------------------#
  img <- as.array(image)
  dim(img) <- c(nrow(img), ncol(img))
  img <- replicate(3, img)
  #----------------------------------------------------------------------------#
  # Next we can run the super pixel segmentation
  #----------------------------------------------------------------------------#
  init <- Image_Segmentation$new()
  spx_km <- init$spixel_segmentation(input_image = image,
    superpixel = 600,
    AP_data = TRUE,
    use_median = TRUE,
    sim_wL = 3,
    sim_wA = 10,
    sim_wB = 10,
    sim_color_radius = box,
    kmeans_method = "kmeans",
    kmeans_initializer = "kmeans++",
    kmeans_num_init = col_depth,
    kmeans_max_iters = 100,
    verbose = FALSE)
  return(spx_km)
}

sis <- function(image) {
  #----------------------------------------------------------------------------#
  # First we need to convert cimg array to simple array
  #----------------------------------------------------------------------------#
  img <- as.array(image)
  dim(img) <- c(nrow(img), ncol(img))
  img <- replicate(3, img)
  #----------------------------------------------------------------------------#
  # Next we can run the super pixel segmentation
  #----------------------------------------------------------------------------#
  init <- Image_Segmentation$new()
  spx <- init$spixel_segmentation(input_image = image,
    superpixel = 600,
    AP_data = TRUE,
    use_median = TRUE,
    sim_wL = 3,
    sim_wA = 10,
    sim_wB = 10,
    sim_color_radius = 10,
    verbose = FALSE)
   return(spx)
}

#' isolating territories from segmented Vesalius images
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,tile) - Segmentation must have been applied beforehand!
#' @param method character describing barcode pooling method.
#' Currently, only "distance" availble
#' @param captureRadius numeric - proportion of maximum distance between
#' barcodes that will be used to pool barcodes together (range 0 - 1).
#' @param global logical - If TRUE, territories will be numbered across all
#' colour segments. If FALSE, territories will be numbered within each colour
#' segment.
#' @param minBar integer - minimum number of barcodes allowed in each territory
#' @param verbose logical - progress message output.
#' @details Segmented images in the form of a Vesalius formatted data frame are
#' further separated into territories. This is accomplished by pooling barcodes
#' that are associated with a colour cluster into territories based on the
#' distance between each barcode.
#'
#' First, \code{isolateTerritories.array} considers the maximum distance
#' between all beads. The captureRadius will define which proportion of this
#' distance should be considered.
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
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile","cluster","territory".
#' "cluster" represents the colour segment the pixel belongs to and "territory"
#' describe the tissue territory after pooling.
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
#' }


#### Requireds refactoring!!!!! Sanity check functions 
isolate_territories <- function(vesalius,
  method = c("distance"),
  trial = "last",
  capture_radius = 0.025,
  global = TRUE,
  minBar = 10,
  verbose = TRUE){

    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # Get stuff out as usual
    # Not super happy with this check
    # it's a bit messy - we might need to consider to do a whole sanity check
    # of inout data and see if that makes sense - this will include checking log
    #--------------------------------------------------------------------------#
    if (is(vesalius)[1L] == "vesaliusObject" && 
      !is.null(vesalius@territories)) {
      if (any(grepl(x= colnames(vesalius@territories),pattern = "Segment")) &&
        trial == "last") {
        #----------------------------------------------------------------------#
        # Making sure you get the last colour segment computed
        # not to be confused with territories
        #----------------------------------------------------------------------#
        trial <- grep(x = colnames(vesalius@territories),
          pattern = "Segment", value = TRUE)

        ter <- vesalius@territories[,c("barcodes","x","y",trial[length(trial)])]
        colnames(ter) <- c("barcodes","x","y","segment")

      } else if(any(grepl(x= colnames(vesalius@territories),pattern = "Segment")) &
                trial != "last") {
        if(length(grep(x = colnames(vesalius@territories),pattern = trial))==0){
            stop(paste(deparse(substitute(trial)),"is not in territory data frame"))
        }
        ter <- vesalius@territories[,c("barcodes","x","y",trial)]
        colnames(ter) <- c("barcodes","x","y","segment")

      } else {
        stop("No image segments have been computed!")
      }

    } else if(is(vesalius)[1L] == "vesaliusObject" & is.null(vesalius@territories)) {
        stop("No image segments have been computed!")
    } else {
        stop("Unsupported format to isolateTerritories function")
    }


    #--------------------------------------------------------------------------#
    # Compute real capture Radius
    #--------------------------------------------------------------------------#

    if(method[1L] == "distance"){
      captureRadius <- sqrt((max(ter$x)-min(ter$x))^2 +
                            (max(ter$y)-min(ter$y))^2) *
                            captureRadius
    }


    #--------------------------------------------------------------------------#
    # Creating new trial column name and adding it to input data
    # The input data here is a subset of the full territory df
    # we at least make sure that we are using the right input
    #--------------------------------------------------------------------------#
    if(any(grepl(x = colnames(vesalius@territories),pattern = "Territory"))){
      previous <- grep(x=colnames(vesalius@territories),
                                    pattern = "Territory",
                                    value = TRUE)
      m <- gregexpr('[0-9]+', previous)
      last <- max(as.numeric(unlist(regmatches(previous,m))))
      ter$trial <- 0
      newTrial <- paste0("Territory_Trial_",last+1)
      #colnames(ter) <- c("barcodes","x","y","segment",paste0("Territory_Trial_",last+1))

    } else {
      ter$trial <- 0
      newTrial <- "Territory_Trial_1"
      #colnames(ter) <- c("barcodes","x","y","segment","Territory_Trial_1")
    }
    #--------------------------------------------------------------------------#
    # Now we can dispatch the necessary information
    # and run the pooling algorithm
    #--------------------------------------------------------------------------#
    segment <- unique(ter$segment)

    for(i in seq_along(segment)){
      .terPool(i,verbose)
      #------------------------------------------------------------------------#
      # We don't need all data in this case only clusters and locations
      # We can just rebuild everything afterwards
      # At least we don't don't need to compute anything unnecessarily
      ## Note other method not in use for now
      # Might not be worthwile to implement them
      # Argument could be removed
      #------------------------------------------------------------------------#

      tmp <- ter[ter$segment == segment[i],]
      tmp <- switch(method[1L],
                    "distance" = .distancePooling.array(tmp,captureRadius,
                                                        minBar))
                  #  "neighbor" = .neighborPooling.array(tmp,captureRadius),
                  #  "watershed" = .watershedPooling.array(tmpImg))
      #------------------------------------------------------------------------#
      # Skipping if does not meet min cell requirements
      # filtering can be done later
      #------------------------------------------------------------------------#
      if(is.null(tmp))next()
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
    if(global){
        ter <- .globalise_territories(ter)
    }
    #--------------------------------------------------------------------------#
    #Now we can add it to the full territory df and update vesalius
    #--------------------------------------------------------------------------#

    colnames(ter) <- c(colnames(ter)[seq(1,length(colnames(ter))-1)], newTrial)
    vesalius <- .updateVesalius(vesalius=vesalius,
                                data=ter,
                                slot="territories",
                                commit = as.list(match.call()),
                                defaults = as.list(args(isolateTerritories)),
                                append=TRUE)
    cat("\n")
    .simpleBar(verbose)
    return(vesalius)
}




.distancePooling.array <- function(img,captureRadius,minBar){
    #--------------------------------------------------------------------------#
    # Select center point of each tile for only one channel
    # Dont need to run it for all channels
    #--------------------------------------------------------------------------#

    imgCopy <- img  %>% distinct(barcodes,.keep_all = TRUE)
    if(nrow(imgCopy)<1){return(NULL)}


    #--------------------------------------------------------------------------#
    # Compute distances
    #--------------------------------------------------------------------------#
    idx <- seq_len(nrow(imgCopy))
    distanceMatrix <- lapply(idx, function(idx,mat){
                            xo <- mat$x[idx]
                            yo <- mat$y[idx]
                            xp <- mat$x
                            yp <- mat$y
                            distance <- sqrt(((abs(xp-xo))^2 + (abs(yp-yo))^2))
                            return(distance)
    }, imgCopy)


    #--------------------------------------------------------------------------#
    # Buildan actual matrix
    #--------------------------------------------------------------------------#

    distanceMatrix <- do.call("rbind",distanceMatrix)
    colnames(distanceMatrix) <- imgCopy$barcodes
    rownames(distanceMatrix) <- imgCopy$barcodes

    #--------------------------------------------------------------------------#
    # If there is only one point
    # In this case you only have one territory as well
    # Don't need to any pooling
    #--------------------------------------------------------------------------#

    if(sum(dim(distanceMatrix))==2){
        return(1)
    }


    #--------------------------------------------------------------------------#
    # Pooling points together
    #--------------------------------------------------------------------------#
    barcodes <- imgCopy$barcodes
    territories <- list()
    count <- 1

    while(length(barcodes) >0){
         #---------------------------------------------------------------------#
         # First lets select a random barcode in the colour segment
         # And create a pool of barcodes to select based on capture radius
         # This first while loop checks if there are still barcodes
         # left in the colour segment
         #---------------------------------------------------------------------#
          tmp <- distanceMatrix[,sample(barcodes,1)]
          pool <- names(tmp)[tmp <= captureRadius]
          inter <- pool
          converge <- FALSE

          while(!converge){
            #------------------------------------------------------------------#
            # This while loop checks if all possible barcodes have been pooled
            # into the current territory
            #------------------------------------------------------------------#
              if(length(inter)==1){
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
                  newPool <- distanceMatrix[,inter]

                  newPool <- unique(unlist(lapply(seq_len(ncol(newPool)),
                                           function(idx,np,captureRadius){

                                            res <- rownames(np)[np[,idx] <=
                                                                captureRadius]
                                            return(res)
                                          },newPool,captureRadius)))
                  #------------------------------------------------------------#
                  # check which barcodes in the new pool overlap with
                  # the ones in the full pool
                  # If there is a perfect overlap then there are no new bacodes
                  # to pool into a territory
                  #------------------------------------------------------------#
                  overlap <- newPool %in% pool

                  if(sum(overlap) != length(newPool)){
                        #------------------------------------------------------#
                        # There are still some new barcodes to pool
                        # lets do some more looping then
                        #------------------------------------------------------#
                        pool <- unique(c(pool,newPool[!overlap]))
                        inter <- unique(newPool[!overlap])
                  } else {
                        #------------------------------------------------------#
                        # it is done ! no more barcodes for this territory
                        #------------------------------------------------------#
                        territories[[count]] <- pool
                        count <- count +1
                        barcodes <- barcodes[!barcodes %in% pool]

                        converge <- TRUE

                  }
              }
          }
      }
    #--------------------------------------------------------------------------#
    # Clean up drop outs
    #--------------------------------------------------------------------------#

    allTers <- img$trial

    for(ter in seq_along(territories)){
        loc <- img$barcodes %in% territories[[ter]]
        if(length(territories[[ter]]) <= minBar){
            allTers[loc] <- "isolated"
        } else {
            allTers[loc] <- ter
        }
    }
      return(allTers)
}