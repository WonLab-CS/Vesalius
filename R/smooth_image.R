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
