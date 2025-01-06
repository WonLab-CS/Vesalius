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



