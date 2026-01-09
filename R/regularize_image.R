
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
    # Param validity pre-future dispatch
    # error handling with future will throw and error and a waring 
    # which makes it messy in unit testing scenarios
    #--------------------------------------------------------------------------#
    if (lambda <= 0) {
      stop("Lambda needs to be a positive numeric")
    }
    if (niter < 1) {
      stop("niter needs to be at least one")
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


