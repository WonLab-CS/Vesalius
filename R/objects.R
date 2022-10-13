################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#

#------------------------------------------------------------------------------#
# Main vesalius object
# This object is the high level object that contain each assay
# We are exclusively working with spatial omics for now
# but this approach means we can add other assys later if needed
#------------------------------------------------------------------------------#
#' The vesaliusObject class
#' 
#' The vesaliusObject class is a high level obect that stores vesalius assay 
#' objects and provides an overview of all assays present. The purpose 
#' of this class is to provide a convenient way to store and process 
#' assays together. This object handles any type of spatial assay and alows
#' to run analysis on multiple assays simultaneously.
#' 
#' @slot assays list of \code{\link{vesalius_assay}}
#' @slot history list containing overview of all analysis trials attempted
#' on each assay present. 
#' @name vesaliusObject-class
#' @rdname vesaliusObject-class
#' @exportClass vesaliusObject
setClass("vesaliusObject",
    slots = list(assays = "list",
        history = "list"),
    validity = function(object) {
        if (!is(object@assays, "list")) {
            stop("Unsupported assay Format")
        }
        if (!is(object@history, "list")) {
            stop("Unsupported log Format")
        }
        return(TRUE)
    }

)

#' build a vesalius object
#' 
#' build a vesalius object from coordinate files and counts matrices
#' 
#' @param coordinates data.frame or list of data.frame in the format 
#' barcodes, x, y (see details).
#' @param counts matrix, sparse matrix or list of matrices or sparse 
#' matrices containing counts (see details).
#' @param assays character vector containing names of the assays 
#' (see details).
#' @param adjust_coordinates character of one of the following
#' "origin" or "norm" (see details).
#' @param verbose logical indicating if progress message should be
#' outputed or not. 
#' @details You can build a vesaliusObject by either supplying 
#' a single matrix/ coordinate combination or by providing a list
#' containing multiple matrices and coordinates. Please note that 
#' both list should be the same length. 
#' 
#' Along side this input data, you can provide a list of names 
#' to your assays. If none are provided or if not enough names
#' are provided, Vesalius will generate a set of names based on 
#' the default assay name "spatial_omics". Please ensure that 
#' the name you provide are unique to each assay.
#' 
#' You can decide if you want to adjust the coordinates to the 
#' origin i.e remove unnecessary space or normalise the coordinates.
#' Norm is not recommened at the moment. 
#' @return A vesaliusObject
#' @export 
#' @examples 
#' data(vesalius)
#' # Single assay object
#' ves <- build_vesalius_object(coordinates, counts)
#' # Building multi assay object
#' coord <- list(coordinates, coordinates)
#' co <- list(counts, counts)
#' assay_names <- c("my_assay", "my_second_assay")
#' ves <- build_vesalius_object(coord, co, assay_names)

build_vesalius_object <- function(coordinates,
    counts,
    assays = "spatial_omics",
    adjust_coordinates = "origin",
    verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # first we check if we have the same length lists 
    # if we only have a single object parsed  we put this into a list
    # in this case it should be equal to 1 anyway
    # We also check if assays have the same length 
    #--------------------------------------------------------------------------#
    if (!is(coordinates, "list")) {
        coordinates <- list(coordinates)
    }
    if (!is(counts, "list")) {
        counts <- list(counts)
    } 
    if (length(counts) != length(coordinates)) {
        stop("Number of count matrices and coordinates do not match")
    }
    #--------------------------------------------------------------------------#
    # since we already know that there are the same number of counts and
    # coordinates we can use only one length value here
    #--------------------------------------------------------------------------#
    
    assays <- check_assays(assays, length(counts), verbose)
    #--------------------------------------------------------------------------#
    # check input and compare them to each other.
    #--------------------------------------------------------------------------#
    input <- vector("list", length(counts))
    names(input) <- assays
    for (i in seq_along(counts)){
        input[[i]] <- build_vesalius_assay(
            coordinates = coordinates[[i]],
            counts =  counts[[i]],
            assay = assays[[i]],
            adjust_coordinates = adjust_coordinates,
            verbose = verbose)
    }
    
   
    vesalius <- new("vesaliusObject",
        assays = input)
    log_summary <- get_log_summary(vesalius)
    vesalius@history <- log_summary
    simple_bar(verbose)
    return(vesalius)
}

## TO DO add show method




#------------------------------------------------------------------------------#
# Object class used for each assay
# here every assay has the same format
#------------------------------------------------------------------------------#
setClassUnion("territory", c("data.frame", NULL))
setClassUnion("embeddings", c("list", NULL))
setClassUnion("DEG", c("list", NULL))


#' The vesalius_assay class 
#' 
#' The vesalius_assay class is the functional unit of vesalius. Each assay is 
#' stored within this class and it contains all the required information to 
#' run analysis on your assay of choice. In this object, you can find spatial
#' tiles, image embeddings, spatial territories, differentially expressed genes
#' (DEG), count matrices (raw and normalised), microscopy images (if present) and
#' a functional log that lets you see what had been run on this object.
#' 
#' We recommend using vesaliusObjects instead of vesalius_arrays but both will
#' work with all exported functions from the vesalius package. 
#' 
#' @slot tiles data.frame containing spatial coordinates and pixels tiles once 
#' they have been computed
#' @slot embeddings list containing latent space embeddings in the form of 
#' data.frames.
#' @slot territories data.frame containing spatial color segments, spatial
#' territories, or layers.
#' @slot DEG list of data.frame for each differentially gene expression trial
#' @slot image list containing associated microscopy images (NOT implemented)
#' @slot log list containing analysis history of the object. 
#' 
#' @name vesalius_assay-class
#' @rdname vesalius_assay-class
#' @exportClass vesalius_assay

setClass("vesalius_assay",
    slots = list(tiles = "data.frame",
        embeddings = "embeddings",
        territories = "territory",
        DEG = "DEG",
        counts  = "list",
        image = "list",
        log = "list"),
    validity = function(object) {
        if (!is(object@tiles, "data.frame")) {
            stop("Unsupported Tile Format")
        }
        if (!is(object@embeddings, "embeddings")) {
            stop("Unsupported Embeddings Format")
        }
        if (!is(object@territories, "territory")) {
            stop("Unsupported territories Format")
        }
        if (!is(object@DEG, "DEG")) {
            stop("Unsupported DEG Format")
        }
        if (!is(object@counts, "list")) {
            stop("Unsupported Count Format")
        }
        if (!is(object@image, "list")) {
            stop("Unsupported image Format")
        }
        if (!is(object@log, "list")) {
            stop("Unsupported log Format")
        }
        return(TRUE)
    }

)

#' build vesalius assay object
#' 
#' build a simple vesalius assay object from single count matrix and spatial
#' coordinate pair.
#' 
#' @param coordinates data.frame in the format 
#' barcodes, x, y.
#' @param counts matrix, sparse matrix containing counts .
#' @param assay character vector containing names of the assays 
#' (see details).
#' @param adjust_coordinates character of one of the following
#' "origin" or "norm" (see details).
#' @param verbose logical indicating if progress message should be
#' outputed or not. 
#' @details 
#' 
#' Along side this input data, you can provide a name 
#' to your assay. If none are provided, 
#' Vesalius will generate a set of names based on 
#' the default assay name "spatial_omics".
#' 
#' You can decide if you want to adjust the coordinates to the 
#' origin i.e remove unnecessary space or normalise the coordinates.
#' Norm is not recommened at the moment. 
#' 
#' We highly recommend using \code{\link{vesaliusObjects}} instead! 
#' @return A vesalius_assay objecy
#' @export 
#' @examples 
#' data(vesalius)
#' # Single assay object
#' ves <- build_vesalius_assay(coordinates, counts)

build_vesalius_assay <- function(coordinates,
    counts,
    assay = "spatial_omics",
    adjust_coordinates = "origin",
    verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # checking inputs
    #--------------------------------------------------------------------------#
    input <- check_inputs(counts,
        coordinates,
        assay,
        adjust_coordinates,
        verbose)
    #--------------------------------------------------------------------------#
    # add assay log
    #--------------------------------------------------------------------------#
    
    vesalius_assay <- new("vesalius_assay",
        counts = input$counts,
        tiles = input$coordinates
        )
   
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(build_vesalius_assay),
      assay = input$assay)
    vesalius_assay <- update_vesalius_assay_log(vesalius_assay,
        commit = commit,
        assay = assay[[1L]])
    return(vesalius_assay)
}

## TO DO add show method