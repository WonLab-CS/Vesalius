################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#
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
#' (DEG), count matrices (raw and normalised), microscopy images (if present)
#' and a functional log that lets you see what had been run on this object.
#' @slot assay character assay name 
#' @slot tiles data.frame containing spatial coordinates and pixels tiles once
#' they have been computed
#' @slot embeddings list containing latent space embeddings in the form of
#' data.frames.
#' @slot active matrix containing active embedding data
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
    slots = list(assay = "character",
        tiles = "data.frame",
        embeddings = "embeddings",
        active = "matrix",
        territories = "territory",
        DEG = "DEG",
        counts  = "list",
        image = "list",
        log = "list"),
    validity = function(object) {
        if (!is(object@assay, "character")) {
            stop("Unsupported assay Format")
        }
        if (!is(object@tiles, "data.frame")) {
            stop("Unsupported Tile Format")
        }
        if (!is(object@embeddings, "embeddings")) {
            stop("Unsupported Embeddings Format")
        }
        if (!is(object@active, "matrix")) {
            stop("Unsupported Active Format")
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

setMethod("show",
    signature = "vesalius_assay",
    definition = function(object) {
        simple_bar(TRUE)
        cat(paste(object@assay,"as vesalius assay containing:\n\n"))
        tiles  <- get_tiles(object)
        if (any(colnames(tiles) == "origin")){
            n_indices <- sum(tiles$origin)
            cat(paste(n_indices,
            "spatial indices used to form pixel tiles \n"))
        } else {
            n_indices <- nrow(tiles)
            cat(paste(n_indices, "spatial indices\n"))
        }
        counts <- get_counts(object)
        cat(paste(nrow(counts), "observations in the count matrix\n"))

        simple_bar(TRUE)
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
#' Norm is not recommened at the moment
#'
#' @return A vesalius_assay objecy
#' @export
#' @examples
#' data(vesalius)
#' # Single assay object
#' ves <- build_vesalius_assay(coordinates, counts)
#' @importFrom methods new

build_vesalius_assay <- function(coordinates,
    counts,
    assay = "spatial_omics",
    adjust_coordinates = "origin",
    verbose = TRUE) {
    simple_bar(verbose)
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
        assay = input$assay,
        counts = input$counts,
        tiles = input$coordinates
        )
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(build_vesalius_assay))
    vesalius_assay <- commit_log(vesalius_assay,
        commit = commit,
        assay = input$assay)
    simple_bar(verbose)
    return(vesalius_assay)
}

## TO DO add show method