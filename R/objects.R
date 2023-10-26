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
#' @slot counts list that containing count matrices. Raw and normalised will
#' be stored here and named by the normalisation method used. 
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
        meta = "list",
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
        if (!is(object@meta, "list")) {
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
        #---------------------------------------------------------------------#
        # initialize the show
        # not that this will need to be updated with the empty object feature
        #---------------------------------------------------------------------#
        simple_bar(TRUE)
        cat(paste(object@assay, "as vesalius assay containing:\n\n"))
        #---------------------------------------------------------------------#
        # get image 
        #---------------------------------------------------------------------#
        img <- length(object@image)
        if (img > 0) {
            cat(paste(img, "image(s)\n"))
        } 
        #---------------------------------------------------------------------#
        # Get log modifications
        #---------------------------------------------------------------------#
        log <- length(object@log)
        cat(paste(log, "modifications applied to base object. \n"))
        #---------------------------------------------------------------------#
        # check and show tiles
        #---------------------------------------------------------------------#
        tiles  <- get_tiles(object)
        if (!is.null(tiles)) {
            cat("\n")
            if (any(colnames(tiles) == "origin")) {
                n_indices <- sum(tiles$origin)
                cat(paste(n_indices,
                    "spatial indices used to form pixel tiles. \n"))
            } else {
                n_indices <- nrow(tiles)
                cat(paste(n_indices, "spatial indices. \n"))
            }
        }
        #---------------------------------------------------------------------#
        # check and show counts
        #---------------------------------------------------------------------#
        counts <- get_counts(object, type = "all")
        if (comment(counts) != "empty" && comment(counts) != "joint") {
            cat("\n")
            active <- comment(counts)
            counts <- counts[[active]]
            cat(paste(nrow(counts), "observations in the",
                active,
                "count matrix. \n"))
        } else if (comment(counts) == "joint"){
            cat("\n")
            cat(paste(length(counts), "count matrices in object \n"))
            cat(paste(paste(names(counts), collapse = " "), "\n"))
        }

        #---------------------------------------------------------------------#
        # check embeddings
        #---------------------------------------------------------------------#
        all_embeds <- object@embeddings
        if (length(all_embeds) > 0 & length(all_embeds) < 4) {
            cat("\n")
            cat(paste(paste(names(all_embeds), sep = " ", collapse = ", "),
                "as embeddings. \n"))
            cat(paste("with", comment(all_embeds), "as active embedding. \n"))
        } else if (length(all_embeds) > 3) {
            cat("\n")
            cat(paste(length(all_embeds), "total embedding trials\n"))
            cat(paste("with", comment(all_embeds), "as active embedding. \n"))
        }
        #---------------------------------------------------------------------#
        # check territories
        #---------------------------------------------------------------------#
        ter <- summarise_territories(object, as_log = FALSE)
        if (any(unlist(ter) > 0)) {
            cat("\n")
            for (i in seq_along(ter)) {
                cat(paste(ter[i], names(ter)[i], "trials. \n"))
            }
        }

        #---------------------------------------------------------------------#
        # check for cell types
        #---------------------------------------------------------------------#
        cells <- search_log(object, "add_cells")
        if (length(cells) > 0) {
            cat("\n")
            if (is.null(tail(cells,1)[[1]]$add_name)){
                cell_col <- get_territories(object)
                n_cells <- length(unique(cell_col[,
                    tail(grep("Cells", colnames(cell_col)))]))
                cat(paste(n_cells, "cell types assigned. \n"))
            } else {
                add_name <- as.character(tail(cells,1)[[1]]$add_name)
                cell_col <- get_territories(object)
                n_cells <- length(unique(cell_col[,
                    tail(grep(add_name, colnames(cell_col)))]))
                cat(paste(n_cells, "cell types assigned. \n"))
            }
        }
        #---------------------------------------------------------------------#
        #check for DEGs
        #---------------------------------------------------------------------#
        deg <- object@DEG
        if (length(deg) != 0) {
            cat("\n")
            all_deg <- sapply(deg, nrow)
            cat(paste(sum(all_deg), "DEGs found across",
                length(all_deg)), "trials. \n")
        }
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
#' Default is NULL. See details.
#' @param counts matrix, sparse matrix containing counts.
#' Default is NULL. See details.
#' @param image connection string or image array
#' @param assay character vector containing names of the assays
#' (see details).
#' @param verbose logical indicating if progress message should be
#' outputed or not.
#' @details
#' The vesalius_assay constructor allows you to create
#' partial or complete vesalius_assay objects.
#' 
#' Partial objects contain only the coordinates.
#' 
#' Complete objects contain both the counts and the coordinates.
#'
#' The main purpose of using partial objects (or empty objects) is 
#' for you to be able to provide your own count matrix. 
#' This will be useful if you want to normalise your data in a 
#' way that is not provided by vesalius.
#' 
#' This approach of using your own data can also apply to embeddings.
#' If you have generated a set of latent space embeddings that you 
#' wish to use instead of those provided by vesalius, you can 
#' use the \code{\link{add_embeddings}} function. 
#' 
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
    counts = NULL,
    image = NULL,
    assay = "spatial_omics",
    scale = "auto",
    unit = "um",
    layer = 1,
    verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # checking inputs
    # sanity check to make sure that what is being parsed is what we expect
    # note that if counts is NULL - this function will return an empty list
    # but commented (tagged as "empty") for input$counts
    #--------------------------------------------------------------------------#
    input <- check_inputs(counts,
        coordinates,
        image,
        assay,
        layer,
        verbose)
    #--------------------------------------------------------------------------#
    # get scale from coordinates
    #--------------------------------------------------------------------------#
    if (scale == "auto") {
        message_switch("scale", verbose)
        scale <- calculate_scale(input$coordinates)
    }
    meta <- list("scale" = list("scale" = scale), "unit" = list("unit" = unit))
    #--------------------------------------------------------------------------#
    # add assay log
    #--------------------------------------------------------------------------#
    vesalius_assay <- new("vesalius_assay",
        assay = assay,
        counts = input$counts,
        tiles = input$coordinates,
        image = input$image,
        meta = meta)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(build_vesalius_assay))
    vesalius_assay <- commit_log(vesalius_assay,
        commit = commit,
        assay = assay)
    simple_bar(verbose)
    return(vesalius_assay)
}

#----------------------------/Vesalius Objects/--------------------------------#
#------------------------------------------------------------------------------#
# INTERNAL OBJECTS 
#------------------------------------------------------------------------------#

# setClassUnion("cost", c("data.frame", "matrix", NULL))
# setClassUnion("signal", c("matrix","Matrix", NULL))
# setClass("map_assay",
#     slot = list(
#         cost_matrix = "cost",
#         seed_coord = "data.frame",
#         query_coord = "data.frame",
#         seed_signal = "signal",
#         query_signal = "signal",
#         mapping = "data.frame",
#         scale = "numeric",
#         neighborhood = "character"))

# build_map_assay <- function(cost_matrix,
#     seed_coord,
#     query_coord,
#     seed_signal,
#     query_signal,
#     mapping,
#     scale,
#     neighborhood) {
#     map_assay <- new("map_assay",
#         cost_matrix = cost_matrix,
#         seed_coord = seed_coord,
#         query_coord = query_coord,
#         seed_signal = seed_signal,
#         query_signal = query_signal,
#         mapping = mapping,
#         scale = scale,
#         neighborhood = neighborhood)
#     return(map_assay)
# }