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


build_vesalius_object <- function(coordinates,
    counts,
    assay = "spatial_omics",
    adjust_coordinates = "origin",
    verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # check input and compare them to each other.
    #--------------------------------------------------------------------------#
    input <- check_inputs(counts,
        coordinates,
        assay,
        adjust_coordinates,
        verbose)
   
    vesalius <- new("vesaliusObject",
        assays = input)
    commit <- create_commit_log(vesalius = vesalius,
      match = as.list(match.call()),
      default = formals(build_vesalius_object),
      assay = names(input))
    vesalius <- commit_log_to_vesalius_assay(vesalius,
      commit,
      names(input))
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

## at the moment internal exclusive but leaving space to expand
build_vesalius_assay <- function(counts,
    coordinates,
    assay = NULL,
    verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # no check at the monent as i won' export this function for now
    # in the future we can easily add the checks here as well.
    #--------------------------------------------------------------------------#
    counts <- list(counts)
    #--------------------------------------------------------------------------#
    # adding some tags - the first one give the name of the count matrix
    # here we only have raw counts so raw
    # the comment section add a general comment attribute that will tell
    # us which is the active count matrix i.e. whcih one we should use
    #--------------------------------------------------------------------------#
    names(counts) <- "raw"
    comment(counts) <- "raw"
    #--------------------------------------------------------------------------#
    # add assay log TODO
    #--------------------------------------------------------------------------#
    
    return(new("vesalius_assay",
        counts = counts,
        tiles = coordinates
        ))
}

## TO DO add show method