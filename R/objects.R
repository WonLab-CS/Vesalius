################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#

#------------------------------------------------------------------------------#
# Vesalius assay log - keepstrack of what is done to each individual assay
# Specific view method should be developped.
#------------------------------------------------------------------------------#

setClass("assay_log",
    slots = list(assay = "character",
        tiles = "list",
        embeddings = "list",
        active = "list",
        territories = "list",
        DEG = "list",
        counts  = "list",
        meta_info = "list")
)
## TO DO add show method 

#------------------------------------------------------------------------------#
# vesalius log keeps track of what is done overall
# this can also include meta information taken from summaries over each assay
#------------------------------------------------------------------------------#
setClass("log",
    slots = list(assays = "list",
        meta_info = "list",
        log = "assay_log")
)

## TO DO add show method 
#------------------------------------------------------------------------------#
# Main vesalius object
# This object is the high level object that contain each assay
# We are exclusively working with spatial omics for now
# but this approach means we can add other assys later if needed
#------------------------------------------------------------------------------#

setClass("vesaliusObject",
    slots = list(assays = "list",
        meta_info = "list",
        history = "log"),
    validity = function(object) {
        if (!is(object@assays, "list")) {
            stop("Unsupported assay Format")
        }
        if (!is(object@meta_info, "list")) {
            stop("Unsupported meta info Format")
        }
        if (!is(object@history, "log")) {
            stop("Unsupported log Format")
        }
        return(TRUE)
    }

)


build_vesalius_object <- function(coordinates,
    counts,
    assay = "spatial_omics",
    adjust_coordinates = c("origin", "norm"),
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
    log <- new("log", assays = as.list(names(input)))
    ves <- new("vesaliusObject",
        assays = input,
        history = log)
    simple_bar(verbose)
    return(ves)
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
        embeddings = "list",
        territories = "territory",
        DEG = "list",
        counts  = "list",
        image = "list",
        log = "assay_log"),
    validity = function(object) {
        if (!is(object@tiles, "data.frame")) {
            stop("Unsupported Tile Format")
        }
        if (!is(object@embeddings, "list")) {
            stop("Unsupported Embeddings Format")
        }
        if (!is(object@territories, "territory")) {
            stop("Unsupported territories Format")
        }
        if (!is(object@DEG, "list")) {
            stop("Unsupported DEG Format")
        }
        if (!is(object@counts, "list")) {
            stop("Unsupported Count Format")
        }
        if (!is(object@image, "list")) {
            stop("Unsupported image Format")
        }
        if (!is(object@log, "assay_log")) {
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
    # add assay log
    #--------------------------------------------------------------------------#
    log <- new("assay_log", assay = assay)
    return(new("vesalius_assay",
        counts = counts,
        tiles = coordinates,
        log = log
        ))
}

## TO DO add show method