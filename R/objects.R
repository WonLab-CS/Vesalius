################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#

#------------------------------------------------------------------------------#
# Vesalius Log object
# Keep track of what has been done
#------------------------------------------------------------------------------#

setClass("log",
    slots = list(assay = "list",
        tiles = "list",
        embeddings = "list",
        activeEmbeddings = "list",
        territories = "list",
        DEG = "list",
        counts  = "list")
)

setMethod("show",
    signature = "log",
        definition = function(object) {
            simple_bar(TRUE)
            #------------------------------------------------------------------#
            # General log information
            #------------------------------------------------------------------#
            if (sum(unlist(slot_ppply(object, length))) == 0) {
               cat("No commits to log tree\n")
            } else {
                #--------------------------------------------------------------#
                # Last update
                #--------------------------------------------------------------#
                cat("Log tree\n")
                last <- slot_apply(object, length)
                for (i in seq_along(last)) {
                    if (last[[i]] == 0) {
                        shape <- "/"
                    } else {
                        shape <- rep("*", times = last[[i]])
                    }
                    cat("|", shape, names(last)[i], "\n")
                }
            }
            simple_bar(TRUE)
        }
)

#------------------------------------------------------------------------------#
# class unions
#------------------------------------------------------------------------------#

setClassUnion("ter", c("data.frame", "NULL"))




### Need to clean these up
### Increase the robustness of this
#------------------------------------------------------------------------------#
# Main vesalius object
#------------------------------------------------------------------------------#
setClass("vesaliusObject",
    slots = list(tiles = "data.frame",
        embeddings = "list",
        activeEmbeddings = "list",
        territories = "ter",
        DEG = "list",
        counts  = "list",
        log = "log"),
    prototype = c(tiles = data.frame(),
                     embeddings = list()),
    validity = function(object) {
        if (is(object@tiles, "data.rame")) {
            stop("Unsupported Tile Format")
        }
        if (is(object@embeddings ,"list")) {
            stop("Unsupported Embeddings Format")
        }
        if (is(object@activeEmbeddings, "list")) {
            stop("Unsupported Active Embeddings Format")
        }
        if (is(object@territories, "ter")) {
            stop("Unsupported territories Format")
        }
        if (is(object@DEG, "list")) {
            stop("Unsupported DEG Format")
        }
        if (is(object@counts, "list")) {
            stop("Unsupported Count Format")
        }
        if (is(object@log, "log")) {
            stop("Unsupported log Format")
        }
        return(TRUE)
    }

)

#------------------------------------------------------------------------------#
# Create new object function
# This one is on hold atm
#------------------------------------------------------------------------------#


build_vesalius_object <- function(coordinates,
    counts,
    assay = "ST",
    adjust_coordinates = c("origin", "norm")) {
    #--------------------------------------------------------------------------#
    # First we will create tiles
    # Definitely need some check here
    #--------------------------------------------------------------------------#
    coordinates <- check_coordinates(coordinates, adjust_coordinates)
    counts <- list(check_counts(counts))
    names(counts) <- "raw"
    assay <- list(assay)
    names(assay) <- "1"
    log <- new("log",
               assay = assay)
    ves <- new("vesaliusObject",
               tiles = coordinates,
               counts = counts,
               log = log)
    return(ves)
}



#------------------------------------------------------------------------------#
# vesalius objecy show method
#------------------------------------------------------------------------------#

setMethod("show",
    signature = "vesaliusObject",
    definition = function(object) {
        simple_bar(TRUE)
        #----------------------------------------------------------------------#
        ### To be modified later
        #----------------------------------------------------------------------#
        cat("Vesalius Object containing:\n\n")
        ntiles <- length(unique(object@tiles$barcodes))
        cat(ntiles, "tiles \n\n")
        #----------------------------------------------------------------------#
        # Showing assays
        #----------------------------------------------------------------------#
        if ((assay <- length(object@log@assay)) > 0) {
            n <- unlist(object@log@assay)
            cat(assay, "assays (", n, ")\n")
            #------------------------------------------------------------------#
            # Adding active embedding
            #------------------------------------------------------------------#
            n <- object@log@assay[[length(object@log@assay)]]
            cat(n, "as active assay\n\n")
        }
        #----------------------------------------------------------------------#
        # Showing norm
        #----------------------------------------------------------------------#
        if ((counts <- length(object@counts)) > 1) {
            n <- names(object@counts)[seq(2, counts)]
            cat(counts - 1, "normalization methods (", n, ")\n")
            #------------------------------------------------------------------#
            # Adding active embedding
            #------------------------------------------------------------------#
            n <- names(object@counts)[counts]
            cat(n, "as active normalized data\n\n")
        }
        #----------------------------------------------------------------------#
        # showing embeds
        #----------------------------------------------------------------------#
        if ((embeds <- length(object@embeddings)) > 0) {
            n <- names(object@embeddings)
            cat(embeds, "embeddings (", n, ")\n")
            #------------------------------------------------------------------#
            # Adding active embedding
            #------------------------------------------------------------------#
            n <- names(object@activeEmbeddings)
            cat(n, "as active Embedding\n\n")
        }
        #----------------------------------------------------------------------#
        #If there are territories
        #----------------------------------------------------------------------#
        if (!is.null(object@territories)) {
            ter <- object@territories[, ncol(object@territories)]
            cat(length(unique(ter)),
                "territories in",
                colnames(object@territories)[ncol(object@territories)],
                "\n\n")
        }
        #----------------------------------------------------------------------#
        #If there are DEG
        #To update based on how we call the base groups
        #----------------------------------------------------------------------#
        if (length(object@DEG) > 0) {
            deg <- object@DEG[[length(object@DEG)]]
            ters <- length(unique(c(deg$seedTerritory, deg$queryTerritory)))
            cat(nrow(deg),
            "Differentially Expressed Genes between",
            ters,
            "territories\n")
        }
        simple_bar(TRUE)
    }

)
