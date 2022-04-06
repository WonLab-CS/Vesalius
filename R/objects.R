################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#

#------------------------------------------------------------------------------#
# Vesalius Log object
# Keep track of what has been done
#------------------------------------------------------------------------------#

setClass("log",
         slots=list(assay = "list",
                    tiles="list",
                    embeddings = "list",
                    activeEmbeddings = "list",
                    territories = "list",
                    DEG = "list",
                    counts  = "list"))

setMethod("show",
          signature = "log",
          definition = function(object){
            .simpleBar(TRUE)
            #------------------------------------------------------------------#
            # General log information
            #------------------------------------------------------------------#

            if(sum(unlist(.slotApply(object, length))) == 0){
               cat("No commits to log tree\n")
            } else {

              #----------------------------------------------------------------#
              # Last update
              #----------------------------------------------------------------#
              cat("Log tree\n")
              last <- .slotApply(object,length)
              for(i in seq_along(last)){
                  if(last[[i]] ==0){
                      shape <- "/"
                  } else {
                      shape <- rep("*", times = last[[i]])
                  }
                  cat("|",shape,names(last)[i],"\n")
              }


            }

            .simpleBar(TRUE)

          })

#------------------------------------------------------------------------------#
# class unions
#------------------------------------------------------------------------------#

setClassUnion("ter",c("data.frame","NULL"))
setClassUnion("DEG",c("data.frame","NULL"))



### Need to clean these up
### Increase the robustness of this
#------------------------------------------------------------------------------#
# Main vesalius object
#------------------------------------------------------------------------------#
setClass("vesaliusObject",
         slots=list(tiles="data.frame",
                    embeddings = "list",
                    activeEmbeddings = "list",
                    territories = "ter",
                    DEG = "DEG",
                    counts  = "list",
                    log = "log"),
         prototype=c(tiles = data.frame(),
                     embeddings = list()),
        validity= function(object){
        if(is(object@tiles)[1L] != "data.frame"){
            stop("Unsupported Tile Format")
        }
        if(is(object@embeddings)[1L]!= "list"){
            stop("Unsupported Embeddings Format")
        }
        if(is(object@activeEmbeddings)[1L]!= "list"){
            stop("Unsupported Active Embeddings Format")
        }
        if(is(object@territories)[1L]!= "data.frame" & is(object@territories)[1L]!= "NULL"){
            stop("Unsupported territories Format")
        }
        if(is(object@DEG)[1L]!= "data.frame" & is(object@DEG)[1L]!= "NULL"){
            stop("Unsupported DEG Format")
        }
        if(is(object@counts)[1L]!= "list"){
            stop("Unsupported Count Format")
        }
        if(is(object@log)[1L]!= "log"){
            stop("Unsupported log Format")
        }
        return(TRUE)
    }

)

#------------------------------------------------------------------------------#
# Create new object function
# This one is on hold atm
#------------------------------------------------------------------------------#


buildVesaliusObject <- function(coordinates,counts,assay = "ST"){
    #--------------------------------------------------------------------------#
    # First we will create tiles
    # Definitely need some check here
    #--------------------------------------------------------------------------#
    coordinates <- .checkCoordinates(coordinates)
    counts <- list(.checkCounts(counts))
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


.checkCounts <- function(counts){

    if(class(counts) == "data.frame"){
      counts <- as(as.matrix(counts),"dgCMatrix")
    } else if(class(counts) == "matrix"){
      counts <- as(counts,"dgCMatrix")
    } else if(class(counts) == "dgCMatrix"){
      counts <- counts
    } else {
      stop("Unsupported count format!")
    }
    return(counts)
}


.checkCoordinates <- function(coordinates){
    #--------------------------------------------------------------------------#
    # Check coordinate input type
    # for now let's put slide seq
    # we can add some more loading functions later and santise from there
    #--------------------------------------------------------------------------#
    if(all(c("barcodes","xcoord","ycoord") %in% colnames(coordinates))){
        coordinates <- coordinates[,c("barcodes","xcoord","ycoord")]
        colnames(coordinates) <- c("barcodes","x","y")
    } else if(all(c("barcodes","x","y") %in% colnames(coordinates))) {
        coordinates <- coordinates[,c("barcodes","x","y")]
    } else {
        stop("Unknown column names")
    }
    return(coordinates)
}


#------------------------------------------------------------------------------#
# vesalius objecy show method
#------------------------------------------------------------------------------#

setMethod("show",
    signature = "vesaliusObject",
    definition = function(object){
      .simpleBar(TRUE)
      #------------------------------------------------------------------------#
      ### To be modified later
      #------------------------------------------------------------------------#
      cat("Vesalius Object containing:\n\n")
      ntiles <- length(unique(object@tiles$barcodes))
      cat(ntiles,"tiles \n\n")
      #------------------------------------------------------------------------#
      # Showing norm
      #------------------------------------------------------------------------#
      if((counts <- length(object@counts)) > 1){
        n <- names(object@counts)[seq(2,counts)]
        cat(counts-1, "normalization methods (",n,")\n")
        #----------------------------------------------------------------------#
        # Adding active embedding
        #----------------------------------------------------------------------#

        n <- names(object@counts)[counts]
        cat(n ,"as active normalized data\n\n")
      }
      #------------------------------------------------------------------------#
      # showing embeds
      #------------------------------------------------------------------------#
      if((embeds <- length(object@embeddings)) > 0){
        n <- names(object@embeddings)
        cat(embeds, "embeddings (",n,")\n")
        #----------------------------------------------------------------------#
        # Adding active embedding
        #----------------------------------------------------------------------#

        n <- names(object@activeEmbeddings)
        cat(n ,"as active Embedding\n\n")
      }

      #------------------------------------------------------------------------#
      # There will always a log for Vesalius objects unless you create a dummy
      # object
      #------------------------------------------------------------------------#
      #branches <- length(object@log)
      #leaves <- sum(sapply(object@log, length))
      #cat(leaves, "leaves between", branches,"branches in log tree\n")
      #------------------------------------------------------------------------#
      #If there are territories
      #------------------------------------------------------------------------#
      if(!is.null(object@territories)){
          ter <-object@territories[,ncol(object@territories)]
          cat(length(unique(ter)),
          "territories in",colnames(object@territories)[ncol(object@territories)],"\n")
      }
      #------------------------------------------------------------------------#
      #If there are DEG
      #To update based on how we call the base groups
      #------------------------------------------------------------------------#
      if(!is.null(object@DEG)){
          deg <- nrow(object@DEG)
          ters <- length(unique(c(object@DEG$seedTerritory,object@DEG$queryTerritory)))
          cat(deg, "Differentially Expressed Genes between", ters,"territories\n")
      }


      .simpleBar(TRUE)
    }

)
