################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Vesalius Objects/--------------------------------#



setClassUnion("mat",c("matrix","dgCMatrix"))
setClassUnion("ter",c("data.frame","NULL"))


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
                    DEG = "data.frame",
                    counts  = "mat",
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
        if(is(object@DEG)[1L] != "data.frame"){
            stop("Unsupported DEG Format")
        }
        if(is(object@counts)[1L]!= "matrix" & is(object@counts)[1L]!= "dgCMatrix"){
            stop("Unsupported Count Format")
        }
        if(is(object@log)[1L]!= "list"){
            stop("Unsupported log Format")
        }
        return(TRUE)
    }

)

#------------------------------------------------------------------------------#
# Create new object function
#------------------------------------------------------------------------------#

.vesaliusObject <- function(tiles = NULL,
                           embeddings = NULL,
                           activeEmbeddings = NULL,
                           territories = NULL,
                           DEG = NULL,
                           counts = NULL
                           log = NULL){

    ves <- new("vesaliusObject",
               tiles = tiles,
               embeddings = embeddings,
               activeEmbeddings = activeEmbeddings,
               DEG = DEG,
               territories = territories,
               counts = counts,
               log = log)
    return(ves)
}

buildVesaliusObject <- function(coordinates,counts){
    #--------------------------------------------------------------------------#
    # First we will create tiles
    #--------------------------------------------------------------------------#
    coordinates <- .checkCoordinates(coordinates)
    counts <- .checkCounts(counts)
    ves <- .vesaliusObject(tiles = coordinates,
                          counts = counts)
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
      cat("Vesalius Object containing:\n")
      ntiles <- length(unique(object@tiles$barcodes))
      cat(ntiles,"tiles \n")
      #------------------------------------------------------------------------#
      # This will change if we add multiple embedding type
      #------------------------------------------------------------------------#
      embeds <- length(object@embeddings)
      n <- names(object@embeddings)
      cat(embeds, "embeddings (",n,")\n")
      #------------------------------------------------------------------------#
      # Adding active embedding
      #------------------------------------------------------------------------#
      n <- names(object@activeEmbeddings)
      cat(n ,"as current active Embedding")
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
          ter <-length(unique(object@territories$territory))
          cat(ter, "territories\n")
      }
      #------------------------------------------------------------------------#
      #If there are DEG
      #To update based on how we call the base groups
      #------------------------------------------------------------------------#
      if(!is.null(object@DEG)){
          deg <- nrow(object@DEG)
          ters <- length(unique(c(object@DEG$seedTerritory,object@DEG$queryTerritory))
          cat(deg, "Differentially Expressed Genes between", ters,"territories\n")
      }


      .simpleBar(TRUE)
    }

)

#------------------------------------------------------------------------------#
# Vesalius Log object
# Keep track of what has been done
#------------------------------------------------------------------------------#

setClass("log",
         slots=list(tiles="list",
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
            leaves <- sum(.slotApply(object, length))
            cat(leaves,"commits to vesalius Log tree")
            #------------------------------------------------------------------#
            # Last update
            #------------------------------------------------------------------#
            last <- .slotApply(object,names)
            lcom <- which.max(sapply(last, function(x){
                    return(max(as.numeric(x)))
            }))
            cat("Last Log commit on",slotNames(object)[lcom],"branche")
            #------------------------------------------------------------------#
            # Show last commit param
            #------------------------------------------------------------------#
            lastCommit <- .slotApply(object,function(x){
                  return(x[[length(x)]])
            })
            for(i in seq_along(lastCommit)){
                cat(slotNames(object)[i],":\n",lastCommit[[i]],"\n")
            }
            .simpleBar(TRUE)

          })
