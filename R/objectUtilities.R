################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Object Utilities/--------------------------------#


.checkVesalius <- function(vesalius,init=FALSE){
    #--------------------------------------------------------------------------#
    # Essentially we want to check what there is in this object
    # to avoid any unnecessary computations
    # Might need to be extended
    #--------------------------------------------------------------------------#
    if(init){
      #------------------------------------------------------------------------#
      # Check to see if there is only a single log entry
      # if so then that means we have a fresh object
      # if not that means we can skip some of the processing
      #------------------------------------------------------------------------#
      log <- getLastLog(vesalius)
      if(length(log) == 1 & names(log)[1L] == "assay"){
          return(TRUE)
      }else{
          return(FALSE)
      }
    } else {
        return(FALSE)
    }


}
.updateVesalius <- function(vesalius,data,slot,commit,defaults,append = TRUE){
    #--------------------------------------------------------------------------#
    # You take a vesalius object, the data you want add to it and which slot
    # it should be added to.
    # First step is to check if the slot we want to update is empty
    # This is not the same as log commit
    # Lets keep those seperate
    #--------------------------------------------------------------------------#

    if(append){
            slot(vesalius,slot) <- c(slot(vesalius,slot),data)
            vesalius <- .commitLog(vesalius,commit,defaults,slot,append)
    } else {
            slot(vesalius,slot) <- data
            vesalius <- .commitLog(vesalius,commit,defaults,slot,append)
    }
    return(vesalius)
}



.commitLog <- function(vesalius,commit,defaults,slot, append = TRUE){
    #--------------------------------------------------------------------------#
    # First let's get all the non symbol arguments
    # symbol argument are argument that require an input i.e no default
    # Trim the last entry as it is always NULL... or is it????
    #--------------------------------------------------------------------------#

    defaults <- defaults[-length(defaults)]
    NonSym <- sapply(defaults,function(x){
                  return(!any(is(x) == "refObject"))
    })
    defaults <- defaults[NonSym]
    #--------------------------------------------------------------------------#
    # Initialise log
    #--------------------------------------------------------------------------#
    defs <- sapply(defaults, function(x){
        if(class(x)=="call"){
            x <- as.character(x)[2L]
        }
        return(unlist(x))
    })

    logdf <- data.frame("Argument" = c("Function",names(defs)),
                        "Value" = c(as.character(commit[[1]]),unname(defs)))

    #--------------------------------------------------------------------------#
    # Update the default df based on mathc call output
    #--------------------------------------------------------------------------#
    if(length(commit)>1){
      for(i in seq(2,length(commit))){

          logdf$Value[logdf$Argument == names(commit)[i]] <- as.character(commit[[i]])
      }
    }
    #--------------------------------------------------------------------------#
    # Check for prior logs
    #--------------------------------------------------------------------------#

    last <- .getLastCommit(vesalius)

    logdf <- list(logdf)
    names(logdf) <- last

    if(append){
        slot(vesalius@log,slot) <- c(slot(vesalius@log,slot),logdf)
    } else {
        slot(vesalius@log,slot) <- logdf
    }


    return(vesalius)

}

.getLastCommit <- function(log){
    log <-.slotApply(log@log,names)
    log <- sapply(lapply(log,as.numeric),max)
    log <- max(log) + 1
    return(log)
}


.slotApply <- function(x,FUN,...){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),...)
    }
    return(result)
}



.assignLogSlot <- function(log,commit,logdf){
    log <- switch(commit,
                  "buildVesaliusEmbeddings" =.slotAssign(log,
                    c("tiles","embeddings","activeEmbeddings"),
                    logdf))
    return(log)
}

.slotAssign <- function(object,sl,input){

    for(i in sl){
      if(length(slot(object,i))!=0){
          slot(object,i) <- c(slot(object,i),input)
      } else {
          slot(object,i) <- input
      }
    }
    return(object)
}

getLastLog <- function(vesalius){
  cl <- class(vesalius@log)
  log <- list()
  for(i in slotNames(cl)){
      tmp <- slot(vesalius@log,i)
      if(length(tmp)>0){
        log[[i]] <- tmp[[length(tmp)]]
      } else {
        log[[i]] <- NULL
      }

  }
  log <- log[!is.null(log)]
  return(log)
}

getCounts <- function(vesalius, type = "raw"){
    counts <- vesalius@counts[[type]]
    if(is.null(counts)){
        stop("No count matrix has been added to Vesalius")
    }
    return(counts)

}

viewLogTree <- function(vesalius){
    return(vesalius@log)
}
