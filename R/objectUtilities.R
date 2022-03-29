
.slotApply <- function(x,FUN,...){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),...)
    }
    return(result)
}

.commitLog <- function(vesalius,commit,defaults){
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
        return(paste(unlist(x),"(Default)"))
    })

    logdf <- data.frame("Argument" = c("Function",names(defs)),
                        "Value" = c(as.character(commit[[1]]),sapply(defs,"[[",1)))
    #--------------------------------------------------------------------------#
    # Update the default df based on mathc call output
    #--------------------------------------------------------------------------#
    if(length(commit)>1){
      for(i in seq(2,length(commit))){
          logdf$Value[logdf$Argument == names(commit)[i]] <- commit[i]
      }
    }
    #--------------------------------------------------------------------------#
    # Check for prior logs
    #--------------------------------------------------------------------------#

    if(sum(unlist(.slotApply(vesalius@log,length)))>0){
      last <- max(sapply(.slotApply(vesalius@log,length),max)) + 1
    }else{
      last <- 1
    }
    logdf <- list(logdf)
    names(logdf) <- last
    vesalius@log <- .assignLogSlot(vesalius@log,as.character(commit[[1]]),logdf)
    return(vesalius)

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
          slot(object,i) <- c(slot(object,i),logdf)
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
