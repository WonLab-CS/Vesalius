
.slotApply <- function(x,FUN,...){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),...)
    }
    return(result)
}

.commitLog <- function(log,commit,defaults){
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
    logdf <- data.frame("Argument" = c("Function",names(defaults)),
                        "Value" = c(commit[[1]],paste(unlist(defaults),"(Default)"))
    #--------------------------------------------------------------------------#
    # Update the default df based on mathc call output
    #--------------------------------------------------------------------------#
    if(length(log)>1){
      for(i in seq(2,length(log))){
          logdf$Value[logdf$Argument == names(commit)[i]] <- commit[i]
      }
    }
    #--------------------------------------------------------------------------#
    # Check for prior logs
    #--------------------------------------------------------------------------#
    if(!is.null(log@log)){
      last <- max(sapply(.slotApply(object,length),max)) + 1
    }else{
      last <- 1
    }
    logdf <- list(logdf)
    names(logdf) <- last
    log <- .assignLogSlot(log,commit[[1]],logdf)
    return(log)

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
      if(length(slot(object,sl))!=0){
          slot(object,sl) <- c(slot(object,sl),logdf)
      } else {
          slot(object,sl) <- input
      }
    }
    return(object)
}
