################################################################################
###############################   Vesalius      ################################
################################################################################

#------------------------------/Sanity Checks /--------------------------------#

check_vesalius <- function(vesalius, init = FALSE) {
    #--------------------------------------------------------------------------#
    # Essentially we want to check what there is in this object
    # to avoid any unnecessary computations
    # Might need to be extended
    #--------------------------------------------------------------------------#
    if (init) {
      #------------------------------------------------------------------------#
      # Check to see if there is only a single log entry
      # if so then that means we have a fresh object
      # if not that means we can skip some of the processing
      #------------------------------------------------------------------------#
      log <- get_last_log(vesalius)
      if (length(log) == 1 && names(log)[1L] == "assay") {
          return(TRUE)
      }else {
          return(FALSE)
      }
    } else {
        return(FALSE)
    }


}

check_embedding <- function(vesalius, embed, dims) {
    #--------------------------------------------------------------------------#
    # Lets check if we have all the right embeddings
    #--------------------------------------------------------------------------#
    embeddings <-
    if (embed == "last") {
      embeddings <- vesalius@activeEmbeddings[[1L]]
    } else {
        in_col <- grep(pattern = embed, x = names(vesalius@embeddings))
        if (length(in_col) == 0) {
            stop(paste(deparse(substitute(embed)),
                ": Unknown embedding selected!"))
        } else if (length(in_col) > 1) {
            warning(paste("More than 1", deparse(substitute(embed)), "embedding.
            Vesalius will use the latest entry - See Vesalius Log"))
            embeddings <- vesalius@embeddings[[in_col[length(in_col)]]]
        } else {
            embeddings <- vesalius@embeddings[[in_col]]
        }
    }
    #--------------------------------------------------------------------------#
    # Let's check if we have the right number of dims
    #--------------------------------------------------------------------------#
    if (length(dims) > ncol(embeddings)) {
      stop(paste0("To many dimesnions supplied! Only",
        ncol(embeddings),
        " present"))
    }
    return(embeddings)
}

check_tiles <- function(vesalius) {
    tiles <- vesalius@tiles
    return(tiles)
}

check_territories <- function(vesalius, trial) {
    territories <- vesalius@territories %||%
        stop("No territories have been computed!")
    if (trial == "last") {
      trial <- colnames(territories)[ncol(territories)]
    } else if (length(grep(x = colnames(vesalius@territories),
        pattern = trial)) == 0) {
      stop(paste(deparse(substitute(trial)), "is not in territory data frame"))
    }
    territories <- territories[, c("x", "y", trial)]
    colnames(territories) <- c("x", "y", "territory")
    territories$territory <- as.factor(territories$territory)
    return(territories)
}

check_segments <- function(vesalius) {

}
