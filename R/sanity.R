################################################################################
###############################   Vesalius      ################################
################################################################################

#------------------------------/Sanity Checks /--------------------------------#

check_vesalius <- function(vesalius, init = FALSE) {
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format to isolate_territories function")
    }
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
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format to isolate_territories function")
    }
    tiles <- vesalius@tiles
    return(tiles)
}

check_territories <- function(vesalius, trial) {
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format to isolate_territories function")
    }
    territories <- vesalius@territories %||%
        stop("No territories have been computed!")
    if (trial == "last") {
      trial <- colnames(territories)[ncol(territories)]
    } else if (length(grep(x = colnames(vesalius@territories),
        pattern = trial)) == 0) {
      stop(paste(deparse(substitute(trial)), "is not in territory data frame"))
    }
    territories <- territories[, c("barcodes", "x", "y", trial)]
    colnames(territories) <- c("barcodes", "x", "y", "trial")
    return(territories)
}

check_segments <- function(vesalius, trial = "last") {
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format to isolate_territories function")
    }
    territories <- vesalius@territories %||%
        stop("No image segments have been computed yet!")
    if (!any(grepl(x = colnames(territories), attern = "Segment"))) {
        stop("No image segments have been computed yet!")
    }
    if (trial == "last") {
        trial <- grep("Segment", colnames(territories), value = TRUE) %>%
            tail(1)
    } else if (length(grep(trial, colnames(territories)) == 0)) {
        stop(
            paste(deparse(substitute(trial)), "is not in territory data frame")
        )
    }
    territories <- territories[, c("barcodes", "x", "y", trial)]
    colnames(territories) <- c("barcodes", "x", "y", "segment")
    return(territories)
}

check_norm <- function(vesalius, norm_method) {
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format to isolate_territories function")
    }
    counts <- vesalius@counts %||%
        stop("Cannot find any counts in vesalius object!")
    if (norm_method == "last") {
        norm_method <- length(counts)
    } else if (length(grep(norm_method, names(counts))) == 0) {
        stop(
            paste0(deparse(substitute(norm_method)), "is not in count list!")
        )
    }
    counts <- as.matrix(counts[[norm_method]])
    return(counts)
}

territory_dispatch <- function(territories, ter_1, ter_2, cells) {
    if (is.null(ter_1) && is.null(ter_2)) {
        territories <- select(territories, c("barcodes", "x", "y", "trial"))
    }else if (!is.null(ter_1) && is.null(ter_2)) {
        territories$trial[!territories$trial %in% ter_1] <- "other"
    }else if (is.null(ter_1) & !is.null(ter_2)) {
        territories$trial[!territories$trial %in% ter_2] <- "other"
    }else {
        territories$trial[!territories$trial %in% c(ter_1, ter_2)] <- "other"
    }
    if (!is.null(cells)) {
        territories$trial[territories$trial %in% ter_1 &&
            territories$barcodes %in% cells] <- paste0(ter_1, collapse = " ")
        territories$trial[territories$trial %in% ter_2 &&
            territories$barcodes %in% cells] <- paste0(ter_2, collapse = " ")
    }
    return(territories)
}