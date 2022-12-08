################################################################################
###############################   Vesalius      ################################
################################################################################

#------------------------------/Sanity Checks /--------------------------------#



#' check input checks the validity of input data to build_vesalius_assay
#' @param counts count matrix
#' @param coordinates coordinate file
#' @param assay string with assay name
#' @param adjust_coordinates string describing how coordinates should be
#' adjusted (origin or norm)
#' @param verbose logical if progress message should be outputed.
#' @return list containing checked counts, coordinates and assay
check_inputs <- function(counts,
    coordinates,
    assay,
    adjust_coordinates,
    verbose) {
    #--------------------------------------------------------------------------#
    # we can check the validity of the coordinates
    #--------------------------------------------------------------------------#
    coordinates <- check_coordinates(coordinates,
        assay,
        adjust_coordinates,
        verbose)
    #--------------------------------------------------------------------------#
    # next let's check counts if they are present
    # if they are we also check if the barcodes between the counts and the 
    # coordinates match
    # We will filter out any barcode that doesn't line up
    # essentially only use the intersection between barcodes 
    #--------------------------------------------------------------------------#
    if (!is.null(counts)) {
        counts <- check_counts(counts, assay, verbose)
        loc <- check_barcodes(colnames(counts), coordinates$barcodes)
        counts <- counts[, loc]
        coordinates <- coordinates[coordinates$barcodes %in% loc,]
        counts <- list(counts)
        names(counts) <- "raw"
        comment(counts) <- "raw"
    } else {
        counts <- list()
        comment(counts) <- "empty"
    }
    
    #--------------------------------------------------------------------------#
    # Finally we return the objects hat will populate the vesalius_assay
    #--------------------------------------------------------------------------#
    return(list("counts" = counts, "coordinates" = coordinates))
}

#' checking overlap between barcodes in counts and coordinates
#' @param counts character vector containing barcode names in count matrix
#' @param coordinates character vector containing barcode names in tile
#' data frame
#' @details Will throw in a warning if the overlap is not perfect.
#' @return overlapping location between counts and coordinates
check_barcodes <- function(counts, coordinates) {
    if (sum(duplicated(counts)) > 0) {
        stop("Duplicated colnames in count matrix!")
    }
    loc <- intersect(counts, coordinates)
    if (length(loc) == 0) {
        stop("Barcodes in count matrix and coordinates do no match")
    }
    if (length(loc) != length(counts)) {
        # Might want to remove this warning 
        # useful for me but could be annoying ofr the user
        warning("Unshared barcodes between counts and coordinates \n
            Using barcode intersection!")
    }
    return(loc)

}

#' check if counts are of the correct type and format
#' @param counts count matrix
#' @param assay string - assay name
#' @param verbose logical if progress message should be printed
#' @return count matrix or error
#' @importFrom methods is as
check_counts <- function(counts, assay, verbose) {
    message_switch("check_counts", verbose, assay = assay)
    if (is(counts, "data.frame")) {
      counts <- as(as.matrix(counts), "dgCMatrix")
    } else if (is(counts, "matrix")) {
      counts <- as(counts, "dgCMatrix")
    } else if (is(counts, "dgCMatrix")) {
      counts <- counts
    } else {
      stop(paste("Unsupported count format in", assay))
    }
    return(counts)
}

#' check if coordinates are of the correct type and format and
#' adjust coordinate value to remove white edge space
#' @param coordinates coordinate data 
#' @param assay string - assay name
#' @param verbose logical if progress message should be printed
#' @param adjust_coordinates string - how should coordinates be adjusted
#' @details adjusts coordinates by either snapping coordinates to origin 
#' or min max normalisation of coordinates. Might add polar for future tests.
#' @return coordinates data.frame or error
#' @importFrom methods is as
check_coordinates <- function(coordinates,
    assay,
    adjust_coordinates = c("origin", "norm"),
    verbose) {
    message_switch("check_coord",verbose, assay = assay)
    if (is(coordinates, "matrix")) {
        coordinates <- as.data.frame(coordinates)
    } else if (is(coordinates, "data.frame")) {
        coordinates <- coordinates
    } else {
       stop(paste("Unsupported coordinate format in", assay))
    }
    #--------------------------------------------------------------------------#
    # Check coordinate input type
    # for now let's put slide seq
    # THIS will need to be update if we use VISIUM and images!!!!
    #--------------------------------------------------------------------------#
    if (all(c("barcodes", "xcoord", "ycoord") %in% colnames(coordinates))) {
        coordinates <- coordinates[, c("barcodes", "xcoord", "ycoord")]
        colnames(coordinates) <- c("barcodes",  "x", "y")
    } else if (all(c("barcodes", "x", "y") %in% colnames(coordinates))) {
        coordinates <- coordinates[, c("barcodes", "x", "y")]
    } else {
        stop("Unknown column names")
    }
    #--------------------------------------------------------------------------#
    # ensuring that we have the right type in each column
    #--------------------------------------------------------------------------#
    coordinates$barcodes <- as.character(coordinates$barcodes)
    coordinates$x <- as.numeric(coordinates$x)
    coordinates$y <- as.numeric(coordinates$y)
    if (adjust_coordinates[1L] == "origin") {
        coordinates$x_orig <- coordinates$x
        coordinates$y_orig <- coordinates$y
        coordinates$x <- (coordinates$x - min(coordinates$x)) + 1
        coordinates$y <- (coordinates$y - min(coordinates$y)) + 1
    } else if (adjust_coordinates[1L] == "norm") {
        coordinates$x_orig <- coordinates$x
        coordinates$y_orig <- coordinates$y
        coordinates$x <- min_max((coordinates$x))
        coordinates$y <- min_max((coordinates$y))
    } else {
        stop("Woops - not sure how you want me to adjust coordinates")
    }
    return(coordinates)
}




#' check embeddings
#' checking user provided embedding matrix 
#' @param embeddings a matrix 
#' @details If a user provides their own embedding matrix, we need 
#' to check if it fits the required format and if not throw, errors,
#' warnings or adjust if possible. 
#' @return a formated matrix
check_embeddings <- function(embeddings) {
    if (is.null(rownames(embeddings))) {
        stop("No rownames provided to embedding matrix! \n")
    }
    barcodes <- rownames(embeddings)
    if (is(embeddings, "data.frame")) {
        embeddings <- as.matrix(embeddings)
        rownames(embeddings) <- barcodes
    }
    return(embeddings)
}

#' check if norm method is present in availble options
#' @param norm string - norm method parsed to function
#' @param use_counts string - which count matrix to use
#' @return norm method string or error
check_norm_methods <- function(norm, use_counts = NULL) {
    if (any(!norm %in% c("log_norm", "SCTransform", "TFIDF", "raw"))) {
        stop("Normalisation method provided do not match available options \n
            Select from: \n
            log_norm, SCTransform, TFIDF, or raw")
    }
    #-------------------------------------------------------------------------#
    # here we are assuming that if the user parses anything else but raw 
    # to use counts they want to use their own normalised matrix when creating
    # the embeddings. we can skip the normalisation method by setting the norm
    # method to raw. 
    #--------------------------------------------------------------------------#
    if (!is.null(use_counts) && use_counts != "raw") {
        norm <- "raw"
    }
    return(norm)
}

#' check if embedding method is present in availble options
#' @param embed embedding method selected by user
#' @return embedding method string or error
check_embed_methods <- function(embed) {
    if (any(!embed %in% c("PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP"))) {
        stop("Embedding method provided do not match available options \n
            Select from: \n
            PCA, PCA_L, UMAP, LSI, LSI_UMAP")
    } else {
        return(embed)
    }
}



#' check embedding selection
#' @param vesalius_assay a vesalius_assay
#' @param embed string embedding selection choice 
#' @param dims integer vector containing embedding dimension to extract 
#' @details we want to check if the embedding that the user requests 
#' is present in the assay. If not return error. If more than one 
#' with that name return last entry of that name and warning. Default is last 
#' that will just take the last embedding created which should be stored in
#' the active slot in the vesalius_Assay object.
#' We also want to be able to select which dimensions we want to return.
#' We make sure that those dimensions can be extracted from the embedding 
#' data frame.
#' @return embedding data frame
check_embedding_selection <- function(vesalius_assay, embed, dims) {
    #--------------------------------------------------------------------------#
    # Lets check if we have all the right embeddings
    #--------------------------------------------------------------------------#
    if (embed == "last") {
        embeddings <- vesalius_assay@active
    } else {
        in_col <- grep(pattern = embed, x = names(vesalius_assay@embeddings))
        if (length(in_col) == 0) {
            stop(paste(deparse(substitute(embed)),
                ": Unknown embedding selected!"))
        } else if (length(in_col) > 1) {
            trial <- grep(x = names(vesalius_assay@embeddings),
                pattern = embed,
                value = TRUE)
            warning(paste("More than one trial contains that name: \n",
                paste(trial, collapse = " ", sep = " "),
                "\nUsing first trial"))
            embeddings <- vesalius_assay@embeddings[[trial[1L]]]
        } else {
            embeddings <- vesalius_assay@embeddings[[in_col]]
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

#' check tile integretiy 
#' @param vesalius_assay a vesalius_assay 
#' @details place holder for in depth checks at the moment only return tiles 
#' @return tile data frame 
check_tiles <- function(vesalius_assay) {
    tiles <- vesalius_assay@tiles
    return(tiles)
}


#' check if territory selection is a valid option
#' @param vesalius_assay a vesalius_assay object
#' @param trial string - trial selection parse by user
#' @details check if the trial selection exists in territory slot
#' Default is last that will take the last entry. This function
#' will also reformat to only include the necessay information.
#' @return data frame contain selected trial
#' @importFrom infix %||%
check_territories <- function(vesalius_assay, trial) {
    territories <- vesalius_assay@territories %||%
        stop("No territories have been computed!")
    if (trial == "last") {
      trial <- colnames(territories)[ncol(territories)]
    } else if (length(grep(x = colnames(vesalius_assay@territories),
        pattern = trial)) == 0) {
      stop(paste(deparse(substitute(trial)), "is not in territory data frame"))
    } else if (length(grep(x = colnames(vesalius_assay@territories),
        pattern = trial)) > 1) {
        trial <- grep(x = colnames(vesalius_assay@territories),
            pattern = trial,
            value = TRUE)
        warning(paste("More than one trial contains that name: \n",
            paste(trial, collapse = " ", sep = " "),
            "\nUsing first trial"))
        trial <- trial[1L]
    }
    territories <- territories[, c("barcodes", "x", "y", trial)]
    colnames(territories) <- c("barcodes", "x", "y", "trial")
    return(territories)
}

#' check if segment selection is a valid option
#' @param vesalius_assay a vesalius_assay object
#' @param trial string - trial selection parse by user
#' @details check if the trial selection exists in territory slot
#' Default is last that will take the last entry. This function
#' will also reformat to only include the necessay information.
#' @return data frame contain selected trial
#' @importFrom infix %||%
#' @importFrom utils tail
check_segments <- function(vesalius_assay, trial = "last") {
    territories <- vesalius_assay@territories %||%
        stop("No image segments have been computed yet!")
    if (!any(grepl(x = colnames(territories), pattern = "Segment"))) {
        stop("No image segments have been computed yet!")
    }
    if (trial == "last") {
        trial <- grep("Segment", colnames(territories), value = TRUE) %>%
            tail(1)
    } else if (length(grep(trial, colnames(territories))) == 0) {
        stop(
            paste(deparse(substitute(trial)), "is not in territory data frame")
        )
    } else if (length(grep(x = colnames(vesalius_assay@territories),
        pattern = trial)) > 1) {
        trial <- grep(x = colnames(vesalius_assay@territories),
            pattern = trial,
            value = TRUE)
        warning(paste("More than one trial contains that name: \n",
            paste(trial, collapse = " ", sep = " "),
            "\nUsing first trial"))
        trial <- trial[1L]
    }
    territories <- territories[, c("barcodes", "x", "y", trial)]
    colnames(territories) <- c("barcodes", "x", "y", "segment")
    return(territories)
}

#' check if select norm method is an available option 
#' @param vesalius_assay a vesalius_assay
#' @param norm_method string - selected normalisation parsed by user
#' @param method DEG method 
#' @param verbose logical if progress message should be ouputed
#' @details Here we check if the normalisation method has been computed in
#' the count slot. We also check if this is used in the context of
#' differential gene expression analysis. DESeq and edgeR require raw counts
#' so if the user parses a valid norm method but is running DEG with one
#' those methods, we ignore request and parse raw count instead.
#' @return count matrix
#' @importFrom infix %||%
check_norm <- function(vesalius_assay,
    norm_method,
    method = NULL,
    verbose = TRUE) {
    counts <- vesalius_assay@counts
    if (length(counts) == 0) {
        stop("Cannot find any counts in vesalius assay!")
    }
    if (norm_method == "last") {
        norm_method <- get_active_count_tag(vesalius_assay)
    } else if (length(grep(norm_method, names(counts))) == 0) {
        stop(
            paste0(deparse(substitute(norm_method)), "is not in count list!")
        )
    } else if (length(grep(norm_method, names(counts))) > 1) {
        trial <- grep(norm_method, names(counts), value = TRUE)
        warning(paste("More than one trial contains that name: \n",
            paste(trial, collapse = " ", sep = " "),
            "Using first trial"))
        trial <- trial[1L]
    }
    if (!is.null(method[1L]) && method[1L] %in% c("DEseq2", "QLF", "LRT")) {
        message_switch("norm_check", verbose, method = method)
        if (norm_method == "empty") {
            stop("No raw counts provided! 
                Using user provided count matrix.
                Results will be sub-optimal!")
        }
        counts <- counts[["raw"]]
        counts <- as.matrix(counts)
    } else {
        counts <- as.matrix(counts[[norm_method]])
    }
    return(counts)
}


#' check if provided cells are contained withing provided territory 
#' @param territory_barcodes character vector containing spatial barcodes
#' of territories
#' @param ter string selected territories in the form of chacater string
#' @param cell_barcodes character vector containing barcodes of cells
#' of interest
#' @param verbose logical if progress message should be outputed
#' @return charcater vector of common barcodes between territory 
#' barcodes and cell barcodes.
check_cells <- function(territory_barcodes, ter, cell_barcodes, verbose) {
    common <- intersect(territory_barcodes, cell_barcodes)
    if (length(common) == 0) {
        message_switch("no_cell", verbose, ter = ter)
        common <- NULL
    }
    return(common)
}

#' check number of spatial indices present
#' @param group count matrix 
#' @param min_spatial_index numeric - min number of spatial indices 
#' that need to be present in group
#' @param id string - id of group 
check_min_spatial_index <- function(group, min_spatial_index, id) {
    if (is.null(dim(group)) || dim(group)[2L] < min_spatial_index) {
        stop(paste0("Territory ", id, " does not contain enough cells.\n
        Territory will be skipped"), call. = FALSE)
        return(NULL)
    } else {
        return(group)
    }
}

#' check if requested territory is in territory list
#' @param territories data frame containing territories
#' @param group vector indicating which territories should be selected
#' from territory data frame.
#' @return a vector of territories that are present in territory data frame
check_group_value <- function(territories, group) {
    present <- unique(territories$trial)
    group_sub <- group[group %in% present]
    if (length(group_sub) == 0) {
        stop(paste("Territory", group, "is not present in the selected trial!"))
    } else if (length(group_sub) < length(group)) {
        warning(paste("Only group territory(ies)",
            paste(group_sub, collaspe = " "),
            "is (are) present in the selected trial. Only those will be used"))
        return(group_sub)
    } else {
        return(group_sub)
    }
}
