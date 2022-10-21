################################################################################
###############################   Vesalius      ################################
################################################################################

#------------------------------/Sanity Checks /--------------------------------#




check_inputs <- function(counts,
    coordinates,
    assay,
    adjust_coordinates,
    verbose) {
    #--------------------------------------------------------------------------#
    # next let's check counts
    # create a list that we name and comment.
    # the comment tells us which count matrix is the default one 
    #--------------------------------------------------------------------------#
    counts <- list(check_counts(counts, assay, verbose))
    names(counts) <- "raw"
    comment(counts) <- "raw"
    #--------------------------------------------------------------------------#
    # now we can check the validity of the coordinates
    #--------------------------------------------------------------------------#
    coordinates <- check_coordinates(coordinates,
        assay,
        adjust_coordinates,
        verbose)
    #--------------------------------------------------------------------------#
    # Finally we can rebuild an assay list filled with vesalius_assay objects 
    #--------------------------------------------------------------------------#
    return(list("counts" = counts,
        "coordinates" = coordinates,
        "assay" = assay))
}



check_counts <- function(counts, assay, verbose) {
    message_switch("check_counts",verbose, assay = assay)
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

check_norm_methods <- function(norm) {
    if (any(!norm %in% c("log_norm", "SCTransform", "TFIDF", "raw"))) {
        stop("Normalisation method provided do not match available options \n
            Select from: \n
            log_norm, SCTransform, TFIDF, or raw")
    } else {
        return(norm)
    }
}
check_embed_methods <- function(embed, n_counts) {
    if (any(!embed %in% c("PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP"))) {
        stop("Embedding method provided do not match available options \n
            Select from: \n
            PCA, PCA_L, UMAP, LSI, LSI_UMAP")
    } else {
        return(embed)
    }
}




check_embedding <- function(vesalius_assay, embed, dims) {
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
            warning(paste("More than 1", deparse(substitute(embed)), "embedding.
            Vesalius will use the latest entry - See Vesalius Log"))
            embeddings <- vesalius_assay@embeddings[[in_col[length(in_col)]]]
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

check_tiles <- function(vesalius_asssay) {
    tiles <- vesalius_asssay@tiles
    return(tiles)
}


check_territories <- function(vesalius_assay, trial) {
    territories <- vesalius_assay@territories %||%
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

check_segments <- function(vesalius_asssay, trial = "last") {
    territories <- vesalius_asssay@territories %||%
        stop("No image segments have been computed yet!")
    if (!any(grepl(x = colnames(territories), pattern = "Segment"))) {
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

check_norm <- function(vesalius_asssay, norm_method, method = NULL, verbose = TRUE) {
    counts <- vesalius_asssay@counts %||%
        stop("Cannot find any counts in vesalius assay!")
    if (norm_method == "last") {
        norm_method <- length(counts)
    } else if (length(grep(norm_method, names(counts))) == 0) {
        stop(
            paste0(deparse(substitute(norm_method)), "is not in count list!")
        )
    }
    if (!is.null(method) && method %in% c("DEseq2", "edgeR")) {
        message_switch("norm_check", verbose, method = method)
        counts <- as.matrix(counts[["raw"]])
    } else {
        counts <- as.matrix(counts[[norm_method]])
    }
    return(counts)
}



check_cells <- function(territory_barcodes, ter, cell_barcodes, verbose) {
    common <- intersect(territory_barcodes, cell_barcodes)
    if (length(common) == 0) {
        message_switch("no_cell", verbose, ter = ter)
        common <- NULL
    }
    return(common)
}

check_min_spatial_index <- function(group, min_spatial_index, id) {
    if (is.null(dim(group)) || dim(group)[2L] < min_spatial_index) {
        warning(paste0("Territory ", id, " does not contain enough cells.\n
        Territory will be skipped"), call. = FALSE)
        return(NULL)
    } else {
        return(group)
    }
}
