################################################################################
###############################   Vesalius      ################################
################################################################################

#------------------------------/Sanity Checks /--------------------------------#



#' check input checks the validity of input data to build_vesalius_assay
#' @param counts count matrix
#' @param coordinates coordinate file
#' @param image connection to image or image array
#' @param assay string with assay name
#' @param adjust_coordinates string describing how coordinates should be
#' adjusted (origin or norm)
#' @param verbose logical if progress message should be outputed.
#' @return list containing checked counts, coordinates and assay
check_inputs <- function(counts,
    coordinates,
    image,
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
        coordinates <- coordinates[coordinates$barcodes %in% loc, ]
        counts <- list(counts)
        names(counts) <- "raw"
        comment(counts) <- "raw"
    } else {
        counts <- list()
        comment(counts) <- "empty"
    }
    #--------------------------------------------------------------------------#
    # checking image 
    #--------------------------------------------------------------------------#
    image <- check_image_input(image)
    #--------------------------------------------------------------------------#
    # Finally we return the objects hat will populate the vesalius_assay
    #--------------------------------------------------------------------------#
    return(list("counts" = counts,
        "coordinates" = coordinates,
        "image" = image))
}

#' check image input 
#' check if input image is in correct format
#' @param image input
#' @return list of images
#' @importFrom imager load.image
check_image_input <- function(image) {
    if (!is.null(image) && is(image, "character")) {
        image <- list(imager::load.image(image))
    } else if (!is.null(image) && is(image, "array")) {
        dims <- length(dim(image))
        if (dims == 0 || dims > 4) {
            stop("Image array dimesnions should be between 2 and 4 \n
            2 = matrix | 3 = array | 4 = cimg array like")
        }
        # Might not need to to do that if Im using a different package
        image <- list(suppressWarnings(as.cimg(image)))
    } else if (!is.null(image) && is(image, "list")) {
        image <- lapply(image, check_image)
    } else if (!is.null(image) && is(image, "cimg")) {
        image <- list(image)
    } else {
        image <- list()
    }
    return(image)
}

#' checking overlap between barcodes in counts and coordinates
#' @param mat_barcodes character vector containing barcode names in matrix
#' (count matrix or embedding matrix)
#' @param coordinates character vector containing barcode names in tile
#' data frame
#' @param throw logical if warning should be thrown or not
#' @details Will throw in a warning if the overlap is not perfect.
#' @return overlapping location between counts and coordinates
check_barcodes <- function(mat_barcodes, coordinates, throw = TRUE) {
    if (sum(duplicated(mat_barcodes)) > 0) {
        stop("Duplicated colnames in matrix!")
    }
    coordinates <- unlist(strsplit(coordinates, "_et_"))
    loc <- intersect(mat_barcodes, coordinates)
    if (length(loc) == 0) {
        stop("Barcodes in matrix and coordinates do no match!")
    }
    if ((length(loc) != length(mat_barcodes) ||
        length(loc) != length(coordinates)) && throw) {
        # Might want to remove this warning
        # useful for me but could be annoying ofr the user
        warning("Unshared barcodes between matrix and coordinates \n
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

#' check if image is present
#' @param vesalius_assay vesalius assay object
#' @return microscopy image
check_image <- function(vesalius_assay) {
    image <- vesalius@image
    if (length(image) == 0) {
        stop("No image have been added to the vesalius_assay")
    }
    image <- image[[1L]][, , , 1:3]
    return(image)
}

#' check if norm method is present in availble options
#' @param norm string - norm method parsed to function
#' @param use_counts string - which count matrix to use
#' @return norm method string or error
check_norm_methods <- function(norm, use_counts = NULL) {
    if (any(!norm %in% c("log_norm", "SCTransform", "TFIDF", "raw"))) {
        stop("Normalisation method provided does not match available options \n
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
        stop("Embedding method provided does not match available options \n
            Select from: \n
            PCA, PCA_L, UMAP, LSI, LSI_UMAP")
    } else {
        return(embed)
    }
}

#' check if smoothing method is present in availble options
#' @param method smoothing method selected by user
#' @return smoothing method string or error
check_smoothing_kernel <- function(method) {
    if (any(!method %in% c("iso", "box", "median"))) {
        stop("Smoothing method provided does not match available options \n
            Select from: \n
            iso, box, median or a combination of these 3 options")
    } else {
        return(method)
    }
}

#' check if across level options are valid for smoothing
#' @param method string - across level option
#' @return acroos level option or error
check_smoothing_level <- function(method) {
    if (any(!method %in% c("min", "max", "mean"))) {
        stop("Smoothing across levels method provided does not match available options \n
            Select from: \n
            min, max, mean")
    } else {
        return(method)
    }
}

#' check if eq method are avilable 
#' @param method string - eq method
#' @return eq method
check_eq_method <- function(method) {
    if (any(!method %in% c("EqualizePiecewise",
        "BalanceSimplest",
        "SPE",
        "EqualizeDP",
        "EqualizeADP",
        "ECDF"))) {
        stop("Equalization method provided does not match available options \n
            Select from: \n
            EqualizePiecewise, BalanceSimplest, SPE, 
            EqualizeDP, EqualizeADP, or ECDF")
    } else {
        return(method)
    }
}

#' check if segmentation options are valid for segmentation
#' @param method string - segmentation method
#' @return segmentation method
check_segmentation_method <- function(method) {
    if (any(!method %in% c("kmeans", "louvain", "leiden"))) {
        stop("Segmentation method provided does not match available options \n
            Select from: \n
            kmeans, louvain, leiden")
    } else {
        return(method)
    }
}

#' check territory isolation method is available 
#' @param method string - territory isolation method
#' @return territory isolation method
check_isolation_method <- function(method) {
    if (any(!method %in% c("distance"))) {
        stop("Isolation method provided does not match available options \n
            Select from: \n
            distance")
    } else {
        return(method)
    }
}

#' check dimension selection method 
#' @param method string - dimensions selecion 
#' @return territory isolation method
check_dim_selection_method <- function(method) {
    if (any(!method %in% c("batty"))) {
        stop("Dimension selection method provided does not match available options \n
            Select from: \n
            batty")
    } else {
        return(method)
    }
}

#' check if DEG method is one of the available options
#' @param method string - DEG method
#' @return DEG method or error
check_deg_method <- function(method) {
    if (any(!method %in% c("wilcox",
        "t.test",
        "chisq",
        "fisher.exact",
        "DESeq2",
        "QLF",
        "LRT",
        "logit"))) {
        stop("Territory identity method provided does not match available options\n
            Select from: \n
            wilcox, t.test, chisq, fisher.exact, DESeq2, QLF, LRT, or logit")
    } else {
        return(method)
    }
}



#' check tensor compression 
#' @param locs rle output of coordinate compression
#' @details Here we just want to check if the output is reasonable or not
#' Essentially if you end up with very little barcodes, this could mean
#' that the user compressed to much.
#' Set warning if compression is set to less 10% of total points
#' @return rle out of coordinate compression
check_tensor_compression <- function(locs) {
    #-------------------------------------------------------------------------#
    # Essentially throw error if you have a single coordinate left
    # keeping it with the < 2 just in case we decide that we want to be 
    # more restrictive and only allow a in of 10 coordinates for example 
    #-------------------------------------------------------------------------# 
    if (length(locs$values) < 2) {
        stop("Tensor compression yielded a single coordinate!
            Consider increasing the tensor_resolution.")
    }
    if (length(locs$values) < 0.1 * sum(locs$length) && 
        !length(locs$values) < 2) {
        warning("Tensor resolution has been reduced as to only retain
        less than 10% of original spatial coordinates")
    }
    
    return(locs)
}

#' check features
#' @param counts a seurat object
#' @details If the user use custom counts, we want to check
#' if variable features have been computed or not. Dim reduction
#' requires a list if features toi be provided or variable features
#' to be computed.
#' @return a list of features
#' @importFrom Seurat VariableFeatures
check_features <- function(counts) {
    if (length(Seurat::VariableFeatures(counts)) == 0) {
        features <- rownames(GetAssayData(counts, slot = "counts"))
    } else {
        features <- Seurat::VariableFeatures(counts)
    }
    return(features)
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
    if (sum(dim(vesalius_assay@active)) == 0) {
            stop("No embeddings have been computed!")
    }
    if (embed == "last") {
        embeddings <- vesalius_assay@active
    } else {
        in_col <- grep(pattern = paste0("^", embed, "$"),
            x = names(vesalius_assay@embeddings))
        if (length(in_col) == 0) {
            stop(paste(deparse(substitute(embed)),
                ": Unknown embedding selected!"))
        } else if (length(in_col) > 1) {
            trial <- grep(x = names(vesalius_assay@embeddings),
                pattern = paste0("^", embed, "$"),
                value = TRUE)
            warning(paste("More than one trial contains that name: \n",
                paste(trial, collapse = " ", sep = " "),
                "\nUsing first trial"))
            embeddings <- vesalius_assay@embeddings[[trial[1L]]]
        } else {
            trial <- grep(x = names(vesalius_assay@embeddings),
                pattern = paste0("^", embed, "$"),
                value = TRUE)
            embeddings <- vesalius_assay@embeddings[[trial]]
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
check_territory_trial <- function(vesalius_assay, trial) {
    if (sum(dim(vesalius_assay@territories)) == 0) {
        stop("No territories have been computed!")
    } else {
        territories <- vesalius_assay@territories
    }
    if (trial == "last") {
      trial <- colnames(territories)[ncol(territories)]
    } else if (length(grep(x = colnames(vesalius_assay@territories),
        pattern = trial)) == 0) {
      stop(paste(deparse(substitute(trial)), "is not in territory data frame"))
    } else if (length(grep(x = colnames(vesalius_assay@territories),
        pattern = paste0("^", trial, "$"))) > 1) {
        trial <- grep(x = colnames(vesalius_assay@territories),
            pattern = paste0("^", trial, "$"),
            value = TRUE)
        warning(paste("More than one trial contains that name: \n",
            paste(trial, collapse = " ", sep = " "),
            "\nUsing first trial"))
        trial <- trial[1L]
    } else {
        trial <- grep(x = colnames(vesalius_assay@territories),
            pattern = paste0("^", trial, "$"),
            value = TRUE)
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
#' @importFrom utils tail
check_segment_trial <- function(vesalius_assay, trial = "last") {
    if (sum(dim(vesalius_assay@territories)) == 0) {
        stop("No segments have been computed yet!")
    } else {
        territories <- vesalius_assay@territories
    }
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
        pattern = paste0("^", trial, "$"))) > 1) {
        trial <- grep(x = colnames(vesalius_assay@territories),
            pattern = paste0("^", trial, "$"),
            value = TRUE)
        warning(paste("More than one trial contains that name: \n",
            paste(trial, collapse = " ", sep = " "),
            "\nUsing first trial"))
        trial <- trial[1L]
    } else {
        trial <- grep(x = colnames(vesalius_assay@territories),
            pattern = paste0("^", trial, "$"),
            value = TRUE)
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
        trial <- get_active_count_tag(vesalius_assay)
    } else if (length(grep(norm_method, names(counts))) == 0) {
        stop(
            paste0(deparse(substitute(norm_method)), "is not in count list!")
        )
    } else if (length(grep(paste0("^", norm_method, "$"),
        names(counts))) > 1) {
        trial <- grep(paste0("^", norm_method, "$"),
            names(counts), value = TRUE)
        warning(paste("More than one trial contains that name: \n",
            paste(trial, collapse = " ", sep = " "),
            "Using first trial"))
        trial <- trial[1L]
    } else {
        trial <-  trial <- grep(paste0("^", norm_method, "$"),
            names(counts), value = TRUE)
    }
    if (!is.null(method[1L]) && method[1L] %in% c("DESeq2", "QLF", "LRT")) {
        message_switch("norm_check", verbose, method = method)
        if (norm_method == "empty") {
            stop("No raw counts provided! 
                Using user provided count matrix.
                Results will be sub-optimal!")
        }
        counts <- counts[["raw"]]
        counts <- as.matrix(counts)
    } else {
        counts <- as.matrix(counts[[trial]])
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
        warning(paste0("No cells of interest are present in ", ter,
            "\n Returning NULL\n"))
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
        warning(paste0("Territory ", id, " does not contain enough cells.\n
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
#' @importFrom infix %||%
check_group_value <- function(territories, group) {
    #-------------------------------------------------------------------------#
    # If NULL is parsed as group - return NULL
    #-------------------------------------------------------------------------#
    if (is.null(group)) {
        return(NULL)
    }
    present <- unique(territories$trial)
    group_sub <- group[group %in% present]
    if (length(group_sub) == 0) {
        stop(paste("Territory", group, "is not present in the selected trial!"))
    } else if (length(group_sub) < length(group)) {
        tmp <- paste(group_sub, collapse = " ")
        warning(paste("Only group territory(ies)", tmp,
            "is (are) present in the selected trial. \n
            Vesalius will discard the other!"))
        return(group_sub)
    } else {
        return(group_sub)
    }
}

#' check binary nature
#' @param buffer list containing output of get_deg_metrics
#' @param test string - test type of scope
#' @details For now we do this check but might be depreciated in the future.
#' @return buffer list
check_binary_nature <- function(buffer) {
    seed <- apply(buffer$seed, 1, table)
    query <- apply(buffer$query, 1, table)
    #-------------------------------------------------------------------------#
    # first we check if there are more than 2 values 
    # don't need to check those value after 
    #-------------------------------------------------------------------------#
    seed <- sapply(seed, function(x) {
        return(length(x) > 2)
    })
    query <- sapply(query, function(x) {
        return(length(x) > 2)
    })
    if (sum(seed) > 0 || sum(query) > 0) {
        stop("Count data is not binary! 
            Cannot create 2 x 2 contingency table")
    }
    return(buffer)
}


#'check_for_zero_counts 
#' check which barcodes contain all 0 counts and filter
#' @param counts count matrix 
#' @return trimmed count matrix
check_for_zero_counts <- function(counts) {
    zeros <- apply(counts, 2, sum) == 0
    if (sum(zeros) > 0) {
        warning(paste("Trimming count matrix!",
        sum(zeros),
        "spatial indices do not contain any counts for selected genes"))
    }
    counts <- counts[, !zeros]
    return(counts)
}


check_template_source <- function(vesalius_assay, from, dimensions) {
    if (from == "last") {
        stop("Please specifiy by name what should 
        be use to build image template.
        Tou can select from:
        territory, segment, embedding, or active embedding")
    }
    type <- regmatches(from,
        regexpr("Segment|Territory", from))
    if (length(type) == 0) {
        from <- ifelse(from == "active", "last", from)
        template <- check_embedding_selection(vesalius_assay,
            embed = from,
            dims = dimensions)[,dimensions]
        type <- "Embedding"
    } else {
        template <- switch(type,
            "Segment" = check_segment_trial(vesalius_assay, from),
            "Territory" = check_territory_trial(vesalius_assay, from))
    }
    comment(template) <- type
    return(template)
}
