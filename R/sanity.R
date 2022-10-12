################################################################################
###############################   Vesalius      ################################
################################################################################

#------------------------------/Sanity Checks /--------------------------------#

check_assays <- function(assays, n_assays, verbose) {
    #--------------------------------------------------------------------------#
    # we want to make sure that they are not unique names
    #--------------------------------------------------------------------------#
    if (length(assays) != length(unique(assays))) {
        message_switch("assay_non_unique_name",verbose)
        assays <- paste0("spatial_omics_", seq(1,n_assays))
    } else if (length(assays) == 1 && assays[1L] == "spatial_omics") {
        message_switch("assay_name", verbose)
        assays <- paste0("spatial_omics_", seq(1,n_assays))
    } else if (!identical(length(assays),n_assays)) {
        stop("Number assays provided do not match number of the number \n
        of counts and coordinates provided")
    }
    return(assays)
}

check_inputs <- function(counts,
    coordinates,
    assays,
    adjust_coordinates,
    verbose) {
    #--------------------------------------------------------------------------#
    # first we check if we have the same length lists 
    # if we only have a single object parsed  we put this into a list
    # in this case it should be equal to 1 anyway
    # We also check if assays have the same length 
    #--------------------------------------------------------------------------#
    if (!is(coordinates, "list")) {
        coordinates <- list(coordinates)
    }
    if (!is(counts, "list")) {
        counts <- list(counts)
    }
    if (length(counts) != length(coordinates)) {
        stop("Number of count matrices and coordinates do not match")
    }
    #--------------------------------------------------------------------------#
    # since we already know that there are the same number of counts and
    # coordinates we can use only one length value here
    #--------------------------------------------------------------------------#
    assays <- check_assays(assays, length(counts), verbose)
    #--------------------------------------------------------------------------#
    # now we can check the validity of each
    #--------------------------------------------------------------------------#
    counts <- mapply(check_counts,
        counts, assays,
        MoreArgs = list(verbose),
        SIMPLIFY = FALSE)
    coordinates <- mapply(check_coordinates,
        coordinates, assays,
        MoreArgs = list(adjust_coordinates,
            verbose),
        SIMPLIFY = FALSE)
    #--------------------------------------------------------------------------#
    # Finally we can rebuild an assay list filled with vesalius_assay objects 
    #--------------------------------------------------------------------------#
    vesalius <- mapply(build_vesalius_assay,
        counts = counts,
        coordinates = coordinates,
        assay = assays,
        MoreArgs = list(verbose = FALSE),
        SIMPLIFY = FALSE)
    names(vesalius) <- assays
    return(vesalius)
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

check_norm_methods <- function(norm, n_counts) {
    if (any(!norm %in% c("log_norm", "SCTransform", "TFIDF", "raw"))) {
        stop("Normalisation method provided do not match available options \n
            Select from: \n
            log_norm, SCTransform, TFIDF, or raw")
    }
    if (length(norm) != 1 && length(norm) != n_counts) {
        stop("Number of normalisation methods provided do not match 
            number of count matrices present.")
    }
    norm <- rep(norm, times = n_counts)
    return(norm)
}
check_embed_methods <- function(embed, n_counts) {
    if (any(!embed %in% c("PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP"))) {
        stop("Embedding method provided do not match available options \n
            Select from: \n
            PCA, PCA_L, UMAP, LSI, LSI_UMAP")
    }
    if (length(embed) != 1 && length(embed) != n_counts) {
        stop("Number of embedding methods provided do not match 
            number of count matrices present.")
    }
    embed <- rep(embed, times = n_counts)
    return(embed)
}

check_embedding <- function(vesalius_assay, embed, dims) {
    #--------------------------------------------------------------------------#
    # Lets check if we have all the right embeddings
    #--------------------------------------------------------------------------#
    if (!is(vesalius_assay, "vesaliusObject")) {
        stop("Unsupported format! Make sure you are parsing a vesalius object")
    }
    if (embed == "last") {
      embeddings <- vesalius_assay@activeEmbeddings[[1L]]
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
        stop("Unsupported format! Make sure you are parsing a vesalius object")
    }
    tiles <- vesalius@tiles
    return(tiles)
}

check_tile_log <- function(tile_log) {
    l_tile_log <- length(tile_log)
    s_tile_log <- sapply(tile_log, length)
    if (all(s_tile_log == 0)) {
        return(0)
    } else if (sum(s_tile_log) < l_tile_log) {
        stop("No all tiles have been computed in vesalius object \n
            Please make sure you are adding processed vesalius objects")
    } else {
        return(NULL)
    }
}

check_territories <- function(vesalius, trial) {
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format! Make sure you are parsing a vesalius object")
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
        stop("Unsupported format! Make sure you are parsing a vesalius object")
    }
    territories <- vesalius@territories %||%
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

check_norm <- function(vesalius, norm_method, method = NULL, verbose = TRUE) {
    if (!is(vesalius, "vesaliusObject")) {
        stop("Unsupported format! Make sure you are parsing a vesalius object")
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
    if (!is.null(method) && method %in% c("DEseq2", "edgeR")) {
        message_switch("norm_check", verbose, method = method)
        counts <- as.matrix(counts[["raw"]])
    } else {
        counts <- as.matrix(counts[[norm_method]])
    }
    return(counts)
}

territory_dispatch <- function(territories, ter_1, ter_2, cells) {
    if (is.null(ter_1) && is.null(ter_2)) {
        territories <- select(territories, c("barcodes", "x", "y", "trial"))
    }else if (!is.null(ter_1) && is.null(ter_2)) {
        territories$trial[!territories$trial %in% ter_1] <- "other"
    }else if (is.null(ter_1) && !is.null(ter_2)) {
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

deg_group_dispatch <- function(ter, seed, query, cells, verbose) {
    if (is.null(seed) && is.null(query)) {
        #----------------------------------------------------------------------#
        # If no territories are provided - then we assume that the user
        # want to look at all territories.
        # This compares each territory to everything else
        #----------------------------------------------------------------------#
        message_switch("deg_dispatch_all_null", verbose)
        seed <- split(ter$barcodes, ter$trial)
        seed_id <- names(seed)

        query <- lapply(seed, function(bar, ter) {
            return(ter$barcodes[!ter$barcodes %in% bar])
        }, ter = ter)
        query_id <- rep("remaining", length(query))
    } else if (!is.null(seed) && is.null(query)) {
        #----------------------------------------------------------------------#
        # if only seed we compare seed to everything else
        # Get initial seed territory
        #----------------------------------------------------------------------#
        seed_id <- paste0(seed, collapse = " ")
        seed <- list(ter[ter$trial %in% seed, "barcodes"])
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        query_id <- "remaining"
        query <- list(ter[!ter$barcodes %in% seed, "barcodes"])
    } else if (is.null(seed) && !is.null(query)) {
        #----------------------------------------------------------------------#
        # if only query we compare query to everything else
        # Get initial query territory
        #----------------------------------------------------------------------#
        query_id <- paste0(query, collapse = " ")
        query <- list(ter[ter$trial %in% query, "barcodes"])
        #----------------------------------------------------------------------#
        # Filter seed based on query
        #----------------------------------------------------------------------#
        seed_id <- "remaning"
        seed <- list(ter[!ter$barcodes %in% query, "barcodes"])
    } else {
        #----------------------------------------------------------------------#
        # if get both filter based on both
        #----------------------------------------------------------------------#
        seed_id <- paste0(seed, collapse = " ")
        seed <- list(ter[ter$trial %in% seed, "barcodes"])
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        query_id <- paste0(query, collapse = " ")
        query <- list(ter[ter$trial %in% query, "barcodes"])
    }
    seed <- mapply(check_cells, territory_barcodes = seed,
        ter = seed_id, MoreArgs = list(cells, verbose), IMPLIFY = FALSE)
    query <- mapply(check_cells, territory_barcodes = query,
        ter = query_id, MoreArgs = list(cells, verbose), SIMPLIFY = FALSE)
    return(list("seed" = seed, "seed_id" = seed_id,
        "query" = query, "query_id" = query_id))
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

