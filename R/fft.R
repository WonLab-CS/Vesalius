################################################################################
################################   Vesalius      ###############################
###############################################################################

#------------------------------/Fourier Transform/----------------------------#

#' @export
integrate_territories <- function(seed_assay,
    query_assay,
    seed_trial = "last",
    query_trial = "last",
    start = "convex",
    partials = 500,
    verbose = TRUE) {
    simple_bar(verbose)
    seed_tiles <- get_tiles(seed_assay)
    query_tiles <- get_tiles(query_assay)
    #-------------------------------------------------------------------------#
    # check and unpack terriotry paths 
    #-------------------------------------------------------------------------#
    seed_assay <- check_territory_trial(seed_assay, seed_trial)
    query_assay <- check_territory_trial(query_assay, query_trial)
    message_switch("unpack_path", verbose)
    seed_path <- unpack_territory_path(trial = seed_assay,
        tiles = seed_tiles,
        start = start,
        verbose = verbose)
    query_path <- unpack_territory_path(trial = query_assay,
        tiles = query_tiles,
        start = start,
        verbose = verbose)
    #-------------------------------------------------------------------------#
    # convert to FS
    #-------------------------------------------------------------------------#
    message_switch("fft", verbose = TRUE)
    ranges_seed <- get_ranges(seed_path)
    ranges_query <- get_ranges(query_path) 
    seed_path <- lapply(seed_path, coordinate_fft,
        partials, ranges_seed)
    query_path <- lapply(query_path, coordinate_fft,
        partials, ranges_query)
    #-------------------------------------------------------------------------#
    # comapring signals - basic correlation for now 
    #-------------------------------------------------------------------------#
    metrics <- signal_similiarity(seed_path, query_path, partials)
    simple_bar(verbose)
    return(list("seed" = seed_path,
        "query" = query_path,
        "cor" = metrics$cor,
        "cov" = metrics$cov))

}

unpack_territory_path <- function(trial,
    tiles,
    start = "convex",
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # First we convert to pixset and detect territory edge 
    #-------------------------------------------------------------------------#
    
    trial_split <- vector("list", length(unique(trial$trial)))
    names(trial_split) <- unique(trial$trial)
    for (i in seq_along(unique(trial$trial))) {
        territory <- unique(trial$trial)[i]
        ter <- right_join(trial, tiles, by = "barcodes") %>%
            filter(trial %in% territory) %>%
            mutate(value = 1) %>%
            select(c("barcodes", "x.y", "y.y", "value", "origin", "trial"))
        colnames(ter) <- c("barcodes", "x", "y", "value", "origin", "trial")
        edge <- extend_boundary(ter, 1) %>%
            detect_edges() %>%
            grow(1) %>%
            as.data.frame()
        edge <- inner_join(edge, ter, by = c("x", "y"))
        trial_split[[i]] <- edge

    }
    #-------------------------------------------------------------------------#
    # next we remove NULLs - this happens when no edge can be deteced 
    #-------------------------------------------------------------------------#
    nulls <- sapply(trial_split, nrow) == 0
    trial_split <- trial_split[!nulls]
    if (sum(nulls) > 0) {
        message_switch("edge_detect", verbose,
            nulls = paste(names(trial_split)[nulls]))
    } else if (length(trial_split) == 0) {
        stop("No edge can be detect in territories! Granularity too high.
         Consider increasing smoothing and/or decreasing segmentationd depth")
    }
    #-------------------------------------------------------------------------#
    # select starting point for path 
    # ATM we convert edge to ordered shape using polar coordinates 
    #-------------------------------------------------------------------------#
    trial <- switch(EXPR = start,
        "convex" = lapply(trial_split, function(trial) {
            ord <- convexify(trial$x,
                trial$y,
                median(trial$x),
                median(trial$y),
                order = TRUE)
            trial <- trial[ord$x, ]
            return(trial)
        }))
    return(trial)
}


get_ranges <- function(path) {
    x_max <- max(sapply(path, function(x) {
        return(max(x$x))
    }))
    x_min <- min(sapply(path, function(x) {
        return(min(x$x))
    }))
    y_max <- max(sapply(path, function(y) {
        return(max(y$y))
    }))
    y_min <- min(sapply(path, function(y) {
        return(min(y$y))
    }))
    return(list("x_max" = x_max,
        "x_min" = x_min,
        "y_max" = y_max,
        "y_min" = y_min))
}

#' converting territory path to FFT 
#' @param path data frame containing ordered x and y coordiates of a territory
#' @param partials int - number of FFT partials to return
#' @param ranges min and max range across coordinates 
coordinate_fft <- function(path, partials, ranges) {
    x <- (ranges$x_max - path$x) / (ranges$x_max - ranges$x_min)
    y <- (ranges$y_max - path$y) / (ranges$y_max - ranges$y_min)
    x <- fft(x)
    y <- fft(y)
    return(list("x" = x, "y" = y))
}

signal_similiarity <- function(seed_path, query_path, partials) {
    

}

circular_xcorr <- function(seed, query) {
    corr <- ccf(seed, query)
}

# n = 100;
# a = rand(1,n); b = rand(1,n);
# time_corr_thresh = .8 * n; freq_corr_thresh = .6 * n;
# time_static = max(a .* b);
# time_shifted = circular_xcorr(a,b);    time_shifted = sum(time_shifted(time_shifted > time_corr_thresh));
# freq_static = max(fft(a) .* fft(b));
# freq_shifted = fft(a .* b);     freq_shifted = sum(freq_shifted(freq_shifted > freq_corr_thresh));
# w1 = 0; w2 = 1; w2 = .7; w3 = 0;
# index = w1 * time_static + w1 * time_shifted + w2 * freq_static + w3 * freq_shifted;