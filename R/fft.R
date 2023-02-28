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
    method = "coherence",
    global = TRUE,
    partials = 500,
    threshold = 0.7,
    weight = rep(0.33, 3),
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
        verbose = verbose) %>%
        normalize_path(global = global)
    query_path <- unpack_territory_path(trial = query_assay,
        tiles = query_tiles,
        start = start,
        verbose = verbose) %>%
        normalize_path(global = global)

    #-------------------------------------------------------------------------#
    # get signal similarity using path signal as in
    #-------------------------------------------------------------------------#
    similarity <- signal_similiarity(seed_path = seed_path,
        query_path = query_path,
        domain = method,
        partials = partials,
        threshold = threshold,
        weight = weight)
    simple_bar(verbose)
    return(list("sim" = similarity,
        "seed" = seed_path,
        "query" = query_path))
}

#' retrieve the points contained in the edge of each territory
#' @param trial name of territory trial that should be selected
#' @param tiles vesalius tiles 
#' @param start string - which point should be used as starting point
#' @param verbose logical - should progress message be print to console
#' @details Here we are using convex as start point. Essentially, we 
#' order the coordinates based on their polar coordinates using the 
#' median coordinate as the center point. 
#' @returns a data frame containing edge of each territory.
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
        edge <- inner_join(edge, ter, by = c("x", "y")) %>%
            select("barcodes") %>% unique()
        edge <- tiles %>% filter(barcodes %in% edge$barcodes & origin == 1)
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
        }),
        "connected" = lapply(trial_split, connected_points))
    return(trial)
}

#' create path from neighboring points 
#' @param trail data frame containign x and y path
#' @return order x y coordinates
#' @importFrom RANN nn2
#' 
connected_points <- function(trial) {
    knn <- RANN::nn2(trial[, c("x","y")], k = nrow(trial))$nn.idx
    ord <- rep(0, nrow(knn))
    ord[1] <- 1
    for (i in seq(2, nrow(knn))) {
        if (i == 2) {
            ord[2] <- knn[(i - 1), 2]
        } else {
            ord[i] <- knn[ord[(i - 1)], min(which(!knn[ord[(i - 1)], ] %in% ord))]
        }
    }
    return(trial[ord, ])
}


#' normalise coordinates of territory edge
#' @param path data frame of coordinates corresponding to 
#' the x and y coordinates of the edge of each territory.
#' @param global logical - using global mion max values or local.
#' @returns list of nromalised x and y coordinates
normalize_path <- function(path, global = TRUE) {
    if (global) {
        ranges <- get_ranges(path)
    } else {
        ranges <- NULL
    }
    normed <- lapply(path, function(path, ranges = NULL) {
        if (is.null(ranges)) {
            ranges <- list()
            ranges$x_max <- max(path$x)
            ranges$x_min <- min(path$x)
            ranges$y_max <- max(path$y)
            ranges$y_min <- min(path$y)
        }
        x <- (ranges$x_max - path$x) / (ranges$x_max - ranges$x_min)
        y <- (ranges$y_max - path$y) / (ranges$y_max - ranges$y_min)
        return(list("x" = x, "y" = y))
    }, ranges = ranges)
    return(normed)
}

#' get coordinate ranges
#' @param path list of x and y coordinates for each terrirtory
#' @return list containing ranges in x dim and y dim 
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


#' compute signal similarity between 2 sets of territories
#' @param seed_path list of x and y coordinates of seed territories
#' @param query_path list of x and y coordinates of query territories
#' @param domain similarity method to use (index or dtw)
#' @param partials int - number of FFT partials to use
#' @param threshold numeric time/freq threshold
#' @param weight vector of 3 numeric weights - weight used in index calcultion
#' @returns a matrix containing similarity scores between each territory in each 
#' data set. 
signal_similiarity <- function(seed_path,
    query_path,
    domain = "index",
    partials = 200,
    threshold = 0.7,
    weight = rep(0.33, 3)) {
    sim <- switch(EXPR = domain,
        "index" = index(seed_path,
            query_path,
            partials,
            threshold,
            weight),
        "dtw" = dynamic_time_warp(seed_path,
            query_path),
        "coherence" = spectral_coherence(seed_path, query_path))
    return(sim)
}

#' circular cross correlation
#' @param seed seed signal as vector 
#' @param query query signal as vector
#' @returns circularised cross correlation
#' @importFrom stats ccf
circular_xcorr <- function(seed, query) {
    corr <- ccf(seed, query, plot = FALSE)$acf
    mid <- length(corr) / 2
    corr <- corr[seq(1, mid)] + corr[seq(mid + 1, length(corr))]
    return(corr)

}

#' calculate weighted index of similarity between 2 signals
#' @param seed_path list of x and y coordinates of seed territories
#' @param query_path list of x and y coordinates of query territories
#' @param partials int - number of FFT partials to use
#' @param threshold numeric time/freq threshold
#' @param weight vector of 3 numeric weights - weight used in index calcultion
#' @returns matrix of similarity indices between each territory
#' @importFrom stats fft
index <- function(seed_path,
    query_path,
    partials = 500,
    threshold = 0.7,
    weight = rep(0.3, 3)) {
    indexed_x <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    indexed_y <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    for (seed in seq_along(seed_path)){
        for (query in seq_along(query_path)){
            #-----------------------------------------------------------------#
            # Get partial max and subset paths
            #-----------------------------------------------------------------#
            seed_sub <- seed_path[[seed]]
            query_sub <- query_path[[query]]
            sub_path <- min(c(length(seed_sub$x), length(query_sub$x)))
            sub_path <- sort(sample(seq(1, sub_path),
                sub_path,
                replace = FALSE))
            seed_sub_x <- seed_sub$x[sub_path]
            query_sub_x <- query_sub$x[sub_path]
            seed_sub_y <- seed_sub$y[sub_path]
            query_sub_y <- query_sub$y[sub_path]
            #-----------------------------------------------------------------#
            # compute signal similarity
            #-----------------------------------------------------------------#
            partial_min <- min(c(length(sub_path), partials))
            thresh <- threshold * partial_min
           
            time_stat_x <- max(seed_sub_x * query_sub_x)
            time_stat_y <- max(seed_sub_y * query_sub_y)

            time_shifted_x <- circular_xcorr(seed_sub_x, query_sub_x)
            time_shifted_x <- sum(time_shifted_x[time_shifted_x > thresh])
            time_shifted_y <- circular_xcorr(seed_sub_y, query_sub_y)
            time_shifted_y <- sum(time_shifted_y[time_shifted_y > thresh])

            freq_static_x <- max(Re(fft(seed_sub_x))[seq(1, partial_min)] *
                Re(fft(query_sub_x))[seq(1, partial_min)])
            freq_static_y <- max(Re(fft(seed_sub_y))[seq(1, partial_min)] *
                Re(fft(query_sub_y))[seq(1, partial_min)])

            freq_shifted_x <- Re(fft(seed_sub_x * query_sub_x))
            freq_shifted_x <- sum(freq_shifted_x[freq_shifted_x > thresh])
            freq_shifted_y <- Re(fft(seed_sub_y * query_sub_y))
            freq_shifted_y <- sum(freq_shifted_y[freq_shifted_y > thresh])
            indexed_x[query, seed] <- weight[1] * time_stat_x +
                weight[1] * time_shifted_x +
                weight[2] * freq_static_x +
                weight[3] * freq_shifted_x
            indexed_y[query, seed] <- weight[1] * time_stat_y +
                weight[1] * time_shifted_y +
                weight[2] * freq_static_y +
                weight[3] * freq_shifted_y
        }
    }
    colnames(indexed_x) <- names(seed_path)
    rownames(indexed_x) <- names(query_path)
    colnames(indexed_y) <- names(seed_path)
    rownames(indexed_y) <- names(query_path)
    return(list("x" = indexed_x, "y" = indexed_y))
}

# #' dynamic time warping for signal similarity
# #' @param seed_path list of x and y coordinates of seed territories
# #' @param query_path list of x and y coordinates of query territories
# #' #' @returns matrix of dynamic time warping distance between each territory
# #' @importFrom dtw dtw
# dynamic_time_warp <- function(seed_path, query_path) {
#     indexed_x <- matrix(0,
#         ncol = length(seed_path),
#         nrow = length(query_path))
#     indexed_y <- matrix(0,
#         ncol = length(seed_path),
#         nrow = length(query_path))
#     for (seed in seq_along(seed_path)){
#         for (query in seq_along(query_path)){
#             indexed_x[query, seed] <- dtw::dtw(seed_path$x,
#                 query_path$x)$distance
#             indexed_y[query, seed] <- dtw::dtw(seed_path$y,
#                 query_path$y)$distance
#         }
#     }
#     colnames(indexed_x) <- names(seed)
#     rownames(indexed_x) <- names(query)
#     colnames(indexed_y) <- names(seed)
#     rownames(indexed_y) <- names(query)
#     return(list("x" = indexed_x, "y" = indexed_y))
# }


spectral_coherence <- function(seed_path, query_path) {
     indexed_x <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    indexed_y <- matrix(0,
        ncol = length(seed_path),
        nrow = length(query_path))
    for (seed in seq_along(seed_path)){
        for (query in seq_along(query_path)){
            #-----------------------------------------------------------------#
            # Get partial max and subset paths
            #-----------------------------------------------------------------#
            seed_sub <- seed_path[[seed]]
            query_sub <- query_path[[query]]
            sub_path <- min(c(length(seed_sub$x), length(query_sub$x)))
            sub_path <- sort(sample(seq(1, sub_path),
                sub_path,
                replace = FALSE))
            x <- cbind(seed_sub$x[sub_path], query_sub$x[sub_path])
            y <- cbind(seed_sub$y[sub_path], query_sub$y[sub_path])
            #-----------------------------------------------------------------#
            # coherence 
            #-----------------------------------------------------------------#
            indexed_x[query, seed] <- max(gsignal::mscohere(x)$coh)
            indexed_y[query, seed] <- max(gsignal::mscohere(y)$coh)
        }
    }
    colnames(indexed_x) <- names(seed_path)
    rownames(indexed_x) <- names(query_path)
    colnames(indexed_y) <- names(seed_path)
    rownames(indexed_y) <- names(query_path)
    return(list("x" = indexed_x, "y" = indexed_y))
}