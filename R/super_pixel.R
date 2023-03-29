################################################################################
################################   Vesalius      ###############################
################################################################################

#------------------------------/Super pixels/----------------------------------#

#' select inital indices
#' @param coordinates data frame containing spatial coordinates of beads
#' @param embeddings matrix containing embedding values - full pixel image
#' @param n_centers numeric number of beads to select as super pixel centers
#' @param max_iter numeric number of iteration before returning result if 
#' no coveregnce.
#' @return barcodes of starting coordinates
#' @importFrom RANN nn2
select_initial_indices <- function(coordinates,
    embeddings,
    type = "bubble",
    n_centers = 500,
    max_iter = 500) {
    indices <- switch(type,
        "linear" = uniform_sampling(coordinates,
            embeddings,
            n_centers),
        "bubble" = bubble_stack(coordinates,
            embeddings,
            n_centers,
            max_iter))
}

linear_sampling <- function(coordinates,
    embeddings,
    n_centers = 500) {
    coordinates <- coordinates[order(coordinates$x, decreasing = FALSE),]
    coordinates <- split(coordinates, coordinates$x) %>%
        lapply(., function(df) {return(df[order(df$y, decreasing = FALSE),])})
    coordinates <- do.call("rbind", coordinates)
    coordinates <- coordinates[seq(1, nrow(coordinates), l = n_centers),]

    in_image <- paste0(embeddings[, "x"], "_", embeddings[, "y"])
    in_background <- paste0(coordinates$x, "_", coordinates$y)
    in_image <- which(in_image %in% in_background)
    return(in_image)
}

bubble_stack <- function(coordinates,
    embeddings,
    n_centers = 500,
    max_iter = 500) {
    coordinates <- coordinates %>% filter(origin == 1)
    #-------------------------------------------------------------------------#
    # intialise background grid, active grod and get nearest neighbors
    # we assume that you will have at least 4 inital points
    #-------------------------------------------------------------------------#
    knn <- RANN::nn2(coordinates[, c("x", "y")], k = nrow(coordinates))
    radius <- max(knn$nn.dist)
    #-------------------------------------------------------------------------#
    # iterated a reduced space until you get the desired numbers of circles  
    #-------------------------------------------------------------------------#
    convergence <- FALSE
    iter <- 1
    while (!convergence) {
        #---------------------------------------------------------------------#
        # intialize everything we need 
        #---------------------------------------------------------------------#
        background_grid <- c()
        active <- seq_len(nrow(coordinates))
        nn_index <- knn$nn.idx
        nn_dist <- knn$nn.dist
        #---------------------------------------------------------------------#
        # as long as their points present we repeat the process of selecting 
        # point that are outisde the selection radius
        #---------------------------------------------------------------------#
        while (length(active) > 0) {
            random_start <- active[sample(x =
                seq(1, l = length(active)),
                size = 1)]
            if (random_start %in% background_grid) browser()
            background_grid <- c(background_grid, random_start)
            removing <- nn_index[random_start,
                nn_dist[random_start, ] <= radius]
            active <- active[!active %in% unique(removing)]
        }
        if (length(background_grid) == n_centers) {
            convergence <- TRUE
        }
        if (iter >= max_iter) {
            convergence <- TRUE
            warning("Max Iteration with not convergence!
            Returning Approximation")
        }
        #---------------------------------------------------------------------#
        # we change the radius based on how many points there were
        # we assume that if there are 2 few barcodes then the radius is to 
        # large 
        # if there are too many barcodes then the radius is 2 small. 
        #---------------------------------------------------------------------#
        if (length(background_grid) < n_centers) {
            radius <- radius - (radius / 2)
        } else {
            radius <- radius + (radius / 2)
        }
        iter <- iter + 1
    }
    #-------------------------------------------------------------------------#
    # Next we need to find what are the indices of these barcodes 
    # within the full pixel image 
    # could use right_join and the likes but no 
    #-------------------------------------------------------------------------#
    in_image <- paste0(embeddings[, "x"], "_", embeddings[, "y"])
    in_background <- paste0(coordinates$x[background_grid],
        "_",
        coordinates$y[background_grid])
    in_image <- which(in_image %in% in_background)
    return(in_image)
}








# #' @importFrom imager imappend imsplit nPix spectrum
# #' @importFrom purrr map_dbl
# slic_segmentation <- function(vesalius_assay,
#     dimensions,
#     col_resolution,
#     embedding,
#     compactness = 1,
#     verbose) {
#     coord <- get_tiles(vesalius_assay) %>%
#     filter(origin == 1)
#     embeddings <- check_embedding_selection(vesalius_assay,
#     embedding,
#     dimensions)
    
#     sc_spat <- max(c(max(coord$x), max(coord$y))  * 0.28)
#     sc_col <-  max(apply(embeddings, 2, sd))

#     #Scaling ratio for pixel values
#     ratio <- (sc_spat / sc_col) / (compactness * 10)
#     embeddings <- embeddings * ratio
#     #Generate initial centers from a grid
#     ind <- round(seq(1, nrow(coord), l = col_resolution))
#     #Run k-means
#     km <- kmeans(embeddings, embeddings[ind, ], iter.max = 100, nstart = 10)

#     clusters <- km$cluster
#     kcenters <- km$centers
#     match_loc <- match(coord$barcodes, names(clusters))
#     clusters <- data.frame(coord, "Segment" = clusters[match_loc])
#     active <- assign_centers(vesalius_assay,
#         clusters,
#         kcenters,
#         dimensions,
#         ratio)
#     new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
#         "Segment") %>%
#         tail(1)
#     colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
#     return(list("segments" = clusters, "active" = active))

# }


#' @importFrom imager imappend imsplit spectrum
#' @importFrom purrr map_dbl map
slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    index_selection = "bubble",
    compactness = 1,
    scaling = 0.3,
    threshold = 0.9,
    max_iter = 5000,
    verbose) {
    if (is(vesalius_assay, "vesalius_assay")) {
      assay <- get_assay_names(vesalius_assay)
      images <- format_ves_to_c(vesalius_assay = vesalius_assay,
        embed = embedding,
        dims = dimensions,
        verbose = FALSE) %>% imappend("cc")
    } else {
      stop("Unsupported format to regularise_image function")
    }
    #-------------------------------------------------------------------------#
    # first we get tiles and get images 
    # then we compute the a scaling metric - this is normally based on 
    # empirical data (0.28 is one way) but it essentially a measure of the
    # scaling value between pixel distance and color distance 
    # we measure the max spread of colors as well
    # all this serves as a way to scale the "spread" of the color values 
    # to match the "spread" of the spatial coordinates
    #-------------------------------------------------------------------------#
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    tiles <- get_tiles(vesalius_assay)
    # max spatial distance 
    sc_spat <- max(dim(images)[1:2]  * scaling)
    # max color distance between two pixels??
    sc_col <-  imsplit(images, "cc") %>% map_dbl(sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(images * ratio, wide = "c") %>%
        right_join(tiles, by = c("x", "y")) %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    #-------------------------------------------------------------------------#
    # Generate initial centers from a grid
    # This will need to be updated to make sure
    # At the moment this will take barcodes within the data set to use as 
    # inital centers.
    # It would be worth looking into using pixels instead. Mainly for low
    # res data sets - we might want to have pixels between spatial indices 
    #-------------------------------------------------------------------------#
    index <- select_initial_indices(get_tiles(vesalius_assay),
        embeddings,
        type = index_selection,
        n_centers = col_resolution,
        max_iter = max_iter)
    #-------------------------------------------------------------------------#
    # Run k-means - you should have read the mnual pages you silly goose
    #-------------------------------------------------------------------------#
    km <- suppressWarnings(kmeans(embeddings,
        embeddings[index, ],
        iter.max = 100,
        nstart = 10))
    centroids <- map(1:spectrum(images), ~ km$centers[km$cluster, 2 + .]) %>%
        do.call("cbind", .)
    centroids <- centroids / ratio
    clusters <- cbind(embeddings[, c("x", "y")], km$cluster) %>%
        as.data.frame() %>%
        right_join(coord, by = c("x", "y"))
    colnames(clusters) <- c("x", "y", "value", "barcodes", "origin")
    match_loc <- match(coord$barcodes, clusters$barcodes)
    clusters <- data.frame(coord, "Segment" = clusters[match_loc, "value"])
    #-------------------------------------------------------------------------#
    # rescaling back to original values
    # rebuilding everything 
    # ON HOLD: combining super pixels into large scale segments
    #-------------------------------------------------------------------------#
    embeddings[, seq(3, ncol(embeddings))] <- centroids
    embeddings <- lapply(seq_len(ncol(centroids)), function(i, embed) {
        ret <- as.data.frame(embed[, c("x", "y", paste0("c.", i))])
        colnames(ret) <- c("x", "y", "value")
        return(ret)
    }, embed = embeddings)
    #clusters <- select_similar(embeddings, coordinates = coord)
    embeddings <- format_c_to_ves(embeddings,
      vesalius_assay,
      dimensions,
      embed = embedding,
      verbose = FALSE)
    #clusters <- connected_pixels(clusters, embeddings, k, threshold, verbose)
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
        "Segment") %>%
    tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = embeddings))
}



generate_slic_graph <- function(spix,
    signals,
    k = "auto",
    scoring_method = "pearson") {
    #-------------------------------------------------------------------------#
    # first we estimate the pixel centers 
    # we will use this as an estimate for nearest neighbor calculation 
    #-------------------------------------------------------------------------#
    center_pixels <- sort(unique(spix$segment))
    centers <- future_lapply(center_pixels, function(center, segments){
        x <- median(segments$x[segments$segment == center])
        y <- median(segments$y[segments$segment == center])
        df <- data.frame("x" = x, "y" = y)
        rownames(df) <- center
        return(df)
    }, segments = spix) %>% do.call("rbind", .)
    #-------------------------------------------------------------------------#
    # Next we get the nearest neighbors
    # If we use the auto option we compute the nearest neighbors based on 
    # delauney traingulation - not that this will results in an uneven number 
    # of neighbors per point. 
    #-------------------------------------------------------------------------#
    if (k != "auto") {
        k <- min(c(k, nrow(centers)))
        knn <- RANN::nn2(centers, k = k)$nn.idx
        rownames(knn) <- rownames(centers)
        graph <- populate_graph(knn)
    } else {
        voronoi <- deldir::deldir(x = as.numeric(centers$x),
            y = as.numeric(centers$y))$delsgs
        center <- seq_len(nrow(centers))
        graph <- lapply(center, function(idx, voronoi){
            tri <- voronoi %>% filter(ind1 == idx | ind2 == idx)
            tri <- unique(c(tri$ind1, tri$ind2))
            graph <- data.frame("e1" = rep(idx, length(tri)),
                "e2" = tri)
            return(graph)
        }, voronoi = voronoi) %>% do.call("rbind",.)
    }
    
    #-------------------------------------------------------------------------#
    # Next we score the graph to see whcih neighbors are similar in color
    #-------------------------------------------------------------------------#
    graph <- score_graph(graph$e1,
        graph$e2,
        signal = signals,
        centers = spix,
        scoring_method = scoring_method)
    return(graph)
}


score_graph <- function(g1,
    g2,
    signal,
    centers,
    scoring_method = "pearson") {
    #-------------------------------------------------------------------------#
    # assuming that if the input to signal is a list 
    # we are comparing 2 data sets
    #-------------------------------------------------------------------------#
    if (is(signal, "list") && length(signal) == 2) {
        seed_signal <- signal[[1]]
        query_signal <- signal[[2]]
    } else {
        seed_signal <- query_signal <- signal
    }
    if (is(centers, "list") && length(centers) == 2) {
        seed_centers <- centers[[1]]
        query_centers <- centers[[2]]
    } else {
        seed_centers <- query_centers <- centers
    }
    #-------------------------------------------------------------------------#
    # Create a score graph 
    # and compute correlation between super pixels
    #-------------------------------------------------------------------------#
    graph <- data.frame("g1" = g1,
        "g2" = g2,
        "cor" = rep(0, length(g1)))
    for (i in seq_len(nrow(graph))) {
        c1 <- seed_signal[, seed_centers$barcodes[
            seed_centers$segment == graph$g1[i]]]
        if (!is.null(nrow(c1))) {
            c1 <- apply(c1, 1, mean)
        }
        c2 <- query_signal[, query_centers$barcodes[
            query_centers$segment == graph$g2[i]]]
        if (!is.null(nrow(c2))) {
            c2 <- apply(c2, 1, mean)
        }
        graph$cor[i] <- cor(c1, c2, method = scoring_method)
    }
    return(graph)
}