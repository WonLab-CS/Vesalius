#' Generate super pixels from ST
#'
#' 
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions numeric vector of latent space dimensions to use.
#' @param embedding character string describing which embedding should
#' be used.
#' @param method character string for which method should be used for
#' segmentation. Select from "slic", 
#' "leiden_slic","louvain_slic"
#' @param col_resolution numeric colour resolution used for segmentation. 
#' (see details)
#' @param compactness numeric - factor defining super pixel compaction.
#' @param scaling numeric - scaling image ration during super pixel 
#' segmentation.
#' @param k numeric - number of closest super pixel neighbors to consider
#' when generating segments from super pixels
#' @param threshold numeric [0,1] - correlation threshold between 
#' nearest neighbors when generating segments from super pixels.
#' @param verbose logical - progress message output.
#
#' @return a vesalius_assay object

#' @export

generate_spix <- function(vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "kmeans",
  col_resolution = 10,
  compactness = 1,
  scaling = 0.5,
  threshold = 0.9,
  index_selection = "bubble",
  verbose = TRUE) {
  simple_bar(verbose)
  message_switch("seg", verbose, method = method)
  #----------------------------------------------------------------------------#
  # Parsing vesalius object so we can recontruct it internally and not need to
  # rebuild intermediates and shift between formats...
  #----------------------------------------------------------------------------#
  method <- check_segmentation_method(method)
  segments <- switch(method,
    "slic" = slic_segmentation(vesalius_assay,
      dimensions = dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      index_selection = index_selection,
      compactness = compactness,
      scaling = scaling,
      threshold = threshold,
      verbose = verbose),
    "louvain_slic" = louvain_slic_segmentation(vesalius_assay,
        dimensions = dimensions,
        col_resolution = col_resolution,
        embedding = embedding,
        compactness = compactness,
        scaling = scaling,
        verbose = verbose),
    "leiden_slic" = leiden_slic_segmentation(vesalius_assay,
        dimensions = dimensions,
        col_resolution = col_resolution,
        embedding = embedding,
        compactness = compactness,
        scaling = scaling,
        verbose = verbose))

  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = segments$active,
    slot = "active",
    append = FALSE)
  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = segments$spix,
    slot = "territories",
    append = TRUE)
  vesalius_assay <- add_active_embedding_tag(vesalius_assay, embedding)
  commit <- create_commit_log(arg_match = as.list(match.call()),
    default = formals(segment_image))
  vesalius_assay <- commit_log(vesalius_assay,
    commit,
    get_assay_names(vesalius_assay))
  simple_bar(verbose)
  return(vesalius_assay)
}

#' @importFrom imager imappend imsplit spectrum
louvain_slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    compactness = 1,
    scaling = 0.3,
    verbose) {
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
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    loc <- match(coord$barcodes, rownames(embeddings))
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    embeddings$barcodes <- barcodes <- rownames(embeddings)
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    rownames(embeddings) <- barcodes
    #-------------------------------------------------------------------------#
    # Run louvain
    #-------------------------------------------------------------------------#
    graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
    clusters <- igraph::cluster_louvain(graph, resolution = col_resolution)
    cluster <- data.frame("cluster" = clusters$membership,
      "barcodes" = clusters$names)


    match_loc <- match(coord$barcodes, cluster$barcodes)
    clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
    active <- create_pseudo_centroids(vesalius_assay,
      clusters,
      dimensions)
    active <- active / ratio
    
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
      "SPIX") %>%
      tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("spix" = clusters, "active" = active))
}


#' @importFrom imager imappend imsplit spectrum
leiden_slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    compactness = 1,
    scaling = 0.3,
    verbose) {
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
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    loc <- match(coord$barcodes, rownames(embeddings))
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    embeddings$barcodes <- barcodes <- rownames(embeddings)
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    rownames(embeddings) <- barcodes
    #-------------------------------------------------------------------------#
    # Run louvain
    #-------------------------------------------------------------------------#
    graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
    clusters <- igraph::cluster_leiden(graph,
      resolution_parameter = col_resolution)
    cluster <- data.frame("cluster" = clusters$membership,
      "barcodes" = clusters$names)


    match_loc <- match(coord$barcodes, cluster$barcodes)
    clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
    active <- create_pseudo_centroids(vesalius_assay,
      clusters,
      dimensions)
    active <- active / ratio
    
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
      "SPIX") %>%
      tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("spix" = clusters, "active" = active))
}



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
    max_iter = 1000,
    verbose) {
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
    
    #tiles <- get_tiles(vesalius_assay)
    embeddings <- check_embedding_selection(vesalius_assay,
      embedding,
      dimensions)[, dimensions]
    # max spatial distance 
    sc_spat <- max(c(max(coord$x), max(coord$y))  * scaling)
    # max color distance between two pixels??
    sc_col <-  apply(embeddings, 2, sd) %>% max()
    #-------------------------------------------------------------------------#
    # Scaling ratio for pixel values - this is what is going to effectively
    # scale color value to match the scale of spatial value 
    # the compactness will define how important the spatial component should 
    # be. Now the way this works here is that we are modulating "color" value
    # and not the pixel value. The higher the compactness the "higher" the 
    # color values will be scaled up to match the spatial coordinates
    #-------------------------------------------------------------------------# 
    ratio <- (sc_spat / sc_col) / (compactness)
    embeddings <- as.data.frame(embeddings * ratio)
    embeddings$barcodes <- rownames(embeddings)
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    #-------------------------------------------------------------------------#
    # Generate initial centers from a grid
    # This will need to be updated to make sure
    # At the moment this will take barcodes within the data set to use as 
    # inital centers.
    # It would be worth looking into using pixels instead. Mainly for low
    # res data sets - we might want to have pixels between spatial indices 
    #-------------------------------------------------------------------------#
    index <- select_initial_indices(coord,
        #embeddings,
        type = index_selection,
        n_centers = col_resolution,
        max_iter = max_iter)
    
    #-------------------------------------------------------------------------#
    # Run k-means - you should have read the mnual pages you silly goose
    #-------------------------------------------------------------------------#
    km <- suppressWarnings(kmeans(embeddings,
        embeddings[index, ],
        iter.max = max_iter))
    centroids <- purrr::map(seq(1, l = ncol(embeddings) - 2),
        ~ km$centers[km$cluster, .]) %>%
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
    
    embeddings[, seq(1, ncol(embeddings) - 2)] <- centroids

    embeddings <- lapply(seq_len(ncol(centroids)), function(i, embed) {
        ret <- as.data.frame(embed[, c("x", "y", paste0("dim_", i))])
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
        "SPIX") %>%
    tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("spix" = clusters, "active" = embeddings))
}




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
        # "linear" = uniform_sampling(coordinates,
        #     embeddings,
        #     n_centers),
        "bubble" = bubble_stack(coordinates,
            #embeddings,
            n_centers,
            max_iter),
        # "pixel" = pixel_sampling(coordinates,
        #     embeddings,
        #     n_centers),
        "random" = random_sampling(coordinates,
          n_centers = n_centers),
        "hex" = hex_grid(coordinates,
          n_centers) )
    return(indices)
}

# linear_sampling <- function(coordinates,
#     embeddings,
#     n_centers = 500) {
#     coordinates <- coordinates[order(coordinates$x, decreasing = FALSE),]
#     coordinates <- split(coordinates, coordinates$x) %>%
#         lapply(., function(df) {return(df[order(df$y, decreasing = FALSE),])})
#     coordinates <- do.call("rbind", coordinates)
#     coordinates <- coordinates[seq(1, nrow(coordinates), l = n_centers), ]

#     in_image <- paste0(embeddings[, "x"], "_", embeddings[, "y"])
#     in_background <- paste0(coordinates$x, "_", coordinates$y)
#     in_image <- which(in_image %in% in_background)
#     return(in_image)
# }

random_sampling <- function(coordinates, n_centers) {
    #idx <- which(coordinates$origin == 1)
    return(sample(seq_len(nrow(coordinates)),
      size = n_centers, replace = FALSE))
}

#' importFrom imager nPix
# pixel_sampling <- function(coordinates, embeddings, n_centers) {
#   images <- coordinates[, c("x", "y")]
#   images$value <- coordinates$origin
#   images <- suppressWarnings(as.cimg(images))
#   ind <- round(seq(1, nPix(images) / spectrum(images), l = n_centers))
#   images <- as.data.frame(images)[ind, ]
#   in_background <- paste0(images[, "x"], "_", images[, "y"])
#   in_image <- paste0(embeddings[,"x"], "_", embeddings[,"y"])
#   in_image <- which(in_image %in% in_background)
#   return(in_image)
# } 

bubble_stack <- function(coordinates,
    n_centers = 500,
    max_iter = 500) {
    coordinates <- coordinates %>% filter(origin == 1)
    #-------------------------------------------------------------------------#
    # intialise background grid, active grod and get nearest neighbors
    # we assume that you will have at least 4 inital points
    # This needs to be changed later - I dont think using knn is the best 
    # option here
    #-------------------------------------------------------------------------#
    knn <- RANN::nn2(coordinates[, c("x", "y")],
      k = floor(nrow(coordinates) * 0.2))
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
            warning("Max Iteration with no convergence!
            Returning Approximation")
        }
        #---------------------------------------------------------------------#
        # we change the radius based on how many points there 
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
    
    return(background_grid)
}


hex_grid <- function(coordinates, n_centers, return_index = TRUE) {
    dims <- ceiling(sqrt(n_centers))
    range_y <- max(coordinates$y)
    range_x <- max(coordinates$x)
    y <- seq(0, range_y + 1, l = dims)
    x_short <- seq(0, range_x, l = dims)
    shift <- (x_short[2] - x_short[1]) / 2
    x_long <- seq(-shift, max(x_short) + shift, l = dims + 1)
    c_x <- c()
    c_y <- c()
    for (i in seq_along(y)){
        if (i %% 2 == 0) {
            c_x <- c(c_x, x_long)
            c_y <- c(c_y, rep(y[i], times = length(x_long)))
        } else {
            c_x <- c(c_x, x_short)
            c_y <- c(c_y, rep(y[i], times = length(x_short)))
        }
    }
    coord <- data.frame("x" = c_x, "y" = c_y)
    if (return_index) {
      background_grid <- RANN::nn2(data = coordinates[, c("x", "y")],
      query = coord,
      k = 1)
      return(background_grid$nn.idx[, 1])
    } else {
      return(coord)
    }
    
}
