################################################################################
################################   Vesalius      ###############################
################################################################################

#------------------------------/Super pixels/----------------------------------#

#' select inital indices
#' @param coordinates data frame containing spatial coordinates
#' @param k numeric - number of intial coordinates
#' @return barcodes of starting coordinates
#' @importFrom future.apply future_apply
select_initial_indices <- function(coordinates,
    k = 500) {
    #---------------------------------------------------------------------#
    # First we get mid point of each box 
    #---------------------------------------------------------------------#
    box_size <- sqrt(k)
    x <- seq(min(coordinates$x),
        max(coordinates$x),
        length.out = (2 * floor(box_size)) - 1)
    x <- x[seq(2, length(x), by = 2)]
    y <- seq(min(coordinates$y),
        max(coordinates$y),
        length.out = (2 * floor(box_size)) - 1)
    y <- y[seq(2, length(y), by = 2)]
    grid <- cbind(rep(x, each = length(y)),
        rep(y, times = length(x)),
        seq(1, length(y) * length(x)))
    grid <- future_apply(coordinates,
        MARGIN = 1,
        closet_grid_point,
        grid)
    
    grid <- do.call("rbind", grid)
    return(grid)
}

#' closest grid point
#' @param grid grid matrix containing coordinate paiur of grid point centers
#' @param coordinates data frame with spatial indices
#' @return closest coordinates to grid mid point
closet_grid_point <- function(coordinates, grid) {
    distance <- sqrt(((grid[, 1] + as.numeric(coordinates[2]))^2 +
        (grid[, 2] + as.numeric(coordinates[3]))^2))
    distance <- which(distance == min(distance))
    coordinates <- data.frame("barcodes" = coordinates[1],
        "xcoord" = coordinates[2],
        "ycoord" = coordinates[3],
        "x_cent" = grid[distance, 1],
        "y_cent" = grid[distance, 2],
        "center" = grid[distance, 3])
    return(coordinates)
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


#' @importFrom imager imappend imsplit nPix spectrum
#' @importFrom purrr map_dbl map
slic_segmentation <- function(vesalius_assay,
    dimensions,
    col_resolution,
    embedding,
    compactness = 1,
    scaling = 0.5,
    k = 6,
    threshold = 0.9,
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
    coord <- get_tiles(vesalius_assay) %>% filter(origin == 1)
    sc_spat <- max(dim(images)[1:2]  * scaling)
    sc_col <-  imsplit(images, "cc") %>% map_dbl(sd) %>% max()

    #Scaling ratio for pixel values
    ratio <- (sc_spat / sc_col) / (compactness * 10)
    #browser()
    embeddings <- as.data.frame(images * ratio, wide = "c") %>% as.matrix
    #Generate initial centers from a grid
    # This will need to be updated to make sure 
    # we only take center within the data set
    index <- round(seq(1, nPix(images) / spectrum(images), l = col_resolution))
    #Run k-means - you should have read the mnual pages you silly goose
    km <- suppressWarnings(kmeans(embeddings,
        embeddings[index, ],
        iter.max = 100,
        nstart = 10))
    embeddings <- map(1:spectrum(images), ~ km$centers[km$cluster, 2 + .]) %>%
        do.call(c, .) %>%
        as.cimg(dim = dim(images))
    embeddings <- embeddings / ratio
    #clusters <- select_similar(embeddings, coordinates = coord)

    #browser()
    embeddings <- format_c_to_ves(imsplit(embeddings, "cc"),
      vesalius_assay,
      dimensions,
      embed = embedding,
      verbose = FALSE)

    clusters <- as.cimg(km$cluster, dim = c(dim(images)[1:2], 1, 1)) %>%
        as.data.frame() %>%
        right_join(coord, by = c("x", "y"))
    match_loc <- match(coord$barcodes, clusters$barcodes)
    clusters <- data.frame(coord, "Segment" = clusters[match_loc, "value"])
    clusters <- connected_pixels(clusters, embeddings, k, threshold, verbose)
    

    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
        "Segment") %>%
    tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = embeddings))
}
