#' segment image
#'
#' segment vesalius images to find initial territories
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions numeric vector of latent space dimensions to use.
#' @param embedding character string describing which embedding should
#' be used.
#' @param method character string for which method should be used for
#' segmentation. Select from "kmeans", "louvain", "leiden", "slic", 
#' "leiden_slic","louvain_slic","som"
#' @param col_resolution numeric colour resolution used for segmentation. 
#' (see details)
#' @param compactness numeric - factor defining super pixel compaction.
#' @param scaling numeric - scaling image ration during super pixel 
#' segmentation.
#' @param threshold numeric [0,1] - correlation threshold between
#' nearest neighbors when generating segments from super pixels.
#' @param index_selection character - method for selecting initial indices
#' @param verbose logical - progress message output.
#' @details Applying image segmentation ensures a reduction in colour
#' complexity.
#'
#' Vesalius provides 7 different methods for clustering colours and
#' reducing color complexity: **Kmeans**, **Louvain**, **Leiden**,
#' **slic**, **leiden_slic**, **louvain_slic**, and **som**
#'
#' In the case of kmeans clustering the \code{col_resolution} argument
#' shows the number of colours that the images should be reduced to.
#' In this case, \code{col_resolution} should be an integer and
#' we suggest first looking at values between 3 and 20.
#' 
#' In the case of **leiden** and **louvain** clustering, the
#' \code{col_resolution} is the granularity of the clustering.
#' In this case, we suggest using values between 0.01 and 1 to start with.
#' We recommned uisng **louvain** clustering over **leiden** in
#' this context.
#' 
#' In the case of slic, the col_resolution define the number of starting
#' points used to generate super pixels. Depending on the number of
#' points there are in the assay, we suggested using 10% of the total 
#' number of points as starting point. 
#' For example, if you have  1000 spatial indices, you can set 
#' col_resolution to 100.
#'
#' The optimal \code{col_resolution} will depend on your interest and 
#' biological question at hand. You might be interested in more or less
#' granular territories. Along with smoothing, the number of segments is
#' one way to control this granularity.
#
#' @return a vesalius_assay object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple smoothing
#' ves <- smooth_image(ves, dimensions = seq(1, 30))
#' 
#' # quick segmentation
#' ves <- segment_image(ves, dimensions = seq(1, 30))
#'}
#' @export

segment_image <- function(vesalius_assay,
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
    "kmeans" = kmeans_segmentation(vesalius_assay,
      col_resolution = col_resolution,
      dimensions = dimensions,
      embedding = embedding,
      verbose = verbose),
    "leiden" = leiden_segmentation(vesalius_assay,
      dimensions = dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      verbose = verbose),
    "louvain" = louvain_segmentation(vesalius_assay,
      dimensions = dimensions,
      col_resolution = col_resolution,
      embedding = embedding,
      verbose = verbose),
    "som" = som_segmentation(vesalius_assay,
      dimensions= dimensions,
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
    data = segments$segments,
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

#' kmeans segmentation function
#' @param vesalius_assay a vesalius_assay
#' @param col_resolution integer or vector of positive integers.
#' Colour depth used for segmentation.
#' @inheritParams segment_image
#' @details Run an interaive kmeans segmentation with the possibility 
#' to run multiple rounds of smoothing. 
#' @return list containing segmented image as an active embedding and
#' territory cluster for all barcodes.
#' @importFrom dplyr right_join
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom stats kmeans
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom imager as.cimg
#' @importFrom utils tail
kmeans_segmentation <- function(vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 10,
  embedding = "last",
  verbose = TRUE) {
  coord <- get_tiles(vesalius_assay) %>%
    filter(origin == 1)
  embeddings <- check_embedding_selection(vesalius_assay,
    embedding,
    dimensions)

  #--------------------------------------------------------------------------#
  # Now lets cluster colours
  # Remember that here colours are your new active embedding values
  #--------------------------------------------------------------------------#
  km <- suppressWarnings(kmeans(embeddings,
    col_resolution,
    iter.max = 10000,
    nstart = 10))
  clusters <- km$cluster
  kcenters <- km$centers
  match_loc <- match(coord$barcodes, names(clusters))
  clusters <- data.frame(coord, "Segment" = clusters[match_loc])
  active <- assign_centers(vesalius_assay, clusters, kcenters, dimensions)
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(list("segments" = clusters, "active" = active))
}


#' leiden segmentation
#'
#' using leiden clustering to cluster colors
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions embedding dimensions used for clustering
#' @param col_resolution clustering resolution used for leiden
#' @param embedding embedding type used for clustering
#' @param verbose logical if progress message should outputed
#' @returns list with updated segmented embedding values 
#' and segment territories.
#' @importFrom igraph cluster_leiden
#' @importFrom dplyr %>%
leiden_segmentation <- function(vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 0.01,
  embedding = "last",
  verbose = TRUE) {
  coord <- get_tiles(vesalius_assay) %>%
    filter(origin == 1)
  embeddings <- check_embedding_selection(vesalius_assay,
    embedding,
    dimensions)[,dimensions]
  
  loc <- match(coord$barcodes, rownames(embeddings))
  embeddings <- cbind(coord[, c("x", "y")], embeddings[loc, ])
  graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
  clusters <- igraph::cluster_leiden(graph,
    resolution = col_resolution)
  cluster <- data.frame("cluster" = clusters$membership,
    "barcodes" = clusters$names)
  match_loc <- match(coord$barcodes, cluster$barcodes)
  clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
  active <- create_pseudo_centroids(vesalius_assay,
    clusters,
    dimensions)
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(list("segments" = cluster, "active" = active))
}



#' louvain segmentation
#'
#' using leiden clustering to cluster colors
#' @param vesalius_assay a vesalius_assay object
#' @param dimensions embedding dimensions used for clustering
#' @param col_resolution clustering resolution used for leiden
#' @param embedding embedding type used for clustering
#' @param verbose logical if progress message should outputed
#' @returns list with updated segmented embedding values 
#' and segment territories.
#' @importFrom  igraph cluster_louvain
#' @importFrom dplyr %>%
louvain_segmentation <- function(vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 0.01,
  embedding = "last",
  verbose = TRUE) {
  coord <- get_tiles(vesalius_assay) %>%
    filter(origin == 1)
  embeddings <- check_embedding_selection(vesalius_assay,
    embedding,
    dimensions)[, dimensions]
  loc <- match(coord$barcodes, rownames(embeddings))
  embeddings <- cbind(coord[, c("x", "y")], embeddings[loc, ])
  graph <- compute_nearest_neighbor_graph(embeddings = embeddings)
  clusters <- igraph::cluster_louvain(graph, resolution = col_resolution)
  cluster <- data.frame("cluster" = clusters$membership,
    "barcodes" = clusters$names)


  match_loc <- match(coord$barcodes, cluster$barcodes)
  clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
  active <- create_pseudo_centroids(vesalius_assay,
    clusters,
    dimensions)
  new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
    "Segment") %>%
    tail(1)
  colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
  return(list("segments" = cluster, "active" = active))
}



#' @importFrom kohonen som somgrid
som_segmentation <- function(vesalius_assay,
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
    bar <- rownames(embeddings)
    embeddings$barcodes <- bar
    embeddings <- embeddings %>%
        right_join(coord, by = "barcodes") %>%
        select(-c("barcodes", "origin")) %>%
        as.matrix()
    colnames(embeddings) <- c(paste0("dim_", dimensions), "x", "y")
    #-------------------------------------------------------------------------#
    # create SOM map of data 
    #-------------------------------------------------------------------------#
    som <- som(as.matrix(embeddings),
      grid = somgrid(xdim = floor(sqrt(col_resolution)),
        ydim = ceiling(sqrt(col_resolution))))
    
    cluster <- data.frame("cluster" = som$unit.classif,
      "barcodes" = bar)
    match_loc <- match(coord$barcodes, cluster$barcodes)
    clusters <- data.frame(coord, "Segment" = cluster$cluster[match_loc])
    active <- create_pseudo_centroids(vesalius_assay,
      clusters,
      dimensions)
    active <- active / ratio
    
    new_trial <- create_trial_tag(colnames(vesalius_assay@territories),
      "Segment") %>%
      tail(1)
    colnames(clusters) <- gsub("Segment", new_trial, colnames(clusters))
    return(list("segments" = clusters, "active" = active))

}



#' assign centroid values to active embedding 
#' @param vesalius_assay a vesalius assy object 
#' @param clusters data.frame containing cluster values
#' @param kcenters matrix containing centroid values for each dimension
#' @param dimensions vector (nummeric / int) describin which latent space
#' @param ratio if used in the context of super pixel - spatial ration
#' dimensiuons shouls be used. 
#' @return matrix for the active embedding usiong color segementation 

assign_centers <- function(vesalius_assay,
  clusters,
  kcenters,
  dimensions,
  ratio = NULL) {
  #sp <- map(1:spectrum(images),~ km$centers[km$cluster,2+.]) %>% do.call(c,.) %>% as.cimg(dim=dim(images))
  active <- vesalius_assay@active
  for (d in seq_along(dimensions)) {
      tmp <- kcenters[clusters$Segment, dimensions[d]]
      if (!is.null(ratio)) {
        tmp <- tmp / ratio 
      }
      active[, dimensions[d]] <- tmp
  }
  return(active)
}


#' compute and greate nearest neighbor graph
#' @param embeddings embedding matrix
#' @param k numeric describing number of nearest neighbors
#' @return igraph object
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_data_frame
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
compute_nearest_neighbor_graph <- function(embeddings, k = 20) {
    knn <- RANN::nn2(embeddings, k = k + 1)$nn.idx
    rownames(knn) <- rownames(embeddings)
    graph <- populate_graph(knn)
    graph <- igraph::graph_from_data_frame(graph, directed = FALSE)
}



#' populate graph with network connections
#' @param chunk data frame with nearest neighbors
#' @returns data frame with network connections
populate_graph <- function(chunk) {
  barcodes <- rownames(chunk)
  nxk <- ncol(chunk)
  template <- data.frame("e1" = rep(barcodes, each = nxk),
    "e2" = barcodes[c(t(chunk))])
  return(template)
}


#' create centroid vlues for louvain and leiden 
#' @param vesalius_assay  vesalius assay object
#' @param clusters data frame containing clusters/segments 
#' @param dimensions vector (nummeric / int) describin which latent space
#' dimensiuons shouls be used.
#' @return matrix for the active embedding usiong color segementation 
create_pseudo_centroids <- function(vesalius_assay, clusters, dimensions) {
  active <- vesalius_assay@active
  for (d in dimensions) {
    for (clust in unique(clusters$Segment)) {
        loc <- match(clusters$barcodes[clusters$Segment == clust],
          rownames(active))
        active[loc, d] <- mean(active[loc, d])
    }
  }
  return(active)
}