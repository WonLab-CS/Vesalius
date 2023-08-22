################################################################################
############################   ST package        ###############################
################################################################################



#' globalise territories
#' converts territories number from territory values per segment
#' to territory values overall the entire spatial omics assay.
#' @param img a vesalius territory data frame containing segments
#' and new territory trial
#' @return the same data frame but with adjusted terriroty values
#' @importFrom dplyr filter %>%
globalise_territories <- function(img) {
  img_tmp <- img %>% filter(trial != "isolated")
  ter <- paste0(img_tmp$Segment, "_", img_tmp$trial)
  all_ter <- unique(ter)
  ter <- seq_along(all_ter)[match(ter, all_ter)]
  img$trial[img$trial != "isolated"] <- ter
  return(img)
}



#' detect_edges
#' detect territory edges in black and white images with sobel filter
#' @param img a cimg image
#' @return a pix set containing deteced edges
#' @importFrom dplyr %>%
#' @importFrom imager imgradient enorm add
detect_edges <- function(img) {
  img <- img %>%
    imager::imgradient("xy") %>%
    imager::enorm() %>%
    imager::add() %>%
    sqrt()
  return(img)
}





#------------------------/ Normalising Embeds /--------------------------------#

#' pixel normalisation dispatch function
#' @param embeds a embedding vector 
#' @param type string "minmax" or "quantileNorm"
#' @details how pixels should be normalised 
#' At the moment only miman is used. Quantile needs to be tested.
norm_pixel <- function(embeds, type = c("minmax", "quantile_norm", "z_norm")) {
    #--------------------------------------------------------------------------#
    # Normalise pixels values
    # at the moment i am only using min max but just in case
    # context: Dr Hong suggest imnvestigating the influence of dat normalisation
    #--------------------------------------------------------------------------#
    embeds <- switch(type[1L],
                     "minmax" = min_max(embeds),
                     "quantile_norm" = quantile_norm(embeds),
                     "z_norm" = z_norm(embeds))
    return(embeds)
}


#' min max normalisation
#' @param x numeric vector
#' @return min max nornalised vector
min_max <- function(x) {
  if (length(table(x)) == 1) {
    return(x)
    warning("Cannot minmax normalise - all values are equal!")
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}

z_norm <- function(x) {
  if (length(table(x)) == 1) {
    return(x)
    warning("Cannot minmax normalise - all values are equal!")
  } else {
    return((x - mean(x)) / sd(x))
  }
}

scale_data_spatial <- function(data, compactness, scale) {
  scale_spat <- max(c(max(data$x), max(data$x))  * scale)
  scale_dat <- apply(data[, !colnames(data) %in% c("x", "y")], 2, sd) %>%
    max()
  ratio <- (scale_spat / scale_dat) / compactness
  data[, !colnames(data) %in% c("x", "y")] <-
    data[, !colnames(data) %in% c("x", "y")] * ratio
  return(data)
}

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <- length(union(a, b))
    return(intersection / union)
}

arrange_knn_matrix <- function(knn) {
   for (i in seq_len(nrow(knn$nn.idx))){
      knn$nn.dists[i, ] <- knn$nn.dists[i, order(knn$nn.idx[i, ])]
   }
   return(knn$nn.dists)
}



chunk <- function(x, n) {
  chunk <- mapply(function(a, b) {
      return(x[a:b])},
    seq.int(from = 1, to = length(x), by = n),
    pmin(seq.int(from = 1, to = length(x), by = n) + (n - 1), length(x)),
    SIMPLIFY = FALSE)
  return(chunk)
}


#-------------------------/ Aligning Assays /--------------------------------#

polar_angle <- function(coord_x, coord_y, center_x, center_y) {
  x <- coord_x - center_x
  y <- coord_y - center_y
  angle <- mapply(function(x, y) {
    if (x >= 0 & y >= 0) angle <- atan(abs(y) / abs(x)) * (180 / pi)
    if (x < 0 & y >= 0) angle <- 180 - (atan(abs(y) / abs(x)) * (180 / pi))
    if (x < 0 & y < 0) angle <- 180 + (atan(abs(y) / abs(x)) * (180 / pi))
    if (x >= 0 & y < 0) angle <- 360 - (atan(abs(y) / abs(x)) * (180 / pi))
    return(angle)
  }, x = x, y = y, SIMPLIFY = TRUE)
  if (is.na(angle)) {
      angle <- 0
  }
  return(angle)
}



cart_coord <- function(d, a) {
  if (a < pi / 2) {
    delta_x <- d * cos(a)
    delta_y <- d * sin(a)
  } else if (a < pi) {
    delta_x <- -d * sin(a - pi/ 2)
    delta_y <- d * cos(a - pi/ 2)
  } else if (a < 3 * pi / 2) {
    delta_x <- -d * cos(a - pi)
    delta_y <- -d * sin(a - pi)
  } else {
    delta_x <- d * sin(a - 3 * pi / 2)
    delta_y <- -d * cos(a - 3 * pi / 2)
  }
  return(list("delta_x" = delta_x, "delta_y" = delta_y))
}

#-----------------------------/ Scaling  /----------------------------------#
#' calculate scale of assay
#' @param coordinates Spatial coordinates as data frame
#' @details Calculate the average distance between spots/beads/indeces
#' @return Single float 
#' @importFrom RANN nn2
calculate_scale <- function(coordinates, q = 0.999) {
    scale <- RANN::nn2(data = coordinates[, c("x", "y")], k = 2)
    scale <- unname(quantile(scale$nn.dist[, 2], q))
    return(scale)
}

#' Calculate the total area of a territory
#' @param vesalius_assay a vesalius_assay object
#' @param territory numeric/character/vector describing wich terriory(ies)
#' that will be use to compute area.
#' @param trial which territory trial should be used
#' @param use_rescaled logical - which value should be use for 
#' @details Compute the total area of a territory 
#' @return territory area
compute_territory_area <- function(vesalius_assay,
  territory = NULL,
  trial = "last",
  use_rescaled = FALSE,
  verbose = TRUE) {
  #---------------------------------------------------------------------------#
  # First lets check the territory selection
  #--------------------------------------------------------------------------#
  territory <- territory %||%
      stop("No specified territory! Cannot compute Area.")
  ter <- check_territory_trial(vesalius_assay, trial)
  territory <- check_group_value(ter, territory)
  ter <- filter(ter, trial %in% territory)
  #--------------------------------------------------------------------------#
  # Next we can get the scale and use this as the capture radius
  #--------------------------------------------------------------------------#
  if (use_rescaled) {
    scale <- vesalius_assay@meta$scale$rescale
  } else {
    scale <- vesalius_assay@meta$scale$scale
  }
  #--------------------------------------------------------------------------#
  # run distance pooling to seperate potential patches
  #--------------------------------------------------------------------------#
  patches <- distance_pooling(ter,
    capture_radius = scale,
    min_spatial_index = 3)
  ter$trial <- patches
  patches <- unique(ter$trial)
  patches <- patches[!patches %in% "isolated"]
  area <- rep(0,length(patches))
  message_switch("area_comp", verbose, patches = length(patches))
  for (i in seq_along(patches)) {
      area[i] <- deldir::deldir(
        as.numeric(ter$x[ter$trial == patches[i]]),
        as.numeric(ter$y[ter$trial == patches[i]]))$del.area
  }
  area <- sum(area)
  return(area)
}