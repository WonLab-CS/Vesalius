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