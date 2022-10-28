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
  ter <- paste0(img_tmp$segment, "_", img_tmp$trial)
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
norm_pixel <- function(embeds, type = c("minmax", "quantileNorm")) {
    #--------------------------------------------------------------------------#
    # Normalise pixels values
    # at the moment i am only using min max but just in case
    # context: Dr Hong suggest imnvestigating the influence of dat normalisation
    #--------------------------------------------------------------------------#
    embeds <- switch(type[1L],
                     "minmax" = min_max(embeds),
                     "quantileNorm" = quantile_norm(embeds))
    return(embeds)
}


#' min max normalisation
#' @param x numeric vector
#' @return min max nornalised vector
min_max <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

#' quantile normalisation 
#' @param x numeric vector
#' @return quatile normalised vector 
quantile_norm <- function(x) {
  x_rank <- apply(x, 2, rank, ties.method = "min")
  x <- data.frame(apply(x, 2, sort))
  x_mean <- apply(x, 1, mean)
  x_final <- apply(x_rank, 2, function(idx, m) {
    return(m[idx])
  }, x_mean)
  rownames(x_final) <- rownames(x)
  return(x_final)
}
