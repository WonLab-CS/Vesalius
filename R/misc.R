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

# z_norm <- function(x) {
#   if (length(table(x)) == 1) {
#     return(x)
#     warning("Cannot minmax normalise - all values are equal!")
#   } else {
#     return((x - mean(x)) / sd(x))
#   }
# }





chunk <- function(x, n, l = NULL) {
    chunk <- mapply(function(a, b) {
      return(x[a:b])},
    seq.int(from = 1, to = length(x), by = n),
    pmin(seq.int(from = 1, to = length(x), by = n) + (n - 1), length(x)),
    SIMPLIFY = FALSE)
    if (!is.null(l) && length(chunk) > l){
        chunk[[l]] <- unlist(chunk[l:length(chunk)])
        chunk <- chunk[seq(1,l)]
    }
    return(chunk)
}




#-------------------------/ Aligning Assays /--------------------------------#

graph_from_voronoi <- function(centers) {
    voronoi <- deldir::deldir(x = as.numeric(centers$x),
        y = as.numeric(centers$y))$delsgs
    center <- seq_len(nrow(centers))
    graph <- lapply(center, function(idx, voronoi){
        tri <- voronoi %>% filter(ind2 == idx)
        tri <- c(tri$ind1, idx)
        graph <- data.frame("from" = rep(idx, length(tri)),
            "to" = tri)
        return(graph)
    }, voronoi = voronoi) %>%
    do.call("rbind", .)
    return(graph)
}

#'@importFrom igraph graph_from_data_frame distances
graph_path_length <- function(graph) {
    gr <- igraph::graph_from_data_frame(graph, directed = FALSE)
    path_length <- igraph::distances(gr)
    return(path_length)
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
