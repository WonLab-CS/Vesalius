################################################################################
############################   ST package        ###############################
################################################################################

# Misc functions function that may or may not be dropped in the future.
# Who knows? I certainly don't. But you designed and wrote this packages.
# That means diddly squat!! Since when do programmers know what they are doing?

#' subSetTerritories gets all barcodes associated with a territory from a
#' Seurat object.
#' @param territories vector of barcode values (generally character strings)
#' @param seurat Seurat object containing all barcode values from ST assay
#' @details Essentially a wrapper function to the Seurat \code{subset} function.
#' The function serves mainly the purpose of a place holder function for
#' future iterations of Vesalius.
#' @return Seurat object only containing the desired barcodes.
#' @examples
#' \dontrun{

#' }

subSetTerritories <- function(territories, seurat) {
  #--------------------------------------------------------------------------#
  # Simplified version for now
  # It might be worth while getting away from seurat later
  # essentially this is a template function
  # LEGACY - might drop after refactoring
  #--------------------------------------------------------------------------#
  seurat <- subset(seurat, cells = territories)
  return(seurat)
}


#' getSeuratCoordinates get barcode coordinates from a seurat object
#' @param seurat Seurat object containing all barcode values from ST assay
#' @details Essentially a wrapper function to the Seurat
#' \code{GetTissueCoordinates} function.
#' Mainly serve as a way if generalising output format.
#' @return Seurat object only containing the desired barcodes.
#' @examples
#' \dontrun{
#' data(Vesalius)
#' }
getSeuratCoordinates <- function(seurat) {
  # Legacy Code - might drop after refactoring!
  ret <- GetTissueCoordinates(seurat)
  if(sum(colnames(ret) %in% c("imagerow", "imagecol")) == 2){
      colnames(ret) <- c("y", "x")
  }
  return(ret)
}



# Used to convert territories per cluster to territories across the whole
# ST array
globalise_territories <- function(img) {
  img_tmp <- img %>% filter(trial != "isolated")
  ter <- paste0(img_tmp$segment, "_", img_tmp$trial)
  all_ter <- unique(ter)
  ter <- seq_along(all_ter)[match(ter, all_ter)]
  img$trial[img$trial != "isolated"] <- ter
  return(img)
}

#---------------------------/Edge Functions/-----------------------------#
# This might need to be moved somewhere else
# Also these functions will need to be cleaned and updated
# can still be used - need to think about what i am going to do with this.
select_similar <- function(img, cols, segment, threshold = "auto") {
  img <- img %>% 
    { . - imfill(dim = dim(img), val = cols[segment, 2:4])} %>%
    imsplit("c") %>%
    enorm()
  img <- !threshold(img, threshold)
  img <- as.cimg(img)
  return(img)
}
detect_edges <- function(img) {
  img <- img %>% 
    imgradient("xy") %>%
    enorm() %>%
    add() %>%
    sqrt()
  return(img)
}

pmap <- function(edges, method = c("inverse", "none")) {
    pmap <- switch(method[1],
      "inverse" = 1 / (1 + edges),
      "none" = edges)
    return(pmap)
}


#------------------------/ Normalising Embeds /--------------------------------#
norm_pixel <- function(embeds, type = c("minmax", "quantileNorm")) {
    #--------------------------------------------------------------------------#
    # Normalise pixels values
    #--------------------------------------------------------------------------#
    embeds <- switch(type[1L],
                     "mixmax" = min_max(embeds),
                     "quantileNorm" = quantile_norm(embeds))
}


min_max <- function(x){
    return((x - min(x)) / (max(x) - min(x)))
}


quantile_norm <- function(x){
  x_rank <- apply(x, 2, rank, ties.method = "min")
  x <- data.frame(apply(x, 2, sort))
  x_mean <- apply(x, 1, mean)
  x_final <- apply(x_rank, 2, function(idx, m) {
    return(m[idx])
  }, x_mean)
  rownames(x_final) <- rownames(x)
  return(x_final)
}
