#' compute fast pearson correlation between matrices
#' @param seed matrix/sparse matrix of seed cells
#' @param query matrix/sparse matrix of query cells
#' @return correlation matrix between all cells in seed and query
pearson_approx <- function(seed, query) { 
  n <- nrow(seed)
  sums <- outer(colSums(query), colSums(seed)) 
  stds <- outer(apply(query, 2, sd), apply(seed, 2, sd)) 
  correlation <- (t(query) %*% seed - sums / n) / stds / n
  return(correlation)
}



#' convert niche composition into matrix format for jaccard compute
#' @param niches list of niches with the cell composition of each niche
#' @return matrix with columns as cells and rows as cell types
make_composition_matrix <- function(niches){
    max_elements <- max(lengths(niches))
    barcodes <- names(niches)
    # Add NA padding 
    niches <- lapply(niches, function(niches, max_elements) {
        niches <- c(niches, rep(NA, times = length(niches) - max_elements))
        return(niches)
    }, max_elements = max_elements)
    niches <- do.call("cbind", niches)
    colnames(niches) <- barcodes
    return(niches)
}