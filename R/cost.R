#' compute fast pearson correlation between matrices
#' @param seed (r x c) seed matrix
#' @param query (r x c') query matrix
#' @return correlation matrix between each columns (c' x c)
pearson_approx <- function(seed, query) { 
  n <- nrow(seed) 
  sums <- outer(colSums(query), colSums(seed)) 
  stds <- outer(apply(query, 2, sd), apply(seed, 2, sd)) 
  correlation <- (t(query) %*% seed - sums / n) / stds / n
  return(correlation)
}



#' create a composition matrix from cell type labels in niches
#' @param niches list of cell labels for each niche in a data set
#' @details The aim is to reformat the niches to make the scoring faster
#' @return a matrix with rows as cell labels and columns as each cell
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