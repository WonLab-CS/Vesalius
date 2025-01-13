#' compute fast pearson correlation between matrices
#'
pearson_approx <- function(seed, query) { 
  n <- nrow(seed) 
  sums <- outer(colSums(query), colSums(seed)) 
  stds <- outer(apply(query, 2, sd), apply(seed, 2, sd)) 
  correlation <- (t(query) %*% seed - sums / n) / stds / n
  return(correlation)
}




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