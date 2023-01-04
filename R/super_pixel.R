################################################################################
################################   Vesalius      ###############################
################################################################################

#------------------------------/Super pixels/----------------------------------#

#' select inital indices
#' @param coordinates data frame containing spatial coordinates
#' @param k numeric - number of intial coordinates
#' @param random logial defining if initial starting points should be
#' selected at random - if false - approximate number
#' @return barcodes of starting coordinates
#' importFrom future.apply, future_lapply
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
    x2 <- mean(tail(x, length(x) - 1) - head(x, length(x) - 1))
    y <- seq(min(coordinates$y),
        max(coordinates$y),
        length.out = (2 * floor(box_size)) - 1)
    y <- y[seq(2, length(y), by = 2)]
    y2 <- mean(tail(y, length(y) - 1) - head(y, length(y) - 1))
    grid <- cbind(rep(x, each = length(y)),
        rep(y, times = length(x)))
    grid <- future_apply(grid,
        MARGIN = 1,
        closet_grid_point,
        coordinates,
        x2,
        y2)
    indices <- do.call("rbind", grid)
    return(indices)
}

#' closest grid point
#' @param grid grid matrix containing coordinate paiur of grid point centers
#' @param coordinates data frame with spatial indices
#' @return closest coordinates to grid mid point
closet_grid_point <- function(grid, coordinates, x2, y2) {
    dist <- sqrt(((coordinates$x - grid[1])^2 + (coordinates$y - grid[2])^2))
    coord <- coordinates[order(dist, decreasing = FALSE), ]
    coord$dist <- dist
    return(coord)
}

