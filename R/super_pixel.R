################################################################################
################################   Vesalius      ###############################
################################################################################

#------------------------------/Super pixels/----------------------------------#

#' select inital indices
#' @param coordinates data frame containing spatial coordinates
#' @param k numeric - number of intial coordinates
#' @return barcodes of starting coordinates
#' @importFrom future.apply future_apply
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
   
    y <- seq(min(coordinates$y),
        max(coordinates$y),
        length.out = (2 * floor(box_size)) - 1)
    y <- y[seq(2, length(y), by = 2)]
    grid <- cbind(rep(x, each = length(y)),
        rep(y, times = length(x)),
        seq(1, length(y) * length(x)))
    grid <- future_apply(coordinates,
        MARGIN = 1,
        closet_grid_point,
        grid)
    
    grid <- do.call("rbind", grid)
    return(grid)
}

#' closest grid point
#' @param grid grid matrix containing coordinate paiur of grid point centers
#' @param coordinates data frame with spatial indices
#' @return closest coordinates to grid mid point
closet_grid_point <- function(coordinates, grid) {
    distance <- sqrt(((grid[, 1] + as.numeric(coordinates[2]))^2 +
        (grid[, 2] + as.numeric(coordinates[3]))^2))
    distance <- which(distance == min(distance))
    coordinates <- data.frame("barcodes" = coordinates[1],
        "xcoord" = coordinates[2],
        "ycoord" = coordinates[3],
        "x_cent" = grid[distance, 1],
        "y_cent" = grid[distance, 2],
        "center" = grid[distance, 3])
    return(coordinates)
}