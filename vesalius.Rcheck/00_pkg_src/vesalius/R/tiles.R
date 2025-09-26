
#' generate tiles
#'
#' generate pixel tiles from punctual spatial assay coordinates
#' @param vesalius_assay a vesalius_assay object
#' @param tensor_resolution numeric (range 0 - 1) describing the compression
#' ratio to be applied to the final image. Default = 1
#' @param filter_grid numeric (range 0 - 1) size of the grid used when filtering
#' outlier beads. Defined as a proportion of total image size. Default = 0.1
#' @param filter_threshold numeric (range 0 -1) describing the quantile
#' @param verbose logical describing if progress message should be outputed
#' @details This function converts punctual coordinates into pixel tiles.
#' The first step is to filter outlier beads. Essentially removing any bead 
#' that lies outside of the main assay (seed \code{\link{filter_grid}}). If 
#' reqested the resolution of the image can be reduced by redudcing the 
#' relative distance between each coordinate pair 
#' (see \code{\link{reduce_tensor_resolution}}). Note that if the 
#' resolution is reduce all spatial indices are retained and counts
#' will be merged together (see \code{\link{adjust_counts}}). 
#' 
#' To create tiles, each coordinates is expanded by using Voronoi 
#' tesselation and each tile is raterised i.e "filled with pixel" 
#' (see \code{\link{rasterise}}). Prior to rasterisation, we
#' remove all tiles that are above an area threshold. These tiles
#' are artefacts of the tesselation algortihm that creates a "box"
#' around all points and uses this as a boudary to generate each tile.
#' They can also be generate by stray beads or holes in the data. 
#' 
#' @returns a data.frame containing barcodes, x and you coordinates
#' of each pixel as well as the original x/y coordinates
#' @importFrom deldir deldir
#' @importFrom dplyr %>% filter
#' @export
generate_tiles <- function(vesalius_assay,
  tensor_resolution = 1,
  filter_grid = 0.01,
  filter_threshold = 0.995,
  verbose = TRUE) {
  env <- identical(.GlobalEnv, parent.frame()) && verbose
  simple_bar(env)
  #----------------------------------------------------------------------#
  # Checking if tiles have been computed 
  #----------------------------------------------------------------------#
  status <- any(search_log(vesalius_assay,
    "generate_tiles",
    return_assay = FALSE))
  
  #----------------------------------------------------------------------#
  # we do the tile check here - if we don't the tiles will be computed on
  # pixels and that will take ages and will lead to nonsense. 
  # I am not sure if holding the orginal data before filtering is worth it
  # While I figure that out, we just return object and skip tlling
  #----------------------------------------------------------------------#
  if (status) {
    warning("Tiles have already been computed!
      Returning vesalius_assay. For a new set of tiles,
      Create a new vesalius_assay
      then 'generate_tiles' or 'generate _embeddings'")
    return(vesalius_assay)
  }
  assay <- get_assay_names(vesalius_assay)
  coordinates <- get_tiles(vesalius_assay)
  message_switch("in_assay",
    verbose = verbose,
    assay = assay,
    comp_type = "Generating Tiles")
  vesalius_assay@meta$orig_coord <- coordinates
  #----------------------------------------------------------------------#
  # Filter outlier beads
  #----------------------------------------------------------------------#
  if (filter_grid > 0 && filter_grid < 1) {
    message_switch("distance_beads", verbose)
    coordinates <- filter_grid(coordinates = coordinates,
      filter_grid = filter_grid)
  }
  #----------------------------------------------------------------------#
  # Reduce resolution
  # Might be worth creating dedicated sanity check functions for this
  # rescaling scaling with scaling factor
  # TO DO add unit check to change unit if needed?
  #----------------------------------------------------------------------#
  if (tensor_resolution < 1 && tensor_resolution > 0) {
    message_switch("tensor_res", verbose)
    coordinates <- reduce_tensor_resolution(coordinates = coordinates,
        tensor_resolution = tensor_resolution)
    rescale <- calculate_scale(coordinates)
    scale_factor <- rescale / vesalius_assay@meta$scale$scale
    vesalius_assay@meta$scale <- list("scale" = vesalius_assay@meta$scale$scale,
      "rescale" = rescale,
      "scale_factor" = scale_factor)
  } 
  #----------------------------------------------------------------------#
  # TESSELATION TIME!
  #----------------------------------------------------------------------#
  message_switch("tess", verbose)
  tesselation <- deldir::deldir(x = as.numeric(coordinates$x),
    y = as.numeric(coordinates$y))
  #----------------------------------------------------------------------#
  # Filtering tiles
  #----------------------------------------------------------------------#
  message_switch("filter_tiles", verbose)
  filtered <- filter_tiles(tesselation, coordinates, filter_threshold)

  #----------------------------------------------------------------------#
  # Resterise tiles
  #----------------------------------------------------------------------#
  message_switch("raster", verbose)
  tiles <- rasterise(filtered)
  #------------------------------------------------------------------------#
  # adjusted counts if necessary
  # essentially merging counts when barcodes overlap
  #------------------------------------------------------------------------#
  message_switch("adj_counts", verbose)
  counts <- get_counts(vesalius_assay, type = "all")
  if (length(counts) > 0) {
    for (i in seq_along(counts)) {
      counts[[i]] <- adjust_counts(tiles,
        counts[[i]],
        throw = FALSE,
        verbose = FALSE)
    }
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = counts,
      slot = "counts",
      append = FALSE)
  }
  #----------------------------------------------------------------------#
  # check if anything has already been added to territories slot
  # happens if you use add_cells - then update with new coord 
  # Add new tiles as well 
  #----------------------------------------------------------------------#
  if (!all(dim(slot(vesalius_assay, "territories")) == c(0, 0))){
    territory_coord <- filter(tiles, origin == 1) %>%
        select(c("barcodes","x","y"))
    tmp <- vesalius_assay@territories
    territory_coord <- left_join(territory_coord, tmp, by = "barcodes")
    territory_coord <- territory_coord[, 
        !colnames(territory_coord) %in% c("x.y","y.y")]
    colnames(territory_coord) <- colnames(tmp)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
        data = territory_coord,
        slot = "territories",
        append = FALSE)
  }
  #----------------------------------------------------------------------#
  # check if there are already some embeddings and then filter
  # This can happen when integrating or when adding custom embeddings
  #----------------------------------------------------------------------#
  if (!all(dim(slot(vesalius_assay,"active")) == c(0, 0))) {
        # filter active first
        tmp_bar <- unique(tiles$barcodes)
        active <- vesalius_assay@active[rownames(vesalius_assay) %in% tmp_bar, ]
        embeds <- lapply(vesalius_assay@embeddings, function(x, bar) {
            return(x[rownames(x) %in% bar, ])
        }, bar = tmp_bar)
        vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
            data = active,
            slot = "active",
            append = FALSE)
        vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
            data = embeds,
            slot = "embeddings",
            append = FALSE)
  }
  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = tiles,
    slot = "tiles",
    append = FALSE)
  commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(generate_tiles))
  vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
  simple_bar(env)
  return(vesalius_assay)
}




#------------------------/ Filtering  Beads /----------------------------------#
#' filter grid
#' 
#' filtered stray barcodes/spatial indices.
#' @param coordinates data frame containing barcodes, / y coordinates of 
#' each barcode.
#' @param filter_grid numeric describing size of the grid to use as proporiton 
#' of array size. 
#' @details we create a grid over the data an check which grid section contain
#' a number of barcodes lower than a quantile threshold. Those barcodes will 
#' be removed. As it stands, I am not satisfied with this function. I think 
#' it is too restrictive and will most likely not be applicable to highly 
#' standardised assay. 
#' @returns a data.frame containing barcodes, x and y coordinates. 
filter_grid <- function(coordinates, filter_grid) {
  #----------------------------------------------------------------------------#
  # Essentially create a grid where each barcode is pooled into a grid space
  # If there are too little barcodes in that grid section then remove
  #----------------------------------------------------------------------------#
  grid_x <- round(coordinates$x * filter_grid)
  grid_y <- round(coordinates$y * filter_grid)
  grid_coord <- paste0(grid_x, "_", grid_y)
  grid <- table(grid_coord)
  grid <- grid[which(grid <= quantile(grid, 0.01))]
  grid_coord <- which(grid_coord %in% names(grid))
  coordinates <- coordinates[-grid_coord, ]
  return(coordinates)
}


#------------------------/ Reducing Resolution /-------------------------------#
### might want to adjust this and use knn instead?
### This is one approach but i guess that we could also use knn or 
### we could used reduced image size at a different point. The point is we 
### don't need full resolution images. It slows things down and for territory
### isolation very finer details are not retained after smoothing anyway. 
### The main issue here is how to account for changes of coordinates 
### once we have reduced the resolution.

#' reduce tensor resoltuon 
#' 
#' reduce the size of the image tensor by merging barcodes together after
#' compressing coordinates values
#' @param coordinates data frame with barcode coordinates 
#' @param tensor_resolution numeric (range 0 - 1) describing the compression
#' ratio to be applied to the final image. Default = 1
#' @details While we compress the coordinates, we retain all barcodes. 
#' Each barcode that overlap with each other are marged together. Their 
#' respective counts will also be merged together. This allows us to 
#' retain all barcodes for downstream analysis. 
#' @return a data frame with barcodes, x and coordinates
#' @importFrom dplyr %>% distinct
reduce_tensor_resolution <- function(coordinates,
  tensor_resolution = 1) {
  #----------------------------------------------------------------------------#
  # we will reduce number of points this way
  # this should keep all barcodes - with overlapping coordinates
  #----------------------------------------------------------------------------#
  coordinates$x <- round(coordinates$x * tensor_resolution) + 1
  coordinates$y <- round(coordinates$y * tensor_resolution) + 1
  #----------------------------------------------------------------------------#
  # Now we get coordinate tags - we use this to find all the merge locations
  # sorting and using rle to ensure that we actually merge them
  #----------------------------------------------------------------------------#
  tag <- paste0(coordinates$x, "_", coordinates$y)

  locs <- rle(sort(tag))
  locs <- check_tensor_compression(locs)
  dup <- locs$values[which(locs$length > 1)]
  #----------------------------------------------------------------------------#
  # creating new merged labels
  #----------------------------------------------------------------------------#
  dup_tags <- lapply(dup, function(dup, tag, barcodes) {
    tmp_locs <- which(tag == dup)
    barcodes <- paste0(barcodes[tmp_locs], sep = "_et_", collapse = "")
    barcodes <- rep(barcodes, times = length(tmp_locs))
    return(barcodes)
  }, tag = tag, barcodes = coordinates$barcodes)
  #----------------------------------------------------------------------------#
  # assigning new merged barcodes
  # for some reason match doesn't seem to work here - it only returns one of the
  # values
  #----------------------------------------------------------------------------#
  locs <- unlist(lapply(dup, function(dup, tag) {
      return(which(tag == dup))
  }, tag = tag))
  coordinates$barcodes[locs] <- unlist(dup_tags)
  coordinates <- coordinates %>% distinct(barcodes, .keep_all = TRUE)
  return(coordinates)
}






#------------------------/ Creating pixel tiles /------------------------------#

#' filter tiles
#' 
#' filter tiles based on tile area and if they share an edge with the 
#' tesselation box 
#' @param tesselation data.frame output from the deldir function
#' @param coordinates data frame with original coordinates
#' @param filter_threshold numeric describing the quantile threshold value
#' to use for area filtering
#' @details Here we want to filter based on the size of the tile
#' under the assumption that very large tiles are probably due to 
#' unpexted space. The issue is that if you don't apply this threshold,
#' subtle pixel patterns are lost wuhtin the sea of pixels. 
#' If the ueser really dont want to filter anything out then you don't
#' We will set a filter threshold pretty high though to make sure 
#' that it does get filter out as default.
#' @return a list with 2 data frame. 1 with filtered tesselation results
#' 2 filtered coordinate file.
#' @importFrom stats quantile
filter_tiles <- function(tesselation, coordinates, filter_threshold) {
  if (filter_threshold == 1) {
    coordinates$ind <- tesselation$ind.orig
    return(list("tess_v" = tesselation$dirsgs,
      "coordinates" = coordinates))
  } else {
    max_area <- quantile(tesselation$summary$dir.area, filter_threshold)
    idx <- which(tesselation$summary$dir.area >= max_area)
    tess_v <- tesselation$dirsgs
    points_to_keep <- tess_v$ind1 %in% idx | tess_v$ind2 %in% idx
    tess_v <- tess_v[!points_to_keep, ]
    coordinates <- coordinates[-idx, ]
    coordinates$ind <- tesselation$ind.orig[-idx]
    return(list("tess_v" = tess_v, "coordinates" = coordinates))
  }
}

#' rasterise tiles
#' 
#' fill tiles with pixel - rasterisation 
#' @param filtered data.frame with voronoi tile coordinates
#' @details Here, we take our tile cooridnates and fill them 
#' with pixels. Essentially, each voronoi tile can be discretised 
#' into a series of pixels and we achieve this by reconstructing a 
#' polygon from the tesselation coordinates and finding all discrete 
#' value that fall within this polygon. 
#' 
#' Note that the polygon coordinates need to be "convexified". 
#' Essentially, the order of the coordinates maters here so we
#' order the coordinates before recontructing the polygon (see
#' \code{\link{convexify}})
#' @return a data frame barcodes and their associated pixels.
#' @importFrom future.apply future_lapply
#' @importFrom sp point.in.polygon
rasterise <- function(filtered) {
    idx <- seq_len(nrow(filtered$coordinates))
    tiles <- future_lapply(idx, function(idx, filtered) {
        #----------------------------------------------------------------------#
        # get indecies from original data
        #----------------------------------------------------------------------#
        x_orig <- filtered$coordinates$x_orig[idx]
        y_orig <- filtered$coordinates$y_orig[idx]
        ind <- filtered$coordinates$ind[idx]
        indx <- filtered$coordinates$x[idx]
        indy <- filtered$coordinates$y[idx]
        tess_v <- filtered$tess_v %>% filter(ind1 == ind | ind2 == ind)
        if (nrow(tess_v) == 0) {
            return(NULL)
        }
        #----------------------------------------------------------------------#
        # create unique set of coordiantes that define tile boundaries
        #----------------------------------------------------------------------#
        coord <- paste0(c(tess_v$x1, tess_v$x2), "_", c(tess_v$y1, tess_v$y2))
        x <- as.numeric(sapply(strsplit(coord[!duplicated(coord)], "_"),
          "[[", 1)
          )
        y <- as.numeric(sapply(strsplit(coord[!duplicated(coord)], "_"),
          "[[", 2)
        )
        #----------------------------------------------------------------------#
        # Convert coordinates into a convex shape for cleaner
        # rasterising.
        #----------------------------------------------------------------------#
        convex <- convexify(x, y, indx, indy)
        x <- ceiling(convex$x)
        y <- ceiling(convex$y)
        #----------------------------------------------------------------------#
        # define max polygon containing all pixels
        #----------------------------------------------------------------------#
        lpx <- min(x) - 1
        hpx <- max(x) + 1
        lpy <- min(y) - 1
        hpy <- max(y) + 1

        max_polygon_x <- rep(seq(lpx, hpx), times = hpy - lpy + 1)
        max_polygon_y <- rep(seq(lpy, hpy), each = hpx - lpx + 1)

        #----------------------------------------------------------------------#
        # Fill triangles with all point in that space
        # And filter out negetive coordinates (associated with tesselation 
        # boundary extension) 
        #----------------------------------------------------------------------#
        cell <- sp::point.in.polygon(max_polygon_x, max_polygon_y, x, y)
        max_x <- max_polygon_x[cell %in% c(1, 2, 3) &
          max_polygon_x > 0 &
          max_polygon_y > 0]
        max_y <- max_polygon_y[cell %in% c(1, 2, 3) &
          max_polygon_y > 0 &
          max_polygon_x > 0]
        #----------------------------------------------------------------------#
        # get original point location 
        # in some cases this is not possible because of the strange shape
        # of the tile
        # in this case we select the mid point of the ordered coordinates
        #----------------------------------------------------------------------#
        cent <- which(max_x == round(indx) & max_y == round(indy))
        centers <- rep(0, length(max_x))
        orig_x <- rep(0, length(max_x))
        orig_y <- rep(0, length(max_x))
        if (length(centers) == 1) {
            cent <- 1
        } else if (length(cent) == 0) {
            cent <- as.integer(length(centers) / 2)
        }
        centers[cent] <- 1
        orig_x[cent] <- x_orig
        orig_y[cent] <- y_orig
        tile <- data.frame("barcodes" = rep(filtered$coordinates$barcodes[idx],
            times = length(max_x)),
          "x" = max_x,
          "y" = max_y,
          "origin" = centers)
        return(tile)
    }, filtered = filtered, future.seed = TRUE)
    tiles <- do.call("rbind", tiles)
    return(tiles)
}

#' convexify
#' 
#' order coordinates based on their angle around a central point
#' @param xside vector of x coordinates 
#' @param yside vector of y coordinates
#' @param indx central point x coordinate 
#' @param indy central point y coordinate
#' @param order logical - if TRUE return order only and not the
#' coordinates 
#' @details For rasterisation, the shape of the polygon must be as
#' convex as possible. To ensure that all points are in some sense 
#' in a convex form we order them based on their polar coordinate 
#' angle. 
#' @return a data frame of ordered x and y coordinates.
convexify <- function(xside, yside, indx, indy, order = FALSE) {
  #----------------------------------------------------------------------------#
  # Converting everything to an angle - from there we can just go clock wise
  # and order the point based on angle
  #----------------------------------------------------------------------------#
  x <- xside - indx
  y <- yside - indy
  angle <- mapply(function(x, y) {
    if (x >= 0 & y >= 0) angle <- atan(abs(y) / abs(x)) * (180 / pi)
    if (x < 0 & y >= 0) angle <- 180 - (atan(abs(y) / abs(x)) * (180 / pi))
    if (x < 0 & y < 0) angle <- 180 + (atan(abs(y) / abs(x)) * (180 / pi))
    if (x >= 0 & y < 0) angle <- 360 - (atan(abs(y) / abs(x)) * (180 / pi))
    return(angle)
  }, x = x, y = y, SIMPLIFY = TRUE)
  if (order) {
    convex <- data.frame("x" = order(angle, decreasing = FALSE),
    "y" = order(angle, decreasing = FALSE))
  } else {
    convex <- data.frame("x" = xside[order(angle, decreasing = FALSE)],
    "y" = yside[order(angle, decreasing = FALSE)])
  }
  
  return(convex)
}
