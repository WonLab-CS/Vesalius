################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Latent space embeddings/--------------------------------#


#' Generate embeddings.
#'
#' Generate image embddings from spatial omics data.
#' @param vesalius_assay vesalius_assay object.
#' @param dim_reduction string describing which dimensionality
#' reduction method should be used. One of the following:
#' "PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP".
#' @param normalisation string describing which normalisation 
#' method to use. One of the following "log_norm", "SCT", "TFIDF", "raw".
#' @param use_count string describing which counts should be used for the 
#' generating emebddings. Default = "raw". See details.
#' @param dimensions numeric describing the number of Principle Components or
#' Latent space dimension to use. Default dimensions = 30
#' @param tensor_resolution numeric (range 0 - 1) describing the compression
#' ratio to be applied to the final image. Default = 1
#' @param filter_grid numeric (range 0 - 1) size of the grid used when filtering
#' outlier beads. Defined as a proportion of total image size. Default = 0.1
#' @param filter_threshold numeric (range 0 -1) describing the quantile
#' threshold at which tiles should be retained (seed details)
#' @param nfeatures numeric describing the number of variable features to use.
#' @param features custom set of features to generate embeddings
#' @param min_cutoff only used when dimensionality reduction method is
#' LSI or LSI_UMAP
#' cutoff for feature to be included in the VariableFeatures for the object.
#' @param remove_lsi_1 logical only used when dimensionality reduction
#' method is LSI or LSI_UMAP
#' indicating if the first LSI component should be removed from further analysis
#' as it usually captures sequencing depth (technical variation)
#' @param verbose logical output progress message or not. Default = TRUE
#' @details The core principle behind vesalius is to convert spatial omics
#' data into an image. \code{generate_embeddings} allows you to convert
#' your omics data into a stack of gray scale images.
#' The stack hight will be equal to the number of dimenions you
#' selected excluding UMAP type as only 3 dimesnions are availbale.
#'
#' If tiles have not yet been generated (see \code{\link{generate_tiles}}),
#' pixel will be generated from the spatial coordinates provided. 
#' If tiles are already present, `generate_embeddings` will skip the
#' tile creation step. 
#'
#' The algorithm is broadly applicable to many spatial omics assays.
#' Vesalius provides a 3 nornalization methods log_norm, SCTransform,
#' and TF-IDF.
#' 
#' The `use_count` argument specifies which count matrix should be used
#' for normalization. This argument is only necessary if you use a custom
#' normalised count matrix. In this case, set this argument to the name
#' you gave your count matrix (see \code{\link{add_counts}}) and
#' `generate_embeddings` will skip the normalization and use your custom
#' count matrix to generate image embeddings. 
#' 
#' Note that it is also possible to add custom embeddings by using the 
#' \code{\link{add_embeddings}} function. 
#'
#'
#'
#' @return vesalius_assay
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run 
#' ves <- generate_embeddings(ves)
#' # maybe we want to try a different method
#' # both will be stored in the object
#' ves <- generate_embeddings(ves, dim_reduction = "UMAP")
#' 
#'}
#' @export


generate_embeddings <- function(vesalius_assay,
  dim_reduction = "PCA",
  normalisation = "log_norm",
  use_count = "raw",
  dimensions = 30,
  tensor_resolution = 1,
  filter_grid = 1,
  filter_threshold = 1,
  nfeatures = 2000,
  features = NULL,
  min_cutoff = "q5",
  remove_lsi_1 = TRUE,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # Check Status of object
    # getting out coordinates and counts from vesalius objects
    # In this case - we get extract either raw counts
    # we will extract the assay based on name (should have been filtered)
    # will always return list even if you parse a single assay
    #--------------------------------------------------------------------------#
    assay <- get_assay_names(vesalius_assay)
    normalisation <- check_norm_methods(normalisation, use_count)
    dim_reduction <- check_embed_methods(dim_reduction)
    status <- search_log(vesalius_assay,
      "generate_tiles",
      return_assay = FALSE)
    #------------------------------------------------------------------------#
    # if there are no tiles present we compute them 
    # otherwise we skip this step - no need to recompute tiles if they are
    # already there
    # NOTE : Maybe it would worth removing this step?
    # Originally it was to avoid re computing tiles everytime the user wants
    # to run a new embedding. But what if the user want to only update
    # the tile? e.g changing the tensor resolution
    #------------------------------------------------------------------------#
    if (!any(status)) {
      #----------------------------------------------------------------------#
      # generate tiles, reduce resoluation and filter out tiles and beads
      #----------------------------------------------------------------------#
      vesalius_assay <- generate_tiles(vesalius_assay,
        tensor_resolution = tensor_resolution,
        filter_grid = filter_grid,
        filter_threshold = filter_threshold,
        verbose = verbose)
      counts <- get_counts(vesalius_assay, type = use_count)
    } else {
      counts <- get_counts(vesalius_assay, type = use_count)
    }
    #--------------------------------------------------------------------------#
    # Now we can start creating colour embeddings
    # This section can be run multiple times
    # for now we dont want to have multiple "tiles" options
    # Once you compute your tiles for an assay you stuck with that 
    # First we pre-process count data -> selecting variable features
    # and normalisation.
    # if use_count is anythin else than raw we consider that the user wants 
    # to use their custom count matrix
    #--------------------------------------------------------------------------#
    counts <- process_counts(counts,
      assay = assay,
      method = normalisation,
      use_count = use_count,
      nfeatures = nfeatures,
      min_cutoff = min_cutoff,
      verbose = verbose)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = counts$norm,
      slot = "counts",
      append = TRUE)
    vesalius_assay <- add_active_count_tag(vesalius_assay,
      norm = ifelse(use_count == "raw", normalisation, use_count))
    
    if (is.null(features)) {
      features <- check_features(counts$SO)
      feature_seclection <- "variable_features"
    } else {
      feature_seclection <- "custom_features"
    }
    #--------------------------------------------------------------------------#
    # Embeddings - get embedding method and convert latent space
    # to color space.
    #--------------------------------------------------------------------------#
    embeds <- embed_latent_space(counts$SO,
      assay = assay,
      dim_reduction,
      dimensions = dimensions,
      features = features,
      remove_lsi_1 = remove_lsi_1,
      verbose = verbose)
    #--------------------------------------------------------------------------#
    # Upodating object
    #--------------------------------------------------------------------------#
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = setNames(list(features), feature_seclection),
      slot = "meta",
      append = TRUE)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeds,
      slot = "embeddings",
      append = TRUE)
    embeds <- embeds[[1L]]
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = embeds,
      slot = "active",
      append = FALSE)
    vesalius_assay <- add_active_embedding_tag(vesalius_assay,
      dim_reduction)
    #----------------------------------------------------------------------#
    # Update objects and add log
    # create log with arguments that have been used in this function 
    # add them to the log list contained in vesalius_assay
    #----------------------------------------------------------------------#
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(generate_embeddings))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      assay)
    simple_bar(verbose)
    return(vesalius_assay)
}

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
  if (filter_grid != 0 && filter_grid != 1) {
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
  # return tiles
  #----------------------------------------------------------------------#
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
#' 
#' 
#' TO DO: replace the tensor resolution reduction with super pixels 
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


#------------------------/ Preprocessing counts /------------------------------#
# NOTE: this might change if we decide do move away from Seurat as dependancy.

#' process counts 
#' 
#' pre-process count matrices
#' @param counts count matrix in the form of a sparse matrix 
#' @param assay character string describing the assay that is being 
#' pre-processed in the vesaliusObject or vesalius_assay
#' @param method character string describing which normalisation method to use.
#' One of the following "log_norm", "SCT", "TFIDF", "raw".
#' @param use_count string describing which counts should be used for the 
#' generating emebddings. Default = "raw".
#' @param nfeatures numeric describing the number of variable features to use.
#' @param min_cutoff only used when dimensionality reduction method is
#' LSI or LSI_UMAP cutoff for feature to be included in the 
#' VariableFeatures for the object.
#' @param verbose logical - progress messages outputed or not
#' @details The `use_count` argument specifies which count matrix should be used
#' for normalization. This argument is only necessary if you use a custom
#' normalised count matrix. In this case, set this argument to the name
#' you gave your count matrix (see \code{\link{add_counts}}) and
#' `generate_embeddings` will skip the normalization and use your custom
#' count matrix to generate image embeddings. 
#' @importFrom Seurat CreateSeuratObject
process_counts <- function(counts,
  assay,
  method = "log_norm",
  use_count = "raw",
  nfeatures = 2000,
  min_cutoff = "q5",
  verbose = TRUE) {
     message_switch("in_assay",
      verbose,
      comp_type = "Pre-processing counts",
      assay = assay)
    #--------------------------------------------------------------------------#
    # We are still going to use Seurat for now
    # rememmber that if we do decide to change things
    # we have to change things in the embbeddings as well
    #--------------------------------------------------------------------------
    counts <- suppressWarnings(Seurat::CreateSeuratObject(counts))
    counts <- switch(method,
                    "log_norm" = log_norm(counts, nfeatures),
                    "SCTransform" = int_sctransform(counts, nfeatures),
                    "TFIDF" = tfidf_norm(counts, min_cutoff = min_cutoff),
                    "raw" = raw_norm(counts, use_count))
    return(counts)
}

#' raw norm
#' 
#' no normalisation applied simply return raw counts
#' @param counts seurat object containing counts
#' @param use_count string describing name that needs to be added to
#' list element. This list will be appended to the count slot in
#' the vesalius_assay. 
#' @details Here, either the user doesn't want to normalise the data or
#' they provide their custom count matrix. In this case, we parse it 
#' as "raw" to avoid writing another function and add the custom name. 
#' @return list with seurat object used later and raw counts to be stored in
#' the vesalius objects 
#' @importFrom Seurat GetAssayData
raw_norm <- function(counts, use_count = "raw") {
    #--------------------------------------------------------------------------#
    # Essentially we want people to be able to parse their matrix
    # If they want to use a different type of norm method that is not present
    # or play around with parameters not provided by vesalius
    # they can do that and just always call norm
    # We are using this just for formating at the moment
    # We have to be a bit hacky with the Seurat object
    #--------------------------------------------------------------------------#
    norm_counts <- list(Seurat::GetAssayData(counts, layer = "counts"))
    counts@assays$RNA@scale.data <- as.matrix(Seurat::GetAssayData(counts,
      slot = "counts"))
    names(norm_counts) <- use_count
    return(list("SO" = counts, "norm" = norm_counts))
}
#' log norm
#' 
#' log normalisation, scaling and top variable features
#' @param counts seurat object containing counts
#' @param nfeatures number of top variable features to select
#' @return list with seurat object used later and normalised counts to be stored
#' in a vesalius object
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat GetAssayData
log_norm <- function(counts, nfeatures) {
  counts <- Seurat::NormalizeData(counts, verbose = FALSE)
  counts <- Seurat::ScaleData(counts, verbose = FALSE)
  counts <- suppressWarnings(Seurat::FindVariableFeatures(counts,
    nfeatures = nfeatures,
    verbose = FALSE))
  norm_counts <- list(Seurat::GetAssayData(counts, layer = "data"))
  names(norm_counts) <- "log_norm"
  return(list("SO" = counts, "norm" = norm_counts))
}

#' SCTransform
#' 
#' SCTransform pre-processing from Seurat
#' @param counts seurat object containing counts
#' @param nfeatures number of top variable features to select
#' @return list with seurat object used later and normalised counts to be stored
#' in a vesalius object
#' @importFrom Seurat SCTransform
#' @importFrom Seurat GetAssayData
int_sctransform <- function(counts, nfeatures) {
    counts <- suppressWarnings(Seurat::SCTransform(counts,
      variable.features.n = nfeatures, verbose = FALSE))
    norm_counts <- list(GetAssayData(counts, layer = "data"))
    names(norm_counts) <- "SCTransform"
    return(list("SO" = counts, "norm" = norm_counts))
}

#' tf idf normalisation
#' 
#' nornalising count using TF IDF
#' @param counts Seurat object containing counts
#' @param min_cutoff min cutoff of for top features
#' list with seurat object used later and normalised counts to be stored
#' in a vesalius object
#' @importFrom Signac RunTFIDF
#' @importFrom Seurat ScaleData
#' @importFrom Signac FindTopFeatures
#' @importFrom Seurat GetAssayData
tfidf_norm <- function(counts, min_cutoff) {
  counts <- Signac::RunTFIDF(counts, verbose = FALSE)
  counts <- Seurat::ScaleData(counts, verbose = FALSE)
  counts <- suppressWarnings(
    Signac::FindTopFeatures(counts, min.cutoff = min_cutoff))
  norm_counts <- list(Seurat::GetAssayData(counts, layer = "data"))
  names(norm_counts) <- "TFIDF"
  return(list("SO" = counts, "norm" = norm_counts))
}



#------------------------/ Color Embeddings /----------------------------------#

#' embed latent space
#' 
#' Embed latent space into grey color scale.
#' @param counts Seurat object containing counts (generally normalised)
#' @param assay charcter string of the assay being used 
#' @param dim_reduction dimensionality reduction method that will be used
#' Select from PCA, PCA_L, UMAP, LSI, LSI_UMAP
#' @param dimensions numeric for number of dimeniosn top retain after 
#' dimensionality reduction
#' @param features custom features used for dim reduction
#' @param remove_lsi_1 logical if first dimension of LSI embedding should be 
#' removed (will soon be depreciated)
#' @param verbose logical if progress messages should be outputed or not
#' @details General method dispatch function for dim reduction methods 
#' @return data frame of normalised embedding values.
embed_latent_space <- function(counts,
  assay,
  dim_reduction,
  dimensions,
  features = NULL,
  remove_lsi_1,
  verbose) {
    message_switch("in_assay",
      verbose,
      comp_type = "Compute Latent Space",
      assay = assay)
    embeds <- switch(dim_reduction,
      "PCA" = embed_pca(counts,
        dimensions = dimensions,
        features = features,
        verbose = verbose),
      "NMF" = embed_nmf(counts,
                        dimensions = dimensions,
                        verbose = verbose),
      "PCA_L" = embed_pcal(counts,
        dimensions = dimensions,
        features = features,
        verbose = verbose),
      "UMAP" = embed_umap(counts,
        dimensions = dimensions,
        features = features,
        verbose),
      "LSI" = embed_lsi(counts,
        dimensions = dimensions,
        features = features,
        remove_lsi_1),
      "LSI_UMAP" = embed_lsi_umap(counts,
        dimensions = dimensions,
        features = features,
        remove_lsi_1))
    return(embeds)
}

#' embed PCA
#' 
#' embed in grey scale using PCA embeddings 
#' @param counts Seurat object containing normalised counts
#' @param dimensions number dimension to retain from PCA
#' @param features custom vector of features
#' @param verbose logical if progress messages should be outputed
#' @return normalised PCA embedding matrix 
#' @importFrom Seurat RunPCA
#' @importFrom Seurat Embeddings
embed_pca <- function(counts,
  dimensions,
  features = NULL,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # First run PCA
    # Parsing feature names just in case no var features have been 
    #--------------------------------------------------------------------------#
    message_switch("pca_tensor", verbose)
    counts <- Seurat::RunPCA(counts,
      npcs = dimensions,
      features = features,
      verbose = FALSE)

    # Progress message embed_rgb_tensor() => Prog.R
    message_switch("pca_rgb_tensor", verbose)
    #--------------------------------------------------------------------------#
    # Here we can just sum and normalise
    # this is going to be much faster
    # transpose at the end so we keep common format
    # ATM only min max norm - could look into using quantile norm as well
    #--------------------------------------------------------------------------#
    pca <- Seurat::Embeddings(counts, reduction = "pca")
    pca <- apply(pca, 2, norm_pixel, "minmax")
    colour_matrix <- list(as.matrix(pca))
    names(colour_matrix) <- "PCA"
    return(colour_matrix)
}

#' embed PCA loading values
#' 
#' embed in grey scale using PCA Loading value 
#' @param counts Seurat object containing normalised counts
#' @param dimensions number dimension to retain from PCA
#' @param features custom vector of features
#' @param verbose logical if progress messages should be outputed
#' @details This approach is a slightly different as it takes 
#' the loading value associted to each gene in a barcode and sums
#' the absolute value of each of those values. 
#' Once all genes in all barcodes have been summed,
#' we normalise the latent space and return the matrix. 
#' @return normalised PCA loading matrix 
#' @importFrom Seurat RunPCA
#' @importFrom Seurat Loadings
#' @importFrom Seurat GetAssayData
#' @importFrom future.apply future_lapply
embed_pcal <- function(counts,
  dimensions,
  features = NULL,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # First run PCA
    #--------------------------------------------------------------------------#
    message_switch("pca_tensor", verbose)
    counts <- Seurat::RunPCA(counts,
      npcs = dimensions,
      features = features,
      verbose = FALSE)

    #--------------------------------------------------------------------------#
    # get laodings and create matrix describing if there are any count values
    #--------------------------------------------------------------------------#
    pca <- Seurat::Loadings(counts, reduction = "pca")
    pca <- apply(pca, 2, function(x) return(abs(x)))
    mat <- as.matrix(Seurat::GetAssayData(counts, slot = "data") > 0)
    colour_matrix <- matrix(0, nrow = ncol(mat), ncol = ncol(pca))
    colnames(colour_matrix) <- colnames(pca)
    rownames(colour_matrix) <- colnames(mat)
    #--------------------------------------------------------------------------#
    # Looping over each PC to get laodings sum and normalising
    #--------------------------------------------------------------------------#
    for (p in seq_len(ncol(pca))) {
      message_switch("pcal_rgb_tensor", verbose, pc = p)
      bars <- future_lapply(seq_len(ncol(mat)), function(idx, mat, pca) {
        c_vec <- as.numeric(mat[names(pca), idx])
        colour <- sum(pca * c_vec)
        return(colour)
      }, mat = mat, pca = pca[, p], future.seed = TRUE)
      bars <- unlist(bars)

      colour_matrix[, p] <- norm_pixel(bars, "minmax")
    }
    colour_matrix <- list(as.matrix(colour_matrix))
    names(colour_matrix) <- "PCA_L"
    return(colour_matrix)
}

#' embed umap
#' 
#' embed in gray scale using UMAP projections
#' @param counts Seurat object containing normalised counts 
#' @param dimensions number of PCs to use for the UMAP projections
#' @param features custom vector of features
#' @param verbose logical if progress messages should be outputed 
#' @details Note that while you can select any number of dimensions
#' the number of UMAP dimensions will always be 3.
#' @return normalised UMAP projection matrix
#' @importFrom Seurat RunPCA
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat FetchData

embed_umap <- function(counts,
  dimensions,
  features = NULL,
  verbose) {
  #----------------------------------------------------------------------------#
  # First run PCA and then UMAP
  # Progress message pca_tensor() => Prog.R
  #----------------------------------------------------------------------------#
  message_switch("pca_tensor", verbose)

  counts <- Seurat::RunPCA(counts,
    npcs = dimensions,
    features = features,
    verbose = FALSE)
  message_switch("umap_rgb_tensor", verbose)
  counts <- suppressWarnings(Seurat::RunUMAP(counts,
    dims = seq_len(dimensions),
    n.components = 3L,
    verbose = FALSE))
  #----------------------------------------------------------------------------#
  # Normalise
  #----------------------------------------------------------------------------#
  counts <- Seurat::FetchData(counts, c("umap_1", "umap_2", "umap_3"))
  counts <- apply(counts, 2, norm_pixel, "minmax")
  counts <- list(as.matrix(counts))
  names(counts) <- "UMAP"
  return(counts)
}

#' embed lsi
#' 
#' embed in grey scale using latent semantic indexing
#' @param counts Seurat object containing normalised counts
#' @param dimensions numeric for number of latent space dimensions to use
#' @param features custom vector of features
#' @param remove_lsi_1 logical if first LSI dimenions should be removed 
#' @param verbose logical if progress messages should be outputed 
#' @returns normalised LSI embedding matrix 
#' @importFrom Signac RunSVD
#' @importFrom Seurat Embeddings
embed_lsi <- function(counts,
  dimensions,
  features = NULL,
  remove_lsi_1,
  verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # Run partial singular value decomposition(SVD) on TF-IDF normalized matrix
  #--------------------------------------------------------------------------#
  message_switch("svd_tensor", verbose)
  
  svd <- Signac::RunSVD(counts,
    n = dimensions + 1,
    features = features,
    verbose = FALSE)

  #--------------------------------------------------------------------------#
  # Getting embedding values and normalize
  #--------------------------------------------------------------------------#
  message_switch("svd_rgb_tensor", verbose)
  if (remove_lsi_1) {
    colour_matrix <- Seurat::Embeddings(svd[["lsi"]])[, -1]
  } else {
    colour_matrix <- Seurat::Embeddings(svd[["lsi"]])[, seq(1, dimensions)]
  }

  colour_matrix <- apply(colour_matrix, 2, norm_pixel, "minmax")
  colour_matrix <- list(as.matrix(colour_matrix))
  names(colour_matrix) <- "LSI"
  return(colour_matrix)
}


#' embed lsi
#' 
#' embed in grey scale using latent semantic indexing 
#' followed by UMAP
#' @param counts Seurat object containing normalised counts
#' @param dimensions numeric for number of latent space dimensions to use
#' @param features custom vector of features
#' @param remove_lsi_1 logical if first LSI dimenions should be removed 
#' @param verbose logical if progress messages should be outputed 
#' @returns normalised LSI embedding matrix 
#' @importFrom Signac RunSVD
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat FetchData

embed_lsi_umap <- function(counts,
  dimensions,
  features = NULL,
  remove_lsi_1,
  verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # Run partial singular value decomposition(SVD) on TF-IDF normalized matrix
  #--------------------------------------------------------------------------#
  message_switch("svd_tensor", verbose)
  
  svd <- Signac::RunSVD(counts,
    n = dimensions + 1,
    features = features,
    verbose = FALSE)

  message_switch("umap_rgb_tensor", verbose)
  if (remove_lsi_1) {
    reduc <-  suppressWarnings(Seurat::RunUMAP(svd,
      reduction = "lsi",
      dims = seq(2, dimensions + 1),
      n.components = 3L,
      verbose = FALSE))
  } else {
    reduc <- suppressWarnings(Seurat::RunUMAP(svd,
        reduction = "lsi",
        dims = seq(1, dimensions),
        n.components = 3L,
        verbose = FALSE))
  }

  #--------------------------------------------------------------------------#
  # Getting embedding values and normalize
  #--------------------------------------------------------------------------#
  colour_matrix <- Seurat::FetchData(reduc, c("UMAP_1", "UMAP_2", "UMAP_3"))

  colour_matrix <- apply(colour_matrix, 2, norm_pixel, "minmax")

  colour_matrix <- list(as.matrix(colour_matrix))
  names(colour_matrix) <- "LSI_UMAP"
  return(colour_matrix)
}
#' embed nmf
#' 
#' embed in grey scale using NMF embeddings 
#' @param counts Seurat object containing normalised counts
#' @param dimensions number dimension to retain from NMF
#' @param verbose logical if progress messages should be outputed
#' @return normalised NMF embedding matrix 
#' @importFrom NMF nmf
#' @importFrom NMF coefficients
#' @importFrom NMF basis
embed_nmf <- function(counts, dimensions, verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # adding this since I don't want to have this package as a dependancy 
  # rather I want it to be selectively loaded via namespace 
  #--------------------------------------------------------------------------#
  if (!isNamespaceLoaded("NMF")) {
    inst <- requireNamespace("NMF", quietly = TRUE)
    if (!inst) {
      stop("NMF is not installed - Please install NMF
      install.packages('NMF')
      https://cran.r-project.org/web/packages/NMF/index.html")
    } else {
      require("NMF")
    }
  }
  #--------------------------------------------------------------------------#
  # Get the normalized count matrix and matrix with variable features 
  # from the Seurat object
  #--------------------------------------------------------------------------#
  vesalius:::message_switch("nmf_tensor", verbose)
  count_matrix <- as.matrix(Seurat::GetAssayData(counts, slot = "data"))
  features <- counts@assays[["RNA"]]@var.features
  count_matrix <- count_matrix[features, ]
  #--------------------------------------------------------------------------#
  # Run NMF
  #--------------------------------------------------------------------------#
  nmf_result <- NMF::nmf(count_matrix,
    rank = dimensions,
    method = "lee",
    .options = "-vp")
  
  #--------------------------------------------------------------------------#
  # Get the NMF projections (W matrix) and normalize
  #--------------------------------------------------------------------------#
  nmf_projections <- t(NMF::coefficients(nmf_result))
  nmf_projections <- apply(nmf_projections, 2, vesalius:::norm_pixel, "minmax")
  nmf_projections <- list(as.matrix(nmf_projections))
  names(nmf_projections) <- "NMF"
  
  return(nmf_projections)
}

