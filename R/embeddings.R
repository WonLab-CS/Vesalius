################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Latent space embeddings/--------------------------------#

#' Create Vesalius Image embeddings from ST data.
#' @param counts count matrix (matrix or sparse matrix) with barcodes as column
#' names. Rows can represent genes (transcriptomics) or bins (DNA).
#' @param coordinates data frame containing barcodes, x coordinates and
#' y coordinates
#' @param dim_reduction character describing which dimensionality reduction
#' method should be used. One of the following:
#' "PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP".
#' @param normalisation character describing which normalisation method to use.
#' One of the following "log", "SCT", "TFIDF", "raw".
#' @param dimensions numeric describing the number of Principle Components or
#' Latent space dimension to use. Default dimensions = 30
#' @param tensor_resolution numeric (range 0 - 1) describing the compression
#' ratio to be applied to the final image. Default = 1
#' @param filter_grid numeric (range 0 - 1) size of the grid used when filtering
#' outlier beads. Defined as a proportion of total image size. Default = 0.1
#' @param filter_threshold numeric (range 0 -1) describing the quantile
#' threshold at which tiles should be retained (seed details)
#' @param nfeatures numeric describing the number of variable features to use.
#' @param min_cutoff only used when dimensionality reduction method is
#' LSI or LSI_UMAP
#' cutoff for feature to be included in the VariableFeatures for the object.
#' @param remove_LSI1 logical only used when dimensionality reduction
#' method is LSI or LSI_UMAP
#' indicating if the first LSI component should be removed from further analysis
#' as it usually captures sequencing depth (technical variation)
#' @param cores numeric number of cores to use. Default = 1
#' @param verbose logical output progress message or not. Default = TRUE
#' @details
#' for more detail please look for Signac
#' "https://satijalab.org/signac/reference/findtopfeatures"
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","slice"
#' Slice represents the PC slice represented.
#' @examples
#' \dontrun{
#' data(Vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' }

build_vesalius_embeddings <- function(vesalius,
  dim_reduction = "PCA",
  normalisation = "log_norm",
  assay = "all",
  dimensions = 30,
  tensor_resolution = 1,
  filter_grid = 0.01,
  filter_threshold = 0.995,
  nfeatures = 2000,
  min_cutoff = "q5",
  remove_lsi_1 = TRUE,
  cores = 1,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # Check Status of object
    # getting out coordinates and counts from vesalius objects
    # In this case - we always get raw counts even if other are present
    # there is a normalisation step in this function
    # we will extract the assay based on name.
    #--------------------------------------------------------------------------#
    coordinates <- get_tiles(vesalius, assay)
    counts <- get_counts(vesalius, type = "raw", assay)
    assays <- names(counts)
    #--------------------------------------------------------------------------#
    # if there are no tiles present we compute them 
    # otherwise we skip this step - no need to recompute tiles if they are
    # already there
    #--------------------------------------------------------------------------#
    if (!is.null(coordinates)) {
      #------------------------------------------------------------------------#
      # generate tiles, reduce resoluation and filter out tiles and beads
      #------------------------------------------------------------------------#
      tiles <- mapply(generate_tiles,
        coordinates,
        assays = assays,
        MoreArgs = list(
          tensor_resolution = tensor_resolution,
          filter_grid = filter_grid,
          filter_threshold = filter_threshold,
          verbose = verbose,
          cores = cores))
      #------------------------------------------------------------------------#
      # adjusted counts if necessary
      # essentially merging counts when barcodes overlap
      #------------------------------------------------------------------------#
      if (tensor_resolution < 1) {
        message_switch("adj_counts", verbose)
        counts <- mapply(adjust_counts,
          tiles,
          counts,
          MoreArgs = list(cores = cores),
          SIMPLIFY = FALSE)
      }
    } 
    #--------------------------------------------------------------------------#
    # Now we can start creating colour embeddings
    # This section can be run multiple times
    # for now we dont want to have multiple "tiles" options
    # Once you compute your tiles for an assay you stuck with that 
    # NOTE: we might want to get away from Seurat as dependancy!!!
    #--------------------------------------------------------------------------#
    normalisation <- check_norm_methods(normalisation, length(counts))
    counts <- mapply(process_counts,
      counts,
      assays = assays,
      method = normalisation,
      MoreArgs = list(
        nfeatures = nfeatures,
        min_cutoff = min_cutoff,
        verbose = verbose),
      SIMPLIFY = FALSE)
    #--------------------------------------------------------------------------#
    # Embeddings - get embedding method and convert latent space
    # to color space.
    #--------------------------------------------------------------------------#
    dim_reduction <- check_embed_methods(dim_reduction, length(counts))
    embeds <- mapply(embed_latent_space,
      lapply(counts, "[[", 1),
      assays = assays,
      dim_reduction,
      MoreArgs = list(
        dimensions = dimensions,
        cores = cores,
        verbose = verbose),
      SIMPLIFY = FALSE)
    #----------------------------------------------------------------------#
    # Update objects and add log
    # update_vesalius => objectUtilies.R
    # Updating all slots that have been modified
    # we also create a commit list => clean argument list when more than
    # 1 argument is supplied. will be commited to individual assasys
    #----------------------------------------------------------------------#
    commit <- create_commit_log(vesalius = vesalius,
      match = as.list(match.call()),
      default = as.list(args(build_vesalius_embeddings)),
      assay = assay,
      normalisation = normalisation,
      dim_reduction = dim_reduction)
    vesalius <- update_vesalius(vesalius = vesalius,
      data = tiles,
      slot = "tiles",
      commit = commit,
      append = FALSE)
    vesalius <- update_vesalius(vesalius = vesalius,
      data = lapply(counts, "[[", 2),
      slot = "counts",
      commit = commit,
      append = TRUE)
    vesalius <- update_vesalius(vesalius = vesalius,
      data = embeds,
      slot = "embeddings",
      commit = commit,
      append = TRUE)
    # Progress message simpleBar => Prog.R
    simple_bar(verbose)
    return(vesalius)
}

generate_tiles <- function(coordinates,
  assays,
  tensor_resolution = 1,
  filter_grid = 0.01,
  filter_threshold = 0.995,
  verbose = TRUE,
  cores = 1) {
  message_switch("in_assay",
    verbose = verbose,
    assay = assays,
    comp_type = "Generating Tiles")
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
  # could be updated to use KNN
  #----------------------------------------------------------------------#
  if (tensor_resolution < 1) {
    message_switch("tensor_res", verbose)
    coordinates <- reduce_tensor_resolution(coordinates = coordinates,
      tensor_resolution = tensor_resolution)
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
  tiles <- rasterise(filtered, cores)
  #----------------------------------------------------------------------#
  # return tiles and adjusted counts
  #----------------------------------------------------------------------#
  return(tiles)
}




#------------------------/ Filtering  Beads /----------------------------------#
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
### maybe that would be better -> aggregate points together so it's "fair"
reduce_tensor_resolution <- function(coordinates, tensor_resolution = 1) {
  #----------------------------------------------------------------------------#
  # we will reduce number of points this way
  # this should keep all barcodes - with overlapping coordinates
  #----------------------------------------------------------------------------#
  coordinates$x <- round(coordinates$x * tensor_resolution)
  coordinates$y <- round(coordinates$y * tensor_resolution)
  #----------------------------------------------------------------------------#
  # Now we get coordinate tags - we use this to find all the merge locations
  # sorting and using rle to ensure that we actually merge them
  #----------------------------------------------------------------------------#
  tag <- paste0(coordinates$x, "_", coordinates$y)

  locs <- rle(sort(tag))
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


adjust_counts <- function(coordinates, counts, cores = 1) {
    #--------------------------------------------------------------------------#
    # First get all barcode names and compare which ones are missing
    #--------------------------------------------------------------------------#
    coord_bar <- unique(coordinates$barcodes)
    coord_bar <- coord_bar[sapply(strsplit(coord_bar, "_et_"), length) > 1]
    if (length(coord_bar) == 0) {
       return(counts)
    }

    #--------------------------------------------------------------------------#
    # next we merge counts together when barcodes have been merged
    #--------------------------------------------------------------------------#
    tmp_bar <- strsplit(coord_bar, "_et_")

    empty <- parallel::mclapply(tmp_bar, function(coord, count) {
      tmp <- rowSums(count[, coord])
      return(tmp)
    }, count = counts, mc.cores = cores)
    empty <- do.call("cbind", empty)
    if (is.null(dim(empty)) && length(empty) != 0) {
        empty <- Matrix(empty, ncol = 1)
    }
    colnames(empty) <- coord_bar
    merged <- cbind(counts[, !colnames(counts) %in% unlist(unique(tmp_bar))],
      empty)
    #--------------------------------------------------------------------------#
    # next we remove any barcodes that were dropped during filtering
    #--------------------------------------------------------------------------#
    merged <- merged[, colnames(merged) %in% coordinates$barcodes]
    return(merged)
}

#------------------------/ Creating pixel tiles /------------------------------#
filter_tiles <- function(tesselation, coordinates, filter_threshold) {
  max_area <- quantile(tesselation$summary$dir.area, filter_threshold)
  idx <- which(tesselation$summary$dir.area >= max_area)
  tess_v <- tesselation$dirsgs
  points_to_remove <- tess_v$ind1 %in% idx | tess_v$ind2 %in% idx
  tess_v <- tess_v[!points_to_remove, ]
  coordinates$ind <- seq_len(nrow(coordinates))
  coordinates <- coordinates[-idx, ]
  return(list("tess_v" = tess_v, "coordinates" = coordinates))
}

rasterise <- function(filtered, cores = 1) {
    idx <- seq_len(nrow(filtered$coordinates))
    tiles <- parallel::mclapply(idx, function(idx, filtered) {

        #----------------------------------------------------------------------#
        # get indecies from original data
        #----------------------------------------------------------------------#
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
        x <- convex$x
        y <- convex$y
        #----------------------------------------------------------------------#
        # define max polygon containing all pixels
        #----------------------------------------------------------------------#
        lpx <- round(min(x)) - 1
        hpx <- round(max(x)) + 1
        lpy <- round(min(y)) - 1
        hpy <- round(max(y)) + 1

        max_polygon_x <- rep(seq(lpx, hpx), times = hpy - lpy + 1)
        max_polygon_y <- rep(seq(lpy, hpy), each = hpx - lpx + 1)

        #----------------------------------------------------------------------#
        # Fill triangles with all point in that space
        #----------------------------------------------------------------------#
        cell <- sp::point.in.polygon(max_polygon_x, max_polygon_y, x, y)
        max_x <- max_polygon_x[cell %in% c(1, 2, 3)]
        max_y <- max_polygon_y[cell %in% c(1, 2, 3)]
        cent <- which(max_x == round(indx) & max_y == round(indy))
        centers <- rep(0, length(max_x))
        centers[cent] <- 1
        tile <- data.frame("barcodes" = rep(filtered$coordinates$barcodes[idx],
            times = length(max_x)),
          "x" = max_x,
          "y" = max_y,
          "origin" = centers)
        return(tile)
    }, filtered = filtered, mc.cores = cores)
    tiles <- do.call("rbind", tiles)
    tiles <- tiles %>% filter(x > 1 & y > 1)
    return(tiles)
}

convexify <- function(xside, yside, indx, indy) {
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
  convex <- data.frame("x" = xside[order(angle, decreasing = FALSE)],
    "y" = yside[order(angle, decreasing = FALSE)])
  return(convex)
}


#------------------------/ Preprocessing counts /------------------------------#
# NOTE: this might change if we decide do move away from Seurat as dependancy.

process_counts <- function(counts,
  assays,
  method = "log",
  nfeatures = 2000,
  min_cutoff = "q5",
  verbose = TRUE) {
     message_switch("in_assay",
      verbose,
      comp_type = "Pre-processing counts",
      assay = assays)
    #--------------------------------------------------------------------------#
    # We are still going to use Seurat for now
    # rememmber that if we do decide to change things
    # we have to change things in the embbeddings as well
    #--------------------------------------------------------------------------
    counts <- CreateSeuratObject(counts[["raw"]], assay = "Spatial")
    counts <- switch(method[1L],
                    "log_norm" = log_norm(counts, nfeatures),
                    "SCTransform" = int_sctransform(counts, assay = "Spatial",
                      nfeatures = nfeatures),
                    "TFIDF" = tfidf_norm(counts, min_cutoff = min_cutoff),
                    "raw" = raw_norm(counts))
    return(counts)
}

raw_norm <- function(counts) {
    #--------------------------------------------------------------------------#
    # Essentially we want people to be able to parse their matrix
    # If they want to use a different type of norm method that is not present
    # or play around with parameters not provided by vesalius
    # they can do that and just always call norm
    # We are using this just for formating at the moment
    #--------------------------------------------------------------------------#
    norm_counts <- list(Seurat::GetAssayData(counts, slot = "counts"))
    names(norm_counts) <- "raw"
    return(list("SO" = counts, "norm" = norm_counts))
}

log_norm <- function(counts, nfeatures) {
  counts <- Seurat::NormalizeData(counts, verbose = FALSE)
  counts <- Seurat::ScaleData(counts, verbose = FALSE)
  counts <- Seurat::FindVariableFeatures(counts,
    nfeatures = nfeatures,
    verbose = FALSE)
  norm_counts <- list(Seurat::GetAssayData(counts, slot = "data"))
  names(norm_counts) <- "log_norm"
  return(list("SO" = counts, "norm" = norm_counts))
}

int_sctransform <- function(counts, assay = "Spatial", nfeatures) {
    counts <- Searat::SCTransform(counts, assay = "Spatial",
      variable.features.n = nfeatures, verbose = FALSE)
    norm_counts <- list(GetAssayData(counts, slot = "data"))
    names(norm_counts) <- "SCT"
    return(list("SO" = counts, "norm" = norm_counts))
}

tfidf_norm <- function(counts, min_cutoff) {
  counts <- Signac::RunTFIDF(counts)
  counts <- Signac::FindTopFeatures(counts, min.cutoff = min_cutoff)
  norm_counts <- list(Seurat::GetAssayData(counts, slot = "data"))
  names(norm_counts) <- "TFIDF"
  return(list("SO" = counts, "norm" = norm_counts))
}



#------------------------/ Color Embeddings /----------------------------------#
embed_latent_space <- function(counts,
  assays,
  dim_reduction,
  dimensions,
  cores,
  verbose) {
    message_switch("in_assay",
      verbose,
      comp_type = "Compute Latent Space",
      assay = assays)
    embeds <- switch(dim_reduction,
      "PCA" = embed_pca(counts,
        dimensions = dimensions,
        verbose = verbose),
      "PCA_L" = embed_pcal(counts,
        dimensions = dimensions,
        cores = cores,
        verbose = verbose),
      "UMAP" = embed_umap(counts,
        dimensions = dimensions,
        verbose),
      "LSI" = embed_lsi(counts,
        dimensions = dimensions,
        remove_LSI1),
      "LSI_UMAP" = embed_lsi_umap(counts,
        dimensions = dimensions,
        remove_LSI1))
    return(embeds)
}
embed_pca <- function(counts,
  dimensions,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # First run PCA
    # Progress message pca_tensor() => Prog.R
    #--------------------------------------------------------------------------#
    message_switch("pca_tensor", verbose)
    counts <- RunPCA(counts, npcs = dimensions, verbose = FALSE)

    # Progress message embed_rgb_tensor() => Prog.R
    message_switch("pca_rgb_tensor", verbose)
    #--------------------------------------------------------------------------#
    # Here we can just sum and normalise
    # this is going to be much faster
    # transpose at the end so we keep common format
    # ATM only min max norm - could look into using quantile norm as well
    #--------------------------------------------------------------------------#
    pca <- Embeddings(counts, reduction = "pca")
    pca <- norm_pixel(pca, "minmax")
    colour_matrix <- list(as.matrix(pca))
    names(colour_matrix) <- "PCA"
    return(colour_matrix)
}

embed_pcal <- function(counts,
  dimensions,
  cores = 1,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # First run PCA
    # Progress message pca_tensor() => Prog.R
    #--------------------------------------------------------------------------#
    message_switch("pca_tensor", verbose)
    counts <- RunPCA(counts, npcs = dimensions, verbose = FALSE)

    #--------------------------------------------------------------------------#
    # get laodings and create matrix describing if there are any count values
    #--------------------------------------------------------------------------#
    pca <- Loadings(counts, reduction = "pca")
    pca <- apply(pca, 2, function(x) return(abs(x)))
    mat <- as.matrix(GetAssayData(counts, slot = "data") > 0)
    colour_matrix <- matrix(0, nrow = ncol(mat), ncol = ncol(pca))
    colnames(colour_matrix) <- colnames(pca)
    rownames(colour_matrix) <- colnames(mat)
    #--------------------------------------------------------------------------#
    # Looping over each PC to get laodings sum and normalising
    #--------------------------------------------------------------------------#
    for (p in seq_len(ncol(pca))) {
      message_switch("pcal_rgb_tensor", verbose, pc = p)
      bars <- parallel::mclapply(seq_len(ncol(mat)), function(idx, mat, pca) {
        c_vec <- as.numeric(mat[names(pca), idx])
        colour <- sum(pca * c_vec)
        return(colour)
      }, mat = mat, pca = pca[, p], mc.cores = cores)
      bars <- unlist(bars)

      colour_matrix[, p] <- norm_pixel(bars,"minmax")
    }
    colour_matrix <- list(as.matrix(colour_matrix))
    names(colour_matrix) <- "PCA_L"
    return(colour_matrix)
}



embed_umap <- function(counts, dimensions, verbose) {
  #----------------------------------------------------------------------------#
  # First run PCA and then UMAP
  # Progress message pca_tensor() => Prog.R
  #----------------------------------------------------------------------------#
  message_switch("pca_tensor", verbose)
  counts <- RunPCA(counts, npcs = dimensions, verbose = FALSE)
  message_switch("umap_rgb_tensor", verbose)
  counts <- RunUMAP(counts,
    dims = seq_len(dimensions),
    n.components = 3L,
    verbose = FALSE)
  #----------------------------------------------------------------------------#
  # Normalise
  #----------------------------------------------------------------------------#
  counts <- FetchData(counts, c("UMAP_1", "UMAP_2", "UMAP_3"))
  counts <- apply(counts, 2, norm_pixel, "minmax")
  counts <- list(as.matrix(counts))
  names(counts) <- "UMAP"
  return(counts)
}


embed_lsi <- function(counts,
  dimensions = dimensions,
  remove_lsi_1,
  verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # Run partial singular value decomposition(SVD) on TF-IDF normalized matrix
  #--------------------------------------------------------------------------#
  message_switch("svd_tensor", verbose)
  svd <- Signac::RunSVD(counts, n = dimensions + 1, verbose = FALSE)

  #--------------------------------------------------------------------------#
  # Getting embedding values and normalize
  #--------------------------------------------------------------------------#
  message_switch("svd_rgb_tensor", verbose)
  if (remove_lsi_1) {
    colour_matrix <- Embeddings(svd[["lsi"]])[, -1]
  } else {
    colour_matrix <- Embeddings(svd[["lsi"]])[, seq(1, dimensions + 1)]
  }

  colour_matrix <- apply(colour_matrix, 2, norm_pixel, "minmax")
  colour_matrix <- list(as.matrix(colour_matrix))
  names(colour_matrix) <- "LSI"
  return(colour_matrix)
}




embed_lsi_umap <- function(counts,
  dimensions,
  remove_lsi_1,
  verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # Run partial singular value decomposition(SVD) on TF-IDF normalized matrix
  #--------------------------------------------------------------------------#
  message_switch("svd_tensor", verbose)
  svd <- Signac::RunSVD(counts, n = dimensions + 1, verbose = FALSE)

  message_switch("umap_rgb_tensor", verbose)
  if (remove_lsi_1) {
    reduc <-  RunUMAP(svd,
      reduction = "lsi",
      dims = seq(2, dimensions + 1),
      n.components = 3L,
      verbose = FALSE)
  } else {
    reduc <- RunUMAP(svd,
        reduction = "lsi",
        dims = seq(1, dimensions),
        n.components = 3L,
        verbose = FALSE)
  }

  #--------------------------------------------------------------------------#
  # Getting embedding values and normalize
  #--------------------------------------------------------------------------#
  colour_matrix <- FetchData(reduc, c("UMAP_1", "UMAP_2", "UMAP_3"))

  colour_matrix <- apply(colour_matrix, 2, norm_pixel, "minmax")

  colour_matrix <- list(as.matrix(colour_matrix))
  names(colour_matrix) <- "LSI_UMAP"
  return(colour_matrix)
}
