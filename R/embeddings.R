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
#' @param normalization string describing which normalisation 
#' method to use. One of the following "log_norm", "SCT", "TFIDF", "none".
#' @param use_counts string describing which counts should be used for the 
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
#' selected excluding UMAP embeddings as only 3 dimesnions are availbale.
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
#' The `use_counts` argument specifies which count matrix should be used
#' for normalization. This argument is only necessary if you use a custom
#' normalised count matrix. In this case, set this argument to the name
#' you gave your count matrix (see \code{\link{add_counts}}) and
#' `generate_embeddings` will skip the normalization and use your custom
#' count matrix to generate image embeddings. 
#' 
#' Note that it is also possible to add custom embeddings by using the 
#' \code{\link{add_embeddings}} function. 
#'
#' Embeddings can be created using a custom set of genes. These genes 
#' can be provided to the `features` argument in the form of a character
#' vector. Note that this will not filter the count matrix hence you 
#' will still have access to the whole matrix for downsteam analysis. 
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
  normalization = "log_norm",
  use_counts = "raw",
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
    normalization <- check_norm_methods(normalization, use_counts)
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
      counts <- get_counts(vesalius_assay, type = use_counts)
    } else {
      counts <- get_counts(vesalius_assay, type = use_counts)
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
      method = normalization,
      use_count = use_counts,
      nfeatures = nfeatures,
      min_cutoff = min_cutoff,
      verbose = verbose)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
      data = counts$norm,
      slot = "counts",
      append = TRUE)
    vesalius_assay <- add_active_count_tag(vesalius_assay,
      norm = ifelse(use_counts == "raw", normalization, use_counts))
    
    
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
#' @importFrom stats setNames
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
    mat <- as.matrix(Seurat::GetAssayData(counts, layer = "data") > 0)
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
  colour_matrix <- Seurat::FetchData(reduc, c("umap_1", "umap_2", "umap_3"))

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
#' @importFrom NMF nmf coefficients
#' @return normalised NMF embedding matrix 
embed_nmf <- function(counts, dimensions, verbose = TRUE) {
#   #--------------------------------------------------------------------------#
#   # adding this since I don't want to have this package as a dependancy 
#   # The problem is that those functions are exported to name space
#   # but NMF only works if the package is attached.
#   #--------------------------------------------------------------------------#
#   inst <- requireNamespace("NMF", quietly = TRUE)
#   if (!inst) {
#       stop("NMF is not installed - Please install NMF
#       install.packages('NMF')
#       https://cran.r-project.org/web/packages/NMF/index.html")
    
#   } else {
#       library("NMF")
#   }


  #--------------------------------------------------------------------------#
  # Get the normalized count matrix and matrix with variable features 
  # from the Seurat object
  #--------------------------------------------------------------------------#
  message_switch("nmf_tensor", verbose)
  features <- check_features(counts)
  count_matrix <- as.matrix(Seurat::GetAssayData(counts, layer = "data"))
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
  nmf_projections <- apply(nmf_projections, 2, norm_pixel, "minmax")
  nmf_projections <- list(as.matrix(nmf_projections))
  names(nmf_projections) <- "NMF"
  
  return(nmf_projections)
}

