
#---------------------/Format conversion Functions/----------------------------#

#' convert vesalius_assay to cimg images 
#' @param vesalius_assay a vesalius_Assay object
#' @param dims integer vector indicating the number of dimensions to select from
#' embeddings.
#' @param embed character indicating which embedding should be selected.
#' Default uses last embedding produced
#' @param verbose logical if progress message should be outputed
#' @return list of cimg images
#' @importFrom dplyr right_join select
#' @importFrom stats median na.exclude
#' @importFrom imager as.cimg
#' @importFrom future.apply future_lapply
format_ves_to_c <- function(vesalius_assay,
  dims,
  embed = "last",
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # This whole thing is too slow
    # for large images this section is always the slowest 
    # the whole function I mean not the sanity checks
    #--------------------------------------------------------------------------#
    embeddings <- check_embedding_selection(vesalius_assay, embed, dims)
    tiles <- check_tiles(vesalius_assay)
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    #--------------------------------------------------------------------------#
    message_switch("vtc", verbose)
    image_list <- future_lapply(seq_along(dims), future_ves_to_cimg,
      embeddings,
      dims,
      tiles,
      future.seed = TRUE)
    return(image_list)
}

#' convert ves embedding to image 
#' @param i index of embedding to use
#' @param embeddings matrix - embedding matrix
#' @param dims dimensions to use
#' @param tiles tile data frame used to reconstruct images
#' @param full_image logical - should the background be returned as well
#' @details using this as a way to run this section in parallel
#' way to slow otherwise. The back ground represents all pixels that are not 
#' part of the Spatial data but constiture the "rest" of the pixels in the image.
#' This tends to happen when you have a non rectangular assay that needs to be fitted
#' into a n * p or n * p * d array.  
#' @return cimg object of embedding
#' @importFrom dplyr right_join select
#' @importFrom stats median na.exclude
#' @importFrom imager as.cimg index.coord
future_ves_to_cimg <- function(i, embeddings, dims, tiles, full_image = TRUE) {
  embeds <- embeddings[, dims[i]]
  embeds <- data.frame(names(embeds), embeds)
  colnames(embeds) <- c("barcodes", as.character(dims[i]))
  cimg <- right_join(tiles, embeds, by = "barcodes")
  colnames(cimg) <- c("barcodes", "x", "y", "origin", "value")
  cimg <- na.exclude(cimg)
  if (!full_image) {
    return(cimg)
  }
  im <- as.cimg(array(median(cimg$value), c(max(cimg$x), max(cimg$y))))
  ind <- imager::index.coord(im, cimg[, c("x", "y"), drop = FALSE])
  im[ind] <- cimg[["value"]]
  return(im)
}

#' convert cimg list to vesalius object embedding
#' @param cimg cimg image list
#' @param vesalius_assay a vesalius_assay object
#' @param dims integer vector of embedding that need to be updated
#' @param embed character string which embedding should be updated
#' @param verbose logical should progress message be outputed
#' @return a vesalius_Assay obejct
#' @importFrom dplyr left_join filter
#' @importFrom stats na.exclude
format_c_to_ves <- function(cimg,
  vesalius_assay,
  dims,
  embed = "last",
  verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # Get stuff out
  # Only need to push the embed selection warning once. It will have been 
  # thrown in the 1st conversion. Makes unit test throw a little tantrum
  # They don't like when more than one warning is pushed
  # Will need to update this at one point 
  # also make it faster
  #--------------------------------------------------------------------------#
  tiles <- check_tiles(vesalius_assay)
  embeddings <- suppressWarnings(check_embedding_selection(
      vesalius_assay, embed, dims))
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  #--------------------------------------------------------------------------#
  message_switch("ctv", verbose)
  for (i in seq_along(dims)) {
    img <- as.data.frame(cimg[[i]])
    barcodes <- left_join(tiles, img, by = c("x", "y")) %>%
      filter(origin == 1) %>%
      na.exclude()
    locs <- match(rownames(embeddings), barcodes$barcodes)
    embeddings[, dims[i]] <- barcodes$value[locs]

  }
  return(embeddings)
}



#' format counts for DESeq
#' @param seed seed count matrix
#' @param query query count matrix
#' @details Always adding a pseudocount of 1 to avoid issues
#' with 0 counts. We also force coercion to int since DESeq does
#' not handle numerics nor does it do internal coersion. 
#' @return DESeq2 object
#' @importFrom DESeq2 DESeqDataSetFromMatrix
format_counts_for_deseq2 <- function(seed, query) {
  seed_tag <- colnames(seed)
  query_tag <- colnames(query)
  merged <- cbind(seed, query) + 1
  mode(merged) <- "integer"
  seed_query_info <- data.frame(row.names = c(seed_tag, query_tag))
  seed_query_info[seed_tag, "group"] <- "seed"
  seed_query_info[query_tag, "group"] <- "query"
  seed_query_info[, "group"] <- factor(x = seed_query_info[, "group"])
  deseq <- DESeq2::DESeqDataSetFromMatrix(
    countData = merged,
    colData = seed_query_info,
    design = ~ group)
  return(deseq)
}

#' format counts for edgeR
#' @param seed seed count matrix
#' @param query query count matrix
#' @return DGEList object from edgeR
#' @importFrom edgeR DGEList
format_counts_for_edger <- function(seed, query) {
  seed <- check_for_zero_counts(seed)
  query <- check_for_zero_counts(query)
  merged <- cbind(seed, query)
  rownames(merged) <- rownames(seed)
  group <- c(rep("seed", ncol(seed)), rep("query", ncol(query)))
  merged <- suppressWarnings(edgeR::DGEList(counts = merged, group = group))
  return(merged)
}

#' format counts for logistic regression
#' @param seed seed count matrix
#' @param query query count matrix
#' @return list of merged counts and meta data formated for logit
format_counts_for_logit <- function(seed, query) {
  seed_tag <- colnames(seed)
  query_tag <- colnames(query)
  merged <- cbind(seed, query)
  seed_query_info <- data.frame(row.names = c(seed_tag, query_tag))
  seed_query_info[seed_tag, "group"] <- "seed"
  seed_query_info[query_tag, "group"] <- "query"
  seed_query_info[, "group"] <- factor(x = seed_query_info[, "group"])
  return(list("merged" = merged, "seed_query_info" = seed_query_info))
}

#' format function call to list for log update
#' @param call a function call argument list
format_call <- function(call) {
  #---------------------------------------------------------------------------#
  # NOTE: if the user put their arguments in "external variable"
  # Only the name of that variable will be parsed not the values themselves
  # could be overcome by writing some function that unwrap and check
  #---------------------------------------------------------------------------#
  for (el in seq_along(call)) {
    tmp <- as.character(call[[el]])
    if(is(call[[el]], "call")) {
      call[[el]] <- tmp[seq(2, length(tmp))]
    } else {
      call[[el]] <- tmp
    }
  }
  names(call) <- c("fun", names(call)[seq(2, length(call))])
  return(as.list(call))
}
