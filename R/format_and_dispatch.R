################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Format conversion Functions/----------------------------#

#' convert vesalius_assay to cimg images 
#' @param vesalius_assay a vesalius_Assay object
#' @param dims integer vector indicating the number of dimensions to select from
#' embeddings.
#' @param embed character indicating which embedding should be selected.
#' Default uses last embedding produced
#' @param correct_background logical indicating if image background should be
#' corrected to reflect median color value or left as black
#' @param verbose logical if progress message should be outputed
#' @return list of cimg images
#' @importFrom dplyr right_join select
#' @importFrom stats median na.exclude
#' @importFrom imager as.cimg
format_ves_to_c <- function(vesalius_assay,
  dims,
  embed = "last",
  correct_background = TRUE,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # For this function we will not create a "method". There are too many
    # other options that we need to account for.
    # Might update this with medici later
    #--------------------------------------------------------------------------#
    embeddings <- check_embedding(vesalius_assay, embed, dims)
    tiles <- check_tiles(vesalius_assay)
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    #--------------------------------------------------------------------------#
    image_list <- list()
    message_switch("vtc", verbose)
    for (i in seq_along(dims)) {
      #------------------------------------------------------------------------#
      # First we create a cimg data frame template from embedding values
      #------------------------------------------------------------------------#
      embeds <- embeddings[, dims[i]]
      embeds <- data.frame(names(embeds), embeds)
      colnames(embeds) <- c("barcodes", as.character(dims[i]))
      cimg <- right_join(tiles, embeds, by = "barcodes")
      colnames(cimg) <- c("barcodes", "x", "y", "origin", "value")
      cimg <- na.exclude(cimg)
      x <- max(cimg$x)
      y <- max(cimg$y)
      #------------------------------------------------------------------------#
      # now we can convert that data frame to a cimg and correct background if
      # required. Note that we are going back and forth between formats
      # This is done so we can fill in the borders and empty space in the array
      # Also need to be done step by step. Using pipe does not work
      #------------------------------------------------------------------------#
      if (correct_background) {
        cimg_tmp <- cimg %>%
          select(c("x", "y", "value"))
        cimg_tmp <- suppressWarnings(as.cimg(cimg_tmp, x = x, y = y))
        cimg_tmp <- as.data.frame(cimg_tmp)
        non_img <- paste0(cimg_tmp$x, "_", cimg_tmp$y)
        in_img <- paste0(cimg$x, "_", cimg$y)
        # Median value for background - Other metrics??
        cimg_tmp[!non_img %in% in_img, "value"] <- median(cimg$value)
        image_list[[i]] <- suppressWarnings(
          as.cimg(cimg_tmp[, c("x", "y", "value")], x = x, y = y))
      } else {
        image_list[[i]] <- suppressWarnings(
          as.cimg(cimg[, c("x", "y", "value")], x = x, y = y))
      }
  }
  return(image_list)
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
  #--------------------------------------------------------------------------#
  tiles <- check_tiles(vesalius_assay)
  embeddings <- check_embedding(vesalius_assay, embed, dims)
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
#' @return DESeq2 object
#' @importFrom DESeq2 DESeqDataSetFromMatrix
format_counts_for_deseq2 <- function(seed, query) {
  seed_tag <- colnames(seed)
  query_tag <- colnames(query)
  merged <- rbind(seed, query)
  seed_query_info <- data.frame(row.names = c(seed_tag, query_tag))
  seed_query_info[seed_tag, "group"] <- "seed"
  seed_query_info[query_tag, "group"] <- "query"
  seed_query_info[, "group"] <- factor(x = seed_query_info[, "group"])
  deseq <- DESeq2::DESeqDataSetFromMatrix(
    countData = merged,
    colData = seed_query_info,
    design = ~ group
  )
  return(deseq)
}

#' format counts for edgeR
#' @param seed seed count matrix
#' @param query query count matrix
#' @return DGEList object from edgeR
#' @importFrom edgeR DGEList
format_counts_for_edger <- function(seed, query) {
  merged <- cbind(seed, query)
  rownames(merged) <- rownames(seed)
  group <- c(rep("seed", ncol(seed)), rep("query", ncol(query)))
  merged <- edgeR::DGEList(counts = merged, group = group)
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


#' dispatch territories labels territories according to which
#' group they belong to.
#' @param territories territories data frame from territoires slot in 
#' vesalius_assay object
#' @param ter_1 integer vector containing territories in group 1
#' @param ter_2 integer vector containing territories in group 2
#' @param cells cell barcodes
dispatch_territory <- function(territories, ter_1, ter_2, cells) {
    if (is.null(ter_1) && is.null(ter_2)) {
        territories <- select(territories, c("barcodes", "x", "y", "trial"))
    }else if (!is.null(ter_1) && is.null(ter_2)) {
        territories$trial[!territories$trial %in% ter_1] <- "other"
    }else if (is.null(ter_1) && !is.null(ter_2)) {
        territories$trial[!territories$trial %in% ter_2] <- "other"
    }else {
        territories$trial[!territories$trial %in% c(ter_1, ter_2)] <- "other"
    }
    if (!is.null(cells)) {
        territories$trial[territories$trial %in% ter_1 &&
            territories$barcodes %in% cells] <- paste0(ter_1, collapse = " ")
        territories$trial[territories$trial %in% ter_2 &&
            territories$barcodes %in% cells] <- paste0(ter_2, collapse = " ")
    }
    return(territories)
}


#' dispatch barcodes to seed and query groups
#' @param ter territories data frame from territoires slot in
#' vesalius_assay object
#' @param seed interger vector indicating which territories should be included
#' in seed group
#' @param query interger vector indicating which territories should be included
#' in query group
#' @param cells cell barcodes
#' @param verbose logical if progress messages should be outputed
#' @return list with seed group and seed id as well as query group and query id 
dispatch_deg_group <- function(ter, seed, query, cells, verbose) {
    if (is.null(seed) && is.null(query)) {
        #----------------------------------------------------------------------#
        # If no territories are provided - then we assume that the user
        # want to look at all territories.
        # This compares each territory to everything else
        #----------------------------------------------------------------------#
        message_switch("deg_dispatch_all_null", verbose)
        seed <- split(ter$barcodes, ter$trial)
        seed_id <- names(seed)

        query <- lapply(seed, function(bar, ter) {
            return(ter$barcodes[!ter$barcodes %in% bar])
        }, ter = ter)
        query_id <- rep("remaining", length(query))
    } else if (!is.null(seed) && is.null(query)) {
        #----------------------------------------------------------------------#
        # if only seed we compare seed to everything else
        # Get initial seed territory
        #----------------------------------------------------------------------#
        seed_id <- paste0(seed, collapse = " ")
        seed <- list(ter[ter$trial %in% seed, "barcodes"])
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        query_id <- "remaining"
        query <- list(ter[!ter$barcodes %in% seed, "barcodes"])
    } else if (is.null(seed) && !is.null(query)) {
        #----------------------------------------------------------------------#
        # if only query we compare query to everything else
        # Get initial query territory
        #----------------------------------------------------------------------#
        query_id <- paste0(query, collapse = " ")
        query <- list(ter[ter$trial %in% query, "barcodes"])
        #----------------------------------------------------------------------#
        # Filter seed based on query
        #----------------------------------------------------------------------#
        seed_id <- "remaning"
        seed <- list(ter[!ter$barcodes %in% query, "barcodes"])
    } else {
        #----------------------------------------------------------------------#
        # if get both filter based on both
        #----------------------------------------------------------------------#
        seed_id <- paste0(seed, collapse = " ")
        seed <- list(ter[ter$trial %in% seed, "barcodes"])
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        query_id <- paste0(query, collapse = " ")
        query <- list(ter[ter$trial %in% query, "barcodes"])
    }
    if (!is.null(cells)) {
      seed <-  mapply(check_cells, territory_barcodes = seed,
        ter = seed_id, MoreArgs = list(cells, verbose), SIMPLIFY = FALSE)
      query <- mapply(check_cells, territory_barcodes = query,
        ter = query_id, MoreArgs = list(cells, verbose), SIMPLIFY = FALSE)
    }
    return(list("seed" = seed, "seed_id" = seed_id,
        "query" = query, "query_id" = query_id))
}