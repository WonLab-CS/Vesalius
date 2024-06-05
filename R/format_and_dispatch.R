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
        territories$trial[territories$trial %in% ter_1 &
            territories$barcodes %in% cells] <- paste0(ter_1, collapse = " ")
        territories$trial[territories$trial %in% ter_2 &
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
#' @param sample barcodes for each sample type
#' @param verbose logical if progress messages should be outputed
#' @details This function generates groups for DEG analysis. It creates 
#' different groups depending on what is being parse to both seed and query. 
#' This could probably be cleaned up and simplified. Always need to return 
#' a list though.
#' @return list with seed group and seed id as well as query group and query id 
dispatch_deg_group <- function(ter, seed, query, cells, sample, verbose) {
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
        seed <- ter[ter$trial %in% seed, "barcodes"]
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        query_id <- "remaining"
        query <- ter[!ter$barcodes %in% seed, "barcodes"]
        #----------------------------------------------------------------------#
        # reformat 
        #----------------------------------------------------------------------#
        seed <- list(seed)
        query <- list(query)
    } else if (is.null(seed) && !is.null(query)) {
        #----------------------------------------------------------------------#
        # if only query we compare query to everything else
        # Get initial query territory
        #----------------------------------------------------------------------#
        query_id <- paste0(query, collapse = " ")
        query <- ter[ter$trial %in% query, "barcodes"]
        #----------------------------------------------------------------------#
        # Filter seed based on query
        #----------------------------------------------------------------------#
        seed_id <- "remaning"
        seed <- ter[!ter$barcodes %in% query, "barcodes"]
        #----------------------------------------------------------------------#
        # reformat 
        #----------------------------------------------------------------------#
        seed <- list(seed)
        query <- list(query)
    } else {
        #----------------------------------------------------------------------#
        # if get both filter based on both
        #----------------------------------------------------------------------#
        seed_id <- paste0(seed, collapse = " ")
        seed <- ter[ter$trial %in% seed, "barcodes"]
        query_id <- paste0(query, collapse = " ")
        query <- ter[ter$trial %in% query, "barcodes"]
        #----------------------------------------------------------------------#
        # reformat 
        #----------------------------------------------------------------------#
        seed <- list(seed)
        query <- list(query)
    }
    if (!is.null(cells)) {
      seed <-  mapply(check_cells, territory_barcodes = seed,
        ter = seed_id, MoreArgs = list(cells), SIMPLIFY = FALSE)
      query <- mapply(check_cells, territory_barcodes = query,
        ter = query_id, MoreArgs = list(cells), SIMPLIFY = FALSE)
    }
    if (!is.null(sample)) {
        seed <-  mapply(dispatch_sample, territory_barcodes = seed,
            ter = seed_id, MoreArgs = list(sample$matched), SIMPLIFY = FALSE)
        seed_id <- paste0(seed_id,"_matched")
        query <- mapply(dispatch_sample, territory_barcodes = query,
            ter = query_id, MoreArgs = list(sample$reference), SIMPLIFY = FALSE)
        query_id <- paste0(query_id,"_reference")
    }
    return(list("seed" = seed, "seed_id" = seed_id,
        "query" = query, "query_id" = query_id))
}



dispatch_sample <- function(territory_barcodes, ter, sample) {
    common <- intersect(territory_barcodes, sample)
    if (length(common) == 0) {
        warning(paste0("No overlap between sample and territory in ", ter,
            "\n Returning NULL\n"))
        common <- NULL
    }
    return(common)
}

#' convert cost to prob
#' @param cost matrix containing cost 
#' @param n_cost number of custom cost matrices parsed
#' @return a probability matrix 
cost_to_prob <- function(cost, n_cost) {
    cost$score <- (n_cost - cost$score) / n_cost
    return(cost)
}


#' dispatch barcodes to subset cost for match clustering
#' @param vesalius_assay vesalius_assay object post cell mapping
#' @param cell_label character - name of column containing cell names
#' if the clustering is to be done by cell types only
#' @param group_identity character - name of column containing group
#' to be used for cluster (i.e. territories, segments, layers etc)
#' @return character string of barcodes
dispatch_cost_groups <- function(vesalius_assay,
    cost,
    trial = NULL,
    group_identity = NULL,
    ref_cells = NULL,
    query_cells = NULL) {
    #-------------------------------------------------------------------------#
    # Get trial and filter sub categories if needed
    # this could added to the sanity checks 
    #-------------------------------------------------------------------------#
    if (!is.null(trial)) {
        trial <- check_cell_labels(vesalius_assay, trial = trial)
        if(!is.null(group_identity)) {
            group_identity <- check_group_value(trial, group_identity)
        }
        trial <- trial$barcodes[trial$trial %in% group_identity]

    } else {
        trial <- get_coordinates(vesalius_assay)$barcodes
    }
    #-------------------------------------------------------------------------#
    # check ref cell abd subset
    #-------------------------------------------------------------------------#
    if (is.null(ref_cells)) {
        ref_cells <- colnames(cost)
    } else {
        ref_cells <- intersect(ref_cells, colnames(cost))
    }
    if (length(ref_cells) == 0){
        stop("No overlap between reference cells provided and cells present in cost matrices")
    }
    #-------------------------------------------------------------------------#
    # check query cells and subset
    #-------------------------------------------------------------------------#
     trial <- clean_trial(trial)
    if (is.null(query_cells)) {
        query_cells <- trial
    } else {
        query_cells <- intersect(query_cells, trial)
    }
    if (length(query_cells) == 0){
        stop("No overlap between query cells provided and cells present in cost matrices")
    }

    return(list("ref" = ref_cells, "query" = query_cells))

}

clean_trial <- function(trial) {
    trial <- sapply(strsplit(trial, "-"), function(str){
        return(paste0(str[seq(1, length(str) - 1)], collapse = "-"))
    })
    return(trial)
}