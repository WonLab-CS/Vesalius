################################################################################
###############################   Vesalius      ################################
################################################################################


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
            ter = seed_id, MoreArgs = list(sample$reference), SIMPLIFY = FALSE)
        seed_id <- paste0(seed_id,"_reference")
        query <- mapply(dispatch_sample, territory_barcodes = query,
            ter = query_id, MoreArgs = list(sample$matched), SIMPLIFY = FALSE)
        query_id <- paste0(query_id,"_matched")
    }
    return(list("seed" = seed, "seed_id" = seed_id,
        "query" = query, "query_id" = query_id))
}


#' dispatch barcodes associted with a specific territory in a specific samples
#' @param territory_barcodes character vector barcodes associated with territory
#' @param ter chracter territory id  
#' @param sample character vector associeted with sample
#' @return character vector of barcodes at the intersection of territory and sample
dispatch_sample <- function(territory_barcodes, ter, sample) {
    common <- intersect(territory_barcodes, sample)
    if (length(common) == 0) {
        warning(paste0("No overlap between sample and territory in ", ter,
            "\n Returning NULL\n"))
        common <- NULL
    }
    return(common)
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