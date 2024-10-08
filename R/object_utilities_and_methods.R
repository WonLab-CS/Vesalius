################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Object Utilities/--------------------------------#


#' adjust count
#'
#' adjust counts after reducing the resolution of the image tensor
#' or after filtering stray beads
#' @param coordinates data frame containing coordinates after reducing 
#' resolution and compressing cooridnates
#' @param counts count matrix
#' @param throw logical - throwing warning message for unshared barcodes
#' @param verbose logical if progress messaged should be outputed.
#' @details This function will check the coordinate file to 
#' see if any barcodes have been merged together. If so,
#' the counts will be adjusted by taking the average count value accross 
#' all barcodes that have been merged together. 
#' @return a count matrix with adjusted count values 
#' @importFrom Matrix Matrix rowMeans
#' @importFrom future.apply future_lapply
adjust_counts <- function(coordinates, counts, throw = TRUE, verbose = TRUE) {
    message_switch("adj_counts", verbose)
    #--------------------------------------------------------------------------#
    # First get all barcode names and compare which ones are missing
    #--------------------------------------------------------------------------#
    coord_bar_uni <- unique(coordinates$barcodes)
    coord_bar <- grep(pattern = "_et_", x = coord_bar_uni, value = TRUE)
    if (length(coord_bar) == 0) {
      loc <- check_barcodes(colnames(counts), coord_bar_uni, throw)
      return(counts[, loc])
    }

    #--------------------------------------------------------------------------#
    # next we merge counts together when barcodes have been merged
    # and take the average count value in that barcode location.
    #--------------------------------------------------------------------------#
    tmp_bar <- strsplit(coord_bar, "_et_")

    empty <- future_lapply(tmp_bar, function(tmp_bar, counts) {
        tmp_bar <- which(colnames(counts) %in% tmp_bar)
        return(Matrix::rowMeans(counts[, tmp_bar]))
    }, counts = counts, future.seed = TRUE)

    empty <- do.call("cbind", empty)
    if (is.null(dim(empty)) && length(empty) != 0) {
        empty <- Matrix::Matrix(empty, ncol = 1)
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


#' update vesalius assay object
#' @param vesalius_assay a vesalius_assay object
#' @param data data that will be inserted into the vesalius_assay
#' @param slot name of the slot that will be updated
#' @param append logical if the dats hould be appended to already existing data
#' @return a vesalius_assay
#' @importFrom methods slot slot<-
update_vesalius_assay <- function(vesalius_assay,
    data,
    slot,
    append = TRUE) {
    #--------------------------------------------------------------------------#
    # First we do some checks
    #--------------------------------------------------------------------------#
    if (append && slot != "territories" && slot != "map") {
        internal <- slot(vesalius_assay, slot)
        tag <- create_trial_tag(names(internal), names(data))
        data <- c(internal, data)
        names(data) <- tag
        slot(vesalius_assay, slot) <- data
    }else if (append && slot == "territories") {
        if (!all(dim(slot(vesalius_assay, slot)) == c(0, 0))) {
            df <-  slot(vesalius_assay, slot)
            data <- data[match(df$barcodes, data$barcodes), ]
            d <- data.frame(df, data[, ncol(data)])
            colnames(d)  <- c(colnames(df), colnames(data)[ncol(data)])
            slot(vesalius_assay, slot) <- d
        } else {
            slot(vesalius_assay, slot) <- data
        }

    } else if (append && slot == "map") {
        if (!all(dim(slot(vesalius_assay, slot)) == c(0, 0))) {
            df <-  slot(vesalius_assay, slot)
            data <- data[match(df[,unique(df$init)], data[,unique(data$init)]), ]
            d <- data.frame(df, data[, ncol(data)])
            colnames(d)  <- c(colnames(df), colnames(data)[ncol(data)])
            slot(vesalius_assay, slot) <- d
        } else {
            slot(vesalius_assay, slot) <- data
        }
    } else {
        slot(vesalius_assay, slot) <- data
    }
    return(vesalius_assay)
}

#' create function log to be commit to the log slot
#' @param arg_match function argument call 
#' @param default default argument values of function
#' @return list with all arguments used in latest function call
create_commit_log <- function(arg_match, default) {
    #--------------------------------------------------------------------------#
    # First lets get the argument that could have more than one value
    # This can be a bit more limited - might need to extend 
    #--------------------------------------------------------------------------#
    arg_match <- format_call(arg_match)
    #--------------------------------------------------------------------------#
    # create list for each assay present - create a seperate "log" for each
    # in many case they will be identical but we want to account for the edge
    # isn't that the guy from U2? 
    #--------------------------------------------------------------------------#
    template <-  default
    for (m in seq_along(arg_match)) {
        template[[names(arg_match)[m]]] <- arg_match[[m]]
    }
    return(template)

}

#' commits function call to log slot 
#' @param vesalius_assay a vesalius_assay object
#' @param commit list containing arguments used in latest function call
#' @param assay assay name used in latest function call
#' @return vesalius_assay with updated log slot
commit_log <- function(vesalius_assay, commit, assay) {
    fun <- get_func_from_commit(commit)
    log <- vesalius_assay@log
    tag <- paste0(assay, "_", fun)
    tag <- create_trial_tag(names(log), tag)
    log <- c(log, list(commit))
    names(log) <- tag
    vesalius_assay@log <- log
    return(vesalius_assay)
}

#' get function name from commit list
#' @param commit commit list
#' @details Using tail since if you make an explicit function call
#' using pkg::func you get pkg as well. Neat. We only want the function
#' @return function name
#' @importFrom utils tail 
get_func_from_commit <- function(commit) {
    return(tail(commit$fun, 1))
}




# This function will do with some more options 
# it would be neat to be able to go through the log 
# not a priority as it would be better for them to keep track of that 
# in their won code

#' search through log for parameter values or names 
#' @param vesalius_assay a vesalius_assay object
#' @param arg string indicating which parameter value or argument name
#' should be searched for
#' @param return_assay logical indicating if the log list
#' should returned or only the value.
#' @details You may search through the log to see if you have used
#' certain parameters and if so which values did you use when running a
#' certain trial.
#' If `return_assay` is `TRUE` then the entire log call will be returned.
#' This will include all parameter values including defaults parse to 
#' function.
#' @return either a list containing log calls or values found in call
#' @export
#' @rdname search_log
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run 
#' ves <- build_vesalius_embeddings(ves)
#' # maybe we want to try a different method 
#' # both will be stored in the object
#' ves <- build_vesalius_embeddings(ves, dim_reduction = "UMAP")
#' 
#' # search log 
#' search_log(ves, "tensor_resolution")
#'}
search_log <- function(vesalius_assay,
    arg,
    return_assay = TRUE) {
    log <- vesalius_assay@log
    locs <- lapply(log, function(x, arg) {
        arg_names <- grepl(pattern = arg, x = x)
        arg_value <- grepl(pattern = arg, x = names(x))
        return(any((arg_names | arg_value)))
    }, arg = arg)
    if (return_assay) {
        return(log[unlist(locs)])
    } else {
        return(unlist(locs))
    }
}


#' create a trail tag name 
#' @param trials character vector containing names of trials 
#' that have already been computed
#' @param tag character string describing which trial tag to add
create_trial_tag <- function(trials, tag) {
    new_trial <- make.names(c(trials, tag), unique = TRUE)
    return(new_trial)
}



#' get assay names from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @return assay names
#' @rdname get_assay_names
#' @export
get_assay_names <- function(vesalius_assay) {
    return(vesalius_assay@assay)
}



#' get counts from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @param type character string for which count matrix to
#' retrieve.
#' @return count matrix
#' @rdname get_counts
#' @export
#' @importFrom methods slot
#' @importFrom utils tail head
get_counts <- function(vesalius_assay, type = "raw") {
    counts <- slot(vesalius_assay, "counts")
    #--------------------------------------------------------------------------#
    # If the user does not parse any counts - there is not slot named raw 
    # Instead it is an empty list that is tagged with an "empty" comment
    # We want to check if that is the case and throw an error
    #--------------------------------------------------------------------------#
    if (type == "all") {
        return(counts)
    } else {
        loc <- grep(type, names(counts), value = TRUE)
        if (length(loc) == 0) {
            stop(paste(type, "is not present in count matrix list"))
        } else if (length(loc) > 1) {
            warning(paste("More than 1 count matrix called", type,
                "Vesalius will return first instance"))
            counts <- counts[[head(loc, 1)]]
        } else {
            counts <- counts[[loc]]
        }
        return(counts)
    }
}


#' get tiles from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @return tiles data frame
#' @rdname get_tiles
#' @export
#' @importFrom methods slot
get_tiles <- function(vesalius_assay) {
    tiles <- slot(vesalius_assay, "tiles")
    return(tiles)
}


#' get coordinates from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @return coordinate data frame
#' @rdname get_coordinayes
#' @export
#' @importFrom methods slot
get_coordinates <- function(vesalius_assay, original = FALSE) {
    tiles <- vesalius_assay@meta$orig_coord
    if (!is.null(tiles) && original){
        return(tiles)    
    } else {
        tiles <- slot(vesalius_assay, "tiles")
        if ("origin" %in% colnames(tiles)) {
            tiles <- tiles %>% filter(origin == 1)
        } 
         return(tiles)
    }   
}

#' get embeddings from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @param active logical if active embedding should be return 
#' or full embedding list.
#' @return embedding matrix
#' @rdname get_embeddings
#' @export
#' @importFrom methods slot
get_embeddings <- function(vesalius_assay, active = TRUE) {
    if (sum(dim(vesalius_assay@active)) == 0) {
            stop("No embeddings have been computed!")
        }
    if (active) {
        tiles <- slot(vesalius_assay, "active")
    } else {
        tiles <- slot(vesalius_assay, "embeddings")
    }
    return(tiles)
}

#' get territories from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @return territories data frame
#' @rdname get_territories
#' @export
#' @importFrom methods slot
get_territories <- function(vesalius_assay) {
    territories <- slot(vesalius_assay, "territories")
    return(territories)
}

#' get markers from vesalius_assay
#' @param vesalius_assay a vesalius_assay
#' @param trial character string describing name of DEG trial 
#' to retrieve.
#' @return marker data frame
#' @rdname get_markers
#' @export
get_markers <- function(vesalius_assay, trial = "last") {
    if (length(vesalius_assay@DEG) == 0) {
         stop("No DEGs have been computed!")
    }
    deg <- vesalius_assay@DEG
    if (trial == "last") {
        trial <- tail(names(deg), 1)
        return(deg[[trial]])
    } else {
        in_deg <- grep(pattern = paste0("^", trial, "$"), x = names(deg))
        if (length(in_deg) == 0) {
            stop(paste(deparse(substitute(trial)),
                ": Unknown embedding selected!"))
        } else if (length(in_deg) > 1) {
            warning(paste("More than 1", deparse(substitute(trial)), "embedding.
            Vesalius will return DEG list"))
            trial <- names(deg)[in_deg]
            return(deg[trial])
        } else {
            trial <- names(deg)[in_deg]
            return(deg[[trial]])
        }
    }
}

#' @export
get_mapping_scores <- function(vesalius_assay) {
    coord <- get_coordinates(vesalius_assay)
    scores <- vesalius_assay@meta$mapping_scores
    scores <- scores[match(coord$barcodes, scores$from), ]
    mapped <- data.frame("barcodes" = scores$from,
        "ref_barcodes" = scores$to,
        "x" = coord$x,
        "y" = coord$y)
    mapped <- cbind(mapped,
        scores[, seq(grep("cost", colnames(scores)), ncol(scores))])
    return(mapped)
}

#' get active embedding
#' @param vesalius_assay a vesalius_assay object
#' @return character string with name of last embedding used
#' @export
get_active_embedding_tag <- function(vesalius_assay) {
    last <- comment(get_embeddings(vesalius_assay, active = FALSE))
    return(last)
}

#' get last count matrix used
#' @param vesalius_assay a vesalius_assay object
#' @return character string with name of last embedding used
#' @export
get_active_count_tag <- function(vesalius_assay) {
    last <- comment(get_counts(vesalius_assay, type = "all"))
    return(last)
}

#' add active embedding tag
#' @param vesalius_assay a vesalius assay object
#' @param embedding embedding
#' @return commented list with active embedding tag
#' @importFrom dplyr %>%
#' @importFrom utils tail
add_active_embedding_tag <- function(vesalius_assay, embedding) {
   if (embedding != "last") {
        active <- names(vesalius_assay@embeddings)
        active <- grep(embedding, active, value = TRUE) %>%
            tail(1)
        comment(vesalius_assay@embeddings) <- active
   }
    return(vesalius_assay)
}

add_integration_tag <- function(vesalius_assay, integrated) {
    comment(vesalius_assay@active) <- integrated
    return(vesalius_assay)
}

#' add active count tag
#' @param vesalius_assay a vesalius assay object
#' @param norm requested count matrix
#' @return commented list with active embedding tag
#' @importFrom dplyr %>%
#' @importFrom utils tail
add_active_count_tag <- function(vesalius_assay, norm) {
    if (norm != "last" && norm != "joint") {
        active <- names(vesalius_assay@counts)
        active <- grep(norm, active, value = TRUE) %>%
            tail(1)
        comment(vesalius_assay@counts) <- active
    } else if (norm == "joint") {
        comment(vesalius_assay@counts) <- norm
    }
    return(vesalius_assay)
}

#' summarise territories
#' @param vesalius_assay a vesalius assay object
#' @param as_log logical defining if log list should be returned
#' @return a list containing summary of territory transformation
#' and manipulations. 
#' @export
summarise_territories <- function(vesalius_assay, as_log = FALSE) {
    if (as_log) {
        log <- vector("list", 4)
        names(log) <- c("Segmentation",
            "Territory",
            "Morphing",
            "Layering")
        log[["Segmentation"]] <- search_log(vesalius_assay,
            "image_segmentation")
        log[["Territory"]] <- search_log(vesalius_assay,
            "isolate_territories")
        log[["Morphing"]] <- search_log(vesalius_assay,
            "territory_morphing")
        log[["Layering"]] <- search_log(vesalius_assay,
            "layer_territory")
    } else {
        log <- vector("list", 4)
        territories <- get_territories(vesalius_assay)
        names(log) <- c("Segmentation",
            "Territory",
            "Morphing",
            "Layering")
        log[["Segmentation"]] <- length(grep("Segment", colnames(territories)))
        log[["Territory"]] <- length(grep("Territory", colnames(territories)))
        log[["Morphing"]] <- length(grep("Morphology", colnames(territories)))
        log[["Layering"]] <- length(grep("Layer", colnames(territories)))
    }
    return(log)
}



#' @importFrom dplyr %>% filter
#' @export
filter_assay <- function(vesalius_assay,
    cells = NULL,
    territories = NULL,
    trial = "last",
    genes = NULL) {
    if (!is.null(territories)) {
        ter_keep <- check_territory_trial(vesalius_assay, trial)
        ter_keep <- ter_keep$barcodes[ter_keep$trial %in% territories]
    } else {
        ter_keep <- c()
    }
    if (!is.null(cells)) {
        cell_keep <- get_coordinates(vesalius_assay)
        cell_keep <- cell_keep$barcodes[cell_keep$barcodes %in% cells]
    } else {
        cell_keep <- c()
    }
    if (!is.null(genes)) {
        gene_keep <- rownames(get_counts(vesalius_assay))
        gene_keep <- gene_keep[gene_keep %in% genes]
    } else {
        gene_keep <- c()
    }
    ### NOTE this is buggy!!! Will not hold to edge cases
    ### Need to fix for long term use!!!!
    keep_barcodes <- unique(c(ter_keep, cell_keep))
    if (length(keep_barcodes) == 0){
        stop("No Barcodes present under current subset conditions!")
    } else {
        coord <- vesalius_assay@meta$orig_coord
        coord <- coord[coord$barcodes %in% keep_barcodes,
            c("barcodes","x_orig","y_orig")]
        colnames(coord) <- c("barcodes","x" ,"y")
        counts <- get_counts(vesalius_assay, "raw")
        counts <- counts[, colnames(counts) %in% keep_barcodes]
    }
    if (length(gene_keep) == 0 && !is.null(genes)) {
        stop("Requested genes are not in count matrix")
    } else if (length(gene_keep) > 0 && !is.null(genes)){
        counts <- counts[rownames(counts) %in% gene_keep, ]
        
    }
    if (length(vesalius_assay@image)  > 0) {
        image <- vesalius_assay@image[[1]]
        scale <- vesalius_assay@meta$scale$scale
    } else {
        image <- NULL
        scale <- "auto"
    }
    vesalius_assay <- build_vesalius_assay(coord,
        counts,
        image = image,
        scale = scale,
        verbose = FALSE)
    return(vesalius_assay)
}


get_cost <- function(vesalius_assay, use_cost = NULL) {
    cost <- vesalius_assay@cost$cost
    if (length(cost) == 0) {
        stop("No cost has been computed")
    }
    
    if (!is.null(use_cost)) {
        if (length(use_cost) >1 & grepl("cost", use_cost)) {
            stop("Cannot use cost with other metrics for clustering!")
        }
        cost <- cost[use_cost]
        if (length(cost) == 0){
            stop(paste(paste(use_cost, collapse = " "), ": not available in score matrix list"))
        }
    }
    return(cost)
}