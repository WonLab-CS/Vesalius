################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Object Utilities/--------------------------------#


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
    if (append && slot != "territories") {
        data <- c(slot(vesalius_assay, slot), data)
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
    if (is.null(names(log))) {
        tag <- "Not_computed"
    } else {
        tag <- grep(x = names(log), pattern = fun, value = TRUE)
    }
    tag <- paste0(assay, "_", create_trial_tag(tag, fun))
    log <- c(log, list(commit))
    names(log)[length(log)] <- tag
    vesalius_assay@log <- log
    return(vesalius_assay)
}

#' get function name from commit list
#' @param commit commit list
#' @return function name
get_func_from_commit <- function(commit) {
    return(commit$fun)
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

#' apply function to slots in object
#' @param x S4 class object
#' @param func closure - function to be applied to slot
#' @param ... any other argument to parse to func
#' @importFrom methods slot slotNames
slot_apply <- function(x, func, ...) {
    cl <- class(x)
    result <- list()
    for (i in slotNames(cl)) {
        result[[i]] <- func(slot(x, i), ...)
    }
    return(result)
}

#' get specific slot out of object if nested in list
#' @param object list containing S4 class objects
#' @param slot_name string name of slot to extract
#' @importFrom methods slot
slot_get <- function(object, slot_name) {
    out <- list()
    for (i in seq_along(object)) {
        out[[i]] <- slot(object[[i]], slot_name)
    }
    return(out)
}



#' assing value to slot in object if nested in list
#' @param object list containing S4 class objects
#' @param sl string - slot name 
#' @param input data that should be assign to slot in S4 object
#' @importFrom methods slot slot<-
slot_assign <- function(object, sl, input) {
    for (i in sl) {
        if (length(slot(object, i)) != 0) {
          slot(object, i) <- c(slot(object, i), input)
        } else {
          slot(object, i) <- input
        }
    }
    return(object)
}



#' create a trail tag name 
#' @param trials character vector containing names of trials 
#' that have already been computed
#' @param tag character string describing which trial tag to add
create_trial_tag <- function(trials, tag) {
    if (is.null(trials) || !any(grepl(x = trials, pattern = tag))) {
        new_trial <-  paste0(tag, "_Trial_1")
    } else {
        previous <- grep(x = trials, pattern = tag, value = TRUE)
        m <- gregexpr("[0-9]+", previous)
        last <- max(as.numeric(unlist(regmatches(previous, m))))
        new_trial <- paste0(tag, "_Trial_", last + 1)
    }
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
get_counts <- function(vesalius_assay, type = "raw") {
    counts <- slot(vesalius_assay, "counts")[[type]]
    return(counts)
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
    deg <- vesalius_assay@DEG %||%
        stop("No DEGs have been computed!")
    if (trial == "last") {
        trial <- tail(names(deg), 1)
        return(deg[[trial]])
    } else {
        in_deg <- grep(pattern = trial, x = names(deg))
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
