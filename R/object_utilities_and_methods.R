################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Object Utilities/--------------------------------#



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
            df <- data.frame(slot(vesalius_assay, slot), data[, ncol(data)])
            colnames(df)  <- c(colnames(slot(vesalius_assay, slot)),
                colnames(data)[ncol(data)])
            slot(vesalius_assay, slot) <- df
        } else {
            slot(vesalius_assay, slot) <- data
        }
    } else {
        slot(vesalius_assay, slot) <- data
    }
    return(vesalius_assay)
}



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


get_func_from_commit <- function(commit) {
    return(commit$fun)
}


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


# This function will do with some more options 
# it would be neat to be able to go through the log 
# not a priority as it would be better for them to keep track of that 
# in their won code
search_log <- function(vesalius_assay, arg, return_assay = TRUE) {
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


slot_apply <- function(x, func, ...) {
    cl <- class(x)
    result <- list()
    for (i in slotNames(cl)) {
        result[[i]] <- func(slot(x, i), ...)
    }
    return(result)
}

slot_get <- function(object, slot_name) {
    out <- list()
    for (i in seq_along(object)) {
        out[[i]] <- slot(object[[i]], slot_name)
    }
    return(out)
}




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

get_last_log <- function(vesalius) {
  cl <- class(vesalius@log)
  log <- list()
  for (i in slotNames(cl)) {
      tmp <- slot(vesalius@log, i)
      if (length(tmp) > 0) {
        log[[i]] <- tmp[[length(tmp)]]
      } else {
        log[[i]] <- NULL
      }

  }
  log <- log[!is.null(log)]
  return(log)
}



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



get_assay_names <- function(vesalius_assay) {
    return(vesalius_assay@assay)
}




get_counts <- function(vesalius_assay, type = "raw") {
    counts <- slot(vesalius_assay, "counts")[[type]]
    return(counts)
}



get_tiles <- function(vesalius_assay) {
    tiles <- slot(vesalius_assay, "tiles")
    return(tiles)
}



view_trial_summary <- function(vesalius_assay) {
    trials <- lapply(c("Segment", "Territory", "Morphology", "Layer"),
        function(id, vesalius_assay) {
            return(grep(x = colnames(vesalius_assay@territories),
                pattern = id, value = TRUE))
        }, vesalius_assay)
    max_trials <- max(sapply(trials, length))
    if (max_trials == 0) {
        stop("No Territory Trials to be found!")
    } else {
        trials <- lapply(trials,
            function(trial, max) {
                if (max - length(trial) == 0) {
                    return(trial)
                } else {
                    return(c(trial, rep("-", times = max - length(trial))))
                }
            }, max = max_trials)
        trials <- do.call("cbind", trials)
        colnames(trials) <- c("Segment", "Territory", "Morphology", "Layer")
    }
    return(trials)
}
