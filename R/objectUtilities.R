################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Object Utilities/--------------------------------#



update_vesalius_assay <- function(vesalius_assay,
    data,
    slot,
    commit,
    defaults,
    assay,
    append = TRUE) {
    #--------------------------------------------------------------------------#
    # Here we get everything in a list format that we can then just append 
    # we can add some checks here to make sure we are formating everything 
    # correctly
    #--------------------------------------------------------------------------#
   
    if (append && slot != "territories") {
        slot(vesalius_assay, slot) <- c(slot(vesalius_assay, slot), data)
    }else if (append && slot == "territories") {
        if (!is.null(slot(vesalius_assay, slot))) {
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

commit_log_to_vesalius_assay <- function(vesalius, commit, assays) {
    for (ass in assays) {
        vesalius@assays[[ass]] <- update_vesalius_assay_log(
            vesalius@assays[[ass]],
            commit = commit[[ass]],
            assay = ass)
    }
    return(vesalius)
}

update_vesalius_assay_log <- function(vesalius_assay, commit, assay) {
    log <- vesalius_assay@log
    if(length(log) == 0) {
        tag <- paste0(assay,"_object_build")
        log <- list(commit)
        names(log) <- tag
    } else {
        tag <- paste0(assay,"_trial_", length(log) - 1)
        log <- c(log, list(commit))
        names(log) <- c(names(log)[-1],tag)
    }
    
    vesalius_assay@log <- log
    return(vesalius_assay)
}




update_vesalius <- function(vesalius,
    data,
    slot,
    commit,
    defaults,
    append = TRUE) {
    #--------------------------------------------------------------------------#
    # First we know that we will always update the assay slot
    # We have already created lists for each assay 
    #--------------------------------------------------------------------------#
    assays <- names(data)
    for (ass in assays) {
        
        vesalius@assays[[ass]] <- update_vesalius_assay(vesalius@assays[[ass]],
            data = data[[ass]],
            slot = slot,
            commit = commit[[ass]],
            defaults = defaults,
            assay = ass,
            append = append)
    }
    return(vesalius)

}


create_commit_log <- function(vesalius, match, default, assay) {
    #--------------------------------------------------------------------------#
    # First lets get the argument that could have more than one value
    # This can be a bit more limited - might need to extend 
    #--------------------------------------------------------------------------#
    match <- format_call(match, assay)
    #--------------------------------------------------------------------------#
    # create list for each assay present - create a seperate "log" for each
    # in many case they will be identical but we want to account for the edge
    # isn't that the guy from U2? 
    #--------------------------------------------------------------------------#
    commit_list <- vector("list", length(assay))
    names(commit_list) <- assay

    for (arg in seq_along(commit_list)) {
        template <-  default
        for (m in seq_along(match)) {
            template[[names(match)[m]]] <- match[[m]][arg]
        }
        commit_list[[arg]] <- template
    }
    return(commit_list)

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
    if (any(grepl(x = trials, pattern = tag))) {
      previous <- grep(x = trials, pattern = tag, value = TRUE)
      m <- gregexpr("[0-9]+", previous)
      last <- max(as.numeric(unlist(regmatches(previous, m))))
      new_trial <- paste0(tag, "_Trial_", last + 1)
    } else {
      new_trial <-  paste0(tag, "_Trial_1")
    }
    return(new_trial)
}

get_assay_names <- function(vesalius) {
    return(names(vesalius@assays))
}



get_counts <- function(vesalius, type = "raw", assay = "all") {
    if (assay != "all") {
        if (!any(assay %in% get_assay_names(vesalius))) {
            stop("Selected assay not in assay list \n")
        } else {
            counts <- vesalius@assays[assay]
        }
    } else {
        counts <- vesalius@assays
        assay <- get_assay_names(vesalius)
    }
    counts <- slot_get(counts, "counts")
    names(counts) <- assay
    return(counts)
}

get_tiles <- function(vesalius, assay = "all") {
    #--------------------------------------------------------------------------#
    # first we check if tiles have been computed 
    # normnally we should have barcode coodinates here but not tiles
    #--------------------------------------------------------------------------#
    if (assay != "all") {
        if (!any(assay %in% get_assay_names(vesalius))) {
            stop("Selected assay not in assay list \n")
        } else {
            tiles <- vesalius@assays[assay]
        }
    } else {
        tiles <- vesalius@assays
        assay <- get_assay_names(vesalius)
    }
    
    tile_log <- check_log(tiles,"build_vesalius_embeddings")
    if (tile_log) {
        return(NULL)
    } else {
        tiles <- slot_get(tiles, "tiles")
        names(tiles) <- assay
        return(tiles)
    }
}


get_log_summary <- function(vesalius) {
    
}

view_trial_summary <- function(vesalius) {
    trials <- lapply(c("Segment", "Territory", "Morphology", "Layer"),
        function(id, vesalius) {
            return(grep(x = colnames(vesalius@territories),
                pattern = id, value = TRUE))
        }, vesalius)
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
