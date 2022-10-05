################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Object Utilities/--------------------------------#


update_vesalius <- function(vesalius,
    data,
    slot,
    commit,
    defaults,
    append = TRUE) {
    #--------------------------------------------------------------------------#
    # You take a vesalius object, the data you want add to it and which slot
    # it should be added to.
    # First step is to check if the slot we want to update is empty
    # This is not the same as log commit
    # Lets keep those seperate
    # Will need some sanitu checks here at one point
    #--------------------------------------------------------------------------#

    if (append && slot != "territories") {
        slot(vesalius, slot) <- c(slot(vesalius, slot), data)
        vesalius <- commit_log(vesalius, commit, defaults, slot)


    }else if (append && slot == "territories") {
        if (!is.null(slot(vesalius, slot))) {
            df <- data.frame(slot(vesalius, slot), data[, ncol(data)])
            colnames(df)  <- c(colnames(slot(vesalius, slot)),
                colnames(data)[ncol(data)])
            slot(vesalius, slot) <- df
        } else {
            slot(vesalius, slot) <- data
        }
        vesalius <- commit_log(vesalius, commit, defaults, slot)
    } else {
        slot(vesalius, slot) <- data
        vesalius <- commit_log(vesalius, commit, defaults, slot)
    }
    return(vesalius)
}



commit_log <- function(vesalius, commit, defaults, slot, append = TRUE) {
    #--------------------------------------------------------------------------#
    # First let's get all the non symbol arguments
    # symbol argument are argument that require an input i.e no default
    # Trim the last entry as it is always NULL... or is it????
    #--------------------------------------------------------------------------#

    defaults <- defaults[-length(defaults)]
    non_sym <- sapply(defaults, function(x) {
                  return(!any(is(x) == "refObject"))
    })
    defaults <- defaults[non_sym]
    #--------------------------------------------------------------------------#
    # Initialise log
    #--------------------------------------------------------------------------#
    defs <- sapply(defaults, function(x) {
        if (class(x) == "call") {
            x <- as.character(x)[2L]
        } else if (is.null(x)) {
            x <- "NULL"
        }
        return(unlist(x))
    })

    logdf <- data.frame("Argument" = c("Function", names(defs)),
        "Value" = c(as.character(commit[[1]]), unname(defs)))

    #--------------------------------------------------------------------------#
    # Update the default df based on mathc call output
    #--------------------------------------------------------------------------#
    if (length(commit) > 1) {
      for (i in seq(2, length(commit))) {
          if (length(as.character(commit[[i]])) > 1) {
              com <- paste0(as.character(commit[[i]]), collapse = "_")
          } else {
              com <- as.character(commit[[i]])
          }
          logdf$Value[logdf$Argument == names(commit)[i]] <- com
      }
    }
    #--------------------------------------------------------------------------#
    # Check for prior logs
    #--------------------------------------------------------------------------#

    last <- get_last_commit(vesalius)
    logdf <- list(logdf)
    names(logdf) <- last
    #--------------------------------------------------------------------------#
    # why did i set this? looking back wouldn't we want to always apped log?
    #--------------------------------------------------------------------------#
    if (append) {
        slot(vesalius@log, slot) <- c(slot(vesalius@log, slot), logdf)
    } else {
        slot(vesalius@log, slot) <- logdf
    }
    return(vesalius)
}

get_last_commit <- function(log) {
    log <- slot_apply(log@log, names)
    log <- sapply(lapply(log, as.numeric), max)
    log <- max(log) + 1
    return(log)
}


slot_apply <- function(x, func, ...) {
    cl <- class(x)
    result <- list()
    for (i in slotNames(cl)) {
        result[[i]] <- func(slot(x, i), ...)
    }
    return(result)
}


assign_log_slot <- function(log, commit, logdf) {
    log <- switch(commit,
        "buildVesaliusEmbeddings" = slot_assign(log,
            c("tiles", "embeddings", "activeEmbeddings"), logdf))
    return(log)
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
    return("new_trial" = new_trial)
}


get_counts <- function(vesalius, type = "raw") {
    counts <- vesalius@counts[[type]]
    if (is.null(counts)) {
        stop("No count matrix has been added to Vesalius")
    }
    return(counts)

}

view_log_tree <- function(vesalius) {
    return(vesalius@log)
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
