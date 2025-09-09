expect_error_or_warning <- function(object,
                                    regexp = NULL,
                                    class = NULL,
                                    expect_s4 = NULL,
                                    info = NULL,
                                    label = NULL) {
  
  # Capture the expression WITHOUT evaluating it yet
  object_quo <- rlang::enquo(object)
  act_label <- label %||% rlang::as_label(object_quo)
  
  warnings_collected <- list()
  error <- NULL
  value <- NULL
  
  # Store original warning settings
  old_warn <- options(warn = 1)  # Make warnings immediate
  on.exit(options(old_warn), add = TRUE)
  
  # Now evaluate with proper warning/error handling
  withCallingHandlers(
    {
      # Main evaluation with error catching
      tryCatch(
        withCallingHandlers(
          {
            # NOW we evaluate the expression
            value <- rlang::eval_tidy(object_quo)
            # Force any pending warnings to fire
            force(value)
          },
          warning = function(w) {
            warnings_collected <<- append(warnings_collected, list(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          error <<- e
          NULL
        }
      )
    },
    warning = function(w) {
      # Catch any stragglers that might escape the inner handler
      warnings_collected <<- append(warnings_collected, list(w))
      invokeRestart("muffleWarning")
    }
  )
  
  # Check if nothing happened
  if (is.null(error) && length(warnings_collected) == 0) {
    testthat::fail(sprintf("%s did not throw an error or warning.", act_label))
    return(invisible(value))
  }
  
  # Build condition description
  conds <- c()
  if (length(warnings_collected) > 0) conds <- c(conds, "warning")
  if (!is.null(error)) conds <- c(conds, "error")
  
  # Select condition object for regexp/class checking
  # Prefer error over warning if both exist
  cond_obj <- if (!is.null(error)) error else warnings_collected[[1]]
  
  # Check regexp pattern
  if (!is.null(regexp)) {
    testthat::expect_match(conditionMessage(cond_obj), regexp, info = info)
  }
  
  # Check condition class
  if (!is.null(class)) {
    if (methods::is(cond_obj, class)) {
      testthat::succeed(sprintf("Condition object is S4 class '%s'", class))
    } else {
      testthat::expect_s3_class(cond_obj, class, info = info)
    }
  }
  
  # Check S4 class of return value (only if no error occurred)
  if (!is.null(expect_s4) && is.null(error)) {
    testthat::expect_true(
      methods::is(value, expect_s4),
      info = sprintf("Expected an S4 object of class '%s', got '%s'",
                     expect_s4, paste(class(value), collapse = ", "))
    )
  }
  
  # Report success with details
  testthat::succeed(
    sprintf("%s threw %s as expected.%s",
            act_label, 
            paste(conds, collapse = " and "),
            if (length(warnings_collected) > 1) 
              sprintf(" (%d warnings total)", length(warnings_collected)) 
            else "")
  )
  
  invisible(value)
}