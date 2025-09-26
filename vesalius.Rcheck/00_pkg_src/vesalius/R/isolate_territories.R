################################################################################
############################   ST package        ###############################
################################################################################



#' isolating territories from vesalius image segments
#' @param vesalius_assay vesalius_Assay object
#' @param method character describing barcode pooling method.
#' Currently, only "distance" availble
#' @param trial character string describing which segmentation trial
#' to use. Default is "last" which is the last segmentation trial used.
#' @param capture_radius numeric - proportion of maximum distance between
#' barcodes that will be used to pool barcodes together (range 0 - 1).
#' @param global logical - If TRUE, territories will be numbered across all
#' colour segments. If FALSE, territories will be numbered within each colour
#' segment.
#' @param min_spatial_index integer - minimum number of barcodes/spots/beads
#' required in each territory
#' @param verbose logical - progress message output.
#' @details Image segments can be further subdivided into 2D
#' seperated territorires. This is accomplished by pooling barcodes
#' that are associated with a colour cluster into territories based on the
#' distance between each barcode.
#'
#' First, \code{isolate_territories} considers the maximum distance
#' between all beads. The \code{capture_radius} will define which 
#' proportion of this distance should be considered.
#'
#' Second, a seed barcode will be selected and all barcodes that are within the
#' capture distance of the seed barcode with be pooled together. This process
#' is then applied on barcodes that have been selected in this manner. The
#' process is repeated until all barcodes have been pooled into a territory.
#' If there are still barcodes remaining, a new seed barcode is selected and the
#' whole process is repeated. NOTE : Territory isolation is applied to each
#' colour segment independently.
#'
#' If a territory does not contain enough barcodes, it will be pooled into the
#' isolated territory. This territory contains all isolated territories
#' regardless of colour cluster of origin.
#'
#' @return a vesalius_assay object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple smoothing
#' ves <- smooth_image(ves, dimensions = seq(1, 30))
#' 
#' # quick segmentation
#' ves <- segment_image(ves, dimensions = seq(1, 30))
#' 
#' # isolate territories
#' ves <- isolate_territories(ves)
#'}
#' @export
isolate_territories <- function(vesalius_assay,
  method = "distance",
  trial = "last",
  capture_radius = 0.05,
  global = TRUE,
  min_spatial_index = 10,
  verbose = TRUE) {
    simple_bar(verbose)
    #--------------------------------------------------------------------------#
    # Get stuff out as usual
    # Not super happy with this check
    # it's a bit messy - we might need to consider to do a whole sanity check
    # of inout data and see if that makes sense - this will include checking log
    #--------------------------------------------------------------------------#
    ter <- check_segment_trial(vesalius_assay, trial) %>%
      na.exclude()
    #--------------------------------------------------------------------------#
    # Compute real capture Radius
    # Only one method for now so this is not necessary 
    # keep it for later
    #--------------------------------------------------------------------------#
    method <- check_isolation_method(method)
    if (method[1L] == "distance") {
      capture_radius <- sqrt(((max(ter$x) - min(ter$x))^2 +
        (max(ter$y) - min(ter$y))^2)) * capture_radius
    }
    #--------------------------------------------------------------------------#
    # Creating new trial column name and adding it to input data
    # The input data here is a subset of the full territory df
    # we at least make sure that we are using the right input
    #--------------------------------------------------------------------------#
    new_trial <- create_trial_tag(colnames(get_territories(vesalius_assay)),
      "Territory") %>%
      tail(1)
    ter$trial <- 0
    #--------------------------------------------------------------------------#
    # Now we can dispatch the necessary information
    # and run the pooling algorithm
    #--------------------------------------------------------------------------#
    segment <- unique(ter$Segment)
    for (i in seq_along(segment)) {
      message_switch("ter_pool", verbose, ter = i)
      #------------------------------------------------------------------------#
      # We don't need all data in this case only clusters and locations
      # We can just rebuild everything afterwards
      # At least we don't don't need to compute anything unnecessarily
      ## Note other method not in use for now
      # Might not be worthwile to implement them
      # Argument could be removed
      #------------------------------------------------------------------------#

      tmp <- ter[ter$Segment == segment[i], ]
      tmp <- switch(method[1L],
        "distance" = distance_pooling(tmp, capture_radius,
          min_spatial_index))
      #------------------------------------------------------------------------#
      # Skipping if does not meet min cell requirements
      # filtering can be done later
      #------------------------------------------------------------------------#
      if (is.null(tmp)) next()
      #------------------------------------------------------------------------#
      # adding territories
      #------------------------------------------------------------------------#
      ter$trial[ter$Segment == segment[i]] <- tmp

    }

    #--------------------------------------------------------------------------#
    # Globalise territories crate a numbering system for all colour clusters
    # If set to false territories only describe territories for any given colour
    # cluster
    #--------------------------------------------------------------------------#
    if (global) {
        ter <- globalise_territories(ter)
    }
    #--------------------------------------------------------------------------#
    #Now we can add it to the full territory df and update vesalius
    #--------------------------------------------------------------------------#

  colnames(ter) <- gsub("trial", new_trial, colnames(ter))
  vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = ter,
    slot = "territories",
    append = TRUE)
  commit <- create_commit_log(arg_match = as.list(match.call()),
    default = formals(isolate_territories))
  vesalius_assay <- commit_log(vesalius_assay,
    commit,
    get_assay_names(vesalius_assay))
  simple_bar(verbose)
  return(vesalius_assay)
}



#' distance pooling beads of colour segment into seperate territories
#' @param img data frame contain all barcodes of a single sgement
#' @param capture_radius numeric proportion of max distance between beads
#' to use as distance threshold between beads
#' @param min_spatial_index numeric minimum number of beads that should
#' be contained in a territory. 
#' @details Beads that are too far away or bead cluster that are below
#' the minimum number of spatial indices will all be pooled under the
#' isolated label. Note that this label is used across color segments.
#' @importFrom dplyr %>% distinct
distance_pooling <- function(img, capture_radius, min_spatial_index) {
    #--------------------------------------------------------------------------#
    # Select center point of each tile for only one channel
    # Dont need to run it for all channels
    #--------------------------------------------------------------------------#
    img_copy <- img  %>%
      distinct(barcodes, .keep_all = TRUE) 
    if (nrow(img_copy) < 1) { return(NULL) }
    #--------------------------------------------------------------------------#
    # Compute distances
    #--------------------------------------------------------------------------#
    idx <- seq_len(nrow(img_copy))
    distance_matrix <- lapply(idx, function(idx, mat) {
      xo <- mat$x[idx]
      yo <- mat$y[idx]
      xp <- mat$x
      yp <- mat$y
      distance <- sqrt(((abs(xp - xo))^2 + (abs(yp - yo))^2))
      return(distance)
    }, img_copy)


    #--------------------------------------------------------------------------#
    # Buildan actual matrix
    #--------------------------------------------------------------------------#

    distance_matrix <- do.call("rbind", distance_matrix)
    colnames(distance_matrix) <- img_copy$barcodes
    rownames(distance_matrix) <- img_copy$barcodes

    #--------------------------------------------------------------------------#
    # If there is only one point
    # In this case you only have one territory as well
    # Don't need to any pooling
    #--------------------------------------------------------------------------#
    if (sum(dim(distance_matrix)) == 2) {
        return(1)
    }


    #--------------------------------------------------------------------------#
    # Pooling points together
    #--------------------------------------------------------------------------#
    barcodes <- img_copy$barcodes
    territories <- list()
    count <- 1

    while (length(barcodes) > 0) {
         #---------------------------------------------------------------------#
         # First lets select a random barcode in the colour segment
         # And create a pool of barcodes to select based on capture radius
         # This first while loop checks if there are still barcodes
         # left in the colour segment
         #---------------------------------------------------------------------#
          tmp <- distance_matrix[, sample(barcodes, 1)]
          pool <- names(tmp)[tmp <= capture_radius]
          inter <- pool
          converge <- FALSE

          while (!converge) {
            #------------------------------------------------------------------#
            # This while loop checks if all possible barcodes have been pooled
            # into the current territory
            #------------------------------------------------------------------#
            if (length(inter) == 1) {
              #------------------------------------------------------------#
              # when there is only one barcodes
              # remove barcode from pool and move on
              #------------------------------------------------------------#
              territories[[count]] <- pool
              barcodes <- barcodes[!barcodes %in% pool]
              count <- count + 1
              converge <- TRUE
            } else {
              #------------------------------------------------------------#
              # Get a new pool from the distance matrix
              # and check which ones are within capture radius
              #------------------------------------------------------------#
              new_pool <- distance_matrix[, inter]
              new_pool <- unique(unlist(lapply(seq_len(ncol(new_pool)),
                function(idx, np, capture_radius) {
                  res <- rownames(np)[np[, idx] <= capture_radius]
                  return(res)
              }, new_pool, capture_radius)))
              #------------------------------------------------------------#
              # check which barcodes in the new pool overlap with
              # the ones in the full pool
              # If there is a perfect overlap then there are no new bacodes
              # to pool into a territory
              #------------------------------------------------------------#
              overlap <- new_pool %in% pool
              if (sum(overlap) != length(new_pool)) {
                #------------------------------------------------------#
                # There are still some new barcodes to pool
                # lets do some more looping then
                #------------------------------------------------------#
                pool <- unique(c(pool, new_pool[!overlap]))
                inter <- unique(new_pool[!overlap])
                converge <- FALSE
              } else {
                #------------------------------------------------------#
                # it is done ! no more barcodes for this territory
                #------------------------------------------------------#
                territories[[count]] <- pool
                count <- count + 1
                barcodes <- barcodes[!barcodes %in% pool]
                converge <- TRUE
              }
            }
          }
      }
    #--------------------------------------------------------------------------#
    # Clean up drop outs
    #--------------------------------------------------------------------------#

    all_ters <- img$trial

    for (ter in seq_along(territories)) {
        loc <- img$barcodes %in% territories[[ter]]
        if (length(territories[[ter]]) <= min_spatial_index) {
            all_ters[loc] <- "isolated"
        } else {
            all_ters[loc] <- ter
        }
    }
      return(all_ters)
}


#' @importFrom imager imsplit  threshold split_connected where
#' @importFrom imagerExtra ThresholdML
#' @importFrom dplyr inner_join
select_similar <- function(img,
  coordinates,
  threshold = 1) {
  img <- img %>%
    imsplit("cc")

  pos <- lapply(img, function(x, threshold){
      ret <- ThresholdML(x, threshold)
      ret <- split_connected(ret)
      return(ret)
  }, threshold = threshold)
  coordinates$Segment <- 1
  count <- 2
  for (i in seq_along(pos)){
    for (j in seq_along(pos[[i]])) {
      tmp <- where(pos[[i]][[j]]) %>%
        inner_join(coordinates, by = c("x", "y"))
         coordinates$Segment[coordinates$barcodes %in% tmp$barcodes &
          coordinates$Segment == 1] <- count
        count <- count + 1
    }
  }
  all_ter <- unique(coordinates$Segment)
  coordinates$Segment <- seq_along(all_ter)[match(coordinates$Segment, all_ter)]
  return(coordinates)
}

#' @importFrom future.apply future_lapply
connected_pixels <- function(clusters,
  embeddings,
  k = 6,
  threshold = 0.90,
  verbose = TRUE) {
    message_switch("connect_pixel", verbose)
    #-------------------------------------------------------------------------#
    # First we need to get super pixel centers 
    #-------------------------------------------------------------------------#
    center_pixels <- sort(unique(clusters$Segment))
    centers <- future_lapply(center_pixels, function(center, segments){
        x <- median(segments$x[segments$Segment == center])
        y <- median(segments$y[segments$Segment == center])
        df <- data.frame("x" = x, "y" = y)
        rownames(df) <- center
        return(df)
    }, segments = clusters) %>% do.call("rbind", .)
    #-------------------------------------------------------------------------#
    # Next we get the nearest neighbors
    #-------------------------------------------------------------------------#
    knn <- RANN::nn2(centers, k = k + 1)$nn.idx
    rownames(knn) <- rownames(centers)
    #-------------------------------------------------------------------------#
    # Next we intialise a graph and then compute correlation
    #-------------------------------------------------------------------------#
    graph <- populate_graph(knn)
    graph$cor <- 0
    for (i in seq_len(nrow(graph))) {
        c1 <- embeddings[clusters$barcodes[clusters$Segment == graph$e1[i]], ]
        if (!is.null(nrow(c1))) {
            c1 <- apply(c1, 2, mean)
        }
        c2 <- embeddings[clusters$barcodes[clusters$Segment == graph$e2[i]], ]
        if (!is.null(nrow(c2))) {
            c2 <- apply(c2, 2, mean)
        }
        graph$cor[i] <- cor(c1, c2, method = "pearson")
    }
    #-------------------------------------------------------------------------#
    # Then we interatively pool pixels together under the a transitive 
    # correlation assumption i.e if A cor B and B cor with C then A cor C
    #-------------------------------------------------------------------------#
    initial_pixels <- unique(graph$e1)
    total_pool <- c()
    segments <-  list()
    count <- 1
    while (length(initial_pixels > 0)) {
        start_pixel <- sample(initial_pixels, size = 1)
        pool <- graph$e2[graph$e1 == start_pixel & graph$cor >= threshold]
        inter <- pool
        total_pool <- c(total_pool, pool)
        converge <- FALSE
        #browser()
        while (!converge) {
            if (length(inter) == 1) {
                segments[[count]] <- pool
                initial_pixels <- initial_pixels[!initial_pixels %in% pool]
                count <- count + 1
                converge <- TRUE
            } else {
                new_pool <- unique(graph$e2[graph$e1 %in% inter &
                  graph$cor >= threshold])
                overlap <- new_pool %in% pool & !new_pool %in% total_pool
                if (sum(overlap) != length(new_pool)) {
                  pool <- unique(c(pool, new_pool[!overlap]))
                  
                  inter <- unique(new_pool[!overlap])
                  converge <- FALSE
                } else {
                  segments[[count]] <- pool
                  total_pool <- c(total_pool, pool)
                  count <- count + 1
                  initial_pixels <- initial_pixels[!initial_pixels %in% pool]
                  converge <- TRUE
                }
            }
        }
    }
    
    for (seg in seq_along(segments)){
        loc <- clusters$Segment %in% segments[[seg]]
        clusters$Segment[loc] <- seg
    }
    return(clusters)
}