################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Format conversion Functions/----------------------------#


.ves_to_c <- function(vesalius,
  dims,
  embed = "last",
  correct_background = TRUE,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # For this function we will not create a "method". There are too many
    # other options that we need to account for.
    # Might update this with medici later
    #--------------------------------------------------------------------------#
    embeddings <- check_embedding(vesalius, embed, dims)
    tiles <- check_tiles(vesalius)
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
      #------------------------------------------------------------------------#
      if (correct_background) {
        cimg_tmp <- cimg %>%
          select(c("x", "y", "value")) %>%
          suppressWarnings() %>%
          as.cimg() %>%
          as.data.frame()
        non_img <- paste0(cimg$x, "_", cimg$y)
        in_img <- paste0(cimg_tmp$x, "_", cimg_tmp$y)
        # Median value for background - Other metrics??
        cimg_tmp[!in_img %in% non_img, "value"] <- median(cimg$value)
        image_list[[i]] <- suppressWarnings(
          as.cimg(cimg_tmp[, c("x", "y", "value")], x = x, y = y))
      } else {
        image_list[[i]] <- suppressWarnings(
          as.cimg(cimg_tmp[, c("x", "y", "value")], x = x, y = y))
      }
  }
  return(image_list)
}



.c_to_ves <- function(cimg,
  vesalius,
  dims,
  embed = "last",
  verbose = TRUE) {
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
    tiles <- check_tiles(vesalius)
    embeddings <- check_embedding(vesalius, embed, dims)
    embed <- ifelse(embed == "last",
      yes = names(embeddings),
      no = embed)
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
      locs <- match(rownames(embed), barcodes$barcodes)
      embeddings[locs, dims[i]] <- barcodes$value

    }

    embeddings <- list(embeddings)
    names(embeddings) <- embed
    vesalius@activeEmbeddings <- embeddings

    return(vesalius)
}



.ves_to_sis <- function(vesalius,
  dims,
  embed = "last",
  correct_background = TRUE,
  verbose = TRUE) {
    #--------------------------------------------------------------------------#
    # First getting the right embedding and stuff out
    #--------------------------------------------------------------------------#
    tiles <- check_tiles(vesalius)
    embeddings <- check_embedding(vesalius, embed, dims)
    #--------------------------------------------------------------------------#
    # Now we can create the image list
    #--------------------------------------------------------------------------#
    image_list <- list()
    for (i in seq_along(dims)) {
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
      #------------------------------------------------------------------------#
      if (correct_background) {
          cimg_tmp <- cimg %>%
            select(c("x", "y", "value")) %>%
            suppressWarnings() %>%
            as.cimg() %>%
            as.data.frame()
          non_img <- paste0(cimg$x, "_", cimg$y)
          in_img <- paste0(cimg_tmp$x, "_", cimg_tmp$y)
          # Median value for background - Other metrics??
          cimg_tmp[!in_img %in% non_img, "value"] <- median(cimg$value)
          sis <- as.matrix(suppressWarnings(
          as.cimg(cimg_tmp[, c("x", "y", "value")], x = x, y = y)))
      } else {
          sis <- as.matrix(suppressWarnings(
          as.cimg(cimg_tmp[, c("x", "y", "value")], x = x, y = y)))
      }

      #------------------------------------------------------------------------#
      # hopefully this will do the trick - OpenImage require 3d Arrays
      # so we will create 3 instances of the same gray scale array
      #------------------------------------------------------------------------#
      img <- replicate(sis, 3)
      dim(img) <- c(nrow(sis), ncol(sis), 3)
      image_list[[i]] <- img
    }
    return(image_list)
}

.sis_to_ves <- function(image,
  object,
  dims,
  embed = "last") {
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
  tiles <- check_tiles(vesalius)
  embeddings <- check_embedding(vesalius, embed, dims)
  embed <- ifelse(embed == "last",
    yes = names(embeddings),
    no = embed)
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  #--------------------------------------------------------------------------#
  for (i in seq_along(dims)){
    img <- .sis_to_df(image[[i]])
    barcodes <- left_join(tiles, img, by = c("x", "y")) %>%
                 filter(origin == 1) %>%
                 na.exclude()
    locs <- match(barcodes$barcodes, rownames(embed))
    embeddings[locs, dims[i]] <- barcodes$value

  }
  embeddings <- list(embeddings)
  names(embeddings) <- embed
  object@activeEmbeddings <- embeddings

  return(object)
}

.sis_to_df <- function(image, is_cimg = TRUE) {
  image <- image$AP_image_data
  y <- rep(seq(1, ncol(image)), each = nrow(image))
  x <- rep(seq(1, nrow(image)), times = ncol(image))
  value <- as.vector(image[seq_len(nrow(image)), seq_len(ncol(image)), 1])
  df <- data.frame("x" = x, "y" = y, "value" = value)
  return(df)
}
