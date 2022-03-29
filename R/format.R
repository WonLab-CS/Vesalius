

.vesToC <- function(object, dims,correctBackground =TRUE,verbose =TRUE){
    #--------------------------------------------------------------------------#
    # Get relevant slots and prepare for joining
    #--------------------------------------------------------------------------#
    tiles <- object@tiles
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    #--------------------------------------------------------------------------#
    imageList <- list()
    .vtc(verbose)
    for(i in seq_along(dims)){
      embeds <- object@embeddings[,dims[i]]
      embeds <- data.frame(names(embeds),embeds)
      colnames(embeds) <- c("barcodes",as.character(dims[i]))
      cimg <- right_join(tiles,embeds, by= "barcodes")
      colnames(cimg) <- c("barcodes","x","y","value")
      cimg <- na.exclude(cimg)
      if(correctBackground){
          cimgTmp <- cimg %>%
             select(c("x","y","value")) %>%
             suppressWarnings %>%
             as.cimg %>%
             as.data.frame
          nonImg <- paste0(cimg$x,"_",cimg$y)
          inImg <- paste0(cimgTmp$x,"_",cimgTmp$y)
          cimgTmp[!inImg %in% nonImg,"value"] <- median(cimg$value)
          imageList[[i]] <- suppressWarnings(as.cimg(cimgTmp[,c("x","y","value")]))
      } else {
         imageList[[i]] <- suppressWarnings(as.cimg(cimg[,c("x","y","value")]))
      }


    }
    return(imageList)
}

.vesToDF <- function(object, dims,verbose =TRUE){
  #--------------------------------------------------------------------------#
  # Get relevant slots and prepare for joining
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  #--------------------------------------------------------------------------#
  # generate a list of images based on the number of dims
  # Remember that any time you do anything to an image
  # it is always applied to each "channel" separately - it will be so
  # much easier to just consider everything as gray scale for this
  #--------------------------------------------------------------------------#
  imageDF <- list()
  .vtdf(verbose)
  for(i in seq_along(dims)){
    embeds <- object@embeddings[,dims[i]]
    embeds <- data.frame(names(embeds),embeds)
    colnames(embeds) <- c("barcodes",as.character(dims[i]))
    cimg <- right_join(tiles,embeds, by= "barcodes")
    cimg$cc <- i
    colnames(cimg) <- c("barcodes","x","y","value","cc")
    cimg <- na.exclude(cimg)
    cimg <- cimg[,c("barcodes","x","y","cc","value")]
    imageDF[[i]] <- cimg
   }
  imageDF <- do.call("rbind", imageDF)
  return(imageDF)
}

.cToVes <- function(cimg,object,dims){
    #--------------------------------------------------------------------------#
    # Get stuff out
    #--------------------------------------------------------------------------#
    tiles <- object@tiles
    embeds <- object@embeddings
    #--------------------------------------------------------------------------#
    # Always going to be a gray scale image.
    # Colour are only used when during viz
    # This allows any arbitrary number of embeds
    #--------------------------------------------------------------------------#
    for(i in seq_along(dims)){
      img <- as.data.frame(cimg[[i]])
      barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                  distinct(barcodes,value) %>% na.exclude()

      locs <- match(rownames(embeds),barcodes$barcodes)
      embeds[locs,dims[i]] <- barcodes$value

    }
    object@embeddings <- embeds

    return(object)
}

.dfToVes <- function(df, object,dims){
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  embeds <- object@embeddings
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  # this is probably a shit way to do it
  # re think for later
  #--------------------------------------------------------------------------#
  for(i in seq_along(dims)){
    img <- filter(df, cc == i)
    barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                distinct(barcodes,value) %>% na.exclude()

    locs <- match(rownames(embeds),barcodes$barcodes)
    embeds[locs,dims[i]] <- barcodes$value

  }
  object@embeddings <- embeds

  return(object)
}

.vesToSIS <- function(object,dims){
    tiles <- object@tiles
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    # also this is a different image format so we will adpt the function
    # we had above.
    #--------------------------------------------------------------------------#
    imageList <- list()
    for(i in seq_along(dims)){
      embeds <- object@embeddings[,dims[i]]
      embeds <- data.frame(names(embeds),embeds)
      colnames(embeds) <- c("barcodes",as.character(dims[i]))
      sis <- right_join(tiles,embeds, by= "barcodes") %>% na.exclude()
      colnames(sis) <- c("barcodes","x","y","value")
      #------------------------------------------------------------------------#
      # Using the cimg format just because they have optimised format conversion
      #------------------------------------------------------------------------#
      sis <- suppressWarnings(as.cimg(sis[,c("x","y","value")]))
      sis <- as.matrix(sis)
      #------------------------------------------------------------------------#
      # hopefully this will do the trick - OpenImage require 3d Arrays
      # so we will create 3 instances of the same gray scale array
      #------------------------------------------------------------------------#
      img <- c(sis,sis,sis)
      dim(img) <- c(nrow(sis),ncol(sis),3)
      imageList[[i]] <- img
    }
    return(imageList)
}

.SISToVes <- function(image,object,dims){
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  embeds <- object@embeddings
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  #--------------------------------------------------------------------------#
  for(i in seq_along(dims)){
    img <- .SISToDF(image[[i]])
    barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                distinct(barcodes,value)
    locs <- match(barcodes$barcodes,rownames(embeds))
    embeds[locs,dims[i]] <- barcodes$value

  }
  object@embeddings <- embeds

  return(object)
}

.SISToDF <- function(image, is.cimg = TRUE){
    image <- image$AP_image_data
    y <- rep(seq(1,ncol(image)), each = nrow(image))
    x <- rep(seq(1,nrow(image)), times = ncol(image))
    value <- as.vector(image[seq_len(nrow(image)),seq_len(ncol(image)),1])

    df <- data.frame("x" = x,"y" = y, "value" = value)
    return(df)
}
