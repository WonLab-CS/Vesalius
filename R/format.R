################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Format conversion Functions/----------------------------#


.vesToC <- function(object,dims,embed = "last",correctBackground =TRUE,verbose =TRUE){
    #--------------------------------------------------------------------------#
    # For this function we will not create a "method". There are too many
    # other options that we need to account for.
    # TODO if embed is not last then you have to switch the last active embed
    # to the one that is being requested.
    #--------------------------------------------------------------------------#
    tiles <- object@tiles

    if(embed == "last"){
        embeddings <- object@activeEmbeddings[[1L]]
    }else{
       embeddings <- object@embeddings[embed]
       if(length(embeddings)==0){
          stop(paste(deparse(substitute(embed)),":Unknown embedding selected!"))
       } else if(length(embeddings)>1){
          warning(paste("More than 1",deparse(substitute(embed)),"embedding.
          Vesalius will use the latest entry - See Vesalius Log"))
          embeddings <- embeddings[[length(embeddings)]]
       } else {
         embeddings <- embeddings[[1L]]
       }
    }
  
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    #--------------------------------------------------------------------------#
    imageList <- list()
    .vtc(verbose)
    for(i in seq_along(dims)){
      embeds <- embeddings[,dims[i]]
      embeds <- data.frame(names(embeds),embeds)
      colnames(embeds) <- c("barcodes",as.character(dims[i]))
      cimg <- right_join(tiles,embeds, by= "barcodes")
      colnames(cimg) <- c("barcodes","x","y","origin","value")
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



.cToVes <- function(cimg,object,dims,embed = "last"){
    #--------------------------------------------------------------------------#
    # Get stuff out
    #--------------------------------------------------------------------------#
    tiles <- object@tiles
    if(embed == "last"){
        embeddings <- object@activeEmbeddings[[1L]]
        embed <- names(object@activeEmbeddings[[1L]])
    }else{
       embeddings <- object@embeddings[embed]
       if(length(embeddings)==0){
          stop(paste(deparse(substitute(embed)),":Unknown embedding selected!"))
       } else if(length(embeddings)>1){
          warning(paste("More than 1",deparse(substitute(embed)),"embedding.
          Vesalius will use the latest entry - See Vesalius Log"))
          embeddings <- embeddings[[length(embeddings)]]
       } else {
         embeddings <- embeddings[[1L]]
       }
    }


    #--------------------------------------------------------------------------#
    # Always going to be a gray scale image.
    # Colour are only used when during viz
    # This allows any arbitrary number of embeds
    #--------------------------------------------------------------------------#
    for(i in seq_along(dims)){
      img <- as.data.frame(cimg[[i]])
      barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                  filter(origin == 1) %>% na.exclude()

      locs <- match(rownames(embeds),barcodes$barcodes)
      embeddings[locs,dims[i]] <- barcodes$value

    }

    embeddings <- list(embeddings)
    names(embeddings) <- embed
    object@activeEmbeddings <- embeddings

    return(object)
}



.vesToSIS <- function(object,dims,embed = "last",correctBackground =TRUE,verbose =TRUE){
    #--------------------------------------------------------------------------#
    # First getting the right embedding and stuff out
    #--------------------------------------------------------------------------#
    tiles <- object@tiles
    if(embed == "last"){
        embeddings <- object@activeEmbeddings[[1L]]
    }else{
       embeddings <- object@activeEmbeddings[embed]
       if(length(embeddings)==0){
          stop(paste(deparse(substitute(embed)),":Unknown embedding selected!"))
       } else if(length(embeddings)>1){
          warning(paste("More than 1",deparse(substitute(embed)),"embedding.
          Vesalius will use the latest entry - See Vesalius Log"))
          embeddings <- embeddings[[length(embeddings)]]
       } else {
         embeddings <- embeddings[[1L]]
       }
    }
    #--------------------------------------------------------------------------#
    # Now we can create the image list
    #--------------------------------------------------------------------------#
    imageList <- list()
    for(i in seq_along(dims)){
      embeds <- object@embeddings[,dims[i]]
      embeds <- data.frame(names(embeds),embeds)
      colnames(embeds) <- c("barcodes",as.character(dims[i]))
      sis <- right_join(tiles,embeds, by= "barcodes") %>% na.exclude()
      colnames(cimg) <- c("barcodes","x","y","origin","value")
      #------------------------------------------------------------------------#
      # Coreecting backgound with median value
      #------------------------------------------------------------------------#
      # Using the cimg format just because they have optimised format conversion
      #------------------------------------------------------------------------#
      if(correctBackground){
          cimgTmp <- cimg %>%
             select(c("x","y","value")) %>%
             suppressWarnings %>%
             as.cimg %>%
             as.data.frame
          nonImg <- paste0(cimg$x,"_",cimg$y)
          inImg <- paste0(cimgTmp$x,"_",cimgTmp$y)
          cimgTmp[!inImg %in% nonImg,"value"] <- median(cimg$value)

          sis <- suppressWarnings(as.cimg(sis[,c("x","y","value")]))
          sis <- as.matrix(sis)
      } else {
        sis <- suppressWarnings(as.cimg(sis[,c("x","y","value")]))
        sis <- as.matrix(sis)
      }

      #------------------------------------------------------------------------#
      # hopefully this will do the trick - OpenImage require 3d Arrays
      # so we will create 3 instances of the same gray scale array
      #------------------------------------------------------------------------#
      img <- replicate(sis,3)
      dim(img) <- c(nrow(sis),ncol(sis),3)
      imageList[[i]] <- img
    }
    return(imageList)
}

.SISToVes <- function(image,object,dims,embed = "last"){
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  if(embed == "last"){
      embeddings <- object@activeEmbeddings[[1L]]
      embed <- names(object@activeEmbeddings[[1L]])
  }else{
     embeddings <- object@embeddings[embed]
     if(length(embeddings)==0){
        stop(paste(deparse(substitute(embed)),":Unknown embedding selected!"))
     } else if(length(embeddings)>1){
        warning(paste("More than 1",deparse(substitute(embed)),"embedding.
        Vesalius will use the latest entry - See Vesalius Log"))
        embeddings <- embeddings[[length(embeddings)]]
     } else {
       embeddings <- embeddings[[1L]]
     }
  }
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  #--------------------------------------------------------------------------#
  for(i in seq_along(dims)){
    img <- .SISToDF(image[[i]])
    barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                 filter(origin == 1) %>% na.exclude()
    locs <- match(barcodes$barcodes,rownames(embeds))
    embeddings[locs,dims[i]] <- barcodes$value

  }
  embeddings <- list(embeddings)
  names(embeddings) <- embed
  object@activeEmbeddings <- embeddings

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
