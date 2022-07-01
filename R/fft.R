################################################################################
###############################   Vesalius      ################################
################################################################################


runTransform <- function(vesalius,
                         transform = c("FFT","DCT"),
                         dims = seq(1,3),
                         chunkSize = 0.1,
                         embedding = "last",
                         type = "ECDF",
                         method = c("median","iso","box"),
                         iter = 1,
                         sigma = 1,
                         box = 20,
                         threshold=0,
                         neuman=TRUE,
                         gaussian=TRUE,
                         na.rm=FALSE,
                         acrossLevels = "min",
                         verbose = TRUE){
    .simpleBar(verbose)

    if(is(vesalius)[1L] == "vesaliusObject"){
      images <- .vesToC(object = vesalius,embed = embedding, dims = dims)
    } else {
      stop("Unsupported format to runTransform function")
    }

    #--------------------------------------------------------------------------#
    # Smoothing the image prior to transform
    #--------------------------------------------------------------------------#
    images <- lapply(images,.internalSmooth,
                                 method = method,
                                 iter = iter,
                                 sigma = sigma,
                                 box = box,
                                 threshold=threshold,
                                 neuman=neuman,
                                 gaussian=gaussian,
                                 na.rm=na.rm,
                                 acrossLevels = acrossLevels)
    #--------------------------------------------------------------------------#
    # balence histo
    # might need to switch this around with smoothing
    #--------------------------------------------------------------------------#
    images <- switch(type,
           "EqualizePiecewise" = lapply(images,EqualizePiecewise,
             N),
           "BalanceSimplest" = lapply(images,BalanceSimplest,
             sleft,sright,range = c(0,1)),
           "SPE" = lapply(image,SPE,
             lambda),
           "EqualizeDP"  = lapply(images,EqualizeDP,
             down,up),
           "EqualizeADP" = lapply(images,EqualizeADP),
           "ECDF" = lapply(images,.ecdf.eq))
    #--------------------------------------------------------------------------#
    # Chunking image into smaller images for each embedding
    # the idea is that we can do this to look at the image at a smaller scale
    # I think this will need to be cleaned up to make sure that we look at Cells
    # and not just sub sections of the images
    #--------------------------------------------------------------------------#
    if(chunkSize <1){
        .chunk(verbose)
        images <- lapply(images,.chunkImage,chunkSize= chunkSize)
        transformed <- switch(transform[1L],
                              "FFT" = lapply(images, function(img){
                                  return(lapply(img,.fft))
                              }),
                              "DCT" = lapply(images, function(img){
                                  return(lapply(img,.dct))
                              }))
    } else if(chunkSize ==1){
        transformed <- switch(transform[1L],
                          "FFT" = lapply(images,.fft),
                          "DCT" = lapply(images,.dct))
    } else {
        stop("chunkSize to large! Must be between 0 and 1")
    }
    #--------------------------------------------------------------------------#
    # more stuff
    #--------------------------------------------------------------------------#
    .simpleBar(verbose)
    return(transformed)

}


.chunkImage <- function(image,chunkSize =1){

    #--------------------------------------------------------------------------#
    # Get dims and convert to cut index
    #--------------------------------------------------------------------------#
    rows <- round(seq(1,dim(image)[1L],l=(1/chunkSize)+1))
    cols <- round(seq(1,dim(image)[2L],l=(1/chunkSize)+1))
    #--------------------------------------------------------------------------#
    # sub image an create list of subimage
    # Images are represented as y = columns and x = rows for some reason
    #--------------------------------------------------------------------------#
    subImages <- list()
    count <- 1
    for(i in seq(1, length(rows)-1)){
        for(j in seq(1,length(cols)-1)){
            #------------------------------------------------------------------#
            # sub images and add the original coordinates
            # at least we can know where each sub image was taken from
            # i could also use string tags...
            #------------------------------------------------------------------#
            tmp <- imsub(image,
                         x >= rows[j] & x <= rows[j+1],
                         y >= cols[i] & y <= cols[i+1])
            subImages[[count]] <- list("img" = tmp,
                                       "coordinates" = data.frame("x" = c(rows[j],rows[j+1]),
                                                                  "y" = c(cols[i],cols[i+1])))

            count <- count +1
        }
    }
    return(subImages)

}

#------------------------------------------------------------------------------#
# Using a line sampling approach 
#------------------------------------------------------------------------------#

runLineTransform <- function(vesalius,){
    
}
#------------------------------------------------------------------------------#
# General functions 
#------------------------------------------------------------------------------#
.fft <- function(image){
    #--------------------------------------------------------------------------#
    # For now I will create a way of plotting both the FFT imag and the original
    # sub plot
    #--------------------------------------------------------------------------#
    fft <- FFT(image[[1]])
    return(list("img" = image[[1]],
                "real" = fft$real,
                "fft" = fft$imag,
                "coordinates" = image[[2]]))

}

.dct <- function(image){
  #--------------------------------------------------------------------------#
  # For now I will create a way of plotting both the DCT imag and the original
  # sub plot
  #--------------------------------------------------------------------------#
  #fft <- FFT(image[[1]])
  #return(list("img" = image[[1]], "fft" = fft$imag, "coordinates" = image[[2]]))
}
