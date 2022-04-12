################################################################################
############################   ST package        ###############################
################################################################################

#--------------/Isolating Territories & composite iamges/----------------------#

buildMosaic <- function(vesalius,
                        dims = seq(1,30),
                        maskType = c("double","pos","neg"),
                        maskThreshold = 0.85,
                        k=10,
                        lambda = 5,
                        method = "simpleMask",
                        smoothingMethod = c("median","iso","box"),
                        iter =1,
                        sigma =1,
                        box =3,
                        acrossLevels =c("mean","min","max"),
                        cores = 1,
                        verbose = TRUE,
                        lambda = 5,
                        niter = 100,
                        regMethod = "TVL2.FiniteDifference",
                        regularization = T,
                        histEqual = T,
                        elementSize = 5,
                        distThresh = 0,
                        sizeThresh = 0){
    #---------------------------------------------------,-----------------------#
    # First let's get the embeddings
    #--------------------------------------------------------------------------#
    .simpleBar(verbose)
    imgList <- .vesToC(vesalius,dims,verbose)
    tiles <- vesalius@tiles %>% distinct(barcodes, .keep_all =TRUE)

    #--------------------------------------------------------------------------#
    # Next we initialise the mosaic building with various stuff
    #--------------------------------------------------------------------------#
    mosaic <- data.frame("barcodes" = tiles$barcodes,
                         "territory" = rep(0,nrow(tiles)))

    #--------------------------------------------------------------------------#
    # eq histo
    # Temp approach  for eq values
    # At the moment only one eq type. I want to be able to automate this
    # selection process
    #--------------------------------------------------------------------------#
    if (histEqual) {
      .eq(verbose)
      eqThresh <- (1-maskThreshold)
      imgList <- equalizeHistogram(imgList,
                                   sleft = eqThresh,sright=eqThresh,
                                   verbose=FALSE,
                                   cores = cores)
    }
    #--------------------------------------------------------------------------#
    # We smooth all arrays -
    # So this goes against the idea of not keeping arrays in memory -
    # At least for this function we don't need to but at the moment lets just
    # check the proof of concept
    #--------------------------------------------------------------------------#
    if (regularization) {
      .reg(verbose)
      imgList <- parallel::mclapply(imgList,.regularise,
                                    lambda = lambda,
                                    mc.cores = cores,
                                    niter = niter,
                                    method = regMethod)
    }
    .smooth(verbose)
    imgList <- parallel::mclapply(imgList,.internalSmooth,
                                  method = smoothingMethod,
                                  iter = iter,
                                  sigma = sigma,
                                  box = box,
                                  threshold=0,
                                  neuman=TRUE,
                                  gaussian=TRUE,
                                  na.rm=FALSE,
                                  acrossLevels = acrossLevels[1L],
                                  mc.cores = cores)


    #--------------------------------------------------------------------------#
    # Run through the embeddings
    # Initialise territory counts
    #--------------------------------------------------------------------------#
    count <- 1
    noInfoEmbed <- c()
    for(i in seq_along(imgList)){
        #----------------------------------------------------------------------#
        # First lets make some maks for each image
        # Count represent each territory - at every round
        #----------------------------------------------------------------------#
        img <- imgList[[i]]
        .msk(i,verbose)
        mosaic <- .mosaic(img,tiles,mosaic,maskThreshold,maskType,method,k, elementSize, distThresh, sizeThresh)
    }
    vesalius@territories <- mosaic
    cat("\n")
    .simpleBar(verbose)
    return(vesalius)
}



.mosaic <- function(img,
                    tiles,
                    mosaic,
                    maskThreshold = 0.9,
                    maskType="double",
                    method ="simpleMask",
                    k=10,
                    elementSize = 5,
                    distThresh = 0,
                    sizeThresh = 0){
    #--------------------------------------------------------------------------#
    # Set up the threshold and what not
    #--------------------------------------------------------------------------#
    if(maskType[1L] == "double" & method == "simpleMask"){
        pix <- .simpleMask(img,tiles,maskThreshold)
        mosaic <- .stackCheck(pix,mosaic)
        pix <- .simpleMask(img,tiles,1-maskThreshold)
        mosaic <- .stackCheck(pix,mosaic)
    } else if(maskType[1L] == "pos" & method == "simpleMask"){
        pix <- .simpleMask(img,tiles,maskThreshold)
        mosaic <- .stackCheck(pix,mosaic)

    }else if(maskType[1L] == "neg" & method == "simpleMask"){
        pix <- .simpleMask(img,tiles,1-maskThreshold)
        mosaic <- .stackCheck(pix,mosaic)
    } else {
        pix <- switch(method,
                    "SIS" = .SISMask(img, tiles,maskThreshold,maskType),
                    "SISKmeans" = .SISKmeansMask(img,tiles,maskThreshold,maskType),
                    "Kmeans" = .kmeansMask(img,tiles,maskThreshold,maskType,k),
                    "KmeansIter" = .kmeansIterMask(img,tiles,maskThreshold),
                    "thresholdML" = .thresholdMlMask(img, tiles, elementSize, maskType=maskType, distThresh=distThresh, sizeThresh = sizeThresh))
        mosaic <- .stackCheck(pix,mosaic)
    }
  return(mosaic)

}

.simpleMask <- function(img,tiles,maskThreshold = 0.9){
  #----------------------------------------------------------------------------#
  # Making mask as checking which barcodes are part of that mask
  # The quantile part here is a bit of a problem. Mainly because it means that
  # you will get quantiles in from your image. Maybe there is no interesting
  # information
  # TO DO check without and assure that you can handle zero sized pixesets
  #----------------------------------------------------------------------------#
  if(maskThreshold > 0.5){

     pix <- img > maskThreshold
     if(sum(pix) ==0) return(NULL)
     pix <- split_connected(pix)
  } else {

     pix <- img < maskThreshold
     if(sum(pix) ==0) return(NULL)
     pix <- split_connected(pix)
  }

  #----------------------------------------------------------------------------#
  # Getting barcodes in that mask
  #----------------------------------------------------------------------------#
  pix <- lapply(pix,function(img, tiles){
      tmp <- where(img)
      ## Need to check this behaviour
      tmp <- inner_join(tmp,tiles, by = c("x","y"))
      tmp <- tmp$barcodes
      return(tmp)
  },tiles)
  #----------------------------------------------------------------------------#
  # Remove empty
  #----------------------------------------------------------------------------#
  nz <- sapply(pix,length) >0
  pix <- pix[nz]

  return(pix)
}

.SISMask <- function(img, tiles,maskThreshold,maskType="double"){
  #----------------------------------------------------------------------------#
  # First we need to convert cimg array to simple array
  #----------------------------------------------------------------------------#
  img <- as.array(img)
  dim(img) <- c(nrow(img),ncol(img))
  img <- replicate(3,img)
  #----------------------------------------------------------------------------#
  # Next we can run the super pixel segmentation
  #----------------------------------------------------------------------------#
  init <- Image_Segmentation$new()
  spx <- suppressWarnings(init$spixel_segmentation(input_image = img,
                                  superpixel = 600,
                                  AP_data = TRUE,
                                  use_median = TRUE,
                                  sim_wL = 1,
                                  sim_wA = 1,
                                  sim_wB = 1,
                                  sim_color_radius = 10,
                                  verbose = FALSE))
  #---------------------------------------------------------------------------#
  # Getting territories from clusters
  #---------------------------------------------------------------------------#
  pix <- .SISToDF(spx)
  if(maskType[1L]=="double"){
      h <- ifelse(max(pix$value) < maskThreshold,max(pix$value),maskThreshold)
      l <- ifelse(min(pix$value) > (1-maskThreshold),min(pix$value),1-maskThreshold)
      #pix <- filter(pix,value == max(value) | value ==min(value))
      pix <- filter(pix,value >=h | value <= l)
      pix <- split(pix, pix$value)
      pix <- lapply(pix,function(img, tiles){
          img <- inner_join(img,tiles, by = c("x","y"))
          img <- img$barcodes
          return(img)
      },tiles)
      nz <- sapply(pix,length) >0
      pix <- pix[nz]
  } else if(maskType[1L] =="pos"){
      h <- ifelse(max(pix$value) < maskThreshold,max(pix$value),maskThreshold)
      #pix <- filter(pix,value == max(value))
      pix <- filter(pix,value >= h)
      pix <- inner_join(pix,tiles, by = c("x","y"))
      pix <- list(pix$barcodes)
  } else {
      l <- ifelse(min(pix$value) > (1-maskThreshold),min(pix$value),1-maskThreshold)
      #pix <- filter(pix,value == min(value))
      pix <- filter(pix,value <= l)
      pix <- inner_join(pix,tiles, by = c("x","y"))
      pix <- list(pix$barcodes)
  }

  #----------------------------------------------------------------------------#
  # Remove empty
  #----------------------------------------------------------------------------#


  return(pix)

}

.SISKmeansMask <- function(img, tiles,maskThreshold,maskType){
  #----------------------------------------------------------------------------#
  # First we need to convert cimg array to simple array
  # No mask threshold at the moment - i have no idea how to select these
  # segments. Will need to ponder this
  #----------------------------------------------------------------------------#
  img <- as.array(img)
  dim(img) <- c(nrow(img),ncol(img))
  img <- replicate(3,img)
  #----------------------------------------------------------------------------#
  # Super pixel segmentation with kmeans
  #----------------------------------------------------------------------------#
  init <- Image_Segmentation$new()
  spx_km <- suppressWarnings(init$spixel_segmentation(input_image = img,
                                superpixel = 600,
                                AP_data = TRUE,
                                use_median = TRUE,
                                sim_wL = 1,
                                sim_wA = 1,
                                sim_wB = 1,
                                sim_color_radius = 10,
                                kmeans_method = "kmeans",
                                kmeans_initializer = "kmeans++",
                                kmeans_num_init = 3,
                                kmeans_max_iters = 100,
                                verbose = FALSE))

   #---------------------------------------------------------------------------#
   # Getting territories from clusters
   #---------------------------------------------------------------------------#
   pix <- .SISToDF(spx_km)

   if(maskType[1L]=="double"){
       h <- ifelse(max(pix$value) < maskThreshold,max(pix$value),maskThreshold)
       l <- ifelse(min(pix$value) > (1-maskThreshold),min(pix$value),1-maskThreshold)
       #pix <- filter(pix,value == max(value) | value ==min(value))
       pix <- filter(pix,value >=h | value <= l)
       pix <- split(pix, pix$value)
       pix <- lapply(pix,function(img, tiles){
           img <- inner_join(img,tiles, by = c("x","y"))
           img <- img$barcodes
           return(img)
       },tiles)
       nz <- sapply(pix,length) >0
       pix <- pix[nz]
   } else if(maskType[1L] =="pos"){
       h <- ifelse(max(pix$value) < maskThreshold,max(pix$value),maskThreshold)
       #pix <- filter(pix,value == max(value))
       pix <- filter(pix,value >= h)
       pix <- inner_join(pix,tiles, by = c("x","y"))
       pix <- list(pix$barcodes)
   } else {
       l <- ifelse(min(pix$value) > (1-maskThreshold),min(pix$value),1-maskThreshold)
       #pix <- filter(pix,value == min(value))
       pix <- filter(pix,value <= l)
       pix <- inner_join(pix,tiles, by = c("x","y"))
       pix <- list(pix$barcodes)
   }

   return(pix)
}


.kmeansMask <- function(img, tiles,maskThreshold,maskType,k=10){
    #--------------------------------------------------------------------------#
    # let's get all those segments using k means
    #--------------------------------------------------------------------------#

    cluster <- kmeans(as.vector(img),k)
    pix <- cluster$centers[cluster$cluster,]
    dim(pix) <- dim(img)
    pix <- as.data.frame(suppressWarnings(as.cimg(pix)))
    if(maskType[1L]=="double"){
        h <- ifelse(max(pix$value) < maskThreshold,max(pix$value),maskThreshold)
        l <- ifelse(min(pix$value) > (1-maskThreshold),min(pix$value),1-maskThreshold)
        #pix <- filter(pix,value == max(value) | value ==min(value))
        pix <- filter(pix,value >=h | value <= l)
        pix <- split(pix, pix$value)
        pix <- lapply(pix,function(img, tiles){
            img <- inner_join(img,tiles, by = c("x","y"))
            img <- img$barcodes
            return(img)
        },tiles)
        nz <- sapply(pix,length) >0
        pix <- pix[nz]
    } else if(maskType[1L] =="pos"){
        h <- ifelse(max(pix$value) < maskThreshold,max(pix$value),maskThreshold)
        #pix <- filter(pix,value == max(value))
        pix <- filter(pix,value >= h)
        pix <- inner_join(pix,tiles, by = c("x","y"))
        pix <- list(pix$barcodes)
    } else {
        l <- ifelse(min(pix$value) > (1-maskThreshold),min(pix$value),1-maskThreshold)
        #pix <- filter(pix,value == min(value))
        pix <- filter(pix,value <= l)
        pix <- inner_join(pix,tiles, by = c("x","y"))
        pix <- list(pix$barcodes)
    }
    return(pix)

}

#' Segment territories by iterative threshold-based segmentation.
#'
#' @param img a cimg object
#' @param tiles from vesalius object
#' @param elementSize an integer describing the size of the element for erosion and dilation. 
#' @param maskType a string describing mask type to use. "pos" segments territories with pixel values higher than the median, "neg" lower than the median and "double" on both sides. It is recommended to use "double".
#' @param distThresh a float describing the minimum mean pixel value distance of a territory from the image median. Should only be used in combination with quantile normalized embeddings.
#' @param sizeThresh an integer describing the minimum size of a territory in pixels.
#' @return Territories: as a list of a list of bar codes.
.thresholdMlMask <- function(img,
                             tiles,
                             elementSize = 5,
                             maskType="double",
                             distThresh = 0,
                             sizeThresh = 0) {
  # run the iterative image segmentation with erosion and dilation
  segImg <- .thresholdMlSmooth(img,
                               maskType,
                               elementSize,
                               distThresh,
                               sizeThresh)
  
  # convert the labeled image into a list of lists of barcodes
  pix <- left_join(tiles, as.data.frame(segImg), by = c("x", "y")) %>%
    distinct(barcodes, value) %>%
    na.exclude()
  
  pix <- filter(pix, value  != 0)
  pix <- split(pix, pix$value)
  
  pix <- lapply(
    pix,
    FUN = function(x)
      x$barcodes
  )
  return(pix)
}

#' Iterative threshold-based segmentation.
#'
#' @param img a cimg object
#' @param elementSize an integer describing the size of the element for erosion and dilation. 
#' @param maskType a string describing mask type to use. "pos" segments territories with pixel values higher than the median, "neg" lower than the median and "double" on both sides.
#' @param distThresh a float describing the minimum mean pixel value distance of a territory from the image median. Should only be used in combination with quantile normalized embeddings
#' @param sizeThresh an integer describing the minimum size of a territory in pixels.
#' @param dev logical: return segmented image before assigning labels to territories. Only for development purposes.
#' @return a cimg object. A segmented image with territories labeled continuously with positive integers and the background labeled as 0.
.thresholdMlSmooth <- function(img,
                               maskType = "double",
                               elementSize = 5,
                               distThresh = 0,
                               sizeThresh = 0) {
  
  # perform first round of image segmentation with 1 threshold.
  segImg1 <- ThresholdML(img, k = 1)
  
  # perform erosion and dilation to close small gaps and remove small artifacts
  segImg1 <- mclosing_square(segImg1, size = elementSize)
  segImg1 <- mopening_square(segImg1, size = elementSize)
  
  # label background as 0 and foreground (segmented territories) as 1
  # assumption that segmented area is smaller than the background.
  label1 <- ifelse(sum(segImg1) > (prod(dim(segImg1)) / 2), 0, 1)
  if (label1 == 0) {
    segImg1 <- abs(segImg1 - 1)
  }
  segment1 <- img[segImg1 == 1]
  # if (sum(segment1) == 0) {
  #   stop("No territories segmented. Possibly too little smoothing or too large structure element.\n")
  # }
  
  # replace the segmented territories with the image median before the second round of segmentation
  imgMedian <- median(as.matrix(img))
  maskedImg <- img
  maskedImg[segImg1 == 1] <- imgMedian
  
  # second round of segmentation, again with 1 threshold
  segImg2 <- ThresholdML(maskedImg, k = 1)
  
  # erosion and dilation
  segImg2 <- mclosing_square(segImg2, elementSize)
  segImg2 <- mopening_square(segImg2, elementSize)
  
  # check if segmented part is 0 or 1 labeled
  label2 <- ifelse(sum(segImg2) > (prod(dim(segImg2)) / 2), 0, 1)
  if (label2 == 0) {
    segImg2 <- abs(segImg2 - 1)
  }
  segment2 <- img[segImg2 == 1]
  # if (sum(segment2) == 0) {
  #   # warning("No territories segmented. Possibly too little smoothing or too large structure element.\n")
  #   return(segImg1)
  # }
  segImgComb <- segImg1
  # if segmentations are on both sides of the median, combine segmentations.
  # otherwise, only the first segmentation will be used.
  if (sign(mean(segment1) - mean(img)) != sign(mean(segment2) - mean(img))) {
    segImgComb[segImg2 == 1] <- max(segImgComb) + 1
  }
  
  # Give every separate territory it's own label
  segImgCombLab <- label(segImgComb)
  # reset inclusions to background label
  segImgCombLab[segImgComb == 0] <- 0
  
  # label territories continuously
  segImgCombLabCont <- segImgCombLab
  segImgCombLabCont[, ] <- 0
  iterator <- 0
  # iterate over labels
  for (l in seq(max(segImgCombLab))) {
    territoryMask <- segImgCombLab == l
    terMean <- mean(img[territoryMask])
    # in case label was an inclusion, it has already been removed
    if (sum(territoryMask) == 0) {
      next
      # Filter size of territories and distance from image median to obtain high information territories
    } else if (abs(terMean - imgMedian) < distThresh) {
      next
    } else if (sum(territoryMask) < sizeThresh) {
      next
    } else {
      iterator <- iterator + 1
    }
    
    # negative and positive option
    if (maskType == "double") {
      segImgCombLabCont[territoryMask] <- iterator
    } else if (maskType == "neg") {
      if (terMean < imgMedian) {
        segImgCombLab[territoryMask] <- 0
        segImgCombLabCont[territoryMask] <- iterator
      }
    } else if (maskType == "pos") {
      if (terMean > imgMedian) {
        segImgCombLab[territoryMask] <- 0
        segImgCombLabCont[territoryMask] <- iterator
      }
    }
  }
  return(segImgCombLabCont)
}


## Need stack check improvement 
.stackCheck <- function(pix,mosaic){
    #--------------------------------------------------------------------------#
    # create a size vector for the mosaic
    #--------------------------------------------------------------------------#
    mosaicSize <- table(mosaic$territory)
    #--------------------------------------------------------------------------#
    # Right next we need to make a comp between the size of each pixset
    #--------------------------------------------------------------------------#
    if(length(pix) ==1){
      count <- max(mosaic$territory)+1
    } else{
      count <- seq((max(mosaic$territory)+1),length.out = length(pix))
    }

    for(i in seq_along(pix)){
        tmp <- mosaic[mosaic$barcodes %in% pix[[i]],"territory"]
        tmp <- mosaicSize[as.character(tmp)]
        tmp <- tmp > length(pix[[i]])

        loc <- which(mosaic$barcodes %in% pix[[i]])[tmp]
        mosaic[loc,"territory"] <- count[i]
    }
    return(mosaic)

}
