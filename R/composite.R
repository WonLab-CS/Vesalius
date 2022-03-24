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
                        verbose = TRUE){
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
    .eq(verbose)
    eqThresh <- (1-maskThreshold)
    imgList <- equalizeHistogram(imgList,
                                 sleft = eqThresh,sright=eqThresh,
                                 verbose=FALSE,
                                 cores = cores)
    #--------------------------------------------------------------------------#
    # We smooth all arrays -
    # So this goes against the idea of not keeping arrays in memory -
    # At least for this function we don't need to but at the moment lets just
    # check the proof of concept
    #--------------------------------------------------------------------------#
    .reg(verbose)
    imgList <- parallel::mclapply(imgList,.regularise,
                                 lambda = lambda,
                                 mc.cores = cores)
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
        mosaic <- .mosaic(img,tiles,mosaic,maskThreshold,maskType,method,k)
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
                    k=10){
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
                    "KmeansIter" = .kmeansIterMask(img,tiles,maskThreshold))
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
