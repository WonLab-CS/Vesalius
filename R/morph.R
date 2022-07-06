################################################################################
############################   ST package        ###############################
################################################################################
# Territory morphology


#' territoryMorphing applies morphological operators to a set of territoriees
#' @param territories data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,cluster, territory) only containing territories of interest.
#' @param morphologyFactor integer or vector of integers describing growth
#' and/or shrink extent.
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,cluster, territory) containing all territories.
#' @param verbose logical - progress message output.
#' @details Territory morphing can manipulate territories by growing, shrinking,
#' filling, and cleaning territories.
#' Growing = Positive integers - Territory will be dilated by x number of pixels
#' Shrinking = Negative integers - Territory will be contracted by x number of
#' pixels
#' Filling = grow followed by shrink.
#' Cleaning = shrink followed by grow.
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile","cluster", and "territory".
#' This data.frame will contain the new territory after morphological operators
#' and will contain barcodes associated with other territories and colour
#' segments.
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' ter <- image %>% filter(territory == 1)
#' # grow
#' g <- territoryMorphing(ter,morphologyFactor = 5, image = image)
#' # shrink
#' s <- territoryMorphing(ter,morphologyFactor = -5, image = image)
#' # fill
#' f <- territoryMorphing(ter,morphologyFactor = c(5,-5), image = image)
#' # clean
#' c <- territoryMorphing(ter,morphologyFactor = c(-5,5), image = image)
#' }

territoryMorphing <- function(vesalius,
                              territory = NULL,
                              trial = "last",
                              morphologyFactor = 0,
                              verbose=TRUE){
      .simpleBar(verbose)
      #------------------------------------------------------------------------#
      # lets make a few checks
      # for now we will assume that either segments or territories can be used
      # TOADD methods to access parameters associated with each trial run
      #------------------------------------------------------------------------#
      if(is.null(territory)){
          stop("No specified territory for territory morphing!")
      } else {
          adjustedTerritoryValues <- min(territory)
      }
      if(!is.null(vesalius@territories) & trial == "last"){
        trial <- colnames(vesalius@territories)[ncol(vesalius@territories)]
        ter <- vesalius@territories[,c("barcodes","x","y",trial)]
        colnames(ter) <- c("barcodes","x","y","trial")
        buffer <- ter
        ter <- filter(ter,trial %in% territory) %>%
               left_join(vesalius@tiles, by = "barcodes") %>%
               mutate(value = 1) %>%
               select(c("barcodes","x.y","y.y","trial","origin","value"))
        colnames(ter) <- c("barcodes","x","y","trial","origin","value")
      } else if(!is.null(vesalius@territories) & trial != "last") {
        if(sum(grepl(x = colnames(vesalius@territories),pattern = trial))==0){
            stop(paste(deparse(substitute(trial)),"is not in territory data frame"))
        }
        ter <- vesalius@territories[,c("barcodes","x","y",trial)]
        colnames(ter) <- c("barcodes","x","y","trial")
        buffer <- ter
        ter <- filter(ter,trial %in% territory) %>%
               left_join(vesalius@tiles, by = "barcodes") %>%
               mutate(value = 1) %>%
               select(c("barcodes","x.y","y.y","trial","origin","value"))
        colnames(ter) <- c("barcodes","x","y","trial","origin","value")
      } else {
        stop("No territories have been computed!")
      }
      .morph(verbose)
      #------------------------------------------------------------------------#
      # getting last name if any to create new column
      # Later it would be good to write these sections as utility functions
      # they often do the same thing and it make things cleaner
      # not high priority - for now we need to make it work
      #------------------------------------------------------------------------#

      if(any(grepl(x = colnames(vesalius@territories),pattern = "Morphology"))){
        previous <- grep(x=colnames(vesalius@territories),
                                      pattern = "Morphology",
                                      value = TRUE)
        m <- gregexpr('[0-9]+', previous)
        last <- max(as.numeric(unlist(regmatches(previous,m))))

        newTrial <- paste0("Morphology_Trial_",last+1)


      } else {
        newTrial <- "Morphology_Trial_1"

      }
      #------------------------------------------------------------------------#
      # First we define territory limits and add a little on each
      # side - this ensures that we won't be clipping any parts of the
      # territory
      #------------------------------------------------------------------------#
      ymin <- ifelse((min(ter$y) - max(abs(morphologyFactor)) *2) <=0,1,
         min(ter$y) - morphologyFactor *2)
      xmin <- ifelse((min(ter$x) - max(abs(morphologyFactor)) *2) <=0,1,
         min(ter$x) - max(abs(morphologyFactor)) *2)
      ymax <- max(ter$y) + max(abs(morphologyFactor)) * 2
      xmax <- max(ter$x) + max(abs(morphologyFactor)) * 2
      #------------------------------------------------------------------------#
      # Now we convert add buffer boundaries, covert to grey scale a
      #------------------------------------------------------------------------#
      ter <- ter %>% select(c("x","y","value")) %>%
              rbind(.,c(xmin,ymin,1),c(xmax,ymax,1)) %>%
              as.cimg()
      #------------------------------------------------------------------------#
      # Now we can do the morphing - grow, erode, clean and fill
      # For some reason I need to create mf seperately
      # maybe due to the fact that the grow/shrink function
      # is parsing weird stuff to imager/C++ ?
      #------------------------------------------------------------------------#
      
      for(i in seq_along(morphologyFactor)){
          mf <- abs(morphologyFactor[i])
          if(morphologyFactor[i] >=0){
              ter <- grow(ter,mf)
          } else {
              ter <- shrink(ter, mf)
          }
      }
      
      
      #------------------------------------------------------------------------#
      # Next we rebuild the image data frame with dilated
      #------------------------------------------------------------------------#
      ter <- ter %>% as.data.frame() 
      ter <- inner_join(ter,vesalius@tiles,by = c("x","y")) %>%
             filter(origin ==1) 
      
      buffer$trial[buffer$barcodes %in% ter$barcodes] <- adjustedTerritoryValues
      colnames(buffer) <- c(colnames(buffer)[seq(1,ncol(buffer)-1)],newTrial)
      
      vesalius <- .updateVesalius(vesalius=vesalius,
                                  data=buffer,
                                  slot="territories",
                                  commit = as.list(match.call()),
                                  defaults = as.list(args(territoryMorphing)),
                                  append=TRUE)
    .simpleBar(verbose)
     return(vesalius)

}



#' extractIdentity compute differentially expressed genes for each territory
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,cluster, territory)
#' @param seedTerritory integer describing which territory should be layered
#' @param layerDepth integer describing maximum number of layers
#' @param morphologyFactor integer or vector of integers describing growth
#' and/or shrink extent.
#' @param verbose logical - progress message output.
#' @details To divide a territory into layers, one approach is to use edge
#' detection. Terriotries are converted to a black and white image containing
#' only the territories of interest. Sobel edge detection is applied along the
#' x and y axis and all barcodes sharing a pixel with the edge are pooled into
#' a layer. The process is applied until no more barcodes can be pooled into
#' a layer.
#'
#' Morphological operators are applied to the isolated territory prior to
#' layering (see \code{territoryMorphing}).
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile","cluster","territory" and "layer".
#' Layer describes the layer to which a barcode belongs.
#' @examples
#' \dontrun{
#' data("vesalius")
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' layer <- layerTerritory.edge(image, seedTerritory = 1)
#' }

layerTerritory <- function(vesalius,
                           territory = NULL,
                           trial = "last",
                           layerDepth = NULL,
                           morphologyFactor = 0,
                           verbose =TRUE){

    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # This might be a bit messier but i'm just going to make a call to
    # territoryMorphing function
    # this has some un necssary transformations but this is cleaner
    # We check first if we need to run this if no morphological operator
    # If we run this then we always take the last trial run
    # this will update the morphology only if it is run
    #--------------------------------------------------------------------------#
    if(any(morphologyFactor !=0)){
       .morph(verbose)
        vesalius <- territoryMorphing(vesalius,
                                      territory,
                                      trial,
                                      morphologyFactor,
                                      verbose =FALSE)
        trial <- "last"
        territory <- min(territory)
    }

    if(!is.null(vesalius@territories) & trial == "last"){
      trial <- colnames(vesalius@territories)[ncol(vesalius@territories)]
      ter <- vesalius@territories[,c("barcodes","x","y",trial)]
      colnames(ter) <- c("barcodes","x","y","trial")
      buffer <- ter
      ter <- right_join(ter,vesalius@tiles, by = "barcodes") %>%
             filter(trial %in% territory) %>%
             mutate(value = 1)%>%
             select(c("barcodes","x.y","y.y","value","origin","trial"))
      colnames(ter) <- c("barcodes","x","y","value","origin","trial")
      terForLoop <- ter
    } else if(!is.null(vesalius@territories) & trial != "last") {
      if(!any(grepl(x = colnames(vesalius@territories),pattern = trial))){
          stop(paste(deparse(substitute(trial)),"is not in territory data frame"))
      }
      ter <- vesalius@territories[,c("barcodes","x","y",trial)]
      colnames(ter) <- c("barcodes","x","y","trial")
      buffer <- ter
      ter <- filter(ter,trial %in% territory) %>%
             right_join(vesalius@tiles, by = "barcodes") %>%
             mutate(value = 1) %>%
             select(c("barcodes","x.y","y.y","value","origin","trial"))
      colnames(ter) <- c("barcodes","x","y","value","origin","trial")
      terForLoop <- ter
    } else {
      stop("No territories have been computed!")
    }

    #------------------------------------------------------------------------#
    # getting last name if any to create new column
    # This will yield two new columns
    #------------------------------------------------------------------------#

    if(any(grepl(x = colnames(vesalius@territories),pattern = "Layer"))){
      previous <- grep(x=colnames(vesalius@territories),
                                    pattern = "Layer",
                                    value = TRUE)
      m <- gregexpr('[0-9]+', previous)
      last <- max(as.numeric(unlist(regmatches(previous,m))))

      newTrial <- paste0("Layer_Trial_",last+1)


    } else {
      newTrial <- "Layer_Trial_1"

    }

    #------------------------------------------------------------------------#
    # First we define territory limits and add a little on each
    # side - this ensures that we won't be clipping any parts of the
    # territory
    #------------------------------------------------------------------------#
    ymin <- ifelse((min(ter$y) - max(abs(morphologyFactor)) *2) <=0,1,
       min(ter$y) - morphologyFactor *2)
    xmin <- ifelse((min(ter$x) - max(abs(morphologyFactor)) *2) <=0,1,
       min(ter$x) - max(abs(morphologyFactor)) *2)
    ymax <- max(ter$y) + max(abs(morphologyFactor)) * 2
    xmax <- max(ter$x) + max(abs(morphologyFactor)) * 2
    #------------------------------------------------------------------------#
    # Now we convert add buffer boundaries, covert to grey scale a
    #------------------------------------------------------------------------#
    ter <- ter %>% select(c("x","y","value")) %>%
            rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
            as.cimg()
    #--------------------------------------------------------------------------#
    # Now we can get edges of shape and compare this to tiles
    # and pool this "edge" into layers
    #--------------------------------------------------------------------------#

    .layerTer(verbose)
    counter <- 1
    layer <- list()

    while(nrow(terForLoop)>0){

      grad <- ter  %>%
              .detectEdges() %>%
              grow(1) %>%
              as.cimg() %>%
              as.data.frame() %>%
              filter(value >0)


      #------------------------------------------------------------------------#
      # getting barcodes from territory
      #------------------------------------------------------------------------#

      edge <- inner_join(grad,terForLoop, by = c("x","y")) %>%
                        select(c("barcodes"))

      #------------------------------------------------------------------------#
      # Resizing ter - removing barcodes that are part of the edge
      #------------------------------------------------------------------------#

      terForLoop <- filter(terForLoop, !barcodes %in% unique(edge$barcodes))

      #------------------------------------------------------------------------#
      # Rebuilding an image but adding a little extra space
      #------------------------------------------------------------------------#
      if(nrow(terForLoop)>0){
      ymin <- ifelse((min(terForLoop$y) - max(abs(morphologyFactor)) *2) <=0,1,
          min(terForLoop$y) - max(abs(morphologyFactor)) *2)
      xmin <- ifelse((min(terForLoop$x) - max(abs(morphologyFactor)) *2) <=0,1,
          min(terForLoop$x) - max(abs(morphologyFactor)) *2)
      ymax <- max(terForLoop$y) + max(abs(morphologyFactor)) * 2
      xmax <- max(terForLoop$x) + max(abs(morphologyFactor)) * 2

      ter <- terForLoop %>% select(c("x","y","value")) %>%
             rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
             as.cimg()
      }
      #------------------------------------------------------------------------#
      # Adding edge to layer list and counting up
      #------------------------------------------------------------------------#
      layer[[counter]] <- unique(edge$barcodes)
      counter <- counter +1

    }

    #--------------------------------------------------------------------------#
    # Now we can add the layers to the original territory
    # A rename layer if required
    #--------------------------------------------------------------------------#

    buffer$trial <- "out"
    for(lay in seq_along(layer)){

        buffer$trial[buffer$barcodes %in% layer[[lay]]] <- lay
    }

    #--------------------------------------------------------------------------#
    # Finally we can split the different layers if we want to combine
    #--------------------------------------------------------------------------#
    layers <- unique(buffer$trial)
    if(!is.null(layerDepth)){
        if(length(layers) < layerDepth){
            warning("Layer depth exceeds layers in Territory -
                     Using layers in territories", immediate. = TRUE)

        } else {
            idx <- floor(seq(1,length(layers), length.out = layerDepth +1))
            for(i in seq(1,length.out = layerDepth)){
                buffer$trial[buffer$trial %in% seq(idx[i], idx[i+1])] <- i
            }
        }
    }
    #--------------------------------------------------------------------------#
    # rename new column
    colnames(buffer) <- c(colnames(buffer)[seq(1,ncol(buffer)-1)],newTrial)

    vesalius <- .updateVesalius(vesalius=vesalius,
                                data=buffer,
                                slot="territories",
                                commit = as.list(match.call()),
                                defaults = as.list(args(layerTerritory)),
                                append=TRUE)
  .simpleBar(verbose)
   return(vesalius)
}
