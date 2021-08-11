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
#' data(Vesalius)
#' }

territoryMorphing <- function(territories,
                              morphologyFactor,
                              image,
                              verbose=TRUE){


      .morph(verbose)
      #------------------------------------------------------------------------#
      # First we define territory limits and add a little on each
      # side - this ensures that we won't be clipping any parts of the
      # territory
      #------------------------------------------------------------------------#
      seed <- territories %>% mutate(value=1)


      ymin <- ifelse((min(seed$y) - max(abs(morphologyFactor)) *2) <=0,1,
         min(seed$y) - morphologyFactor *2)
      xmin <- ifelse((min(seed$x) - max(abs(morphologyFactor)) *2) <=0,1,
         min(seed$x) - max(abs(morphologyFactor)) *2)
      ymax <- max(seed$y) + max(abs(morphologyFactor)) * 2
      xmax <- max(seed$x) + max(abs(morphologyFactor)) * 2
      #------------------------------------------------------------------------#
      # Now we convert add buffer boundaries, covert to grey scale a
      #------------------------------------------------------------------------#
      seed <- seed %>% select(c("x","y","cc","value")) %>%
              rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
              as.cimg() %>% grayscale()
      #------------------------------------------------------------------------#
      # Now we can do the morphing - grow, erode, clean and fill
      # For some reason I need to create mf seperately
      # maybe due to the fact that the grow/shrink function
      # is parsing weird stuff to imager/C++ ?
      #------------------------------------------------------------------------#
      for(i in seq_along(morphologyFactor)){
          mf <- abs(morphologyFactor[i])
          if(morphologyFactor[i] >=0){
              seed <- grow(seed,mf)
          } else {
              seed <- shrink(seed, mf)
          }
      }
      #------------------------------------------------------------------------#
      # Next we rebuild the image data frame with dilated
      #------------------------------------------------------------------------#
      seed <- seed %>% as.data.frame()
      seed <- inner_join(seed,image,by = c("x","y"))%>%
             select(c("barcodes","x","y","cc.y","z","tile",
                      "cluster","territory"))%>%
             filter(cc.y ==1)
      colnames(seed) <- c("barcodes","x","y","cc","value","tile",
                         "cluster","territory")

     return(seed)

}

#' layerTerritory.concave layers selected territories using concave hulling
#' @param image data.frame - Vesalius formatted data.frame (i.e. barcodes,
#' x,y,cc,value,cluster, territory) 
#' @param seedTerritory integer describing which territory should be layered
#' @param layerDepth integer describing maximum number of layers
#' @param concavity numeric describing concavity of the territory. 1 for highly
#' detailed structures. At higher values the territory will be considered as 
#' a convex hull. 
#' @param length_threshold numeric describing the minimum distance between 
#' points to be considered during hulling. High values = simpler shapes 
#' @param morphologyFactor integer or vector of integers describing growth 
#' and/or shrink extent.
#' @param captureRadius numeric - proportion of maximum distance between 
#' barcodes that will be used to pool barcodes together (range 0 - 1).
#' @param minBar integer - minimum number of barcodes allowed in each territory 
#' @param verbose logical - progress message output. 
#' @details To divide a territory into layers, one approach is to use concave 
#' hulling. This approach considers all barcodes locations (only center pixels)
#' as points on a 2D grid and extracts the "outer" layer of barcodes. 
#' 
#' This outer layer is extracted and defined as the first layer of the 
#' territory. The process is applied until no more barcodes can be pooled into
#' a layer. 
#' 
#' It should be noted that due to the geometrical nature of this approach, 
#' sub-territories are isolated prior to concave hulling to ensure that the 
#' hulling is only applied to contiguous barcodes (defined by captureRadius).
#' 
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","tile","cluster","territory" and "layer".
#' Layer describes the layer to which a barcode belongs. 
#' @examples 
#' \dontrun{
#' data(Vesalius)
#' }

layerTerritory.concave <- function(image,
                                   seedTerritory = NULL,
                                   layerDepth = NULL,
                                   concavity=1,
                                   length_threshold=0, 
                                   morphologyFactor = 3,
                                   captureRadius=0.2,
                                   minBar = 10, 
                                   verbose =TRUE){

    .simpleBar(verbose)

    #--------------------------------------------------------------------------#
    # Get pixels that are part of that territory
    #--------------------------------------------------------------------------#
    .seedSelect(verbose)
    ter <- image %>% filter(territory %in% seedTerritory)
    image <- filter(image,cc == 1)

    #--------------------------------------------------------------------------#
    # Splitting territory if it consists of multiple sub territories
    # Otherwise the concave hulling parameters become tricky to tune
    # not ideal - will need to fix and improve
    #--------------------------------------------------------------------------#

    pooled <- isolateTerritories.array(ter, method = c("distance"),
                                       captureRadius = captureRadius,
                                       global=TRUE,
                                       minBar = minBar,
                                       verbose =FALSE)
    colnames(pooled)[colnames(pooled) == "territory"] <- "subTerritory"

    pooled <- split(pooled, pooled$subTerritory)

    for(te in seq_along(pooled)){
      layer <- list()
      #------------------------------------------------------------------------#
      # Next for convenience we will convert to grey scale and dillate
      # as Default lets set dillation to 3 pixels
      #------------------------------------------------------------------------#
      .dilate(verbose)
      ter <- pooled[[te]]
      ter <- ter %>% select(c("x","y","cc","value")) %>%
             as.cimg %>% grayscale %>% grow(morphologyFactor)
      #------------------------------------------------------------------------#
      # Once we have dilated we can start iterative layayering
      # We will do this on tile centers and that is what we really care about
      # First let's convert back to a data frame and add back the barcodes
      #------------------------------------------------------------------------#
      .rebuildDF(verbose)

      ter <- as.data.frame(ter)
      ter <- right_join(ter,pooled[[te]],by = c("x","y"))
      ter <- inner_join(ter, image,by = c("x","y")) %>%
             select(c("barcodes.x","x","y","cluster.x",
                      "territory","subTerritory","tile.y"))
      colnames(ter) <-c("barcodes","x","y","cluster",
                        "territory","subTerritory","tile")
      ter <- ter %>% filter(tile == 1)

      #------------------------------------------------------------------------#
      # Next we will iteratively pool beads into single layers
      #------------------------------------------------------------------------#
      .layerTer(verbose)
      counter <- 1

      while(nrow(ter)>0){
          edge <- concaveman::concaveman(as.matrix(ter[,c("x","y")]),
                                concavity,length_threshold)
          colnames(edge) <- c("x","y")
          edge <- inner_join(as.data.frame(edge), ter, by = c("x","y"))
          edge$layer <- counter
          ter <- filter(ter,!barcodes %in% edge$barcodes)
          layer[[counter]] <- edge

          counter <- counter +1
      }
      layer <- do.call("rbind",layer)
      #------------------------------------------------------------------------#
      # Now we can pool layers together to equal the desired layer number
      # with a few checks for the number of layers
      #------------------------------------------------------------------------#
      layers <- unique(layer$layer)
      if(!is.null(layerDepth)){
          if(length(layers) < layerDepth){
              warning("Layer depth exceeds layers in Territory -
                       Using layers in territories", immediate. = TRUE)

          } else {
              idx <- floor(seq(1,length(layers), length.out = layerDepth +1))
              for(i in seq(1,length.out = layerDepth)){
                  layer$layer[layer$layer %in% seq(idx[i], idx[i+1])] <- i
              }
          }
      }

      pooled[[te]] <- layer

    }

    pooled <- do.call("rbind", pooled)

    .simpleBar(verbose)
    return(pooled)


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
#' data(Vesalius)
#' }

layerTerritory.edge <- function(image,
                                seedTerritory = NULL,
                                layerDepth = NULL,
                                morphologyFactor = 3,
                                verbose =TRUE){

    .simpleBar(verbose)

    #--------------------------------------------------------------------------#
    # Get pixels that are part of that territory
    # and prepare df for image transformations
    # adding extra layer so we don't loose any information
    # There is a whole lot of empty space but I'll keep for Now
    # Removing it is trivial but I need the original coordinates and barcodes
    # I will use the dilaltion factor as a way of expanding the image
    # Initial image set up
    #--------------------------------------------------------------------------#
    .seedSelect(verbose)
    ter <- image %>% filter(territory %in% seedTerritory) %>% mutate(value=1)

    ymin <- ifelse((min(ter$y) - max(abs(morphologyFactor)) *2) <=0,1,
        min(ter$y) - max(abs(morphologyFactor)) *2)
    xmin <- ifelse((min(ter$x) - max(abs(morphologyFactor)) *2) <=0,1,
        min(ter$x) - max(abs(morphologyFactor)) *2)
    ymax <- max(ter$y) + max(abs(morphologyFactor)) * 2
    xmax <- max(ter$x) + max(abs(morphologyFactor)) * 2


    #--------------------------------------------------------------------------#
    # Dilate territory to ensure that we cover the outer layers as well
    #--------------------------------------------------------------------------#
    .morph(verbose)
    terTmp <- ter %>% select(c("x","y","cc","value")) %>%
              rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
              as.cimg %>% grayscale()
    for(i in seq_along(morphologyFactor)){
        mf <- abs(morphologyFactor[i])
        if(morphologyFactor[i] >=0){
            terTmp <- grow(terTmp,mf)
        } else {
            terTmp <- shrink(terTmp, mf)
        }
    }

    terForLoop <- terTmp %>% as.data.frame()

    terForLoop <- inner_join(terForLoop,image,by = c("x","y"))%>%
                  select(c("barcodes","x","y","cc.y","z",
                           "tile","cluster","territory"))%>%
                  filter(cc.y ==1)
    colnames(terForLoop) <- c("barcodes","x","y","cc","value",
                              "tile","cluster","territory")
    ter <- terForLoop
    colnames(ter) <- c("barcodes","x","y","cc","value","tile",
                       "cluster","territory")


    #--------------------------------------------------------------------------#
    # Now we can get edges of shape and compare this to tiles
    # and pool this "edge" into layers
    #--------------------------------------------------------------------------#

    .layerTer(verbose)
    counter <- 1
    layer <- list()
    while(nrow(terForLoop)>0){
      grad <- terTmp  %>%
              imgradient("xy") %>%
              enorm() %>%
              add() %>%
              sqrt() %>%
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

      terTmp <- terForLoop %>% select(c("x","y","cc","value")) %>%
                rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
                as.cimg
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

    ter$layer <-0
    for(lay in seq_along(layer)){

        ter$layer[ter$barcodes %in% layer[[lay]]] <- lay
    }

    #--------------------------------------------------------------------------#
    # Finally we can split the different layers if we want to combine
    #--------------------------------------------------------------------------#
    layers <- unique(ter$layer)
    if(!is.null(layerDepth)){
        if(length(layers) < layerDepth){
            warning("Layer depth exceeds layers in Territory - 
                     Using layers in territories", immediate. = TRUE)

        } else {
            idx <- floor(seq(1,length(layers), length.out = layerDepth +1))
            for(i in seq(1,length.out = layerDepth)){
                ter$layer[ter$layer %in% seq(idx[i], idx[i+1])] <- i
            }
        }
    }
    .simpleBar(verbose)

    return(ter)
}
