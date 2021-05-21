################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#


#' Smoothing image array - wrapper function for imager smoothing functions
#' @param img 4 dimensional array or cimg object
#' @param method character describing smoothing method to use "median" , "iso"  or "box"
#' @param sigma numeric - standard deviation associated with the isoblur
#' @param box numeric describing box size (centered around center pixel) for smoothing
#' @param threshold numeric - discard pixels that are too low in value (see \code{medianblur})
#' @param neuman logical describing If Neumann boundary conditions should be used,
#' Dirichlet otherwise (default true, Neumann)
#' @param gaussian logical - use gaussian filter
#' @param na.rm logical describing if NA values should be removed
#' @param iter numeric describing number of smoothing rounds to apply to image
#' @details Applying smoothing to image data.
#' @return Array or cimg object after bieng smoothed
smoothArray <- function(img,method = c("median","iso","box"),
                        sigma = 1, box = 20, threshold=0, neuman=TRUE,
                        gaussian=TRUE,na.rm=FALSE, acrossLevels = "min",
                        iter = 1,verbose = TRUE,invert = FALSE){
    #--------------------------------------------------------------------------#
    # Class check
    #--------------------------------------------------------------------------#

    if(!is.cimg(img)){
      imgCopy <- as.cimg(select(img, c("x", "y", "cc", "value")))
    }else {
      imgCopy <- img
    }

    if(invert){
      #------------------------------------------------------------------------#
      ## This is effing clonky
      ## Probably need to re work this
      ## It works but it just sucks that I have to rebuild the whole thing
      ## Unless.... I "hack" the as.cimg function and add an argument myself?
      ## That just sounds like a stupid idea with extra steps.
      #------------------------------------------------------------------------#
      .invertCols(verbose)
      imgCopy <- img %>%
             select(c("x","y","cc","value")) %>%
             as.cimg %>%
             as.data.frame
      nonImg <- paste0(img$x,"_",img$y)
      inImg <- paste0(imgCopy$x,"_",imgCopy$y)
      imgCopy[!inImg %in% nonImg,"value"] <- 1
      imgCopy <- as.cimg(imgCopy)
    }

    #--------------------------------------------------------------------------#
    # Running multiple smoothing itteration
    # TODO change the verbose part to account for multiple signa/box values
    #--------------------------------------------------------------------------#
    for(i in seq_len(iter)){
        .smooth(i,verbose)
        for(j in method){
          imgCopy <- switch(j,
                            "median" = map_il(box,~medianblur(imgCopy,.,threshold)),
                            "iso" = map_il(sigma,~isoblur(imgCopy,.,neuman,gaussian,na.rm)),
                            "box" = map_il(box,~boxblur(imgCopy,.,neuman)))
          imgCopy <- switch(acrossLevels,
                            "min" = parmin(imgCopy),
                            "max" = parmax(imgCopy),
                            "mean" = average(imgCopy))
        }

    }

    imgCopy <- as.data.frame(imgCopy)



    img <- right_join(imgCopy, img, by  = c("x","y","cc")) %>%
           select(c("barcodes","x","y","cc","value.x","tile")) %>% tibble

    colnames(img) <- c("barcodes","x","y","cc","value","tile")
    return(img)
}






equalizeHistogram <- function(image, type = "BalanceSimplest",N=1,smax=1,
                         sleft =1,sright =1,lambda=0.1,up=100, down = 10,
                         range = c(0,1),invert=FALSE,verbose=TRUE){
    .simpleBar(verbose)
    if(invert){
      #------------------------------------------------------------------------#
      ## This is effing clonky
      ## Probably need to re work this
      ## It works but it just sucks that I have to rebuild the whole thing
      ## Unless.... I "hack" the as.cimg function and add an argument myself?
      ## That just sounds like a stupid idea with extra steps.
      #------------------------------------------------------------------------#
      .invertCols()
      img <- image %>%
             select(c("x","y","cc","value")) %>%
             as.cimg %>%
             as.data.frame
      nonImg <- paste0(image$x,"_",image$y)
      inImg <- paste0(img$x,"_",img$y)
      img[!inImg %in% nonImg,"value"] <- 1
      img <- img %>% as.cimg() %>% imsplit("c")
    } else {
      img <-image %>% select(c("x","y","cc","value")) %>% as.cimg %>% imsplit("c")
    }
    .eq(verbose)
    #--------------------------------------------------------------------------#
    # Equalizing histogram - dpending on the image different methods will
    # work differently. It worth keeping in mind that these are made and
    # optimised for real images. BalanceSimplest seems to be the best for our
    # pseudo images...
    #--------------------------------------------------------------------------#
    img <- switch(type,
                    "EqualizePiecewise" = lapply(img,EqualizePiecewise,N) %>% imappend("c"),
                    "BalanceSimplest" = lapply(img,BalanceSimplest,sleft,sright,range) %>% imappend("c"),
                    "SPE" = lapply(img,SPE,lambda) %>% imappend("c"),
                    "EqualizeDP"  = lapply(img,EqualizeDP,down,up) %>% imappend("c"),
                    "EqualizeADP" = lapply(img,EqualizeADP) %>% imappend("c"),
                    "ECDF" = lapply(img,ecdf.eq) %>% imappend("c"))
    .rebuildDF(verbose)
    img <- as.data.frame(img)
    img <- right_join(img, image, by  = c("x","y","cc")) %>%
           select(c("barcodes","x","y","cc","value.x","tile")) %>% tibble

    colnames(img) <- c("barcodes","x","y","cc","value","tile")
    .simpleBar(verbose)
    return(img)
}


ecdf.eq <- function(im){
   return(as.cimg(ecdf(im)(im),dim=dim(im)))
}


#' iterativeSegmentation - collapse and smooth points into colour startingSegements
#' @param img data frame containing RGB code for each bead
#' @param colDepth numeric describing final number of colour segments in image
#' @param segIter numeric describing number of segmentation iterations
#' @param smoothIter numeric describing number of smoothing iterations
#' @param startingSegements numeric describing the number of initial colour segments to consider.
#' @param method character describing the smoothing method to be used (median , iso or box)
#' @param sigma numeric - standard deviation associated with the isoblur
#' @param box numeric describing box size (centered around center pixel) for smoothing
#' @param threshold numeric - discard pixels that are too low in value (see \code{medianblur})
#' @param neuman logical describing If Neumann boundary conditions should be used,
#' Dirichlet otherwise (default true, Neumann)
#' @param gaussian logical - use gaussian filter
#' @param na.rm logical describing if NA values should be removed
#' @param verbose logical indicating if progress messages should be outputed to console.
#' @details For each iteration, the algorithm collapse the colour space into n colour segments.
#' Then, smoothing is applied to nearest neighbours. The process is repeated for each iteration.
#' Note that the number of segments from starting segments to final colour depth follows the following
#' Function : \code{seq(startingSegements,colDepth+1,length.out = iter)}
#' @return a data frame containing RGB code for each bead.

iterativeSegmentation.array <- function(img,colDepth = 10,
                                        smoothIter = 1,
                                        method = c("median","iso","box"),
                                        acrossLevels = "min",
                                        sigma = 1, box = 20, threshold=0, neuman=TRUE,
                                        gaussian=TRUE, useCenter=TRUE, na.rm=FALSE,
                                        invert = FALSE,
                                        verbose = TRUE){

  .simpleBar(verbose)
  #----------------------------------------------------------------------------#
  # Comverting image to data.frame
  # Need to check how things work here
  ## could be cleaner if you just go back and forth between data frame and cimg
  #----------------------------------------------------------------------------#
  if(is.cimg(img)){
      img <- as.data.frame(img)
  }

  #----------------------------------------------------------------------------#
  # Segmenting image by iteratively decreasing colour depth and smoothing
  #----------------------------------------------------------------------------#
  for(j in seq_along(colDepth)){



    #--------------------------------------------------------------------------#
    # Segmentation is done by using kmeans clustering
    # We will only run clustering on 1 pixel of each tile
    # This will ensure faster run time and cleaner results
    #--------------------------------------------------------------------------#
    img <- smoothArray(img,method =method,sigma = sigma, box = box,
                       threshold=threshold, neuman=neuman, gaussian=gaussian,
                       na.rm=FALSE,acrossLevels=acrossLevels,
                       iter = smoothIter,invert = invert,verbose=verbose)


    cat("\n")
    .seg(j,verbose)
    cat("\n")
    if(useCenter){
      tmpImg <- img %>% filter(tile == 1) %>% group_by(cc) %>% distinct(barcodes,.keep_all = TRUE)

      #--------------------------------------------------------------------------#
      # clustering with each colour together
      # It cause issues down the line
      # How do you recombine unique colours backtogether
      # If you can find a new way of doing iot great but in the mean time...
      #--------------------------------------------------------------------------#
      colours <- select(tmpImg, c("cc","value"))
      colours <- data.frame(colours$value[colours$cc ==1],
                            colours$value[colours$cc ==2],
                            colours$value[colours$cc ==3])
      km <- kmeans(colours,colDepth[j],iter.max = 200,nstart = 50)
      cluster <- km$cluster
      Kcenters <- km$centers
      tmpImg$cluster <- 0

      for(i in seq(1,3)){
          tmpImg$cluster[tmpImg$cc == i] <- cluster
          tmpImg$value[tmpImg$cc == i] <- Kcenters[cluster,i]
      }


      #--------------------------------------------------------------------------#
      # Replacing values in original image
      #--------------------------------------------------------------------------#

      img <- inner_join(img,tmpImg, by = c("barcodes","cc")) %>%
             select(c("barcodes","x.x","y.x","cc","value.y","tile.x","cluster"))

      colnames(img) <- c("barcodes","x","y","cc","value","tile","cluster")
    } else {
      colours <- select(img, c("cc","value"))
      colours <- data.frame(colours$value[colours$cc ==1],
                            colours$value[colours$cc ==2],
                            colours$value[colours$cc ==3])
      km <- kmeans(colours,colDepth[j],iter.max = 200,nstart = 10)

      cluster <- km$cluster
      Kcenters <- km$centers
      img$cluster <- 0

      for(i in seq(1,3)){
          img$cluster[img$cc == i] <- cluster
          img$value[img$cc == i] <- Kcenters[cluster,i]
      }
      img <- img %>% group_by(cc,barcodes) %>%
           mutate(value = mean(value),cluster = .top_cluster(cluster)) %>%
           ungroup

    }
  }

  .simpleBar(verbose)
  return(img)
}

.intKmeans <- function(img,depth,cluster = TRUE){
    if(cluster){
      k <- kmeans(img,depth, iter.max = 200)$cluster
    } else {
      k <- kmeans(img,depth, iter.max = 200)
      k <- as.vector(k$centers)[k$cluster]
    }
    return(k)
}

.top_cluster <- function(cluster){
    top <- table(cluster)
    top <- names(top)[order(top, decreasing = T)]
    return(as.numeric(top[1]))
}



#' isolating territories from spatial transcriptomic data
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @param dropRadius numeric describing the proportion of total distance to consider when pooling beads together
#' @details In order to select territories, colours will be collapsed into \code{colDepth} number of
#' dominant colours. For each colour, territories will be isolated based on distance between point.
#' The \code{dropRadius} represents the proportion of total distance (x and y values may differ) to use in order
#' to pool beads/spots together. The algorithm selects a point and pools all neighboring beads into a territory.
#' The process is applied to the nearest neighbors until no more points can be pooled into the territory.
#' If beads remain for a given colour, the process is repeated until all beads/spots are put into a territory.
#' @return a data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster.

isolateTerritories.array <- function(img, method = c("distance"),
                                     captureRadius = 0.025,global=TRUE,
                                     minCell = 10,verbose =TRUE){

    #--------------------------------------------------------------------------#
    # Compute real capture Radius
    #--------------------------------------------------------------------------#
    .simpleBar(verbose)
    if(method[1L] == "distance"){
      captureRadius <- sqrt((max(img$x)-min(img$x))^2 +
                            (max(img$y)-min(img$y))^2) * captureRadius
    }


    #--------------------------------------------------------------------------#
    # IsolatingTerritories dispatching
    #--------------------------------------------------------------------------#
    clust <- unique(img$cluster)

    img$territory <- 0
    for(i in seq_along(clust)){
      .terPool(i,verbose)
      #------------------------------------------------------------------------#
      # We don't need all data in this case only clusters and locations
      # We can just rebuild everything afterwards
      # At least we don't don't need to compute anything unnecessarily
      ## Note other method not in use for now
      # Might not be worthwile to implement them
      # Argument could be removed
      #------------------------------------------------------------------------#

      tmp <- filter(img,cluster == clust[i] & cc == 1)

      ter <- switch(method[1L],
                    "distance" = .distancePooling.array(tmp,captureRadius,minCell),
                    "neighbor" = .neighborPooling.array(tmp,captureRadius),
                    "watershed" = .watershedPooling.array(tmpImg))
      #------------------------------------------------------------------------#
      # Skipping if does not meet min cell requirements
      # filtering can be done later
      #------------------------------------------------------------------------#
      if(is.null(ter))next()
      #------------------------------------------------------------------------#
      # adding territories to each channel
      #------------------------------------------------------------------------#
      for(j in seq(1,3)){
          img$territory[img$cluster == clust[i] & img$cc == j] <- ter
      }
    }

    #--------------------------------------------------------------------------#
    # Filtering out territories that do not meat he minimum number of cell
    # criteria - by default they should have a territory value of 0
    #--------------------------------------------------------------------------#
    img <- img %>% filter(territory!=0)

    #--------------------------------------------------------------------------#
    # Globalise territories crate a numbering system for all colour clusters
    # If set to false territories only describe territories for any given colour
    # cluster
    #--------------------------------------------------------------------------#

    if(global){
        img <- .globaliseTerritories(img)
    }
    cat("\n")
    .simpleBar(verbose)
    return(img)
}

.globaliseTerritories <- function(img,seurat=FALSE){
    if(!seurat){
      imgTmp <- img %>% filter(territory != "isolated")
      ter <- paste0(imgTmp$cluster,"_", imgTmp$territory)
      allTer <- unique(ter)
      ter <- seq_along(allTer)[match(ter,allTer)]
      img$territory[img$territory != "isolated"] <- ter
      return(img)
    } else {
      imgTmp <- img %>% filter(territory != "isolated")
      ter <- paste0(img$seurat_clusters,"_", img$territory)
      allTer <- unique(ter)
      ter <- seq_along(allTer)[match(ter,allTer)]
      img$seurat_clusters[img$territory != "isolated"] <- ter
      return(img)
    }


}


.distancePooling.array <- function(img,captureRadius,minCell){
    #--------------------------------------------------------------------------#
    # Select center point of each tile for only one channel
    # Dont need to run it for all channels
    # !!!!! Hot fix - this will need to chnaged upstream !!!!!
    # it is possible to have a colour cluster that does not have any center beads
    # This is generally an artifact of smoothing and not using center for the
    # segmentation. This mainly occurs when colours are inverted
    #--------------------------------------------------------------------------#

    imgCopy <- img %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)
    if(nrow(imgCopy)<1){return(NULL)}


    #--------------------------------------------------------------------------#
    # Compute distances
    #--------------------------------------------------------------------------#
    idx <- seq_len(nrow(imgCopy))
    distanceMatrix <- lapply(idx, function(idx,mat){
                            xo <- mat$x[idx]
                            yo <- mat$y[idx]
                            xp <- mat$x
                            yp <- mat$y
                            distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)
                            return(distance)
    }, imgCopy)


    #--------------------------------------------------------------------------#
    # Buildan actual matrix
    #--------------------------------------------------------------------------#

    distanceMatrix <- do.call("rbind",distanceMatrix)
    colnames(distanceMatrix) <- imgCopy$barcodes
    rownames(distanceMatrix) <- imgCopy$barcodes

    #--------------------------------------------------------------------------#
    # If there is only one point
    # In this case you only have one territory as well
    # Don't need to any pooling
    #--------------------------------------------------------------------------#

    if(sum(dim(distanceMatrix))==2){
        return(1)
    }


    #--------------------------------------------------------------------------#
    # Pooling points together
    #--------------------------------------------------------------------------#
    barcodes <- imgCopy$barcodes
    territories <- list()
    count <- 1

    while(length(barcodes) >0){
          tmp <- distanceMatrix[,sample(barcodes,1)]
          pool <- names(tmp)[tmp <= captureRadius]
          inter <- pool
          converge <- FALSE

          while(!converge){
              if(length(inter)==1){
                  territories[[count]] <- pool
                  barcodes <- barcodes[!barcodes %in% pool]
                  count <- count + 1
                  converge <- TRUE
              } else {
                  newPool <- distanceMatrix[,inter]

                  newPool <- unique(unlist(lapply(seq_len(ncol(newPool)), function(idx,np,captureRadius){

                                            res <- rownames(np)[np[,idx] <= captureRadius]
                                            return(res)
                                          },newPool,captureRadius)))

                  overlap <- newPool %in% pool

                  if(sum(overlap) != length(newPool)){

                        pool <- unique(c(pool,newPool[!overlap]))
                        inter <- unique(newPool[!overlap])
                  } else {
                        territories[[count]] <- pool
                        count <- count +1
                        barcodes <- barcodes[!barcodes %in% pool]

                        converge <- TRUE

                  }
              }
          }
      }
    #--------------------------------------------------------------------------#
    # Clean up drop outs
    #--------------------------------------------------------------------------#

      allTers <- img$territory

      for(ter in seq_along(territories)){
            loc <- img$barcodes %in% territories[[ter]]
            if(length(territories[[ter]]) <= minCell){
               allTers[loc] <- "isolated"
            } else {
               allTers[loc] <- ter
            }



      }
      return(allTers)
}


### Selecting similar colours from
selectSimilar <- function(image,territory = NULL,threshold = "auto"){
    #--------------------------------------------------------------------------#
    # Getting seed colour fron each territories
    # maybe base r would be faster in this case ??
    # we can run it on more than territory maybe ??
    #--------------------------------------------------------------------------#
    if(is.null(territory)){
        cols <- image %>% group_by(cluster,cc) %>%
              select(value) %>% unique %>%
              pivot_wider(names_from = "cc", values_from = "value")
    } else if(length(territory ==1)){
        cols <- image %>% filter(territory == territory) %>%
                group_by(cc) %>% select(value) %>%
                unique() %>%
                pivot_wider(names_from = "cc", values_from = "value")
    } else {
        cols <- image %>% filter(territory %in% territory) %>%
                group_by(cc) %>% select(value) %>%
                median() %>%
                pivot_wider(names_from = "cc", values_from = "value")
    }
    #--------------------------------------------------------------------------#
    # First convertin image to actual image
    #--------------------------------------------------------------------------#
    img <- as.cimg(image[,c("x","y","cc","value")])

    #--------------------------------------------------------------------------#
    # Next we can loop over each colour palette
    # i.e each territory cluster
    #--------------------------------------------------------------------------#
    byTer <- vector("list", nrow(cols))
    for(i in seq_along(byTer)){
        tmp <- img %>%
               { . - imfill(dim=dim(img),val=cols[i,2:4]) } %>%
                imsplit("c") %>%
                enorm
        tmp <- !threshold(tmp,threshold)
        byTer[[i]] <- tmp
    }

    return(byTer)
}




.fineGrain <- function(territories,minCell,combine =FALSE){
    ### To be removed ? might not be worthwhile to keep this section
    #if(combine){
      #--------------------------------------------------------------------------#
      # Consider what you need to change for min cell and genes?
      #--------------------------------------------------------------------------#

      #allTags <- strsplit(territories$barcodes,"_")
      #allLength <- sapply(allTags, length)

      #territories <- rep(territories$territory, times = allLength)
      #return(data.frame("barcodes" = unlist(allTags),"territory"= territories))

    #} else {
      #--------------------------------------------------------------------------#
      # Split all barcode tags by seperator
      #--------------------------------------------------------------------------#
      #territories<- territories %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)

      allTags <- unlist(strsplit(territories$barcodes,"_"))
      if(length(allTags) <= minCell){
          return(NULL)
      } else {
        return(allTags)
      }

    #}


}



### might be dropped in final version
addTerritories <- function(dat,coordinates = NULL, global = TRUE){

    ters <- names(dat)

    dat <- lapply(seq_along(ters), function(idx,ters,dat){
                  dat[[idx]]$territory <- ters[idx]
                  return(dat[[idx]])
    },ters,dat)

    if(!is.null(coordinates)){
        dat <- mapply(function(dat,coordinates){
                      tmp <- cbind(dat,coordinates[,c("x","y")])
                      return(tmp)
        },dat,coordinates, SIMPLIFY = FALSE)
    }


    dat <- do.call("rbind",dat)
    if(global){
      dat <- .globaliseTerritories(dat, seurat = TRUE)
    }

    return(dat)
}

.getSeuratCoordinates <- function(seurat){
    return(seurat@images$slice1@coordinates)
}




layerTerritory.concave <- function(image,seedTerritory = NULL,layerDepth = NULL,concavity=1,
  length_threshold=0, dilationFactor = 3,captureRadius=0.2,minCell = 10, verbose =TRUE){

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
                                         captureRadius = captureRadius,global=TRUE,
                                         minCell = minCell,verbose =FALSE)
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
             as.cimg %>% grayscale %>% grow(dilationFactor)
      #------------------------------------------------------------------------#
      # Once we have dilated we can start iterative layayering
      # We will do this on tile centers and that is what we really care about
      # First let's convert back to a data frame and add back the barcodes and what
      #not
      #------------------------------------------------------------------------#
      .rebuildDF(verbose)

      ter <- as.data.frame(ter)
      ter <- right_join(ter,pooled[[te]],by = c("x","y"))
      ter <- inner_join(ter, image,by = c("x","y")) %>%
             select(c("barcodes.x","x","y","cluster.x","territory","subTerritory","tile.y"))
      colnames(ter) <-c("barcodes","x","y","cluster","territory","subTerritory","tile")
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
              warning("Layer depth exceeds layers in Territory - using layers in territories", immediate. = TRUE)

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


layerTerritory.edge <- function(image,seedTerritory = NULL,layerDepth = NULL,
  dilationFactor = 3,verbose =TRUE){

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
    dilationFactor <- ifelse(dilationFactor <= 0,1,dilationFactor)
    ymin <- ifelse((min(ter$y) - dilationFactor *2) <=0,1,min(ter$y) - dilationFactor *2)
    xmin <- ifelse((min(ter$x) - dilationFactor *2) <=0,1,min(ter$x) - dilationFactor *2)
    ymax <- max(ter$y) + dilationFactor * 2
    xmax <- max(ter$x) + dilationFactor * 2


    #--------------------------------------------------------------------------#
    # Dilate territory to ensure that we cover the outer layers as well
    #--------------------------------------------------------------------------#
    .dilate(verbose)
    terTmp <- ter %>% select(c("x","y","cc","value")) %>%
              rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
              as.cimg %>% grayscale() %>% grow(dilationFactor)

    terForLoop <- terTmp %>% as.data.frame()

    terForLoop <- inner_join(terForLoop,image,by = c("x","y"))%>%
           select(c("barcodes","x","y","cc.y","z","tile","cluster","territory"))%>%
           filter(cc.y ==1)
    colnames(terForLoop) <- c("barcodes","x","y","cc","value","tile","cluster","territory")
    ter <- terForLoop
    colnames(ter) <- c("barcodes","x","y","cc","value","tile","cluster","territory")


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
      ymin <- ifelse((min(terForLoop$y) - dilationFactor *2) <=0,1,min(terForLoop$y) - dilationFactor *2)
      xmin <- ifelse((min(terForLoop$x) - dilationFactor *2) <=0,1,min(terForLoop$x) - dilationFactor *2)
      ymax <- max(terForLoop$y) + dilationFactor * 2
      xmax <- max(terForLoop$x) + dilationFactor * 2

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
            warning("Layer depth exceeds layers in Territory - using layers in territories", immediate. = TRUE)

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
