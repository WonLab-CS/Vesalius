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
                        gaussian=TRUE,na.rm=FALSE, acrossLevels = "min", iter = 1,verbose = TRUE){
    #--------------------------------------------------------------------------#
    # Class check
    #--------------------------------------------------------------------------#

    if(!is.cimg(img)){
      imgCopy <- as.cimg(select(img, c("x", "y", "cc", "value")))
    }else {
      imgCopy <- img
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






equalizeContrast <- function(image, type = "BalanceSimplest",N=1,smax=1,
                         sleft =1,sright =1,lambda=0.1,up=100, down = 10,range = c(0,1)){
    img <-image %>% select(c("x","y","cc","value")) %>% as.cimg %>% imsplit("c")
    img <- switch(type,
                    "EqualizePiecewise" = lapply(img,EqualizePiecewise,N) %>% imappend("c"),
                    "BalanceSimplest" = lapply(img,BalanceSimplest,sleft,sright,range) %>% imappend("c"),
                    "SPE" = lapply(img,SPE,lambda) %>% imappend("c"),
                    "EqualizeDP"  = lapply(img,EqualizeDP,down,up) %>% imappend("c"),
                    "EqualizeADP" = lapply(img,EqualizeADP) %>% imappend("c"))
    img <- as.data.frame(img)
    img <- right_join(img, image, by  = c("x","y","cc")) %>%
           select(c("barcodes","x","y","cc","value.x","tile")) %>% tibble

    colnames(img) <- c("barcodes","x","y","cc","value","tile")
    return(img)
}

## Have to run them internally



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
                       iter = smoothIter,verbose=verbose)


    cat("\n")
    .seg(j,verbose)
    cat("\n")
    if(useCenter){
      tmpImg <- img %>% filter(tile == 1) %>% group_by(cc) %>% distinct(barcodes,.keep_all = TRUE)
      #tmpImg <- tmpImg %>% group_by(cc) %>%
      #          mutate(cluster = .intKmeans(value,colDepth[j],cluster = TRUE),
      #                 Kcenters = .intKmeans(value,colDepth[j],cluster = FALSE))

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
      #tmpImg <- tmpImg%>% group_by(cc,cluster) %>%
      #          mutate(value = Kcenters)

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
    #--------------------------------------------------------------------------#
    barcodesInitial <- img$barcodes
    imgCopy <- img %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)


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



### Cant remember why I took this approach
### Lets keep it ultra simple for now
### Dont create this intermediate object only if required
extractTerritories <- function(img,seurat,terIdent = NULL,combine = FALSE,
                               minCell = 10, verbose = TRUE,cores = 1){
    #--------------------------------------------------------------------------#
    # if combine is null we consider that we want to have all
    # territories prepared for clustering analysis
    #--------------------------------------------------------------------------#
    if(is.null(terIdent)){
        territories <- img %>% filter(tile == 1) %>% distinct(barcodes, .keep_all = FALSE)
        territories <- split(territories, territories$territory)
        #----------------------------------------------------------------------#
        # Revert back to fine grain if inage array was reduced
        #----------------------------------------------------------------------#
        #territories <- parallel::mclapply(territories,.fineGrain,minCell,mc.cores = cores)
        #territories <- territories[!sapply(territories,is.null)]
        cells <- sapply(territories,nrow) > minCell
        territories <- territorries[cells]
        territories <- lapply(territories,"$",territory)
        territories <- parallel::mclapply(territories,.subSetTerritories,seurat,mc.cores= cores)
        return(territories)
    } else {
        #----------------------------------------------------------------------#
        # Filtering only the barcodes that are in the territories of interest
        # To ADD check for combine - numeric values
        #----------------------------------------------------------------------#
        territories <- img %>% filter(tile==1 & territory %in% terIdent)%>%
                       distinct(barcodes, .keep_all = FALSE)
        territories <- territories$barcodes
        #----------------------------------------------------------------------#
        # Revert back to fine graind and return null if territory does not
        # Conatin enough cells - in this case throw in a warning
        #----------------------------------------------------------------------#
        #territories <- .fineGrain(territories,minCell)

        if(length(territories) < minCell){

            warning("Territory selection does not contain enough cells - NULL returned")
            return(NULL)
        }else{
            territories <- .subSetTerritories(territories,seurat)
            return(territories)
        }
    }


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

.subSetTerritories <- function(territories,seurat){
    #--------------------------------------------------------------------------#
    # Simplified version for now
    # It might be worth while getting away from seurat later
    # essentially this is a template function
    #--------------------------------------------------------------------------#

    seurat <- subset(seurat, cells = territories)
    return(seurat)
}


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




extractLayers <- function(image,territory = NULL,layerDepth = 1,concavity=1,length_threshold=0){
    layer <- list()
    if(is.null(territory)){
        # to add
        # loop over all territories
    } else {
        ter <- image %>% filter(territory %in% territory)
        edge <- concaveman::concaveman(as.matrix(ter[,c("x","y")]),
                            concavity,length_threshold)
        edge <- paste0(edge[,1],edge[,2],sep = "_")
        points <- paste0(ter$x,ter$y,sep = "_")
        edge <- points %in% edge
        #layer

    }
}
