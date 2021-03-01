################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#

#' Applying smoothing filter to beads by smoothing bead to dominant colour of nearest neighbors
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @param nn nearest neighbors list
#' @param k numeric describing number of neighbors to use for smoothing
#' @param iter numeric describing the number of smoothing iterations that should be performed
#' @param angleBias logical describing if angle of nearest neighbor should be account for.
#' @param threshold numeric between 0 and 1 describing which quantile should be used to select dominant colour.
#' @param logical describing if progress message should outputed to console.
#' @details Applying morphological filters on non image data. Morphological filters converts pixels values
#' to the average value of neighboring pixels. The same process is applied here however instead of using
#' pixel, this function uses closest beads/spots in space (according to x/coordinates). Depending on how nearest
#' neighbors were selected, nearest neighbors should be beads/spots that surround the center bead in all directions.
#' @return data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns but with adjusted RGB values.
smoothToDominant <- function(img,nn,method = c("median","exp","linear"), iter = 1,angleBias = FALSE, verbose =TRUE){
    ## Building smoothing kernel


    ### For each smoothing iteration
    ### Note that this only applies to this function and is not the same
    ### As segmentation iteration even if both are used in the segmentation
    ### process
    for(it in seq_len(iter)){
      ###########################
      if(verbose){
        cat(paste("Smoothing Iteration",it, "\r"))
      }
      ###########################

      for(i in seq_along(nn)){
          #tmp <- img[img$barcodes %in% nn[[i]]$barcodes,c("R","G","B") ]
          ##### Applying angle bias
          ### Bug fix -- this is happening somewhere else
          ### but for now
          if(nrow(nn[[i]])==0) next()

          img <- switch(method[1L],
                        "median" = .medianKernel(img, nn[i], angleBias),
                        "exp" = .expKernel(img, nn[i], angleBias),
                        "linear" = .linearKernel(img, nn[i], angleBias))


      }

      img<-img[!is.na(img$R) &
               !is.na(img$G) &
               !is.na(img$B) ,]

    }

    return(img)
}

.medianKernel <- function(img,nn, angleBias=FALSE){
    if(angleBias){
      ### setting up quadrants
      q1 <- c(0,90)
      q2 <- c(90,180)
      q3 <- c(180,270)
      q4 <- c(270,360)
      if(any(colnames(img) == c("cluster"))){
          ### TO ADD
      } else {
          stop("angle Bias smoothing only possible when using iterative segmentation")
      }

    } else {
         loc <- names(nn)
         nn <- nn[[1L]]
         tmp <- img[img$barcodes %in% nn$barcodes,c("R","G","B")]
         tmp <- apply(tmp,2,median)
         img[img$barcodes == names(nn),c("R","G","B")] <- tmp
    }

    return(img)
}

.expKernel <- function(img,nn, angleBias=FALSE){
    if(angleBias){
      ### setting up quadrants
      q1 <- c(0,90)
      q2 <- c(90,180)
      q3 <- c(180,270)
      q4 <- c(270,360)
      if(any(colnames(img) == c("cluster"))){
          ### TO ADD
      } else {
          stop("angle Bias smoothing only possible when using iterative segmentation")
      }

    } else {

         loc <- names(nn)
         nn <- nn[[1L]]

         tmp <- img[match(nn$barcodes,img$barcodes),c("R","G","B")]
         maxDist <- max(nn$distance)
         minDist <- min(nn$distance)
         ### Not sure if this is the best number of increments
         ### But It can get out of hand if you stack them

         increments <- seq(maxDist, minDist, length.out = 6)
         ## First step adding increment values
         tmp$inc <- rep(0,nrow(tmp))

         for(inc in seq(1,length.out = length(increments)-1)){
            tmp[nn$distance < increments[inc] &
                nn$distance >= increments[inc+1],"inc"] <- ceiling(exp(inc))
         }
         ## Replicating points
         r <- c() ; g <- c() ; b <- c()
         for(i in seq_len(nrow(tmp))){
            r <- c(r,rep(tmp$R[i],tmp$inc[i]))
            g <- c(g,rep(tmp$G[i],tmp$inc[i]))
            b <- c(b,rep(tmp$B[i],tmp$inc[i]))
         }
         tmp <- cbind(r,g,b)

         tmp <- apply(tmp,2,median)
         img[img$barcodes == loc,c("R","G","B")] <- tmp
    }

    return(img)
}

.linearKernel <- function(img,nn, angleBias=FALSE){
    if(angleBias){
      ### setting up quadrants
      q1 <- c(0,90)
      q2 <- c(90,180)
      q3 <- c(180,270)
      q4 <- c(270,360)
      if(any(colnames(img) == c("cluster"))){
          ### TO ADD
      } else {
          stop("angle Bias smoothing only possible when using iterative segmentation")
      }

    } else {
         loc <- names(nn)
         nn <- nn[[1L]]
         tmp <- img[match(nn$barcodes,img$barcodes),c("R","G","B")]
         ### Not sure if this is the best number of increments
         ### But It can get out of hand if you stack them

         tmp$inc <- order(nn$distance, decreasing =TRUE)
         ## Replicating points
         r <- c() ; g <- c() ; b <- c()
         for(i in seq_len(nrow(tmp))){
            r <- c(r,rep(tmp$R[i],tmp$inc[i]))
            g <- c(g,rep(tmp$G[i],tmp$inc[i]))
            b <- c(b,rep(tmp$B[i],tmp$inc[i]))
         }
         tmp <- cbind(r,g,b)
         tmp <- apply(tmp,2,median)
         img[img$barcodes == loc,c("R","G","B")] <- tmp
    }

    return(img)
}

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
                        gaussian=TRUE,na.rm=FALSE, iter = 1){
    #--------------------------------------------------------------------------#
    # Class check
    #--------------------------------------------------------------------------#

    if(!is.cimg(img)){
      impCopy <- as.cimg(select(img, c("x", "y", "cc", "value")))
    }else {
      impCopy <- img
    }

    #--------------------------------------------------------------------------#
    # NOTE Export cimg object as data frame and vice versa
    # It might be easier be easier
    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Running multiple smoothing itteration
    #--------------------------------------------------------------------------#
    for(i in seq_len(iter)){
        for(j in method){
          print(j)
          impCopy <- switch(j,
                            "median" = medianblur(impCopy,box, threshold),
                            "iso" = isoblur(impCopy,sigma,neuman,gaussian,na.rm),
                            "box" = boxblur(impCopy,box,neuman))
        }

    }

    imgCopy <- as.data.frame(impCopy)



    img <- right_join(imgCopy, img, by  = c("x","y","cc")) %>%
           select(c("barcodes","x","y","cc","value.x","tile")) %>% tibble

    colnames(img) <- c("barcodes","x","y","cc","value","tile")
    return(img)
}




#' Enhance colour contrast by equalizing the colour histogram
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @details Histogram equalisation ensure that colours are uniformly distributed between all possible colours
#' Conceptually, this is similar to flatening a curve.
#' @return data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns with adjusted RGB values.
enhance <- function(img){
    ## checking classes
    if(is(img) == "data.frame"){
      img[,c("R","G","B")] <- apply(img[,c("R","G","B")],2,
                              .equalizeHist)
    } else if( is(img) == "array"){
        for(i in seq(1,3)){
            tmp <- img[,,,i]
            tmp <- .equalizeHist(as.vector(tmp))
            img[,,,i] <- tmp
        }
    } else {
      stop("Unkown image type")
    }


    return(img)
}

.equalizeHist <- function(img,range = c(0,1),levels = 256) {
      breaks <- seq(range[1L], range[2L], length.out = levels + 1L)
      h <- hist.default(img, breaks = breaks, plot = FALSE)
      cdf <- cumsum(h$counts)
      cdf_min <- min(cdf[cdf>0])

      equalized <- ( (cdf - cdf_min) / (length(img) - cdf_min) * (range[2L] - range[1L]) ) + range[1L]
      bins <- round( (img - range[1L]) / (range[2L] - range[1L]) * (levels-1L) ) + 1L

      res <- equalized[bins]
      return(res)
}

#' iterativeSegmentation - collapse and smooth points into colour startingSegements
#' @param img data frame containing RGB code for each bead
#' @param nn list containing nearest spatial neighbours
#' @param colDepth numeric describing final number of colour segments in image
#' @param segIter numeric describing number of segmentation iterations
#' @param smoothIter numeric describing number of smoothing iterations
#' @param method character describing smoothing method to use (median ,exp or linear)
#' @param startingSegements numeric describing the number of initial colour segments to consider.
#' @param eq logical indicating if colour histogram should be equalised before each iteration.
#' @param angleBias logical describing if angle to nearest neighbor should be considered (NOT IMPLEMENTED)
#' @param verbose logical indicating if progress messages should be outputed to console.
#' @details For each iteration, the algorithm collapse the colour space into n colour segments.
#' Then, smoothing is applied to nearest neighbours. The process is repeated for each iteration.
#' Note that the number of segments from starting segments to final colour depth follows the following
#' Function : \code{seq(startingSegements,colDepth+1,length.out = iter)}
#' @return a data frame containing RGB code for each bead.

iterativeSegmentation.bead <- function(img,nn,colDepth = 9,segIter = 10,
                                  smoothIter = 1,method = c("median","exp","linear"),
                                  startingSegements = 256,
                                  eq = TRUE,angleBias=FALSE, verbose = TRUE){
  ## Check
  if(segIter > (startingSegements - colDepth)){
      warning("Iterations surpass number of startingSegements \n
              Consider increasing startingSegements or reducing iterations")
      segIter <- startingSegements - colDepth
  }
  ## Building colour segments
  segments <- c(ceiling(seq(startingSegements,colDepth+1,length.out = segIter)),colDepth)

  for(j in segments){
      ## clustering colours
      if(verbose){
          cat(paste("Current Number of segments: ", j,"\r"))
      }

      clust <- kmeans(img[,c("R","G","B")],j,iter.max = 200,nstart = 10)
      img$cluster <- clust$cluster
      clusters <- unique(clust$cluster)
      ## collapsing coulours for output

      for(i in seq_along(clusters)){
          tmp <- img[img$cluster == clusters[i],]
          tmp$R <- median(tmp$R)
          tmp$G <- median(tmp$G)
          tmp$B <- median(tmp$B)
          img[img$cluster == clusters[i],] <- tmp
      }

      ### Once we have go through the last segmentation we don't want any more smoothing
      ### Will slightly change the colour values
      if(j != segments[length(segments)]){
          img <- smoothToDominant(img,nn,method =method[1L],iter=smoothIter,angleBias = angleBias,verbose = FALSE)
      }


      if(eq){
        img <- enhance(img)
      }
  }
  return(img)
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
                                        sigma = 1, box = 20, threshold=0, neuman=TRUE,
                                        gaussian=TRUE, useCenter=TRUE, na.rm=FALSE,
                                        verbose = TRUE){

  #----------------------------------------------------------------------------#
  ## Building colour segments
  #----------------------------------------------------------------------------#


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

      if(verbose){
          #cat(paste("Current Number of segments: ", j,"\r"))
      }

    #--------------------------------------------------------------------------#
    # Segmentation is done by using kmeans clustering
    # We will only run clustering on 1 pixel of each tile
    # This will ensure faster run time and cleaner results
    #--------------------------------------------------------------------------#
    img <- smoothArray(img,method =method,sigma = sigma, box = box,
                       threshold=threshold, neuman=neuman, gaussian=gaussian,
                       na.rm=FALSE,iter = smoothIter)




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

isolateTerritories.bead <- function(img,dropRadius = 0.025,colDepth =12){
      if(any(colnames(img)=="cluster")){
        clusters <- unique(as.numeric(as.character(img$cluster)))
        ## Finding territories in each
        clusterCount <- 1
        img$territories <- NA
        for(i in seq_along(clusters)){
            ## all barcodes from cluster i
            barcodes <- img[img$cluster == clusters[i],]
            message(paste0("Pooling cluster ",i))
            ## pooling
            pool <- .distancePooling(barcodes,)

            img[img$cluster == clusters[i],]<- pool

        }

      } else {
        clust <- kmeans(img[,c("R","G","B")],colDepth,iter.max = 200,nstart = 10)
        img$cluster <- clust$cluster
        clusters <- unique(as.numeric(as.character(img$cluster)))
        ## Finding territories in each
        clusterCount <- 1
        img$territories <- NA
        for(i in seq_along(clusters)){
            ## all barcodes from cluster i
            barcodes <- img[img$cluster == clusters[i],]
            message(paste0("Pooling cluster ",i))
            ## pooling
            pool <- .distancePooling(barcodes,)

            img[img$cluster == clusters[i],]<- pool

        }

      }
      return(img)
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

isolateTerritories.array <- function(img, method = c("distance","neighbor","watershed"),
                                     captureRadius = 0.025,local=TRUE){

    #--------------------------------------------------------------------------#
    # Compute real capture Radius
    #--------------------------------------------------------------------------#
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
      #------------------------------------------------------------------------#
      # We don't need all data in this case only clusters and locations
      # We can just rebuild everything afterwards
      # At least we don't don't need to compute anything unnecessarily
      #------------------------------------------------------------------------#
      tmp <- filter(img,cluster == clust[i] & cc == 1)

      ter <- switch(method[1L],
                    "distance" = .distancePooling.array(tmp,captureRadius),
                    "neighbor" = .neighborPooling.array(tmp,captureRadius),
                    "watershed" = .watershedPooling.array(tmpImg))

      for(j in seq(1,3)){
          img$territory[img$cluster == clust[i] & img$cc == j] <- ter
      }
    }

    if(!local){
        img <- .globaliseTerritories(img)
    }

    return(img)
}

.globaliseTerritories <- function(img,seurat=FALSE){
    if(!seurat){
      ter <- paste0(img$cluster,"_", img$territory)
      allTer <- unique(ter)
      ter <- seq_along(allTer)[match(ter,allTer)]
      img$territory <- ter
      return(img)
    } else {
      ter <- paste0(img$seurat_clusters,"_", img$territory)
      allTer <- unique(ter)
      ter <- seq_along(allTer)[match(ter,allTer)]
      img$seurat_clusters <- ter
      return(img)
    }


}


.distancePooling.array <- function(img,captureRadius){
    #--------------------------------------------------------------------------#
    # Select center point of each tile for only one channel
    # Dont need to run it for all channels
    #--------------------------------------------------------------------------#
    barcodesInitial <- img$barcodes
    imgCopy <- img %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)

    cat("Pooling \n")
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
    # Compute real capture Radius
    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Buildan actual matrix
    #--------------------------------------------------------------------------#
    distanceMatrix <- do.call("rbind",distanceMatrix)
    colnames(distanceMatrix) <- imgCopy$barcodes
    rownames(distanceMatrix) <- imgCopy$barcodes



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

            allTers[loc] <- ter
      }




      return(allTers)
}








extractTerritories <- function(img,seurat, combine = FALSE,keepCoarse = FALSE,
                               minCell = 10,cores = 1, verbose = TRUE){

    #--------------------------------------------------------------------------#
    # for combine we will just add a meta data columns
    # Parsing Seurat object instead
    # If needed can change it to simple tables later
    #--------------------------------------------------------------------------#
    if(combine){

      seurat@meta.data$territory <- 0
      territories <- img %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)
      territories <- .fineGrain(img, combine = TRUE)
      loc <- match(territories$barcodes,rownames(seurat@meta.data))
      seurat@meta.data$territory[loc[!is.na(loc)]]<- territories$territory[!is.na(loc)]
      return(seurat)

    } else {
      #--------------------------------------------------------------------------#
      # Splitting between territories
      # Will need to add a tag to see if territories are pooled locally or globally
      # If local , then split based on cluster number
      # with pre clean step
      #--------------------------------------------------------------------------#
      territories <- img %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)
      territories <- split(territories, territories$territory)


      #--------------------------------------------------------------------------#
      # Essentially, do we want to keep the lowered resolution or not
      # High resolution goes back to the original beads locations and data
      #--------------------------------------------------------------------------#
      if(keepCoarse){
              ##
      } else {
          territories <- parallel::mclapply(territories,.fineGrain,minCell,mc.cores = cores)

          territories <- territories[!sapply(territories,is.null)]
          territories <- parallel::mclapply(territories,.subSetTerritories,seurat,mc.cores= cores)
      }
      return(territories)
    }
}

.fineGrain <- function(territories,minCell,combine =FALSE){
    if(combine){
      #--------------------------------------------------------------------------#
      # Consider what you need to change for min cell and genes?
      #--------------------------------------------------------------------------#

      allTags <- strsplit(territories$barcodes,"_")
      allLength <- sapply(allTags, length)

      territories <- rep(territories$territory, times = allLength)
      return(data.frame("barcodes" = unlist(allTags),"territory"= territories))

    } else {
      #--------------------------------------------------------------------------#
      # Split all barcode tags by seperator
      # NOTE this will need to be changed once cluster vs territory is implemented
      #--------------------------------------------------------------------------#
      #territories<- territories %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)

      allTags <- unlist(strsplit(territories$barcodes,"_"))
      if(length(allTags) <= minCell){
          return(NULL)
      } else {
        return(allTags)
      }

    }


}

.subSetTerritories <- function(territories,seurat){
    #--------------------------------------------------------------------------#
    # Simplified version for now
    # It might be worth while getting away from seurat later
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
