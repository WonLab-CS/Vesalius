################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#

#' Applying morphological filter to beads by smoothing bead to dominant colour of nearest neighbors
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

         tmp$inc <- order(nn$distance)
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





#' Enhance colour contrast by equalizing the colour histogram
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @details Histogram equalisation ensure that colours are uniformly distributed between all possible colours
#' Conceptually, this is similar to flatening a curve.
#' @return data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns with adjusted RGB values.
enhance <- function(img){
    img[,c("R","G","B")] <- apply(img[,c("R","G","B")],2,
                            .equalizeHist)

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
#' @param iter numeric describing number of segmentation/smoothing iterations
#' @param startingSegements numeric describing the number of initial colour segments to consider.
#' @param eq logical indicating if colour histogram should be equalised before each iteration.
#' @param verbose logical indicating if progress messages should be outputed to console.
#' @details For each iteration, the algorithm collapse the colour space into n colour segments.
#' Then, smoothing is applied to nearest neighbours. The process is repeated for each iteration.
#' Note that the number of segments from starting segments to final colour depth follows the following
#' Function : \code{seq(startingSegements,colDepth+1,length.out = iter)}
#' @return a data frame containing RGB code for each bead.

iterativeSegmentation <- function(img,nn,colDepth = 9,segIter = 10,
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

isolateTerritories <- function(img,dropRadius = 0.025,colDepth =12){
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

.distancePooling <- function(img,dropRadius= 0.025){
        dropDistance <- sqrt((max(img$xcoord)-min(img$xcoord))^2 +
                         (max(img$ycoord)-min(img$ycoord))^2) * dropRadius


        idx <- seq_len(nrow(img))
        distanceMatrix <- lapply(idx, function(idx,mat){
                            xo <- mat$xcoord[idx]
                            yo <- mat$ycoord[idx]
                            xp <- mat$xcoord
                            yp <- mat$ycoord
                            distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)
                            return(distance)
        }, img)

        distanceMatrix <- do.call("rbind",distanceMatrix)
        colnames(distanceMatrix) <- img$barcodes
        rownames(distanceMatrix) <- img$barcodes

        barcodes <- img$barcodes
        territories <- vector("list", length(barcodes))
        count <- 1

        while(length(barcodes) >0){
            #print(length(barcodes)/length(territories)*100)
            tmp <- distanceMatrix[,barcodes[1]]
            pool <- rownames(distanceMatrix)[tmp < dropDistance]
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

                    newPool <- unique(unlist(lapply(seq_len(ncol(newPool)), function(idx,np,dropDistance){
                                            res <- rownames(np)[np[,idx] < dropDistance]
                                            return(res)
                                          },newPool,dropDistance)))
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



        nulls <- sapply(territories, is.null)
        territories <- territories[!nulls]

        for(ter in seq_along(territories)){
            loc <- match(territories[[ter]],img$barcodes)
            img[ loc[!is.na(loc)],"territories"] <- ter
        }
        return(img)

}

#' findSubCluster - Clustering of territories
#' @param img data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster.
#' @param SO Seurat Object with full spatial assay
#' @param by character describing if subclustering should be caried out on colour cluster or territories.
#' @param varF integer describing number of variable features
#' @param npcs integer describing the number of Principle components to compute
#' @param resolution numeric describing granularity of clustering
#' @details The following function is a wrapper function that performs the standard Seurat analysis.
#' @return A list of lists. The first level represents each cluster of territory. For each cluster or territory,
#' the list element contains a Seurat Object with sub cluster and markers associated to each subcluster.
findSubClusters <- function(img,SO, by = c("cluster","territory"),varF = 2000,npcs = 20, resolution = 0.5){
    if(by[1L] == "cluster"){
        ## subsetting by cluster
        clusters <- unique(img$cluster)
        subCluster <- vector("list", length(clusters))
        names(subCluster) <- as.character(clusters)
        for(i in seq_along(clusters)){
            barcodes <- img[img$cluster == clusters[i],"barcodes"]
            tmpSO <- subset(SO, cells = barcodes)
            tmpSO <- SCTransform(tmpSO, assay = "Spatial")
            tmpSO <- FindVariableFeatures(tmpSO, selection.method = "vst", nfeatures = varF)
            tmpSO <- ScaleData(tmpSO)
            tmpSO <- RunPCA(tmpSO,ncps = npcs)
            tmpSO <- RunUMAP(tmpSO,reduction = "pca", dims = seq(1,npcs))
            tmpSO <- FindNeighbors(tmpSO,reduction="pca", dims = seq(1,npcs))
            tmpSO <- FindClusters(tmpSO,resolution = resolution)
            markers <-  FindAllMarkers(tmpSO)
            subCluster[[i]] <- list("subClusters"= tmpSO,"markers" = markers)

        }
    } else {
        message("coming Soon")
        subCluster <- NULL
    }

    return(subCluster)
}



#' Extracting Territories from image data
#' @param img data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster.
#' @param territories list of vectors (cluster and territory).
#' @details In order to combine territories together, this function takes a list of territories that should be
#' combined into one. The list elements are vectors containing first the cluster number and then the territory
#' number within that cluster.
#'@return  data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster where
#'combined territories will take the value of the first cluster/territory provided.

extractTerritories <- function(img,territories){
    ## first get all section out of img
    clust <- sapply(territories,"[[",1L)
    ter <- sapply(territories,"[[",2L)
    # combination territories
    combTer <- paste0(clust,"_",ter)
    # index territories
    indexTer <- paste0(img$cluster,"_",img$territories)


    ## get index
    subImg <- indexTer %in% combTer

    subImg <- img$barcodes[subImg]


    return(subImg)

}
