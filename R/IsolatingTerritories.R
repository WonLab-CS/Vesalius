################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#

#' Applying morphological filter to beads by smoothing bead to dominant colour of nearest neighbors
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @param nn nearest neighbors list
#' @param k numeric describing number of neighbors to use for smoothing
#' @param threshold numeric between 0 and 1 describing which quantile should be used to select dominant colour.
#' @details Applying morphological filters on non image data. Morphological filters converts pixels values
#' to the average value of neighboring pixels. The same process is applied here however instead of using
#' pixel, this function uses closest beads/spots in space (according to x/coordinates). Depending on how nearest
#' neighbors were selected, nearest neighbors should be beads/spots that surround the center bead in all directions.
#' @return data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns but with adjusted RGB values.
smoothToDominant <- function(img,nn, k =10, threshold = 0.5){
    imgCopy <- img
    for(i in seq_along(nn)){
        cat(paste((i/length(nn))*100),"%","\r")
        tmp <- img[img$barcodes %in% nn[[i]]$barcodes,c("R","G","B") ]

        tmp <- apply(tmp,2,quantile, threshold, na.rm =TRUE)
        imgCopy[imgCopy$barcodes == names(nn)[i],c("R","G","B")] <- tmp
    }
    imgCopy<-imgCopy[!is.na(imgCopy$R) &
                     !is.na(imgCopy$G) &
                     !is.na(imgCopy$B) ,]
    return(imgCopy)
}

#' Enhance colour contrast by equalizing the colour histogram
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @details Histogram equalisation ensure that colours are uniformly distributed between all possible colours
#' Conceptually, this is similar to flatening a curve.
#' @return data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns with adjusted RGB values.
enhance <- function(img){
    img[,c("R","G","B")] <- apply(img[,c("R","G","B")],2,
                            equalizeHist)

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

#' isolating territories from spatial transcriptomic data
#' @param img data frame with barcodes, xcoord, ycoord, R, G, and B as coloumns
#' @param colDepth integer describing the number of colours to extract
#' @param dropRadius numeric describing the proportion of total distance to consider when pooling beads together
#' @details In order to select territories, colours will be collapsed into \code{colDepth} number of
#' dominant colours. For each colour, territories will be isolated based on distance between point.
#' The \code{dropRadius} represents the proportion of total distance (x and y values may differ) to use in order
#' to pool beads/spots together. The algorithm selects a point and pools all neighboring beads into a territory.
#' The process is applied to the nearest neighbors until no more points can be pooled into the territory.
#' If beads remain for a given colour, the process is repeated until all beads/spots are put into a territory.
#' @return a data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster.

isolateTerritories <- function(img,colDepth = 8,dropRadius = 0.025){
      ## clustering colours
      clust <- kmeans(img[,c("R","G","B")],colDepth)
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

      ## Finding territories in each
      clusterCount <- 1
      img$territories <- NA
      for(i in seq_along(clusters)){
          ## all barcodes from cluster i
          barcodes <- img[img$cluster == clusters[i],]
          message(paste0("Pooling cluster ",i))
          ## pooling
          pool <- distancePooling(barcodes,)

          img[img$cluster == clusters[i],]<- pool

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
                    print(sum(overlap))
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
            tmpSO <- NormalizeData(tmpSO)
            tmpSO <- FindVariableFeatures(tmpSO, selection.method = "vst", nfeatures = varF)
            tmpSO <- ScaleData(tmpSO)
            tmpSO <- RunPCA(tmpSO,ncps = npcs)
            tmpSO <- RunUMAP(tmpSO,reduction = "pca", dims = seq(1,ncps))
            tmpSO <- FindNeighbors(tmpSO,reduction="pca", dims = seq(1,ncps))
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



#' Combining Territories
#' @param img data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster.
#' @param territories list of vectors (cluster and territory).
#' @details In order to combine territories together, this function takes a list of territories that should be
#' combined into one. The list elements are vectors containing first the cluster number and then the territory
#' number within that cluster.
#'@return  data frame with barcodes, xcoord, ycoord, R, G, B, cluster number, territory number within that cluster where
#'combined territories will take the value of the first cluster/territory provided.

combineTerritories <- function(img,territories){
    ## first get all section out of img
    clust <- sapply(territories,"[[",1L)
    ter <- sapply(territories,"[[",2L)
    # combination territories
    combTer <- paste0(clust,"_",ter)
    # index territories
    indexTer <- paste0(img$cluster,"_",img$territories)

    ## references territory
    cl1 <- clust[1L]
    ter1 <- ter[1L]
    ## get index
    subImg <- indexTer %in% combTer
    #browser()
    img[subImg,"cluster"] <- cl1
    img[subImg,"territories"] <- ter1

    return(img)

}
