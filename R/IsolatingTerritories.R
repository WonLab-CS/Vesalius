################################################################################
############################   ST package        ###############################
################################################################################

#----------------------/Isolating Territories/---------------------------------#


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


enhance <- function(img, saturation = 0.3){
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


findSubClusters <- function(img,SO, by = c("cluster","territory"),varF = 2000,ncps = 20, resolution = 0.5){
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
            tmpSO <- RunPCA(tmpSO,ncps = ncps)
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



##### The oppostite now
### manuallycollapsing

### this is going to be a bit odd in the approache
## parsing list ?

# cola <- list(list(c(1,2),c(2,3),c(1,5)),list(c(1,2),c(2,3),c(1,5)))
# keep in simple for now

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
