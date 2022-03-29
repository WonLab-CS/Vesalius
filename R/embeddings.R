################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/Latent space embeddings/--------------------------------#

#' Create Vesalius Image embeddings from ST data.
#' @param counts count matrix (matrix or sparse matrix) with barcodes as column
#' names. Rows can represent genes (transcriptomics) or bins (DNA).
#' @param coordinates data frame containing barcodes, x coordinates and
#' y coordinates
#' @param method character describing which dimensionality reduction method
#' should be used (PCA or UMAP)
#' @param pcs numeric describing the number of Principle components to be used.
#' Default pcs = 30
#' @param tensorResolution numeric (range 0 - 1) describing the compression
#' ratio to be applied to the final image. Default = 1
#' @param filterGrid numeric (range 0 - 1) size of the grid used when filtering
#' outlier beads. Defined as a proportion of total image size. Default = 0.1
#' @param filterThreshold numeric (range 0 -1) describing the quantile threshold
#' at which tiles should be retained (seed details)
#' @param nfeatures numeric describing the number of variable features to use.
#' @param loadings logical if loading values should be used instead of
#' embeddings when coverting PCA to RGB. Default = FALSE
#' @param cores numeric number of cores to use. Default = 1
#' @param verbose logical output progress message or not. Default = TRUE
#' @details
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","slice"
#' Slice represents the PC slice represented.
#' @examples
#' \dontrun{
#' data(Vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' }

buildVesaliusEmbeddings <- function(vesalius,
                             method = c("PCA","PCA_L","UMAP","LSI","LSI_UMAP"),
                             pcs = 30,
                             tensorResolution = 1,
                             filterGrid =0.01,
                             filterThreshold = 0.995,
                             nfeatures = 2000,
                             cores = 1,
                             verbose = TRUE){
    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # get coord and counts
    # TODO update this when you have written utility functions
    #--------------------------------------------------------------------------#
    coordinates <- vesalius@tiles
    counts <- vesalius@counts

    #--------------------------------------------------------------------------#
    # Filter outlier points
    # Both 0 and 1 threshold will lead to the similar results i.e no filtering
    # This will need to be tweaked - Visium won't have this issue
    # might need to add a platform tag to auto change the filter
    #--------------------------------------------------------------------------#
    if(filterGrid != 0 & filterGrid != 1){
      .distanceBeads(verbose)
      coordinates <- .filterGrid(coordinates = coordinates,
                                 filterGrid = filterGrid)
    }
    #--------------------------------------------------------------------------#
    # If we want to reduce the number of tiles
    # reduced resolution
    # TO DO look at knn algorithm for this one it might be better
    #--------------------------------------------------------------------------#
    if(tensorResolution < 1){
      .tensorRes(verbose)
      coordinates <- .reduceTensorResolution(coordinates = coordinates,
                                             tensorResolution = tensorResolution)
      .adjCounts(verbose)
      counts <- .adjustCounts(coordinates, counts,cores)
    }
    #--------------------------------------------------------------------------#
    # TESSELATION TIME!
    #--------------------------------------------------------------------------#
    .tess(verbose)
    tesselation <- deldir::deldir(x = as.numeric(coordinates$x),
                                  y = as.numeric(coordinates$y))
    #--------------------------------------------------------------------------#
    # Filtering tiles
    #--------------------------------------------------------------------------#
    .fTiles(verbose)
    filtered <- .filterTiles(tesselation,coordinates,filterThreshold)

    #--------------------------------------------------------------------------#
    # Resterise tiles
    #--------------------------------------------------------------------------#
    .raster(verbose)
    tiles <- .rasterise(filtered, cores)
  
    vesalius@tiles <- tiles
    #--------------------------------------------------------------------------#
    # Now we can start creating colour embeddings
    # start with basic Seurat pre-processing
    # Will Need to adjust this!!! Once LSI comes into play we will need to adpat
    # this section of the code to ensure that we have the right approach
    # It might be time to move away from Seurat....
    #### TO ADD -> pre processing section
    #--------------------------------------------------------------------------#

    #counts <- .checkCounts(counts,verbose)
    .buildSO(verbose)
    counts <- CreateSeuratObject(counts, assay ="Spatial")
    counts <- NormalizeData(counts,verbose = FALSE)
    counts <- ScaleData(counts, verbose = FALSE)
    counts <- FindVariableFeatures(counts,nfeatures = nfeatures,verbose = FALSE)
    vesalius@counts <- GetAssayData(counts,slot = "data")
    #--------------------------------------------------------------------------#
    # Embeddings
    #--------------------------------------------------------------------------#
    if(!method[1L] %in% c("PCA","PCA_L","UMAP","LSI","LSI_UMAP")){
        method <- "none"
    }
    embeds <- switch(method[1L],
                     "PCA" = .embedPCA(counts,pcs = pcs,cores = cores,
                              verbose = verbose),
                     "PCA_L" = .embedPCAL(counts,pcs = pcs,cores = cores,
                                verbose = verbose),
                     "UMAP" = .embedUMAP(counts,pcs = pcs,verbose),
                     "LSI" = .embedLSI(...),
                     "LSI_UMAP" = .embedLSIUMAP(...),
                     "none" = stop("Unsupported embedding type!"))
    vesalius@embeddings <- embeds
    vesalius@activeEmbeddings <- embeds

    #--------------------------------------------------------------------------#
    # create run log
    #--------------------------------------------------------------------------#
    newLog <- as.list(match.call())
    vesalius <- .commitLog(vesalius = vesalius,
                      commit = newLog,
                      defaults = as.list(args(buildVesaliusEmbeddings)))
    .simpleBar(verbose)
    return(vesalius)
}




.filterGrid <- function(coordinates,filterGrid){
  #----------------------------------------------------------------------------#
  # Essentially create a grid where each barcode is pooled into a grid space
  # If there are too little barcodes in that grid section then remove
  #----------------------------------------------------------------------------#
  gridX <- round(coordinates$x * filterGrid)
  gridY <- round(coordinates$y * filterGrid)
  gridCoord <- paste0(gridX,"_",gridY)
  grid <- table(gridCoord)
  grid <- grid[which(grid <= quantile(grid,0.01))]
  gridCoord <- which(gridCoord %in% names(grid))
  coordinates <- coordinates[-gridCoord,]
  return(coordinates)
}


### might want to adjust this and use knn instead?
### maybe that would be better -> aggregate points together so it's "fair"
.reduceTensorResolution <- function(coordinates,tensorResolution = 1){
  #----------------------------------------------------------------------------#
  # we will reduce number of points this way
  # this should keep all barcodes - with overlapping coordinates
  #----------------------------------------------------------------------------#
  coordinates$x <- round(coordinates$x * tensorResolution)
  coordinates$y <- round(coordinates$y * tensorResolution)
  #----------------------------------------------------------------------------#
  # Now we get coordinate tags - we use this to find all the merge locations
  # sorting and using rle to ensure that we actually merge them
  #----------------------------------------------------------------------------#
  tag <- paste0(coordinates$x,"_",coordinates$y)

  locs <- rle(sort(tag))
  dup <- locs$values[which(locs$length >1)]
  #----------------------------------------------------------------------------#
  # creating new merged labels
  #----------------------------------------------------------------------------#
  dupTags <- lapply(dup,function(dup,tag,barcodes){
    tmpLocs <- which(tag == dup)
    barcodes<- paste0(barcodes[tmpLocs],sep ="_et_", collapse ="")
    barcodes <- rep(barcodes,times = length(tmpLocs))
    return(barcodes)
  },tag = tag,barcodes = coordinates$barcodes)
  #----------------------------------------------------------------------------#
  # assigning new merged barcodes
  # for some reason match doesn't seem to work here - it only returns one of the
  # values
  #----------------------------------------------------------------------------#
  locs <- unlist(lapply(dup, function(dup,tag){
      return(which(tag == dup))
  },tag = tag))
  coordinates$barcodes[locs] <- unlist(dupTags)
  coordinates <- coordinates %>% distinct(barcodes,.keep_all = TRUE)
  return(coordinates)
}


.adjustCounts <- function(coordinates, counts,cores = 1){
    #--------------------------------------------------------------------------#
    # First get all barcode names and compare which ones are missing
    #--------------------------------------------------------------------------#
    coordBar <- unique(coordinates$barcodes)
    countBar <- unique(colnames(counts))
    concat <- sapply(strsplit(coordBar,"_et_"),length) > 1
    coordBar <- coordBar[concat]
    countBar <- countBar[!concat]

    #--------------------------------------------------------------------------#
    # next we merge counts together when barcodes have been merged
    #--------------------------------------------------------------------------#
    tmpBar <- strsplit(coordBar,"_et_")

    empty <- parallel::mclapply(tmpBar, function(coord,count){
        tmp <- rowSums(count[,coord])
        return(tmp)
    },count = counts, mc.cores = cores)

    empty <- do.call("cbind",empty)

    colnames(empty) <- coordBar

    merged <- cbind(counts[,countBar],empty)
    return(merged)
}

.filterTiles <- function(tesselation,coordinates,filterThreshold){
  maxArea <- quantile(tesselation$summary$dir.area, filterThreshold)
  idx <- which(tesselation$summary$dir.area >= maxArea)
  tessV <- tesselation$dirsgs
  pointsToRemove <- tessV$ind1 %in% idx | tessV$ind2 %in% idx
  tessV <- tessV[!pointsToRemove,]
  coordinates$ind <- seq_len(nrow(coordinates))
  coordinates <- coordinates[-idx,]
  return(list("tessV" = tessV,"coordinates" = coordinates))
}

.rasterise <- function(filtered,cores = 1){
    idx <- seq_len(nrow(filtered$coordinates))
    tiles <- parallel::mclapply(idx, function(idx,filtered){

        #----------------------------------------------------------------------#
        # get indecies from original data
        #----------------------------------------------------------------------#
        ind <- filtered$coordinates$ind[idx]
        indx <- filtered$coordinates$x[idx]
        indy <- filtered$coordinates$y[idx]
        tessV <- filtered$tessV %>% filter(ind1 == ind | ind2 == ind)
        if(nrow(tessV) == 0){
            return(NULL)
        }
        #----------------------------------------------------------------------#
        # create unique set of coordiantes that define tile boundaries
        #----------------------------------------------------------------------#
        coord <- paste0(c(tessV$x1,tessV$x2),"_",c(tessV$y1,tessV$y2))
        x <- as.numeric(sapply(strsplit(coord[!duplicated(coord)],"_"),"[[",1))
        y <- as.numeric(sapply(strsplit(coord[!duplicated(coord)],"_"),"[[",2))

        convex <- .convexify(x,y,indx,indy)
        x <- convex$x
        y <- convex$y
        #----------------------------------------------------------------------#
        # define max polygon containing all pixels
        #----------------------------------------------------------------------#
        lpx <- round(min(x)) - 1
        hpx <- round(max(x)) + 1
        lpy <- round(min(y)) - 1
        hpy <- round(max(y)) + 1

        maxPolygonX <- rep(seq(lpx,hpx), times = hpy - lpy +1)
        maxPolygonY <- rep(seq(lpy,hpy), each = hpx - lpx +1)

        #----------------------------------------------------------------------#
        # Fill triangles with all point in that space
        #----------------------------------------------------------------------#
        cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
        maxX <- maxPolygonX[cell %in% c(1,2,3)]
        maxY <- maxPolygonY[cell %in% c(1,2,3)]
        tile <- data.frame("barcodes" = rep(filtered$coordinates$barcodes[idx],
                                            times = length(maxX)),
                           "x" = maxX,
                           "y" = maxY)
        return(tile)
    },filtered = filtered, mc.cores = cores)
    tiles <- do.call("rbind",tiles)

    tiles <- tiles %>% filter(x > 1 & y > 1)
    return(tiles)
}

.convexify <- function(xside,yside,indx,indy){
  #----------------------------------------------------------------------------#
  # Converting everything to an angle - from there we can just go clock wise
  # and order the point based on angle
  #----------------------------------------------------------------------------#
  x <- xside - indx ; y <- yside - indy
  angle <- mapply(function(x,y){
    if(x >= 0 & y >= 0) angle <- atan(abs(y)/abs(x))*(180/pi)
    if(x < 0 & y >= 0) angle <- 180 - (atan(abs(y)/abs(x))*(180/pi))
    if(x < 0 & y < 0) angle <- 180 + (atan(abs(y)/abs(x))*(180/pi))
    if(x >= 0 & y < 0) angle <- 360 - (atan(abs(y)/abs(x))*(180/pi))
      return(angle)
  },x=x,y=y, SIMPLIFY = TRUE)
  convex <- data.frame(x=xside[order(angle,decreasing = F)],
                       y=yside[order(angle,decreasing = F)])
  return(convex)
}

.normPix <- function(embeds,type = c("minmax","quantileNorm")){
    #--------------------------------------------------------------------------#
    # Normalise pixels values

    #--------------------------------------------------------------------------#
    embeds <- switch(type[1L],
                     "mixmax" = .minMax(embeds),
                     "quantileNorm" = midFix(embeds))
}


### might move this over to misc
.minMax <- function(embeds){
    embeds <- apply(embeds,2,function(embeds){
        return((embeds - min(embeds)) / (max(embeds) - min(embeds)))
    })
    return(embeds)
}
.quantileNorm <- function(embeds){

    embeds_rank <- apply(embeds,2,rank,ties.method="min")
    embeds <- data.frame(apply(embeds, 2, sort))
    embeds_mean <- apply(embeds, 1, mean)

    index_to_mean <- function(my_index, my_mean){
      return(my_mean[my_index])
    }

    embeds_final <- apply(embeds_rank, 2, function(idx,m){
        return(m[idx])
    },embeds_mean)
    rownames(embeds_final) <- rownames(embeds)
    return(embeds_final)


}

.embedPCA <- function(counts,pcs,loadings = TRUE,cores = 1,verbose = TRUE){
    #--------------------------------------------------------------------------#
    # First run PCA
    #--------------------------------------------------------------------------#
    .pcaTensor(verbose)
    counts <- RunPCA(counts, npcs = pcs, verbose = FALSE)

    .embedRGBTensor(verbose)
    #--------------------------------------------------------------------------#
    # Here we can just sum and normalise
    # this is going to be much faster
    # transpose at the end so we keep common format
    #--------------------------------------------------------------------------#
    pca <- Embeddings(counts, reduction = "pca")
    #pca <- apply(pca,2,function(x)return(abs(x)))
    pca <- apply(pca,2,function(x){
      x <- (x - min(x)) / (max(x) - min(x))
      return(x)
    })
    colourMatrix <- list(as.matrix(pca))
    names(colourMatrix) <- "PCA"
    return(colourMatrix)
}

.embedPCAL <- function(counts,pcs,cores = 1,verbose = TRUE){
    #--------------------------------------------------------------------------#
    # First run PCA
    #--------------------------------------------------------------------------#
    .pcaTensor(verbose)
    counts <- RunPCA(counts, npcs = pcs, verbose = FALSE)

    #--------------------------------------------------------------------------#
    # get laodings and create matrix describing if there are any count values
    #--------------------------------------------------------------------------#
    pca <- Loadings(counts, reduction = "pca")
    pca <- apply(pca,2,function(x)return(abs(x)))
    mat <- as.matrix(GetAssayData(counts, slot = "data") > 0)
    colourMatrix <- matrix(0,nrow = ncol(mat),ncol = ncol(pca))
    colnames(colourMatrix) <- colnames(pca)
    rownames(colourMatrix) <- colnames(mat)
    #--------------------------------------------------------------------------#
    # Looping over each PC to get laodings sum and normalising
    #--------------------------------------------------------------------------#
    for(p in seq_len(ncol(pca))){
      .pcaRGBTensor(verbose,p)
      bars <- parallel::mclapply(seq_len(ncol(mat)), function(idx,mat,pca){
        cVec <- as.numeric(mat[names(pca),idx])
        colour <- sum(pca * cVec)
      },mat = mat, pca = pca[,p],mc.cores = cores)
      bars <- unlist(bars)

      colourMatrix[,p] <- (bars - min(bars)) / (max(bars) - min(bars))
    }
    colourMatrix <- list(as.matrix(pca))
    names(colourMatrix) <- "PCA_L"
    return(colourMatrix)
}



.embedUMAP <- function(counts,pcs,verbose){


  #----------------------------------------------------------------------------#
  # First run PCA and then UMAP
  #----------------------------------------------------------------------------#
  .pcaTensor(verbose)
  counts <- RunPCA(counts, npcs = pcs, verbose = FALSE)
  .umapRGBTensor(verbose)
  counts <- RunUMAP(counts, dims = seq_len(pcs),n.components = 3L,verbose = F)
  #----------------------------------------------------------------------------#
  # Normalise
  #----------------------------------------------------------------------------#
  counts <- FetchData(counts, c("UMAP_1","UMAP_2","UMAP_3"))
  counts <- apply(counts,2,function(x){
                     return((x-min(x))/(max(x) - min(x)))
                   })

  counts <- list(as.matrix(counts))
  names(counts) <- "UMAP"
  return(counts)
}
