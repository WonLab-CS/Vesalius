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
#' @param min.cutoff only used when dimensionality reduction method is LSI or LSI_UMAP
#' cutoff for feature to be included in the VariableFeatures for the object.
#' for more detail please look "https://satijalab.org/signac/reference/findtopfeatures"
#' @param loadings logical if loading values should be used instead of
#' embeddings when coverting PCA to RGB. Default = FALSE
#' @param cores numeric number of cores to use. Default = 1
#' @param remove_LSI1 logical only used when dimensionality reduction method is LSI or LSI_UMAP
#' indicating if the first LSI component should be removed from further analysis
#' as it usually captures sequencing depth (technical variation)
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
                             norm = c("log","SCT","TFIDF","raw"),
                             pcs = 30,
                             tensorResolution = 1,
                             filterGrid =0.01,
                             filterThreshold = 0.995,
                             nfeatures = 2000,
                             min.cutoff = "q5",
                             cores = 1,
                             remove_LSI1 = TRUE,
                             verbose = TRUE){
    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # Check Status of object
    # Leaving space for changes and updates
    #--------------------------------------------------------------------------#
    status <- .checkVesalius(vesalius,init = TRUE)
    if(status){
        #----------------------------------------------------------------------#
        # Get raw counts
        #----------------------------------------------------------------------#
        coordinates <- vesalius@tiles
        counts <- getCounts(vesalius,"raw")
        #----------------------------------------------------------------------#
        # Filter outlier beads
        #----------------------------------------------------------------------#
        if(filterGrid != 0 & filterGrid != 1){
          .distanceBeads(verbose)
          coordinates <- .filterGrid(coordinates = coordinates,
                                     filterGrid = filterGrid)
        }
        #----------------------------------------------------------------------#
        # Reduce resolution
        # could be updated to use KNN
        #----------------------------------------------------------------------#
        if(tensorResolution < 1){
          .tensorRes(verbose)
          coordinates <- .reduceTensorResolution(coordinates = coordinates,
                                            tensorResolution = tensorResolution)
          .adjCounts(verbose)
          counts <- .adjustCounts(coordinates, counts,cores)
        }
        #----------------------------------------------------------------------#
        # TESSELATION TIME!
        #----------------------------------------------------------------------#
        .tess(verbose)
        tesselation <- deldir::deldir(x = as.numeric(coordinates$x),
                                      y = as.numeric(coordinates$y))
        #----------------------------------------------------------------------#
        # Filtering tiles
        #----------------------------------------------------------------------#
        .fTiles(verbose)
        filtered <- .filterTiles(tesselation,coordinates,filterThreshold)

        #----------------------------------------------------------------------#
        # Resterise tiles
        #----------------------------------------------------------------------#
        .raster(verbose)
        tiles <- .rasterise(filtered, cores)

        vesalius <- .updateVesalius(vesalius=vesalius,
                                    data=tiles,
                                    slot="tiles",
                                    commit = as.list(match.call()),
                                    defaults = as.list(args(buildVesaliusEmbeddings)),
                                    append=FALSE)
    } else {
      coordinates <- vesalius@tiles
      counts <- getCounts(vesalius,"raw")
    }
    #--------------------------------------------------------------------------#
    # Now we can start creating colour embeddings
    # This section can be run multiple times
    # for now we dont want to have multiple "tiles" options
    # It's going to make things really m
    #--------------------------------------------------------------------------#

    #counts <- .checkCounts(counts,verbose)
    .buildSO(verbose)
    counts <- .processCounts(counts,
                             vesalius= vesalius,
                             commit = as.list(match.call()),
                             method = norm,
                             nfeatures = nfeatures,
                             min.cutoff = min.cutoff)



    vesalius <- .updateVesalius(vesalius=vesalius,
                                  data=counts$norm,
                                  slot="counts",
                                  commit = as.list(match.call()),
                                  defaults = as.list(args(buildVesaliusEmbeddings)),
                                  append=TRUE)

    #--------------------------------------------------------------------------#
    # Embeddings
    #--------------------------------------------------------------------------#
    if(!method[1L] %in% c("PCA","PCA_L","UMAP","LSI","LSI_UMAP")){
        method <- "none"
    }
    embeds <- switch(method[1L],
                     "PCA" = .embedPCA(counts$SO,pcs = pcs,cores = cores,
                              verbose = verbose),
                     "PCA_L" = .embedPCAL(counts$SO,pcs = pcs,cores = cores,
                                verbose = verbose),
                     "UMAP" = .embedUMAP(counts$SO,pcs = pcs,verbose),
                     "LSI" = .embedLSI(counts$SO,pcs = pcs,remove_LSI1),
                     "LSI_UMAP" = .embedLSIUMAP(counts$SO,pcs = pcs,remove_LSI1),
                     "none" = stop("Unsupported embedding type!"))

    vesalius <- .updateVesalius(vesalius=vesalius,
                                data=embeds,
                                slot="embeddings",
                                commit = as.list(match.call()),
                                defaults = as.list(args(buildVesaliusEmbeddings)),
                                append=TRUE)

    vesalius <- .updateVesalius(vesalius=vesalius,
                                data=embeds,
                                slot="activeEmbeddings",
                                commit = as.list(match.call()),
                                defaults = as.list(args(buildVesaliusEmbeddings)),
                                append=FALSE)


    .simpleBar(verbose)
    return(vesalius)
}




#------------------------/ Filtering  Beads /----------------------------------#
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


#------------------------/ Reducing Resolution /-------------------------------#
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
    coordBar <- coordBar[sapply(strsplit(coordBar,"_et_"),length) > 1]
    if(length(coordBar) ==0){
        return(counts)
    }

    #--------------------------------------------------------------------------#
    # next we merge counts together when barcodes have been merged
    #--------------------------------------------------------------------------#
    tmpBar <- strsplit(coordBar,"_et_")


    empty <- parallel::mclapply(tmpBar, function(coord,count){

        tmp <- rowSums(count[,coord])
        return(tmp)
    },count = counts, mc.cores = cores)

    empty <- do.call("cbind",empty)
    if(is.null(dim(empty)) & length(empty) !=0){
        empty <- Matrix(empty,ncol=1)
    }
    colnames(empty) <- coordBar
    
    merged <- cbind(counts[,!colnames(counts) %in% unlist(unique(tmpBar))],empty)
    #--------------------------------------------------------------------------#
    # next we remove any barcodes that were dropped during filtering
    #--------------------------------------------------------------------------#
    merged <- merged[,colnames(merged) %in% coordinates$barcodes]
    return(merged)
}
#------------------------/ Creating pixel tiles /------------------------------#
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
        cent <- which(maxX == round(indx) &
                      maxY == round(indy))
        centers <- rep(0, length(maxX))
        centers[cent] <- 1
        tile <- data.frame("barcodes" = rep(filtered$coordinates$barcodes[idx],
                                            times = length(maxX)),
                           "x" = maxX,
                           "y" = maxY,
                           "origin" = centers)
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


#------------------------/ Preprocessing counts /------------------------------#
.processCounts <- function(counts,vesalius,commit, method = "log",nfeatures = 2000, min.cutoff = "q5"){
    #--------------------------------------------------------------------------#
    # We are still going to use Seurat for now
    # rememmber that if we do decide to change things
    # we have to change things in the embbeddings as well
    #--------------------------------------------------------------------------#
    counts <- CreateSeuratObject(counts, assay ="Spatial")
    counts <- switch(method[1L],
                    "log" = .logNorm(counts, nfeatures),
                    "SCT" = .SCTransform(counts,assay= "Spatial",
                                         nfeatures = nfeatures),
                    "TFIDF" = .TFIDFNorm(counts,min.cutoff=min.cutoff),
                    "raw" = .rawNorm(counts))


    #if(length(vesalius@log@counts)>0){
    #  last <- vesalius@log@counts
    #  last <- sapply(last, function(x){
    #      return(filter(x,Argument == "norm") %>% select(Value) %>% as.character())
    #  })
    #  if(!is.null(commit[["norm"]])){
    #      new <- ifelse(any(commit[["norm"]] %in% last),FALSE,TRUE)
    #  } else{
    #      new <- FALSE
    #  }
    #}else{
    #  new <- TRUE
    #}

  #  if(new){
  #      update <- list("update")
  #      names(update) <- "update"
  #      counts <- c(counts,update)
  #  } else {
  #    update <- list("noUpdate")
  #    names(update) <- "update"
  #    counts <- c(counts,update)
  #  }
  return(counts)
}

.rawNorm <- function(counts){
    #--------------------------------------------------------------------------#
    # Essentially we want people to be able to parse their matrix
    # If they want to use a different type of norm method that is not present
    # or play around with parameters not provided by vesalius
    # they can do that and just always call norm
    # We are using this just for formating at the moment
    #--------------------------------------------------------------------------#
    normCounts <- list(GetAssayData(counts,slot = "counts"))
    names(normCounts) <- "raw"
    return(list("SO" = counts, "norm" = normCounts))
}

.logNorm <- function(counts, nfeatures){
  counts <- NormalizeData(counts,verbose = FALSE)
  counts <- ScaleData(counts, verbose = FALSE)
  counts <- FindVariableFeatures(counts,nfeatures = nfeatures,verbose = FALSE)
  normCounts <- list(GetAssayData(counts,slot = "data"))
  names(normCounts) <- "log"
  return(list("SO" = counts, "norm" = normCounts))
}

.SCTransform<- function(counts,assay= "Spatial",nfeatures){
    counts <- SCTransform(counts,assay= "Spatial",
              variable.features.n = nfeatures,verbose=FALSE)
    normCounts <- list(GetAssayData(counts,slot = "data"))
    names(normCounts) <- "SCT"
    return(list("SO" = counts, "norm" = normCounts))
}

.TFIDFNorm <- function(counts, min.cutoff) {
  counts <- RunTFIDF(counts)
  counts <- FindTopFeatures(counts, min.cutoff = min.cutoff)
  normCounts <- list(GetAssayData(counts, slot = "data"))
  names(normCounts) <- "TFIDF"
  return(list("SO" = counts, "norm" = normCounts))
}

#------------------------/ Normalising Embeds /--------------------------------#
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

#------------------------/ Color Embeddings /----------------------------------#
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
    colourMatrix <- list(as.matrix(colourMatrix))
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


.embedLSI <- function(counts,pcs = pcs,remove_LSI1){

  #--------------------------------------------------------------------------#
  # Run partial singular value decomposition(SVD) on TF-IDF normalized matrix
  #--------------------------------------------------------------------------#
  svd <- RunSVD(counts, n = pcs + 1, verbose = FALSE)

  #--------------------------------------------------------------------------#
  # Getting embedding values and normalize
  #--------------------------------------------------------------------------#
  if (remove_LSI1 == TRUE) {
    embeddings <-
      Embeddings(svd[["lsi"]])[, -1]
  } else{
    embeddings <-
      Embeddings(svd[["lsi"]])[, 1:30]
  }

  embeddings <- apply(embeddings,2,function(x){
    x <- (x - min(x)) / (max(x) - min(x))
    return(x)
  })

  colourMatrix <- list(as.matrix(embeddings))
  names(colourMatrix) <- "LSI"
  return(colourMatrix)

}




.embedLSIUMAP <- function(counts,pcs = pcs,remove_LSI1){

  #--------------------------------------------------------------------------#
  # Run partial singular value decomposition(SVD) on TF-IDF normalized matrix
  #--------------------------------------------------------------------------#
  svd <- RunSVD(counts, n = pcs + 1, verbose = FALSE)

  if (remove_LSI1 == TRUE) {
    reduc <-
      RunUMAP(
        svd,
        reduction = 'lsi',
        dims = 2:(pcs + 1),
        n.components = 3L,
        verbose = F
      )
  } else{
    reduc <-
      RunUMAP(
        svd,
        reduction = 'lsi',
        dims = 1:pcs,
        n.components = 3L,
        verbose = F
      )
  }

  #--------------------------------------------------------------------------#
  # Getting embedding values and normalize
  #--------------------------------------------------------------------------#
  embeddings <- FetchData(reduc, c("UMAP_1", "UMAP_2", "UMAP_3"))

  embeddings <- apply(embeddings, 2, function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  })

  embeddings <- list(as.matrix(embeddings))
  names(embeddings) <- "LSI_UMAP"
  return(embeddings)

}
