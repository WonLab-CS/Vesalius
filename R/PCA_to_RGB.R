################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/PCA RGB Functions/--------------------------------------#

#' Embed PCA loadings into RGB colour space
#' @param SO Seurat Object after normalisation, scaling and
#' finding variable features.
#' @param slices number of PCA slices to consider (Integer/numeric).
#' Slice 1 will embed PC 1 to PC3, slice 2 will embed PC4 to PC6, etc
#' @param adjusted logical indicating if RGB code should be adjusted to
#' colour contribution per channel.
#' @param rgbWeight logical indicating if RGB channels should be weighted by
#' the variance associated to the PC it used for embedding.
#' @param countWeight logical describing if colour embedding
#' should be weighted by count number.
#' @param conserveSparse logical indicating if sparse matrix format
#' should be conserved.
#' @param verbose logical - progress message output.
#' @details The core concept of Vesalius is to embed PCA loading values into
#' the RGB colour space. This is achieved by:
#'
#' \deqn{ B_{(c,i)} = \sum_{i=1} \mid L_(c,i) \mid }
#'
#' with B_((i,c)), barcode i in color channel c and L, the loading values
#' associated to barcode i in color channel c.
#' Notice that we take the absolute value of PCA loadings.
#' In this case, we are only interested in variance and
#' not direction of variance.
#'
#' Barcode values are then min/max normalised to ensure that all colour
#' codes ranged from 0 to 1.
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


rgbPCA<- function(SO,
                  slices = 1,
                  adjusted = FALSE,
                  rgbWeight=FALSE,
                  countWeight = FALSE,
                  conserveSparse = TRUE,
                  verbose=TRUE){


    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # First we get the count values
    #--------------------------------------------------------------------------#
    if(is(SO) == "Seurat"){
        slide <- GetAssayData(SO, slot = "data")
    } else {
       stop("SO is not a Seurat Object.")
    }
    #--------------------------------------------------------------------------#
    # Sparse matrices are more memory efficient but also slower
    #--------------------------------------------------------------------------#
    if(!conserveSparse){
        .consSparse(verbose)
        slide <- as.matrix(slide)
    }

    #--------------------------------------------------------------------------#
    # PCA from seurat
    #--------------------------------------------------------------------------#
    .pca(verbose,slices)
    pca <- RunPCA(SO,npcs = slices*3,verbose =FALSE)

    image_slice <- vector("list",slices)
    slice_start <- seq(1,slices*3,by=3)
    #--------------------------------------------------------------------------#
    # Going through each slice - This could be changed in the future
    # We could go to grey scale and then use all PC
    #--------------------------------------------------------------------------#
    for(sl in seq_len(slices)){
    #--------------------------------------------------------------------------#
    # Getting those loading values - also computing variance per PC
    # Variance by PC is use to weight the channels
    # Same Principle is used when converting colour to grey scale
    #--------------------------------------------------------------------------#
      loadings<-Loadings(pca[["pca"]])[,seq(slice_start[sl], slice_start[sl]+2)]

      if(rgbWeight){
        .pcadj(verbose,sl)
        varPerPC <- apply(loadings,2,var)
        varPerPC <- varPerPC /sum(varPerPC)
      }


      loadings <- apply(loadings,2,function(x)return(abs(x)))


      # RGB conversion
      rgb <- vector("list", 3)
      names(rgb) <- c("R","G","B")
      .rgb(verbose,sl)
      #------------------------------------------------------------------------#
      # Going through each colour
      #------------------------------------------------------------------------#
      for(i in seq_along(rgb)){
          rgb[[i]] <-rep(0,ncol(slide))
          names(rgb[[i]]) <- colnames(slide)
          for(j in seq_len(ncol(slide))){
             #-----------------------------------------------------------------#
             # Getting all genes and associated loading values
             #-----------------------------------------------------------------#
             genes <- rownames(slide)[slide[,j]!=0]

             cs <- loadings[,i][names(loadings[,i]) %in% genes]

             #-----------------------------------------------------------------#
             # Just in case a barcode does not contain genes present in the
             # loadings
             # This could happen if the selection of variable
             # features is too stringent
             #-----------------------------------------------------------------#
             if(length(cs) ==0){
               rgb[[i]][j] <- 0
             }else {
             #-----------------------------------------------------------------#
             # weighted by count number
             # Generally There is very little difference
             #-----------------------------------------------------------------#
               if(countWeight){

                  gcount <- slide[slide[,j]!=0,j]
                  gcount <- gcount[!is.na(match(genes,names(cs)))]

                  cs <- cs * gcount
                  cs <- tail(cumsum(cs),1)
               } else {
                  cs <- tail(cumsum(cs),1)
               }
               rgb[[i]][j] <- cs
             }

          }
      }

      #------------------------------------------------------------------------#
      ## Normalising RGB channels
      #------------------------------------------------------------------------#
      .norm(verbose)
      rgb <- lapply(rgb, function(x){
                    x <- (x - min(x)) / (max(x) - min(x))
                    return(x)
      })
      #------------------------------------------------------------------------#
      # adjusted RGB colour - NOT THE SAME AS THE WEIGHTED PC
      #------------------------------------------------------------------------#
      if(adjusted){
          .adj(verbose)
          sums <- sapply(seq_along(rgb[[1]]), function(idx,rgb){
                          return(sum(rgb[[1]][idx],
                                     rgb[[2]][idx],
                                     rgb[[3]][idx]))
          },rgb)
          rgb <- lapply(rgb, function(rgb,sums){
                        return(rgb/sums)
          },sums)
      }
      #------------------------------------------------------------------------#
      # Weighted colours based on PC variance
      #------------------------------------------------------------------------#
      if(rgbWeight){

          for(i in seq_along(rgb)){
            rgb[[i]] <- rgb[[i]] * varPerPC[i]
          }
      }

      image_slice[[sl]] <- rgb
    }
    #--------------------------------------------------------------------------#
    # List with three different levels.
    # Level 1 : each list element is a different slice
    # Level 2 : each list element is R G B channel for that slice
    # Level 3 : numeric vectors containing colour code for each barcode
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    # We can use this list to get the coordinates and assign the RGB coulours
    # to their correct location
    #--------------------------------------------------------------------------#

    coordinates <- getSeuratCoordinates(SO)

    image_slice <- lapply(seq_along(image_slice),.assignRGBtoPixel,
                          image_slice,coordinates)

    if(length(image_slice) >1){
       image_slice <- do.call("rbind",image_slice)
    } else {
       image_slice <- image_slice[[1]]
    }


    .simpleBar(verbose)
    return(image_slice)
}



#' Embed 3D UMAP projections into RGB colour space
#' @param SO Seurat Object after normalisation, scaling and
#' finding variable features.
#' @param pcs number of Principle components to use for UMAP projections
#' @param adjusted logical indicating if RGB code should be adjusted to
#' colour contribution per channel.
#' @param conserveSparse logical indicating if sparse matrix format
#' should be conserved.
#' @param verbose logical - progress message output.
#' @details Emedding 3D UMAP projections into the RGB color space is quite
#' straight forward: Each UMAP dimension is extracted and min/max normalised to
#' ensure that each value reanges between 0 and 1. The first dimension is
#' assigned to red color channel, the second dimension to the green color
#' channel and the third dimension assgined to the blue color channel.
#' UMAP projections consider more varience per channel than PCA embeddings but
#' this additional information in further compressed. Using PCA slice provides
#' a way to visualize and explore different anatomical structure in more depth.
#' Choosing between rgbPCA and rgbUMAP will depend on the question at hand and
#' the depth of the analysis.
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","slice"
#' @examples
#' \dontrun{
#' data(Vesalius)
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbUMAP(image,pcs = 30)
#' }


rgbUMAP<- function(SO,
                  pcs = 30,
                  adjusted = FALSE,
                  conserveSparse = TRUE,
                  verbose=TRUE){


    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # First we get the count values
    #--------------------------------------------------------------------------#
    if(is(SO) == "Seurat"){
        slide <- GetAssayData(SO, slot = "data")
    } else {
       stop("SO is not a Seurat Object.")
    }
    #--------------------------------------------------------------------------#
    # Sparse matrices are more memory efficient but also slower
    #--------------------------------------------------------------------------#
    if(!conserveSparse){
        .consSparse(verbose)
        slide <- as.matrix(slide)
    }

    #--------------------------------------------------------------------------#
    # PCA from seurat
    #--------------------------------------------------------------------------#

    reduc <- RunPCA(SO,npcs = pcs,verbose =FALSE)

    reduc <- RunUMAP(reduc,dims =seq_len(pcs),n.components = 3L,verbose=F)


    #--------------------------------------------------------------------------#
    # Getting those loading values - also computing variance per PC
    # Variance by PC is use to weight the channels
    # Same Principle is used when converting colour to grey scale
    #--------------------------------------------------------------------------#
    embeddings <- FetchData(reduc, c("UMAP_1","UMAP_2","UMAP_3"))

    embeddings <- apply(embeddings,2,function(x){
                       return((x-min(x))/(max(x) - min(x)))
                     })

    #------------------------------------------------------------------------#
    # adjusted RGB colour - NOT THE SAME AS THE WEIGHTED PC
    #------------------------------------------------------------------------#
    if(adjusted){
        .adj(verbose)
        sums <- sapply(seq_along(rgb[[1]]), function(idx,rgb){
                          return(sum(rgb[[1]][idx],
                                     rgb[[2]][idx],
                                     rgb[[3]][idx]))
        },rgb)
        rgb <- lapply(rgb, function(rgb,sums){
                        return(rgb/sums)
        },sums)
    }

    #--------------------------------------------------------------------------#
    # We can use this list to get the coordinates and assign the RGB coulours
    # to their correct location
    #--------------------------------------------------------------------------#

    coordinates <- getSeuratCoordinates(SO)

    tmp <- vector("list", 3)
    for(i in seq_len(ncol(embeddings))){
        tmp[[i]] <- data.frame(rownames(coordinates),
                               coordinates[,c("x","y")],
                               rep(i,length(embeddings[,i])),
                               embeddings[,i],
                               rep(1, length(embeddings[,i])))

    }
    tmp <- do.call("rbind", tmp)
    colnames(tmp) <- c("barcodes","x","y","cc","value","slice")
    tmp <- tmp[!is.na(tmp$value),]



    .simpleBar(verbose)
    return(tmp)
}



#Internal to assign colour code back to its location
.assignRGBtoPixel <- function(idx,rgb,coordinates, na.rm = TRUE){
    #--------------------------------------------------------------------------#
    # This is just to create a data frame with colour channel and colour value
    # This follows the cimg format
    # we will also add a slice clause
    #--------------------------------------------------------------------------#

    tmp <- vector("list", 3)
    rgb <- rgb[[idx]]
    for(i in seq_along(rgb)){
        code <- rep(NA,nrow(coordinates))
        locs <- match(rownames(coordinates),names(rgb[[i]]))
        code[!is.na(locs)] <- rgb[[i]][locs[!is.na(locs)]]
        tmp[[i]] <- data.frame(rownames(coordinates),
                               coordinates[,c("x","y")],
                               rep(i,length(code)),
                               code,
                               rep(idx,length(code)))

    }
    tmp <- do.call("rbind", tmp)
    colnames(tmp) <- c("barcodes","x","y","cc","value","slice")
    if(na.rm){
        tmp <- tmp[!is.na(tmp$value),]
    }

    return(tmp)
}





#' buildImageArray creating an Image array from coloured coordinates
#'
#' @param coordinates a Vesalius data.frame with "barcodes","x","y","cc",
#' "value","slice". Output from \code{rgbPCA}.
#' @param sliceID integer - PC slice used in image building. Only one value!
#' @param invert logical describing if colour pattern should be inverted
#' (i.e 1-R,1-G,1-B)
#' @param na.rm logical indicating if NA should be removed from bead list.
#' @param resolution numeric (range 0 - 100) describing a percentage of original
#' image size. Used to reduce image size.
#' @param filterGrid numeric between (range 0 - 1) describing size of grid used
#' to filter stray beads.
#' @param filterThreshold numeric (range 0 - 1) describing the quantile threshold
#' at which barcodes and tiles should be retained (seed details)
#'
#' @param interpolation_type Method of interpolation during image resizing:
#' \describe{
#'    \item{-1}{no interpolation: raw memory resizing.}
#'    \item{0}{no interpolation: additional space is filled according to
#'             boundary_conditions.}
#'    \item{1}{Nearest-neighbor interpolation}
#'    \item{2}{moving average interpolation}
#'    \item{3}{linear interpolation}
#'    \item{4}{grid interpolation}
#'    \item{5}{cubic interpolation}
#'    \item{6}{lanczos interpolation}
#' }
#'
#' @param cores integer describing number of cores used (Default = 1)
#' @param verbose logical - progress message output.
#' @details Vesalius converts PCA loading values into an RGB colour code
#' associated with a coordinates on a Spatial Transcriptomic Assay. The
#' \code{buildImageArray} converts coloured coordinates into an actual image
#' via the use of Voronoi diagrams and Tile rasterisation.
#'
#' Each barcodes coordinate is used as the center point of a voronoi tile and
#' an artificial "box" is created surrounding barcode coordinates. This
#' creates superfluous tiles that can be detrimental to further analysis.
#' To remove excessive tiles, Vesalius takes a two step process. First, Vesalius
#' removes any barcode that is too far away from other barcodes (likely stray
#' barcodes). This is controlled by the filterGrid argument.
#' Second, Vesalius filters out tiles that exceed a certain area threshold.
#' This is controlled by the filterThreshold argument.
#'
#'
#' A filterThreshold of 0.99 means that 99 \% or barcodes and tile triangles
#' will be retained.
#'
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value",and "tile".
#' @examples
#' \dontrun{
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' # Slice ID = PCA slice
#' image <- buildImageArray(image, sliceID = 1)
#' # If rgbUMAP was run - no need to specifiy slice
#' image <- buildImageArray(image)
#' }


buildImageArray <- function(coordinates,
                            sliceID = 1,
                            invert=FALSE,
                            na.rm = TRUE,
                            resolution = 100,
                            filterGrid =0.01,
                            filterThreshold=1,
                            interpolation_type =1,
                            cores=1,
                            verbose = TRUE){
  .simpleBar(verbose)
  #----------------------------------------------------------------------------#
  # check data too see if it contains right slices
  # We will run everything on one slice first and then just rebuild at the end
  # at least we don't need to run tesselation on each slice
  #----------------------------------------------------------------------------#

  #coordinates <- .checkVesalius(coordinates, sliceID,verbose)

  #----------------------------------------------------------------------------#
  # Type changed
  #----------------------------------------------------------------------------#
  .typeChange(verbose)
  coordinates$value <- as.numeric(coordinates$value)
  coordinates$cc <- as.integer(coordinates$cc)
  #----------------------------------------------------------------------------#
  # Inverting colours
  #----------------------------------------------------------------------------#
  if(invert){
      .invertCols(verbose)
      coordinates$value <- 1 - coordinates$value
  }

  #----------------------------------------------------------------------------#
  # Removing all "outer point" - tend to make tesselation messy
  #----------------------------------------------------------------------------#
  if(filterGrid != 0 & filterGrid != 1){
    .distanceBeads(verbose)
    coordinates <- .filterGrid(coordinates = coordinates,
                               filterGrid = filterGrid)
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
  coordinates <- right_join(coordinates,tiles, by = "barcodes")

  coordinates <- coordinates[, c("barcodes","x.y","y.y","cc","value","tile")]

  colnames(coordinates) <- c("barcodes","x","y","cc","value","tile")

  #----------------------------------------------------------------------------#
  # Decreasing reolsution
  #----------------------------------------------------------------------------#
  if(resolution<100){
      .res(verbose)
      coordinates <- .resShift(coordinates, resolution,interpolation_type,na.rm)
  }
  .simpleBar(verbose)
  return(coordinates)

}





.resShift <- function(image,resolution = 50,interpolation_type = 1,na.rm=TRUE){
  #----------------------------------------------------------------------------#
  # Converting image back image for resizing
  # Here we are using the imager function as the intrapolation is already
  # optimised and also provides a few different ways of doing so
  # NOTE testing regularisation - will need to integrate a bit better
  #----------------------------------------------------------------------------#

  img <- as.cimg(select(image, c("x", "y", "cc", "value")))
  img <- resize(img,size_x = -(resolution),
                    size_y = -(resolution),
                    interpolation_type = interpolation_type)

  img <- tibble(as.data.frame(img))
  #----------------------------------------------------------------------------#
  # Gray scale looses the cc column so adding it back in for later
  #----------------------------------------------------------------------------#
  if(!"cc" %in% colnames(img)){

     img <- data.frame(img[,c("x","y")],rep(1,nrow(img)),img$value)
     colnames(img) <- c("x", "y", "cc", "value")
  }
  #----------------------------------------------------------------------------#
  # Scaling original coordinates and joining them the new table
  # This is so dumb - why did I not think of this before...
  #----------------------------------------------------------------------------#
  image$x <- floor((image$x / max(image$x))* max(img$x)) +1
  image$y <- floor((image$y / max(image$y))* max(img$y)) +1
  image <- image %>% distinct()

  img <- right_join(img, image, by  = c("x","y","cc")) %>%
         select(c("barcodes","x","y","cc","value.y","tile")) %>% distinct()
  colnames(img) <- c("barcodes","x","y","cc","value","tile")
  #----------------------------------------------------------------------------#
  # remove NAs
  #----------------------------------------------------------------------------#
  if(na.rm){
      img <- img %>% na.exclude()
  }
  return(img)
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

.rasteriseOld <- function(filtered,cores = 1){
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
                           "tile" = centers)
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




#' exportRGB.csv export RGB coordinates as csv.
#' @param coordinates Vesalius data.frame with barcodes, x, y, cc, value,and
#' slice.
#' @param file character - file name or vector of files names if multiple slices
#' will be exported.
#' @param slice integer - slice to be exported
#' @param split logical - If TRUE, slices will be exported in separately
#' @details Export rgb colour code using slice information. If no file name or
#' file names are provided, Vesalius will use its internal file name
#' (pca_to_rgb.csv).
#'
#' If one or more slice is provided, slice number will be included in file name.
#' If split is true, indvidual files for each selected files will be created.
#' @return csv file at desired location containing rgb coloured coordinates.

exportRGB.csv <- function(coordinates,
                          file = NULL,
                          slice = NULL,
                          split = TRUE){
  if(!is.null(slice)){
      coordinates <- coordinates[coordinates$slice %in% slice]
  }
  if(split){
      coordinates <- split(coordinates,coordinates$slice)
      for(i in seq_along(coordinates)){
          if(is.null(file[1L])){
              file <- paste0("pca_to_rgb_slice",coordinates[[i]]$slice[1L],
                           ".csv")
              write.csv(coordinates[[i]], file = file,row.names=F)
          }else{
              if(length(file) <length(coordinates)){
                 stop("Not enough file names provided")
              }

              write.csv(coordinates[[i]], file = file[i],row.names=F)
          }
      }
  } else {
      if(is.null(file)){
          file <- paste0("pca_to_rgb.csv")
      }
      write.csv(coordinates, file = file , row.names =F)
  }

}
