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
#' @param filterThreshold numeric (range 0 -1) describing the quantile threshold
#' at which barcodes and tiles should be retained (seed details)
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
#' barcodes). Second, Vesalius filters out tiles that are tied to the outer box
#' and that exceed a certain area threshold. These triangles are nearly always
#' related to boundary tiles. Both of these filtering steps are controlled via
#' the \code{filterThreshold} argument.
#'
#' A filterThreshold of 0.99 means that 99 \% or barcodes and tile triangles
#' will be retained.
#'
#'
#' @return Returns a Vesalius data.frame with "barcodes","x","y","cc",
#' "value",and "tile".
#' @examples
#' \dontrun{
#' data(Vesalius)
#' }


buildImageArray <- function(coordinates,
                            sliceID = 1,
                            invert=FALSE,
                            na.rm = TRUE,
                            resolution = 100,
                            filterThreshold=0.999,
                            interpolation_type =1,
                            cores=1,
                            verbose = TRUE){
  .simpleBar(verbose)
  #----------------------------------------------------------------------------#
  # check data too see if it contains right slices
  # We will run everything on one slice first and then just rebuild at the end
  # at least we don't need to run tesselation on each slice
  #----------------------------------------------------------------------------#

  coordinates <- .checkVesalius(coordinates, sliceID,verbose)

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
  if(filterThreshold <1){
    .distanceBeads(verbose)
    #--------------------------------------------------------------------------#
    # We only need one "set" of coordinates - once we have distance metrics
    # we can filter out on all other
    # Let's just assume for now there might be more than one slice
    #--------------------------------------------------------------------------#
    tmp <- coordinates %>% filter(cc ==1)
    idx <- seq_len(nrow(tmp))
    minDist <- parallel::mclapply(idx, function(idx,mat){
                        xo <- mat$x[idx]
                        yo <- mat$y[idx]
                        xp <- mat$x
                        yp <- mat$y
                        distance <- sqrt(((abs(xp-xo))^2 + (abs(yp-yo))^2))
                        distance <- distance[distance !=0]
                        distance <- sort(distance,decreasing = FALSE)[1:25]
                        return(mean(distance))
    }, tmp, mc.cores=cores)
    minDist <- unlist(minDist)
    distanceThrehsold <- quantile(minDist,filterThreshold)
    keeps  <- tmp$barcodes[minDist <= distanceThrehsold]
    tmp <- tmp %>% filter(barcodes %in% keeps)
  } else {
    tmp <- coordinates %>% filter(cc ==1)

  }
  #----------------------------------------------------------------------------#
  # TESSELATION TIME!
  # Only need to tesslated on one slice - all the other will be the same tiles
  # We will do all the replacements later
  #----------------------------------------------------------------------------#
  .tess(verbose)
  tesselation <- deldir::deldir(x = as.numeric(tmp$x),
                                y = as.numeric(tmp$y))

  #----------------------------------------------------------------------------#
  # Extracting Voronoi/dirichlet tessaltion - not intested in Delauney
  #----------------------------------------------------------------------------#

  ## Voronoi tessaltion coordinates
  tessV <- tesselation$dirsgs
  ## If the "point" is on the edge of the box used for tessaltion
  tessV <- tessV[!tessV$bp1 & !tessV$bp2,]
  #tessV <- .filterTesselation(tessV)
  ## Contains original data - just to ensure that that indeces line up
  coord <- cbind(tmp,tesselation$summary)
  #----------------------------------------------------------------------------#
  # This section is equivalent to triangle rasterisation
  #----------------------------------------------------------------------------#
  .raster(verbose)
  allIdx <- parallel::mclapply(seq_len(nrow(tessV)),.fillTesselation,
    tess = tessV,coord= coord,
    ogCoord = coordinates,
    mc.cores = cores)

  #----------------------------------------------------------------------------#
  # Filtering triangle that exceed max tile size
  #----------------------------------------------------------------------------#
  allIdx <- .imageBuild(allIdx,filterThreshold,verbose)

  #----------------------------------------------------------------------------#
  # Rebuilding new cooirdinates based on tesselation borders
  # cooridinates are shifted for both tesselation cooridantes and starting coord
  #----------------------------------------------------------------------------#
  lpx <- round(min(allIdx$x))
  lpy <- round(min(allIdx$y))


  allIdx <- mutate(allIdx, x = x -lpx +1,y=y-lpy +1) %>%
            tibble

  #----------------------------------------------------------------------------#
  # Decreasing reolsution
  #----------------------------------------------------------------------------#
  if(resolution<100){
      .res(verbose)
      allIdx <- .resShift(allIdx, resolution,interpolation_type,na.rm)
  }
  .simpleBar(verbose)
  return(allIdx)

}


.checkVesalius <- function(coordinates, sliceID,verbose){
    #--------------------------------------------------------------------------#
    # First let's check if we have multiple slice in input data
    #--------------------------------------------------------------------------#

    slices <- unique(coordinates$slice)
    inSlice <- slices %in% sliceID

    #--------------------------------------------------------------------------#
    # Let's do some checks
    #--------------------------------------------------------------------------#
    if(sum(inSlice) > 1){
        warning("More than one slice provided!
                 Only lowest slice value will be used",
                immediate. = TRUE)
    } else if(sum(inSlice) <1){
        stop("SliceID is not present in coordinates.")
    }
    .checkVes(min(sliceID),verbose)
    coordinates <- coordinates %>% filter(slice == min(sliceID))

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




.imageBuild <- function(allIdx, filterThreshold,verbose){
    #--------------------------------------------------------------------------#
    # First get area and extract distribution
    #--------------------------------------------------------------------------#
    if(filterThreshold <1){
      .filter(verbose)
      triArea1 <- sapply(allIdx,function(x){
                        return(x[[2]][1])
      })
      triArea2 <- sapply(allIdx,function(x){
                        return(x[[2]][2])
      })

      areaThresh <- quantile(c(triArea1,triArea2),filterThreshold)
    }


    #--------------------------------------------------------------------------#
    # Next filter triangle that exceed limit
    #--------------------------------------------------------------------------#
    tri1 <- lapply(allIdx,function(x){
                  return(x[[1]][[1]])
    })
    if(filterThreshold<1){
      tri1 <- tri1[triArea1 < areaThresh]
    }

    tri1 <- as.data.frame(do.call("rbind",tri1))

    tri2 <- lapply(allIdx,function(x){
                  return(x[[1]][[2]])
    })
    if(filterThreshold <1){
      tri2 <- tri2[triArea2 < areaThresh]
    }

    tri2 <- as.data.frame(do.call("rbind",tri2))

    #--------------------------------------------------------------------------#
    # Rebuild full data frame
    #--------------------------------------------------------------------------#
    tri <- rbind(tri1,tri2)
    return(tri)


}

.fillTesselation <- function(idx,tess,coord,ogCoord){
    #--------------------------------------------------------------------------#
    ## Get data
    #--------------------------------------------------------------------------#
    allEdge <- tess[idx,]

    #--------------------------------------------------------------------------#
    ## First split data between both triangles
    # Why am idoing this? Not sure to be honest
    # could be good to refactor this as well
    # I think it was so I could use multi core at a higher level
    # Yep that was it. I could run both triangles in parallel
    # But it would not decrease comp time that much and we don't want to mess
    # around with weird core scheduling and what not
    #--------------------------------------------------------------------------#
    t1 <- .fillTesselationTriangle(allEdge,coord,ogCoord,tri = 1)
    a1 <- nrow(t1)

    t2 <- .fillTesselationTriangle(allEdge,coord,ogCoord,tri = 2)
    a2 <- nrow(t2)


    #--------------------------------------------------------------------------#
    # Rebuild list with filled triangle for both points
    #--------------------------------------------------------------------------#
    allIn <- list(t1,t2)
    return(list(allIn,c(a1,a2)))

}

.fillTesselationTriangle <- function(allEdge,coord,ogCoord,tri){
  if(tri ==1){

    #--------------------------------------------------------------------------#
    # Boundaries of tesselation
    #--------------------------------------------------------------------------#
    lpx <- round(min(c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"]))) - 1
    hpx <- round(max(c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"]))) + 1
    lpy <- round(min(c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"]))) - 1
    hpy <- round(max(c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"]))) + 1

    maxPolygonX <- rep(seq(lpx,hpx), times = hpy - lpy +1)
    maxPolygonY <- rep(seq(lpy,hpy), each = hpx - lpx +1)

    #--------------------------------------------------------------------------#
    # Creating triangles
    #--------------------------------------------------------------------------#
    x <- round(c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"]))
    y <- round(c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"]))

    #--------------------------------------------------------------------------#
    # Fill triangles with all point in that space
    #--------------------------------------------------------------------------#
    cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
    maxX <- maxPolygonX[cell %in% c(1,2,3)]
    maxY <- maxPolygonY[cell %in% c(1,2,3)]

    ccVals <- ogCoord[ogCoord$barcodes == coord$barcodes[allEdge$ind1],
                      c("cc","value")]
    ccVals <- ccVals[rep(seq_len(nrow(ccVals)),each = length(maxX)),]
    maxX <- rep(maxX, times = 3)
    maxY <- rep(maxY, times = 3)
    barcode <- rep(coord[allEdge$ind1,"barcodes"],length(maxX))
    cent <- which(maxX == round(coord[allEdge$ind1,"x"]) &
                  maxY == round(coord[allEdge$ind1,"y"]))
    centers <- rep(0, length(maxX))
    centers[cent] <- 1

    allIn <- data.frame(barcode,maxX,maxY,ccVals,centers)
    colnames(allIn) <- c("barcodes","x","y","cc","value","tile")
    return(allIn)

  } else {
    #--------------------------------------------------------------------------#
    # Boundaries of tesselation
    #--------------------------------------------------------------------------#
    lpx <- round(min(c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"]))) - 1
    hpx <- round(max(c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"]))) + 1
    lpy <- round(min(c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"]))) - 1
    hpy <- round(max(c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"]))) + 1

    maxPolygonX <- rep(seq(lpx,hpx), times = hpy - lpy +1)
    maxPolygonY <- rep(seq(lpy,hpy), each = hpx - lpx +1)

    #--------------------------------------------------------------------------#
    # Creating triangles
    #--------------------------------------------------------------------------#
    x <- round(c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"]))
    y <- round(c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"]))

    #--------------------------------------------------------------------------#
    # Fill triangles with all point in that space
    #--------------------------------------------------------------------------#
    cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
    maxX <- maxPolygonX[cell %in% c(1,2,3)]
    maxY <- maxPolygonY[cell %in% c(1,2,3)]
    ccVals <- ogCoord[ogCoord$barcodes == coord$barcodes[allEdge$ind2],
                      c("cc","value")]
    ccVals <- ccVals[rep(seq_len(nrow(ccVals)),each =length(maxX)),]
    maxX <- rep(maxX, times = 3)
    maxY <- rep(maxY, times = 3)
    barcode <- rep(coord[allEdge$ind2,"barcodes"],length(maxX))
    cent <- which(maxX == round(coord[allEdge$ind2,"x"]) &
                  maxY == round(coord[allEdge$ind2,"y"]))
    centers <- rep(0, length(maxX))
    centers[cent] <- 1
    allIn <- data.frame(barcode,maxX,maxY,ccVals,centers)
    colnames(allIn) <- c("barcodes","x","y","cc","value","tile")
    return(allIn)
  }

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
