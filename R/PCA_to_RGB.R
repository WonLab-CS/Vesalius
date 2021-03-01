################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/PCA RGB Functions/--------------------------------------#

#' Embed PCA loading to RGB colour space
#' @param slide gene count matrix with barcodes as columns and genes as rows
#' @param SO Seurat Object after normalisation, scaling and finding variable features.
#' The description of how to build Seurat Objects is described in the Seurat Vignettes.
#' @param slices number of PCA slices to consider (Integer/numeric).
#' Slice 1 will embed PC 1 to PC3, slice 2 will embed PC4 to PC6, etc
#' @param adjusted logical indicating if RGB code should be adjusted to rgb (rgb coordinate)
#' @param countWeight logical describing if colour embeding should be weighted by count number.
#' @param conserveSparse logical indicating if sparse matrix format should be conserved.
#' @param trim numeric describing the quantile at which coulour histogram should be trimmed (two tailed).
#' Default set to 0 where no histogram trimming will be applied.
#' @details rgbPCA embeds PCA loadings in the RGB colour space by cumulative summing of PCA loadings for all gene
#' present on every bead. The absolute values of PCA loading are used for the summation.
#' The cumulative sum is then normalise in order to bound scores between 0 and 1. If countWegiht
#' is set to TRUE, loading value for each genes (per bead) are multiplied by the count number of that gene for each bead.
#' It should be noted that this process is only applied to genes present on a bead/spot. There is no consideration for
#' genes that are not present on the bead/spot.
#' @return a list with three different levels.
#' Level 1 : each list element is a different slice
#' Level 2 : each list element is R G B channel for that slice
#' Level 3 : numeric vectors containing colour code for each barcode

rgbPCA<- function(slide,SO,slices = 1,adjusted = FALSE,rgbWeight=FALSE,countWeight = FALSE, conserveSparse = TRUE, trim = 0){



    #--------------------------------------------------------------------------#
    # Sparse matrices are more memory efficient but also slower
    #--------------------------------------------------------------------------#
    if(!conserveSparse){
        slide <- as.matrix(slide)
    }

    #--------------------------------------------------------------------------#
    # PCA from seurat
    #--------------------------------------------------------------------------#
    pca <- RunPCA(SO,npcs = slices*3)

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
      loadings <-Loadings(pca[["pca"]])[,seq(slice_start[sl], slice_start[sl]+2)]

      if(rgbWeight){
        varPerPC <- apply(loadings,2,var)
        varPerPC <- varPerPC /sum(varPerPC)
      }


      loadings <- apply(loadings,2,function(x)return(abs(x)))


      # RGB conversion
      rgb <- vector("list", 3)
      names(rgb) <- c("R","G","B")
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
      ## trimming histogram - unused for now
      #------------------------------------------------------------------------#
      cutoff <- do.call("cbind",lapply(rgb,trimHistogram,trim))
      cutoff <- apply(cutoff,1,sum) ==3

      #------------------------------------------------------------------------#
      ## Normalising RGB channels
      #------------------------------------------------------------------------#
      rgb <- lapply(rgb, function(x,cutoff){
                    x <- x[cutoff]

                    x <- (x - min(x)) / (max(x) - min(x))
                    return(x)
      }, cutoff = cutoff)
      #------------------------------------------------------------------------#
      # adjusted RGB colour - NOT THE SAME AS THE WEIGHTED PC
      #------------------------------------------------------------------------#
      if(adjusted){
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
    return(image_slice)
}


#' Assigning the RGB colour code to its specific barcode location
#'
#' @param rgb list of colour codes in the format describe in the Detail section.  Output of rgbPCA.
#' @param coordinates data frame containing x and y coordinates of barcodes
#' @param na.rm logical describing if NA should be removed
#' @details assignRGBtoPixel remaps the colour code to x/y coordinates. \code{rgb} is the output of \code{rgbPCA}.
#' As such, it should be a list with three different levels.
#' Level 1 : each list element is a different slice
#' Level 2 : each list element is R G B channel for that slice
#' Level 3 : numeric vectors containing colour code for each barcode
#' @return a data frame with five columns: barcodes (barcode names), xcoord (x coordinates),
#' ycoord (y coordinates), R (Red colour channel),G (Green colour channel) ,and B (Blue Colour channel).

assignRGBtoPixel <- function(rgb,coordinates, na.rm = TRUE){
    # for each slice find coresponding barcodes
    for(i in seq_along(rgb)){
        code <- rep(NA,nrow(coordinates))
        locs <- match(coordinates$barcodes,names(rgb[[i]]))
        code[!is.na(locs)] <- rgb[[i]][locs[!is.na(locs)]]
        coordinates <- data.frame(coordinates,code)

    }
    colnames(coordinates) <- c("barcodes","xcoord","ycoord","R","G","B")
    if(na.rm){
        coordinates <- coordinates[!is.na(coordinates$R),]
    }
    # Data frame with barcode name, x/y coordinates and associated RGB colour code
    return(coordinates)
}






assingRGBtoPixelQuickBlock <- function(rgb, coordinates,resolution = 200,drop =TRUE,na.rm =TRUE){
  # for each slide find coresponding cooridinates
  for(i in seq_along(rgb)){
      code <- rep(NA,nrow(coordinates))
      locs <- match(coordinates$barcodes,names(rgb[[i]]))
      code[!is.na(locs)] <- rgb[[i]][locs[!is.na(locs)]]
      coordinates <- data.frame(coordinates,code)

  }
  colnames(coordinates) <- c("barcodes","xcoord","ycoord","R","G","B")
  if(na.rm){
      coordinates <- coordinates[!is.na(coordinates$R),]
  }
  ### this is going to be a pain
  ### Type coersion - has the tendency of converting values to "factors"
  ### Useful but a real pain in this case
  for(cols in seq(2,ncol(coordinates))){
    coordinates[,cols] <- as.numeric(as.character(coordinates[,cols]))
  }

  ### Building template Sparse Matrix for each pixel
  rgbMat <- Matrix(0,ncol=3,nrow = resolution * resolution)
  ## Matrix index in long format
  idx <- rep(seq(1,resolution), each = resolution)
  idy <- rep(seq(1,resolution), times = resolution)
  ## Barcode templates
  barcode <- rep("NONE",resolution * resolution)
  ## Coordinate chunking
  xcor <- seq(min(coordinates$xcoord), max(coordinates$xcoord),length.out = resolution+1) - min(coordinates$xcoord) + 1
  ycor <- seq(min(coordinates$ycoord), max(coordinates$ycoord),length.out = resolution+1) - min(coordinates$ycoord) + 1

  for(cols in seq_len(resolution)){

     for(rows in seq_len(resolution)){
       print(paste(cols,rows))
           ## Finding all barcodes in each "pixel"
          xtmp <- which(coordinates$xcoord >= xcor[cols] & coordinates$xcoord <= xcor[cols+1])
          ytmp <- which(coordinates$ycoord >= ycor[rows] & coordinates$ycoord <= ycor[rows+1])
          inter <- intersect(xtmp,ytmp)
          pos <- which(idx == cols & idy == rows)
          if(length(inter) == 0){
            ## just in case pixels are "empty" in terms of beads
            next()
          }
          ## average colours over all barcodes in pixel
          temp <- apply(coordinates[inter,c("R","G","B")],2,mean)
          # Assigne new colour code to its location in the sparse matrix
          rgbMat[pos,] <- as.matrix(temp, ncol = 3)
          ## adding all barcodes associated to that pixel to the barcode tag
          barcode[pos] <- paste0(barcode[pos],coordinates[inter,"barcodes"],sep = "_",collapse ="_")


     }
  }
  colnames(rgbMat) <- c("R","G","B")

  ## Removing all white/ black pixel -
  ## These represent empty spaces and are not required to rebuild the matrix
  if(drop){
    zeros <- apply(rgbMat,1,sum) != 0

    rgbMat <- rgbMat[zeros,]

    idx <- idx[zeros]
    idy <- idy[zeros]
    barcode <- barcode[zeros]
  }

  ## Return a data frame with the barcode list for each pixel,
  ## x and y coordinates in the matrix and
  ## the colour assigned to that matrix coordinate

  rgbMat <- data.frame(barcode,idx,idy,rgbMat)
  return(rgbMat)

}


#' buildImageArray creating an Image array from coloured coordinates
#'
#' @param coordninates data.frame containing the coordinatesof each spot/bead.
#' If colour coordinates were previosuly computed this data.frame may also contain
#' the RGB colour code associated to each location.
#' @param rgb list containing colour code for each bead. Required if colour
#' is not already contained in \code{coordinates}.
#' @param invert logical describing if colour pattern should be inverted (i.e 1-R,1-G,1-B)
#' @param na.rm logical indicating if NA should be removed from bead list.
#' @param cores numeric describing number of cores used (Default = 1)
#' @param verbose logical describing if progress message should be displayed in consolde
#' @details In order to build an Image type object, it is require to create a 4
#' dimensional array with rows as x coordinates, columns as y coordinates, slices and colour channels.
#' Note that the slice dimension is required by the \code{imager} package.
#' @return a 4 dimensional array - convertable to \code{cimg} objects.

buildImageArray <- function(coordinates,rgb=NULL,invert=FALSE,na.rm = TRUE,resolution = 1,filterThreshold=0.999,cores=1, verbose = TRUE){
  #----------------------------------------------------------------------------#
  # Class checks - this will be removed once S4 classes are functional
  #----------------------------------------------------------------------------#

  if(is.null(rgb)){
      if(!any(colnames(coordinates) %in% c("R","G","B"))){
        stop("Coordinates do not contain colour values")
      }
  } else if(!class(rgb) == "list"){
      stop("RGB - Unrecognised format")
  } else if(!is.null(rgb) & class(rgb)=="list") {
    .assignCoordinates(verbose)
    for(i in seq_along(rgb)){
        code <- rep(NA,nrow(coordinates))
        locs <- match(coordinates$barcodes,names(rgb[[i]]))
        code[!is.na(locs)] <- rgb[[i]][locs[!is.na(locs)]]
        coordinates <- data.frame(coordinates,code)

    }
    colnames(coordinates) <- c("barcodes","xcoord","ycoord","R","G","B")
    if(na.rm){
        coordinates <- coordinates[!is.na(coordinates$R),]
    }
  }

  #----------------------------------------------------------------------------#
  # Type changed
  #----------------------------------------------------------------------------#
  .typeChange(verbose)
  for(cols in seq(2,ncol(coordinates))){
    coordinates[,cols] <- as.numeric(as.character(coordinates[,cols]))
  }



  #----------------------------------------------------------------------------#
  # Inverting colours
  #----------------------------------------------------------------------------#
  if(invert){
      .invertCols(verbose)
      coordinates$R <- 1 - coordinates$R
      coordinates$G <- 1 - coordinates$G
      coordinates$B <- 1 - coordinates$B
  }

  #----------------------------------------------------------------------------#
  # Decreasing reolsution and making bigger tiles
  #----------------------------------------------------------------------------#
  if(resolution <1){
    .res(verbose)
    widthX <- max(coordinates$xcoord) - min(coordinates$xcoord) + 1
    widthY <- max(coordinates$ycoord) - min(coordinates$ycoord) + 1

    xbox <- round(seq(min(coordinates$xcoord),max(coordinates$xcoord), length.out = widthX * resolution))
    ybox <- round(seq(min(coordinates$ycoord),max(coordinates$ycoord), length.out = widthY * resolution))
    ## Testing new version might need to remove the section above
    coordinates <- .resShift(coordinates,xbox,ybox,cores)
  }



  #----------------------------------------------------------------------------#
  # Removing all "outer point" - tend to make tesselation messy
  #----------------------------------------------------------------------------#
  if(filterThreshold <1){
    .distanceBeads(verbose)
    idx <- seq_len(nrow(coordinates))
    minDist <- parallel::mclapply(idx, function(idx,mat){
                        xo <- mat$xcoord[idx]
                        yo <- mat$ycoord[idx]
                        xp <- mat$xcoord
                        yp <- mat$ycoord
                        distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)
                        distance <- distance[distance !=0]
                        distance <- sort(distance,decreasing = FALSE)[1:25]
                        return(mean(distance))
    }, coordinates, mc.cores=cores)
    minDist <- unlist(minDist)
    distanceThrehsold <- quantile(minDist,filterThreshold)
    coordinates <- coordinates[minDist <= distanceThrehsold,]
  }


  #----------------------------------------------------------------------------#
  # TESSELATION TIME!
  #----------------------------------------------------------------------------#
  .tess(verbose)
  tesselation <- deldir::deldir(x = as.numeric(coordinates$xcoord),
                                y = as.numeric(coordinates$ycoord))

  #----------------------------------------------------------------------------#
  # Extracting Voronoi/dirichlet tessaltion - not intested in Delauney
  #----------------------------------------------------------------------------#

  ## Voronoi tessaltion coordinates
  tessV <- tesselation$dirsgs
  ## If the "point" is on the edge of the box used for tessaltion
  tessV <- tessV[!tessV$bp1 & !tessV$bp2,]
  #tessV <- .filterTesselation(tessV)
  ## Contains original data - just to ensure that that indeces line up
  ogCoord <- cbind(coordinates,tesselation$summary)



  #----------------------------------------------------------------------------#
  # This section is equivalent to triangle rasterisation
  #----------------------------------------------------------------------------#
  .raster(verbose)
  allIdx <- parallel::mclapply(seq_len(nrow(tessV)),Vesalius:::.fillTesselation,
    tess = tessV,coord= ogCoord,mc.cores = cores)



  #----------------------------------------------------------------------------#
  # Filtering triangle that exceed max tile size
  #----------------------------------------------------------------------------#

  .filter(verbose)
  allIdx <- .filterTriangle(allIdx,filterThreshold)
  #----------------------------------------------------------------------------#
  # Rebuilding new cooirdinates based on tesselation borders
  # cooridinates are shifted for both tesselation cooridantes and starting coord
  #----------------------------------------------------------------------------#
  lpx <- round(min(allIdx$x))
  lpy <- round(min(allIdx$y))


  allIdx <- mutate(allIdx, x = x -lpx +1,y=y-lpy +1) %>%
            tibble

  return(allIdx)

}



#.resShift <- function(coord){

#}


.resShift <- function(coord, xbox , ybox,cores){
  #----------------------------------------------------------------------------#
  # Building template Sparse Matrix for each pixel
  #----------------------------------------------------------------------------#

  #rgbMat <- Matrix(0,ncol=3,nrow = length(xbox) * length(ybox))

  #----------------------------------------------------------------------------#
  # Subdivision index vector
  # For all possible subdivisions extract all combinations
  # This will hopefully speed things up
  #----------------------------------------------------------------------------#

  xs <- rep(xbox[seq(1,length(xbox)-1)],times = length(ybox)-1)
  xe <- rep(xbox[seq(2,length(xbox))], times = length(ybox)-1)


  ys <- rep(ybox[seq(1,length(ybox)-1)], each = length(xbox)-1)
  ye <- rep(ybox[seq(2,length(ybox))], each = length(xbox)-1)

  idx <- rep(seq(1,length(xbox)-1), times = length(ybox) -1)
  idy <- rep(seq(1,length(ybox)-1), each = length(xbox) -1)


  #----------------------------------------------------------------------------#
  # Looping over idx
  #----------------------------------------------------------------------------#

  rgbMat <- parallel::mcmapply(.pixelate,xs,xe,ys,ye,idx,idy,MoreArgs = list(coord), mc.cores = cores)


  nulls <- sapply(rgbMat,is.null)

  rgbMat <- rgbMat[!nulls]

  rgbMat <- do.call("rbind",rgbMat)
  colnames(rgbMat) <- c("barcodes","xcoord","ycoord","R","G","B")

  return(rgbMat)

}


.pixelate <- function(xs,xe,ys,ye,idx,idy,coord){
    #--------------------------------------------------------------------------#
    # Select code in each box
    #--------------------------------------------------------------------------#

    tmp <- filter(coord,xcoord >= xs & xcoord < xe & ycoord >= ys & ycoord < ye)
    barcodes <- paste0(tmp$barcodes,sep = "", collapse = "_")
    tmp <- select(tmp,c("R","G","B"))

    #--------------------------------------------------------------------------#
    # If nothing in interval box
    #--------------------------------------------------------------------------#
    if(nrow(tmp)<1){
        return(NULL)
    }

    #--------------------------------------------------------------------------#
    # Median colour for interval
    #--------------------------------------------------------------------------#
    tmp <- data.frame(barcodes,idx,idy,matrix(apply(tmp,2,median),ncol = 3))
    colnames(tmp) <- c("barcodes","xcoord","ycoord","R","G","B")
    return(tmp)
}

.findTileCenters <- function(x,y,ogCoord){
    pos <- paste0(x,"_",y)
    posCent <- paste0(round(ogCoord$x),"_",round(ogCoord$y))

    centersLoc <- match(posCent,pos)
    centers <- rep(FALSE, length(pos))
    centers[centersLoc[!is.na(centersLoc)]] <- TRUE
    return(centers)

}


.filterTriangle <- function(allIdx, filterThreshold){
    #--------------------------------------------------------------------------#
    # First get area and extract distribution
    #--------------------------------------------------------------------------#
    triArea1 <- sapply(allIdx,function(x){
                      return(x[[2]][1])
    })
    triArea2 <- sapply(allIdx,function(x){
                      return(x[[2]][2])
    })

    areaThresh <- quantile(c(triArea1,triArea2),filterThreshold)

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

.fillTesselation <- function(idx,tess,coord){
    #--------------------------------------------------------------------------#
    ## Get data
    #--------------------------------------------------------------------------#
    allEdge <- tess[idx,]

    #--------------------------------------------------------------------------#
    ## First split data between both triangles
    #--------------------------------------------------------------------------#
    t1 <- .fillTesselationTriangle(allEdge,coord,tri = 1)
    a1 <- nrow(t1)

    t2 <- .fillTesselationTriangle(allEdge,coord,tri = 2)
    a2 <- nrow(t2)


    #--------------------------------------------------------------------------#
    # Rebuild list with filled triangle for both points
    #--------------------------------------------------------------------------#
    allIn <- list(t1,t2)
    return(list(allIn,c(a1,a2)))

}

.fillTesselationTriangle <- function(allEdge,coord,tri){
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
    x <- c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"])
    y <- c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"])

    #--------------------------------------------------------------------------#
    # Fill triangles with all point in that space
    #--------------------------------------------------------------------------#
    cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
    maxX <- round(maxPolygonX[cell %in% c(1,2,3)])
    maxY <- round(maxPolygonY[cell %in% c(1,2,3)])
    value <- c(rep(coord[allEdge$ind1,"R"],length(maxX)),
               rep(coord[allEdge$ind1,"G"],length(maxX)),
               rep(coord[allEdge$ind1,"B"],length(maxX)))

    cc <- rep(1:3,each = length(maxX))
    maxX <- rep(maxX, times = 3)
    maxY <- rep(maxY, times = 3)
    barcode <- rep(coord[allEdge$ind1,"barcodes"],length(maxX))
    cent <- which(maxX == round(coord[allEdge$ind1,"x"]) & maxY == round(coord[allEdge$ind1,"y"]))
    centers <- rep(0, length(maxX))
    centers[cent] <- 1
    allIn <- data.frame(barcode,maxX,maxY,cc,value,centers)
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
    x <- c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"])
    y <- c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"])

    #--------------------------------------------------------------------------#
    # Fill triangles with all point in that space
    #--------------------------------------------------------------------------#
    cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
    maxX <- round(maxPolygonX[cell %in% c(1,2,3)])
    maxY <- round(maxPolygonY[cell %in% c(1,2,3)])
    value <- c(rep(coord[allEdge$ind2,"R"],length(maxX)),
               rep(coord[allEdge$ind2,"G"],length(maxX)),
               rep(coord[allEdge$ind2,"B"],length(maxX)))
    cc <- rep(1:3,each = length(maxX))
    maxX <- rep(maxX, times = 3)
    maxY <- rep(maxY, times = 3)
    barcode <- rep(coord[allEdge$ind2,"barcodes"],length(maxX))
    cent <- which(maxX == round(coord[allEdge$ind2,"x"]) & maxY == round(coord[allEdge$ind2,"y"]))
    centers <- rep(0, length(maxX))
    centers[cent] <- 1
    allIn <- data.frame(barcode,maxX,maxY,cc,value,centers)
    colnames(allIn) <- c("barcodes","x","y","cc","value","tile")
    return(allIn)
  }

}


.convexify <- function(x,y){
    ### fixing this
    ### Still not creating proper polygons for some reason

    startX <- which(x==min(x))
    if(length(startX)>1){
      start <- startX[y[startX] == min(y[startX])]
    } else {
      start <- startX
    }

    xst <- x[start]
    yst <- y[start]

    xrest <- x[!seq_along(x) %in% start]
    yrest <- y[!seq_along(y) %in% start]

    top <- yrest >=yst
    bottom <- yrest <yst

    topx <- xrest[top]
    topy <- yrest[top]

    bottomx <- xrest[bottom]
    bottomy <- yrest[bottom]

    topfx <- topx[order(topx,decreasing = F)]
    topfy <- topy[order(topx,decreasing = F)]
    bottomfx <- bottomx[order(bottomx,decreasing = T)]
    bottomfy <- bottomy[order(bottomx,decreasing = T)]

    x <- c(xst,topfx,bottomfx)
    y<- c(yst,topfy,bottomfy)
    convex <- cbind(x,y)
    return(convex)

}
