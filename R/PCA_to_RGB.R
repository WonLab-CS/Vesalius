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




#' Assigning the RGB colour code to coordinates in lower resolution.
#'
#' @param rgb list of colour codes in the format describe in the Detail section.  Output of rgbPCA.
#' @param coordinates data frame containing x and y coordinates of barcodes
#' @param resolution size of "image matrix" to be used. n by n matrix where n = resolution
#' @param drop logical describing if empty location should be removed.
#' @param na.rm logical describing if NA should be removed
#' @details assignRGBtoPixel remaps the colour code to x/y coordinates at lower resolution.
#' This means that beads/spots will be assigned to a "pixel" within a matrix then exported as a data frame.
#' The current version does not return assigned colours as array images.
#' \code{rgb} is the output of \code{rgbPCA}.
#' As such, it should be a list with three different levels.
#' Level 1 : each list element is a different slice
#' Level 2 : each list element is R G B channel for that slice
#' Level 3 : numeric vectors containing colour code for each barcode
#' @return a data frame with five columns: barcodes (barcode names), xcoord (column index),
#' ycoord (row index), R (Red colour channel),G (Green colour channel) ,and B (Blue Colour channel).


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



### temp just to be able to load intermediate objects
### This will need to be changed in the final version
buildImageArray <- function(coordinates,rgb=NULL,invert=FALSE,na.rm = TRUE,cores=1){
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
  for(cols in seq(2,ncol(coordinates))){
    coordinates[,cols] <- as.numeric(as.character(coordinates[,cols]))
  }



  #----------------------------------------------------------------------------#
  # Inverting colours
  #----------------------------------------------------------------------------#
  if(invert){
      coordinates$R <- 1 - coordinates$R
      coordinates$G <- 1 - coordinates$G
      coordinates$B <- 1 - coordinates$B
  }

  #----------------------------------------------------------------------------#
  # Removing all "outer point" - tend to make tesselation messy
  #----------------------------------------------------------------------------#
  idx <- seq_len(nrow(coordinates))
  minDist <- sapply(idx, function(idx,mat){
                      xo <- mat$xcoord[idx]
                      yo <- mat$ycoord[idx]
                      xp <- mat$xcoord
                      yp <- mat$ycoord
                      distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)
                      distance <- distance[distance !=0]

                      return(min(distance))
  }, coordinates)

  distanceThrehsold <- quantile(minDist,0.99)
  coordinates <- coordinates[minDist <= distanceThrehsold,]

  #----------------------------------------------------------------------------#
  # TESSELATION TIME!
  #----------------------------------------------------------------------------#
  tesselation <- deldir::deldir(x = as.numeric(coordinates$xcoord),
                                y = as.numeric(coordinates$ycoord))

  #----------------------------------------------------------------------------#
  # Extracting Voronoi/dirichlet tessaltion - not intested in Delauney
  #----------------------------------------------------------------------------#
  #browser()
  ## Voronoi tessaltion coordinates
  tessV <- tesselation$dirsgs
  ## If the "point" is on the edge of the box used for tessaltion
  tessV <- tessV[!tessV$bp1 & !tessV$bp2,]
  #tessV <- .filterTesselation(tessV)
  ## Contains original data - just to ensure that that indeces line up
  ogCoord <- data.frame(coordinates,tesselation$summary)


  #----------------------------------------------------------------------------#
  # Rebuilding new cooirdinates based on tesselation borders
  # cooridinates are shifted for both tesselation cooridantes and starting coord
  #----------------------------------------------------------------------------#
  lpx <- min(c(tessV$x1,tessV$x2))
  hpx <- max(c(tessV$x1,tessV$x2))
  lpy <- min(c(tessV$y1,tessV$y2))
  hpy <- max(c(tessV$y1,tessV$y2))


  #----------------------------------------------------------------------------#
  # This section is equivalent to triangle rasterisation
  #----------------------------------------------------------------------------#
  allIdx <- parallel::mclapply(seq_len(nrow(tessV)),Vesalius:::.fillTesselation,
    tess = tessV,coord= ogCoord,mc.cores = cores)


  #----------------------------------------------------------------------------#
  # Fill in the gaps with the coulours
  #----------------------------------------------------------------------------#

  img <- array(1,dim = c(hpx - lpx +1,hpy -lpy +1,1,3))

  for(i in seq_along(allIdx)){

      cat((i/length(allIdx))*100 , "% \r")
      #------------------------------------------------------------------------#
      # First fill in Triangle 1
      #------------------------------------------------------------------------#
      x1 <- allIdx[[i]][["t1"]][["fill"]]$x - lpx
      y1 <- allIdx[[i]][["t1"]][["fill"]]$y - lpy
      coord <- allIdx[[i]][["t1"]][["ind"]]

      img[x1,y1,1L,1L] <- ogCoord[coord,"R"]
      img[x1,y1,1L,2L] <- ogCoord[coord,"G"]
      img[x1,y1,1L,3L] <- ogCoord[coord,"B"]
      #------------------------------------------------------------------------#
      # Next fill in Triangle 2
      #------------------------------------------------------------------------#
      x2 <- allIdx[[i]][["t2"]][["fill"]]$x - lpx
      y2 <- allIdx[[i]][["t2"]][["fill"]]$y - lpy
      coord <- allIdx[[i]][["t2"]][["ind"]]

      img[x2,y2,1L,1L] <- ogCoord[coord,"R"]
      img[x2,y2,1L,2L] <- ogCoord[coord,"G"]
      img[x2,y2,1L,3L] <- ogCoord[coord,"B"]
  }

  return(img)

}


.fillTesselation <- function(idx,tess,coord){
    #--------------------------------------------------------------------------#
    ## Get data
    #--------------------------------------------------------------------------#
    allEdge <- tess[idx,]

    #--------------------------------------------------------------------------#
    ## First split data between both triangles
    #--------------------------------------------------------------------------#
    t1 <- list("fill" = .fillTesselationTriangle(allEdge,coord,tri = 1), "ind" = allEdge$ind1)
    t2 <- list("fill" = .fillTesselationTriangle(allEdge,coord,tri = 2), "ind" = allEdge$ind2)

    #--------------------------------------------------------------------------#
    # Rebuild list with filled triangle for both points
    #--------------------------------------------------------------------------#
    allIn <- list("t1" = t1 , "t2"=t2)
    return(allIn)

}

.fillTesselationTriangle <- function(allEdge,coord,tri){
  if(tri ==1){
    #------------------------------------------------------------------------#
    # Boundaries of tesselation
    #------------------------------------------------------------------------#
    lpx <- round(min(c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"]))) - 1
    hpx <- round(max(c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"]))) + 1
    lpy <- round(min(c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"]))) - 1
    hpy <- round(max(c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"]))) + 1

    maxPolygonX <- rep(seq(lpx,hpx), times = hpy - lpy +1)
    maxPolygonY <- rep(seq(lpy,hpy), each = hpx - lpx +1)

    #------------------------------------------------------------------------#
    # Creating triangles
    #------------------------------------------------------------------------#
    x <- c(allEdge$x1,allEdge$x2,coord[allEdge$ind1,"x"])
    y <- c(allEdge$y1,allEdge$y2,coord[allEdge$ind1,"y"])

    #------------------------------------------------------------------------#
    # Fill triangles with all point in that space
    #------------------------------------------------------------------------#
    cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
    maxX <- round(maxPolygonX[cell %in% c(1,2,3)])
    maxY <- round(maxPolygonY[cell %in% c(1,2,3)])

    allIn <- cbind(maxX,maxY)
    colnames(allIn) <- c("x","y")
    return(allIn)

  } else {
    #------------------------------------------------------------------------#
    # Boundaries of tesselation
    #------------------------------------------------------------------------#
    lpx <- round(min(c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"]))) - 1
    hpx <- round(max(c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"]))) + 1
    lpy <- round(min(c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"]))) - 1
    hpy <- round(max(c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"]))) + 1

    maxPolygonX <- rep(seq(lpx,hpx), times = hpy - lpy +1)
    maxPolygonY <- rep(seq(lpy,hpy), each = hpx - lpx +1)

    #------------------------------------------------------------------------#
    # Creating triangles
    #------------------------------------------------------------------------#
    x <- c(allEdge$x1,allEdge$x2,coord[allEdge$ind2,"x"])
    y <- c(allEdge$y1,allEdge$y2,coord[allEdge$ind2,"y"])

    #------------------------------------------------------------------------#
    # Fill triangles with all point in that space
    #------------------------------------------------------------------------#
    cell <- point.in.polygon(maxPolygonX,maxPolygonY,x,y)
    maxX <- round(maxPolygonX[cell %in% c(1,2,3)])
    maxY <- round(maxPolygonY[cell %in% c(1,2,3)])

    allIn <- cbind(maxX,maxY)
    colnames(allIn) <- c("x","y")
    return(allIn)
  }

}



.generateMaxDistancePolygon <- function(centerCoordinates,ogCoord,resolution=360){
    #--------------------------------------------------------------------------#
    # Making a reference polygon from bead locations
    # The idea is to trim anything that is not in that polygon
    # Removing all the outer points
    # There will be a loss of some points for sure
    #--------------------------------------------------------------------------#
    allPoints <- .polarConversion.array(ogCoord,centerCoordinates)

    #--------------------------------------------------------------------------#
    # creating sectors = resoltuion
    # Each sector contains point that is the furtherst from center
    #--------------------------------------------------------------------------#
    sectors <- seq(0,360,length.out = resolution+1)
    angle <- allPoints$angle
    quads <- rep(0,length(sectors))

    for(qu in seq(1,length(sectors)-1)){
        tmp <- allPoints$distance[angle >= sectors[qu] &
                                  angle<sectors[qu+1]]
        tmp <- sort(tmp,decreasing =TRUE)

        quads[qu] <- tmp[1L]

    }
    #--------------------------------------------------------------------------#
    # quantiled distance - used a threshold value
    # all pixel that are futher away from center point that this will be removed
    #--------------------------------------------------------------------------#

    return(quantile(quads,0.95))
}



.filterTesselation <- function(tess){
  #----------------------------------------------------------------------------#
  # Finding thge best way to filter out un wanted points
  #----------------------------------------------------------------------------#

  edgeConnected <- tess[tess$bp1 | tess$bp2,c("ind1","ind2")]
  edgeConnected <- unique(as.vector(edgeConnected))
  tess <- tess[tess$ind1 %in% edgeConnected | tess$ind2 %in% edgeConnected]

  return(tess[numTess,])
}



.trimOuter <- function(all,borders){
    #--------------------------------------------------------------------------#
    # For each triangle find which points overlap
    # For many sectors it will return 0 but that is okay as we can just collapse
    # everything and zeroes will just be discarded
    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Finding points in each sector - trimming outer
    #--------------------------------------------------------------------------#

    p <- point.in.polygon(all[,"maxPolygonX"],
                            all[,"maxPolygonY"],
                            borders[,1L],
                            borders[,2L])
    keep <- p %in% c(1,2,3)
    trimmed <- all[keep,]
    return(trimmed)
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
