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

rgbPCA<- function(slide,SO,slices = 1,countWeight = FALSE, conserveSparse = TRUE, trim = 0){

    # Converting back to a matrix (faster but uses more memory)
    if(!conserveSparse){
        slide <- as.matrix(slide)
    }
    # Seurat PCA
    pca <- RunPCA(SO,npcs = slices*3)

    image_slice <- vector("list",slices)
    slice_start <- seq(1,slices*3,by=3)

    for(sl in seq_len(slices)){
      # Extract Loadings and absolute value of loadings
      loadings <-Loadings(pca[["pca"]])[,seq(slice_start[sl], slice_start[sl]+2)]
      loadings <- apply(loadings,2,function(x)return(abs(x)))

      # RGB conversion
      rgb <- vector("list", 3)
      names(rgb) <- c("R","G","B")

      for(i in seq_along(rgb)){
          rgb[[i]] <-rep(0,ncol(slide))
          names(rgb[[i]]) <- colnames(slide)
          for(j in seq_len(ncol(slide))){

             genes <- rownames(slide)[slide[,j]!=0]
             # Extract loadings of genes for barcode j
             cs <- loadings[,i][names(loadings[,i]) %in% genes]

             # Just in case a barcode does not contain genes present in the
             # loadings
             # This could happen if the selection of variable
             # features is too stringent
             if(length(cs) ==0){
               rgb[[i]][j] <- 0
             }else {
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
      ## trimming
      cutoff <- do.call("cbind",lapply(rgb,trimHistogram,trim))
      cutoff <- apply(cutoff,1,sum) ==3


      ## Normalising RGB channels
      rgb <- lapply(rgb, function(x,cutoff){
                    x <- x[cutoff]
                    x <- (x - min(x)) / (max(x) - min(x))
                    return(x)
      }, cutoff = cutoff)

      ##
      image_slice[[sl]] <- rgb
    }




    ## Return image slice list

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
