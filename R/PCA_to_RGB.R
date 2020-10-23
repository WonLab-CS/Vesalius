################################################################################
###############################   Vesalius      ################################
################################################################################

#---------------------/PCA RGB Functions/--------------------------------------#

## Running pca with rgb conversion
## This function uses the seurat PCA functions
# slide = gene count matrix (digital expression)
# SO = Seurat Object - In order to run PCA and find variable features using the
# Seurat package ,a Seurat Object need to be generated. The generation of such
# object is describe in the analysis pipeline below.
# slice = number of PCA slice to compute
#( e.g: 1 slice computes PC1 to PC3, 2 slices computes PC1 to PC6 -) )
# trim = quantile limits for colour histogram trimming (Not used)

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
                    x <- x [cutoff]
                    x <- (x - min(x)) / (max(x) - min(x))
                    return(x)
      }, cutoff = cutoff)

      ##
      image_slice[[sl]] <- rgb
    }




    ## Return image slice list
    ## Level 1 : each list element is a different slice
    ## Level 2 : each list element is R G B channel for that slice
    ## Level 3 : numeric vectors containing colour code for each barcode
    return(image_slice)
}


## Assigning the RGB colour code to its specific barcode location
## rgb = output of rgbPCA_Seurat
## image slice list
## Level 1 : each list element is a different slice
## Level 2 : each list element is R G B channel for that slice
## Level 3 : numeric vectors containing colour code for each barcode

## coordinates = data frame of barcode coordinates (bead_location.csv)
## drop = remove NA rgb colour codes

assignRGBtoPixel <- function(rgb,coordinates, drop = TRUE){
    # for each slice find coresponding barcodes
    for(i in seq_along(rgb)){
        code <- rep(NA,nrow(coordinates))
        locs <- match(coordinates$barcodes,names(rgb[[i]]))
        code[!is.na(locs)] <- rgb[[i]][locs[!is.na(locs)]]
        coordinates <- data.frame(coordinates,code)

    }
    colnames(coordinates) <- c("barcodes","xcoord","ycoord","R","G","B")
    if(drop){
        coordinates <- coordinates[!is.na(coordinates$R),]
    }
    # Data frame with barcode name, x/y coordinates and associated RGB colour code
    return(coordinates)
}




## Generating a pixel matrix with colour codes
## rgb = output of rgbPCA_Seurat
## image slice list
## Level 1 : each list element is a different slice
## Level 2 : each list element is R G B channel for that slice
## Level 3 : numeric vectors containing colour code for each barcode

## coordinates = data frame of barcode coordinates (bead_location.csv)
## Resolution = size of matrix to generate
## drop = remove pixel that do not any colours( this is just to remove white pixel)
## the idea is to reduce the size of the data frame that will be returned
## A data frame is returned in order to keep track of barcode combinations that are
## associated to pixels. Keeping track of barcodes and their location is crucial!
## na.rm = removing NAs if needed
## NOTE : This function is memory greedy in R (gc())

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
