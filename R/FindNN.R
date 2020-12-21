################################################################################
###############################   Vesalius      ################################
################################################################################


#/--------------------------/ nearest Neighbors /------------------------------/

#' Finding Nearest spatial neighbors
#'
#' @param slide a slide containing the barcodes and associated x/y coordinates.
#' @param k the number of nearest neighbors to extract
#' @param box the size of the box used to find nearest neighbors.
#'     Default is set at 0.05. This represents the proportion of total distance.
#' @param sectors the number of angle sectors used to find nearest neighbors.
#'    Default is set at four.
#' @param cores number of cores used to find nearest neighbors
#'
#' @details In order to ensure that the nearest neighbors are distributed in all
#'    directions, this function find all neighbors within a certain distance (box) and
#'    categorises all points based on their angle compared to the center point. These points are
#'    then sorted. The resulting points are the points that are both the closest
#'    but also found in all directions compared to a center point. If direction is not
#'    considered important, then set sectors to 1.
#'
#' @return A list of nearest neighbors. Each element in the list represents a single
#'     barcode used as center point. Each elements in the list is a data frame with neighbor barcodes
#'     x and y coordninates of said neighbors as well as distance and angle from center point.


findNearestNeighbors <- function(slide, k = 5 ,box = 0.05,sectors =4,cores=1){

    ## cell id locations
    cellLoc <- slide[!duplicated(slide$barcodes),c("barcodes","xcoord","ycoord")]
    xcut <- (max(slide$xcoord)-min(slide$xcoord))*box
    ycut <- (max(slide$ycoord)-min(slide$ycoord))*box
    box <- c(xcut,ycut)
    ## lets try local expansion

    localcoor <- mclapply(seq_len(nrow(cellLoc)), .localCoordinates,
                       cellLoc = cellLoc,
                       k=k,totalCells = nrow(cellLoc),
                       box = box,
                       sectors =sectors,
                       getNeighbors = TRUE,
                       mc.cores=cores)
    names(localcoor) <- cellLoc$barcodes
    cat("\n")
    return(localcoor)


}


### Internal function to FindNearestNeighbors - see above
### Function is not exported
### idx - barcodes index
.localCoordinates <- function(idx,cellLoc, k,totalCells, box, sectors = 8,getNeighbors = FALSE){

    ## setting boudaries from coordinates
    xstart <- as.numeric(cellLoc[idx,"xcoord"])
    ystart <- as.numeric(cellLoc[idx,"ycoord"])
    xcor <- c((xstart + (box[1])),(xstart - (box[1])))
    ycor <- c((ystart + (box[2])),(ystart - (box[2])))
    origin <- c(xstart,ystart)
    ## expanding from source
    ## a need to expand from each quadrant

    quads <- vector("list", sectors)
    tops <- ceiling(k/sectors)
    if(tops <1) tops <-1
    sectors <- seq(0,360,length.out = sectors+1)
    cellLoc <- cellLoc[(as.numeric(cellLoc$xcoord) >= min(xcor)) &
                              (as.numeric(cellLoc$xcoord) <= max(xcor)) &
                              (as.numeric(cellLoc$ycoord) >= min(ycor)) &
                              (as.numeric(cellLoc$ycoord) <= max(ycor)),]

    ## transform x;y coordinates into polar coordinates using start cell as origin
    polarTransform <- lapply(seq(1,nrow(cellLoc)),.polarConversion.bead,
                             points = cellLoc, origin =origin)

    polarTransform <- as.data.frame(do.call("rbind",polarTransform))
    polarTransform$distance <- as.numeric(as.character(polarTransform$distance))
    polarTransform$angle <- as.numeric(as.character(polarTransform$angle))
    ## removing nan
    polarTransform <- polarTransform[!polarTransform$distance ==0,]
    polarTransform <- polarTransform[!is.na(polarTransform$distance),]


    angle <- polarTransform[,"angle"]
    for(qu in seq(1,length(sectors)-1)){
        tmp <- polarTransform[which(angle >= sectors[qu] &
                                           angle<sectors[qu+1]),]
        tmp <- tmp[order(tmp[,"distance"]),]
        if(nrow(tmp) == 0) next()
        tmp$sector <- qu
        quads[[qu]] <- tmp
    }

    ### find nearest neighbors
    if(getNeighbors){
        cellCoordinates <- lapply(quads, function(x,t){
                                  if(length(x)!=0){
                                      return(x[seq(1,t),])
                                  } else {
                                    return(x)
                                  }
                                  },tops)
        cellCoordinates <- do.call("rbind", cellCoordinates)
        distance <- cellCoordinates$distance[!is.na(cellCoordinates$distance)]
        angle <- cellCoordinates$angle[!is.na(cellCoordinates$angle)]

        cellCoordinates <- cellLoc[cellLoc$barcodes %in% cellCoordinates$barcodes,]

        cellCoordinates$distance <- distance
        cellCoordinates$angle <- angle
    } else {
        cellCoordinates <- do.call("rbind",quads)
        distance <- cellCoordinates$distance
        angle <- cellCoordinates$angle

        cellCoordinates <- cellLoc[cellLoc$barcodes %in% cellCoordinates$barcodes,]
        cellCoordinates$distance <- distance
        cellCoordinates$angle <- angle
    }

    ## progress



    return(cellCoordinates)
}

.polarConversion.bead <- function(idx,points,origin){
    ## let's consider that points is going to be a two column matrix or df
    local <- as.vector(as.matrix(points[idx,c("xcoord","ycoord")]))

    ## distance from origin
    xo <- origin[1]
    yo <- origin[2]
    xp <- local[1]
    yp <- local[2]
    distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)
    ## quadrant test
    x <- xp-xo ; y <- yp-yo
    if(x >= 0 & y >= 0) angle <- atan(abs(yp-yo)/abs(xp-xo))*(180/pi)
    if(x < 0 & y >= 0) angle <- 180 - (atan(abs(yp-yo)/abs(xp-xo))*(180/pi))
    if(x < 0 & y < 0) angle <- 180 + (atan(abs(yp-yo)/abs(xp-xo))*(180/pi))
    if(x >= 0 & y < 0) angle <- 360 - (atan(abs(yp-yo)/abs(xp-xo))*(180/pi))

    if(is.null(points$cluster)){
      polar <- c(points[idx,"barcodes"],distance,angle)
      names(polar) <- c("barcodes","distance","angle")
    } else {
      polar <- c(points[idx,c("barcodes","cluster")],distance,angle)
      names(polar) <- c("barcodes","cluster","distance","angle")
    }

    return(polar)
}

.polarConversion.array <- function(points,origin){
    #--------------------------------------------------------------------------#
    # generating all neigbors and angles from center
    #--------------------------------------------------------------------------#

    xo <- origin[1]
    yo <- origin[2]
    xp <- points[,"x"]
    yp <- points[,"y"]

    #--------------------------------------------------------------------------#
    # Distance from center point
    #--------------------------------------------------------------------------#
    distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)

    #--------------------------------------------------------------------------#
    # Angle from center point
    #--------------------------------------------------------------------------#

    x <- xp-xo ; y <- yp-yo
    angle <- mapply(function(x,y){
      if(x >= 0 & y >= 0) angle <- atan(abs(y)/abs(x))*(180/pi)
      if(x < 0 & y >= 0) angle <- 180 - (atan(abs(y)/abs(x))*(180/pi))
      if(x < 0 & y < 0) angle <- 180 + (atan(abs(y)/abs(x))*(180/pi))
      if(x >= 0 & y < 0) angle <- 360 - (atan(abs(y)/abs(x))*(180/pi))
        return(angle)
    },x=x,y=y, SIMPLIFY = TRUE)
    
    polar <- data.frame("distance" = as.numeric(distance),"angle" = as.numeric(angle))

    return(polar)
}
