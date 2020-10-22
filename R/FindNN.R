################################################################################
############################   ST package        ###############################
################################################################################


#/-------------------------/ nearest Neighbours /------------------------------/


findNearestNeighborsTrial <- function(slide, k = 5 ,sectors =4,cores=1){
    # computing polar cooddinates for all points
    idx <- seq_len(nrow(slide))
    message("Building Polar Matrix \n")
    polarMatrix <- mclapply(idx, function(idx,mat,sectors){
                      ## priming distances
                      cat(paste0(round(idx/nrow(mat),2),"% Completed"),"\r")
                      xo <- mat$xcoord[idx]
                      yo <- mat$ycoord[idx]
                      xp <- mat$xcoord
                      yp <- mat$ycoord
                      # computing euclidean distance between all points
                      # using point "idx" as center
                      distance <- sqrt((abs(xp-xo))^2 + (abs(yp-yo))^2)
                      names(distance) <- mat$barcodes
                      # Defining x and y for angle calculation
                      x <- xp-xo ; y <- yp-yo
                      ## computing all angle from center
                      angle <- mapply(function(x,y){
                        if(x >= 0 & y >= 0) angle <- atan(abs(yp-yo)/abs(xp-xo))*(180/pi)
                        if(x < 0 & y >= 0) angle <- 180 - (atan(abs(yp-yo)/abs(xp-xo))*(180/pi))
                        if(x < 0 & y < 0) angle <- 180 + (atan(abs(yp-yo)/abs(xp-xo))*(180/pi))
                        if(x >= 0 & y < 0) angle <- 360 - (atan(abs(yp-yo)/abs(xp-xo))*(180/pi))
                        return(angle)
                      },x,y)
                      names(angle) <- mat$barcodes
                      ## Assigning angle to sector
                      sectors <- seq(0,360,length.out = sectors+1)
                      for(sect in seq(1,length(sectors)-1)){
                          angle[which(angle >= sectors[sect] &
                                      angle < sectors[sect+1])] <- sect
                      }

                      return(list("dist" = distance, "angle" = angle))
                    }, mat = slide,
                       sectors = sectors,
                       mc.cores = cores)

    ## Extracting Distances and sectors
    message("Building Distance Matrix \n")
    distanceMatrix <- mclapply(polarMatrix, function(x)return(x$dist),
                               mc.cores=cores)
    cat("Building Sector Matrix \n")
    sectorMatrix <- mclapply(polarMatrix, function(x)return(x$angle),
                             mc.cores=cores)

    ## Split distance by sectors
    distanceMatrix <- mcmapply(function(dist,sector){

                              return(split(dist,sector))
                            },dist = distanceMatrix,
                              sector = sectorMatrix,mc.cores = cores)

    ## Selecting nearest neighbors in each sector
    nn <- ceiling(k/sectors)
    if(nn <1) nn <-1
    cat("Extracting Nearest Neighbours \n")
    nearest <- mclapply(distanceMatrix, function(dist,nn){
                        sectors <- names(nn)
                        nn <- lapply(dist,function(x,nn,sector){
                                     tmp <- order(x,decreasing = FALSE)[seq(1,nn)]
                                     sector <- rep(sector, length(tmp))
                                     nn <- data.frame(names(tmp),tmp,sector)
                                     return(nn)
                        },nn = nn)
                        nn <- do.call("rbind",nn)
                        return(nn)
    })

    return(nearest)

}



findNearestNeighbors <- function(slide, k = 5 ,box = 0.05,sectors =4,cores=1){

    ## cell id locations
    #cellLoc <- slide[!duplicated(slide$cellID),c("cellID","x","y")]
    cellLoc <- slide[!duplicated(slide$barcodes),c("barcodes","xcoord","ycoord")]
    xcut <- (max(slide$xcoord)-min(slide$xcoord))*box
    ycut <- (max(slide$ycoord)-min(slide$ycoord))*box
    box <- c(xcut,ycut)
    ## lets try local expansion

    localcoor <- mclapply(seq_len(nrow(cellLoc)), localCoordinates,
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




## Under the assumption that we will be using apply
## will go row by row with this instead
localCoordinates <- function(idx,cellLoc, k,totalCells, box, sectors = 8,getNeighbors = FALSE){
    ##

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
    polarTransform <- lapply(seq(1,nrow(cellLoc)),.polarConversion,
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
    cat(paste(round(idx/,, "out of all Cells","\r"))


    return(cellCoordinates)
}

.polarConversion <- function(idx,points,origin){
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
