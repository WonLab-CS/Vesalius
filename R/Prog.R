################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Progress Message/--------------------------------#

# Section exclusively dedicated to progress message and pretty print
# I really don't want to include this in the main functions
# That is just messy

.assignCoordinates <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Assigning Colours to Coordinates \n"))

}

.typeChange <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Class coersion to numeric \n"))

}

.invertCols <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Inverting Colours \n"))

}

.res <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Coarse Grain Image Building \n"))

}


.distanceBeads <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Filtering outlier beads \n"))

}

.tess <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Generating Voronoi Tesselation \n"))

}

.raster <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Rasterising Tiles \n"))


}

.filter <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")
    cat( paste(t," Filtering Triangles that exceed area threshold\n"))
    cat(bar,"\n")

}


.simpleBar <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")

}

.degAllProg <- function(ter,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()

    cat( paste(t," Computing DEG for territory",ter,"against all\n"))


}

.degEachProg <- function(seed,query,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()

    cat( paste(t," Computing DEG for territory",seed,"against territory", query,"\r"))

}

.seg <- function(seg,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t,"Segmentation Iteration ", seg,"\r"))
}

.smooth <- function(smooth,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t,"Smoothing Iteration ", smooth,"\r"))
}

.terPool <- function(ter,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t,"Pooling territory ", ter,"\r"))
}
