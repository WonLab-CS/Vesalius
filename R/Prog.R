################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Progress Message/--------------------------------#

# Section exclusively dedicated to progress message and pretty print
# I really don't want to include this in the main functions
# That is just messy

.consSparse <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Converting Sparse Matrix to Matrix \n"))
}
.pca <- function(verbose =TRUE,slices =1){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Running Principle Component Analysis in",slices,"slices", "\n"))
}

.pcaTensor<- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Running Principle Component Analysis \n"))
}
.rgb <- function(verbose =TRUE,slice){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Converting Loading Values to RGB in slice",slice, "\n"))
}

.pcaRGBTensor <- function(verbose =TRUE,pc){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Converting Loading Values to RGB in PC",pc, "\n"))
}

.embedRGBTensor <-function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Converting Embeddings to RGB \n"))
}

.umapRGBTensor <- function(verbose){
  if(!verbose) return(NULL)
  t <- Sys.time()
  cat( paste(t," Converting UMAP Embeddings to RGB \n"))
}
.adj <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Adjusting RGB values \n"))
}

.norm <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Normalising RGB values \n"))
}
.pcadj <- function(verbose =TRUE,slice){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Adujsting RGB values using total variance in each PC - slice",slice,"\n"))
}

.assignCoordinates <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Assigning Colours to Coordinates \n"))
}

.typeChange <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Class coersion to numeric \n"))
}

.invertCols <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Inverting Colours \n"))
}

.res <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Resizing Image \n"))
}

.reg <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Regularising Image \n"))
}


.distanceBeads <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Filtering outlier beads \n"))
}

.tess <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Generating Voronoi Tesselation \n"))
}

.raster <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Rasterising Tiles \n"))
}

.filter <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Filtering Triangles that exceed area threshold\n"))
}
.findCenter <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Finding Tile centers in image\n"))
}


.simpleBar <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    widthDisplay<-round(options()$width)
    bar <- paste0("#", paste0(rep("-",widthDisplay*0.9),collapse=""),"#")
    cat(bar,"\n")

}

.degProg <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Computing Differentially Expressed Genes\n"))
}

.degEachProg <- function(seed,query,verbose =TRUE){
    if(!verbose) return(NULL)
    cat( paste("===>",seed,"VS",query,"<===","\r"))
}


.extractTerProg <- function(seed,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Extracting territory",seed,"\n"))
}

.seg <- function(seg,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Segmentation Iteration ", seg,"\r"))
}

.smooth <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Smoothing Image Arrays \n"))
}

.terPool <- function(ter,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Pooling territory ", ter,"\r"))
}

.checkCounts <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Checking and extracting Counts \n"))
}

.checkCoord <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Checking Coordinates \n"))
}

#.checkVesalius <- function(verbose =TRUE){
#    if(!verbose) return(NULL)
#    t <- Sys.time()
#    cat( paste(t," Checking and Converting Vesalius Object \n"))
#}
.seedSelect <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Extracting Seed territories \n"))
}

.querySelect <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Extracting Query territories \n"))
}
.eq <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Equalizing Histogram \n"))
}
.rebuildDF <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Rebuilding Data Frame from image \n"))
}

.morph <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Converting to pixset and morphing territory \n"))
}

.layerTer <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Generating layers\n"))
}

.checkVes <- function(sliceID,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t," Checking Vesalius Input - Using slice ",sliceID, "\n"))
}

.tensorRes <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Reducing tensor resolution\n"))
}
.fTiles <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Filtering Tiles\n"))
}

.buildSO <- function(verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Pre-processing Count Data\n"))
}
.adjCounts <- function(verbose=TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Adjusting count matrix\n"))
}

.vtc <- function(verbose=TRUE){
  if(!verbose) return(NULL)
  t <- Sys.time()
  cat( paste(t," Converting Vesalius to Image\n"))
}

.msk <- function(em,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t," Applying Masks in dim",em, "\r"))
}
