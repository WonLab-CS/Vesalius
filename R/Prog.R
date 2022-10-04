################################################################################
###############################   Vesalius      ################################
################################################################################

#----------------------------/Progress Message/--------------------------------#

# Section exclusively dedicated to progress message and pretty print
# I really don't want to include this in the main functions
# That is just messy

#---------------------------/Embedding Messages/-------------------------------#

.conserve_sparse <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Converting Sparse Matrix to Matrix \n"))
}
.pca_tensor <- function(verbose =TRUE){
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Running Principle Component Analysis \n"))
}
.pca_rgb_tensor <- function(verbose = TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Converting PCA Embediing Values to RGB", "\n"))
}
.pcal_rgb_tensor <- function(verbose = TRUE, pc) {
    if ( !verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Converting Loading Values to RGB in PC", pc, "\n"))
}

.svd_tensor <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Running Single Value Decomposition \n"))
}
.svd_rgb_tensor <- function(verbose = TRUE){ 
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t, " Converting Embeddings to RGB \n"))
}

.umap_rgb_tensor <- function(verbose) {
  if(!verbose) return(NULL)
  t <- Sys.time()
  cat( paste(t, " Converting UMAP Embeddings to RGB \n"))
}

#---------------------/Image Creation Messages/---------------------------#
.distance_beads <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Filtering outlier beads \n"))
}

.tess <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Generating Voronoi Tesselation \n"))
}
.raster <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t, " Rasterising Tiles \n"))
}
.filter <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Filtering Triangles that exceed area threshold\n"))
}
.tensor_res <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Reducing tensor resolution\n"))
}
.f_tiles <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Filtering Tiles\n"))
}

.build_so <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Pre-processing Count Data\n"))
}
.adj_counts <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Adjusting count matrix\n"))
}

#-------------------/Image processing Messages/----------------------------#
.reg <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Regularising Image \n"))
}
.seg <- function(seg,verbose = TRUE) {
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t," Segmentation Iteration ", seg,"\r"))
}
.smooth <- function(verbose = TRUE) {
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t," Smoothing Image Arrays \n"))
}

.ter_pool <- function(ter,verbose =TRUE) {
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Pooling territory ", ter, "\r"))
}
.eq <- function(verbose =TRUE){
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Equalizing Histogram \n"))
}
.morph <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Converting to pixset and morphing territory \n"))
}
.layer_ter <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Generating layers\n"))
}


#-----------------------------/Prettify/-----------------------------------#
.simple_bar <- function(verbose = TRUE) {
    if(!verbose) return(NULL)
    widthDisplay <- round(options()$width)
    bar <- paste0("#",
        paste0(rep("-", widthDisplay * 0.9), collapse = ""),
        "#")
    cat(bar, "\n")

}


#----------------------------/DEG messages/-------------------------------#
.deg_prog <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Computing Differentially Expressed Genes\n"))
}

.deg_each_prog <- function(seed, query, verbose = TRUE) {
    if(!verbose) return(NULL)
    cat(paste("===>", seed, "VS", query, "<===", "\r"))
}

#------------------------/Object Sanity messages/---------------------------#

.check_coord <- function(verbose = TRUE){
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Checking Coordinates \n"))
}
.seed_select <- function(verbose = TRUE) {
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Extracting Seed territories \n"))
}
.query_select <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Extracting Query territories \n"))
}
.extract_ter_prog <- function(seed, verbose = TRUE) {
    if (!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Extracting territory", seed, "\n"))
}
.rebuild_df <- function(verbose = TRUE) {
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t, " Rebuilding Data Frame from image \n"))
}

#----------------------/Format Conversion messages/------------------------#
.vtc <- function(verbose=TRUE){
  if(!verbose) return(NULL)
  t <- Sys.time()
  cat( paste(t," Converting Vesalius to Image\n"))
}
.vtdf <- function(verbose=TRUE){
  if(!verbose) return(NULL)
  t <- Sys.time()
  cat( paste(t," Converting Vesalius to Data frame\n"))
}



#------------------------------------/TBD/---------------------------------#
.msk <- function(em,verbose =TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat(paste(t," Applying Masks in dim",em, "\r"))
}

.chunk <- function(verbose=TRUE){
    if(!verbose) return(NULL)
    t <- Sys.time()
    cat( paste(t," Chunking images into bite size pieces\n"))
}
