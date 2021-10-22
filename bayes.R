
library(BayesSpace)
library(SingleCellExperiment)
library(Seurat)


# Standard libraries - just in case you are a command line fiend
library(utils)
library(stats)
library(graphics)
library(grDevices)
### Set seed
set.seed(1)
slideTag <- c("Puck_200115_08")

slideBeads <-c("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv")

slideCounts <- c("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz")


#-----------------------/Slide-seq analysis/-----------------------------------#
#------------------------------------------------------------------------------#
# Next we can set a few parameters for the pre-processing
#------------------------------------------------------------------------------#
## number of variable features
nfeatures <- 2000

## FindVariableFeatures method
method <- "vst"

## Number of PCA slices - each slice will cover 3 PCs
slices <- 3

## Plot output directory
plots <- "~/group/slide_seqV2/plots/"

## Invert colour code - F = black dominant & T = White dominant
## Colour inversion is not strickly speaking necessary. Light colours are
## however easier to visualise.
invert <- TRUE

## Image processing Matrices output directory
IP <- "~/group/slide_seqV2/IP/"

#------------------------------------------------------------------------------#
# Bayes Space - custom functions
# will replace the internal function for the purpose of testing
#------------------------------------------------------------------------------#

.find_neighbors <- function(sce, platform) {
    # we ar going to use switch instead
    # it will be easir for our purpose

    neighbors <- switch(platform,
                        "Visium" =.bs(sce,data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                                                     y.offset=c( 0, 0, -1, -1,  1, 1))),
                         "ST" = .bs(sce,data.frame(x.offset=c( 0, 1, 0, -1),
                                               y.offset=c(-1, 0, 1,  0))),
                         "SS" = .ss(sce))



    return(neighbors)
}

## Get array coordinates (and label by index of spot in SCE)


.bs <- function(sce,offset){
  spot.positions <- colData(sce)[, c("col", "row")]
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))

  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset

  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions),
                     as.data.frame(spot.positions),
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)

  ## Shift to zero-indexing for C++
  neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1

  ## Group neighbor indices by spot
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary,
                         neighbors$spot.idx.neighbor), ]


  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)

  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  df_j <- map(df_j, function(nbrs) discard(nbrs, function(x) is.na(x)))

  ## Log number of spots with neighbors
  n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
                                 ncol(sce), " spots.")
  return(df_j)
}

.ss <- function(sce){
  # We will try to use the same format just to reduce error risk
  # I am assuming that the code is optimised for that input
  # Also keeping the col and row thing just because I don't know
  # how many other instance there are of needing col and row later
  # I will make the assumption that 6 neighbors are good enough
  # after some filtering
  message("custom shit is working")
  spot.positions <- colData(sce)[, c("col", "row")]
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  dist <- parallel::mclapply(seq_len(nrow(spot.positions)), function(idx,mat){
                      xo <- mat$col[idx]
                      yo <- mat$row[idx]
                      xp <- mat$col
                      yp <- mat$row
                      distance <- sqrt(((abs(xp-xo))^2 + (abs(yp-yo))^2))
                      distance <- distance[distance !=0]
                      index <- order(distance,decreasing = FALSE)[1:6]
                      index <- index - 1
                      index <- sort(index)
                      return(index)
  }, spot.positions, mc.cores=4)


  return(dist)
}

  #----------------------------------------------------------------------------#
  # Loading coordinates
  #----------------------------------------------------------------------------#
  bead <- ReadSlideSeq(slideBeads)

  message(paste("Loading Count file",slideTag))
  #----------------------------------------------------------------------------#
  # Unconventional loading - however required as some data sets
  # Fail to load  - This ensures all data sets can be loaded
  #----------------------------------------------------------------------------#
  ctmp <- read.table(slideCounts, header = TRUE )
  rownames(ctmp) <- ctmp[,1]
  ctmp <- ctmp[,-1]

  #----------------------------------------------------------------------------#
  # Creating seurat spatial object
  # NOTE this code is taken from the Seurat source code as it does not seem that
  # Slide seq loading function are all exported
  # If this has been updated - this section can be changed accordingly
  #----------------------------------------------------------------------------#
  seu <- CreateSeuratObject(ctmp, assay ="Spatial")
  bead <- bead[Cells(x = seu)]
  DefaultAssay(object = bead) <- "Spatial"
  seu[["slice1"]] <- bead

  #----------------------------------------------------------------------------#
  # Getting top var features
  #----------------------------------------------------------------------------#
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 30)

  var_features <- VariableFeatures(seu)

  raw <- seu@assays$Spatial@counts
  raw <- raw[rownames(raw) %in% var_features,]
  counts <- GetAssayData(seu, slot ="data")
  counts <- counts[rownames(counts) %in% var_features,]

  sce <- SingleCellExperiment(list(counts=raw,logcounts = counts))

  coord <- GetTissueCoordinates(seu)[,c("y","x")]
  colnames(coord) <- c("row","col")
  coord$sizeFactor <- 1
  # rescale coords
  #coord$row <- floor((coord$row/max(coord$row)) * (0.1 * max(coord$row)))

  #coord$col <- floor((coord$col/max(coord$col)) * (0.1 * max(coord$col)))
  colData(sce) <- DataFrame(coord)

  genes <- data.frame(rownames(counts), rownames(counts))
  colnames(genes) <- c("gene_id","gene_name")
  rownames(genes) <- rownames(counts)
  rowData(sce) <- genes

  reducedDim(sce, "PCA") <- as.data.frame(seu[["pca"]]@cell.embeddings)

  meta <- list(sample = "Puck_200115_08", dataset = "SSV2",
          BayesSpace.data = list(platform = "SS",is.enhanced=FALSE))
  metadata(sce) <- meta
  bayes <- spatialCluster(sce,q=30)




  bayes_col <- length(levels(bayesData$spatial.cluster))
  bayes_pal <- colorRampPalette(brewer.pal(8, "Accent"))


  coord_bayes <- ggplot(bayesData, aes(col,row,col = as.factor(spatial.cluster))) +
              geom_point(size = 0.05, alpha = 0.8) +
              theme_void() +
              scale_color_manual(values = bayes_pal(bayes_col)) +
              theme(legend.text = element_text(size = 12),
                    legend.position = "left",
                    plot.title = element_text(size=15)) +
              guides(colour = guide_legend(override.aes = list(size=5)))+
              labs(colour = "Cluster nr.", title = "BayesSpace - Clusters")



### BayesSpace
library(BayesSpace)
library(SingleCellExperiment)
library(Seurat)
set.seed(1)
input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)

time <-vector("list",length(input))
count <- 1
n <-c(7,7,7,7,5,5,5,5,7,7,7,7)
for(i in input){

  s <- Sys.time()
  output <- paste0(i,"/BayesSpace")
  if(!dir.exists(output)){
      dir.create(output)
  }
  tag <- strsplit(i, "/")[[1]]
  tag <- tag[length(tag)]
  ## The file name requirements between Seurat and BayesSpace are different
  ## I'm not going to start renaming everything just for that
  ## SO loading with Seurat and then parsing to BayesSpace

  sec <- Read10X(i)
  img <- Read10X_Image(i)


  sce <- SingleCellExperiment(list(counts=sec))

  coord <- img@coordinates[,c("row","col")]
  coord <- coord[colnames(sec),]
  colData(sce) <- DataFrame(coord)

  genes <- data.frame(rownames(sec), rownames(sec))
  colnames(genes) <- c("gene_id","gene_name")
  rownames(genes) <- rownames(sec)
  rowData(sce) <- genes
  meta <- list(sample =tag, dataset = "DPLFC",
          BayesSpace.data = list(platform = "Visium",is.enhanced=FALSE))
  metadata(sce) <- meta

  #set.seed(149)
  sce <- spatialPreprocess(sce, platform="Visium",skip.PCA=FALSE,
                             n.PCs=15, n.HVGs=2000, log.normalize=TRUE)

  sce <- spatialCluster(sce,q=n[count],d=15,nrep=50000, gamma=3, save.chain=TRUE,platform = "Visium")

  file <- paste0(i,"/BayesSpace/",tag,"_BayesSpace.rda")
  save(sce,file = file)
  e<- Sys.time()
  time[[count]] <- e - s
  count <- count +1

}
save(time,file = "BayesSpace_time.Rda")
