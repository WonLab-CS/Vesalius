
#-----------------------------/Loading all data/-------------------------------#
#------------------------------------------------------------------------------#
# First we load all data sets
# Here we load both the mouse Hippocampus and the mouse embryo
#------------------------------------------------------------------------------#
slideTagBrain <- "Puck_200115_08"
slideBeadsBrain <-"~/group/slide_seqV2/Puck_200115_08_bead_locations.csv"
slideCountsBrain <- "~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz"
#slideIPBrain <- "~/group/slide_seqV2/IP/PCA_to_RGB_log_2000_slice1_Puck_200115_08.csv"

brainCoord <- utils::read.csv(slideBeadsBrain, header=F)
brainCoord <- brainCoord[-1,]
colnames(brainCoord)<- c("barcodes","xcoord","ycoord")

beadBrain <- ReadSlideSeq(slideBeadsBrain)

brainCounts <- read.table(slideCountsBrain, header = TRUE )
rownames(brainCounts) <- brainCounts[,1]
brainCounts <- brainCounts[,-1]

#------------------------------------------------------------------------------#
# We still create a Seurat object
# It's just that Seurat is such a great package to handling single cell data
# Hats off to the Seurat team
#------------------------------------------------------------------------------#
brainCounts <- CreateSeuratObject(brainCounts, assay ="Spatial")
beadBrain <- beadBrain[Cells(x = brainCounts)]
DefaultAssay(object = beadBrain) <- "Spatial"
brainCounts[["slice1"]] <- beadBrain

#------------------------------------------------------------------------------#
# We keep a Raw data set. We will use the raw counts when we subset out the
# beads in isolated territories. We want to avoid using the values resulting
# from overall normalisation.
# We want to avoid to "normalise" already normalised and scaled data.
#------------------------------------------------------------------------------#
brainRaw <- brainCounts



#----------------------------/Seurat Clustering/-------------------------------#
#------------------------------------------------------------------------------#
# Pre-processing data using the same approaches as the one used when converting
# PCA loadings to RGB colour space
#------------------------------------------------------------------------------#
brainCounts <- NormalizeData(brainCounts)
brainCounts <- FindVariableFeatures(brainCounts, nfeatures = 2000)
brainCounts <- ScaleData(brainCounts)

#------------------------------------------------------------------------------#
# Running Seurat pipeline and plotting results
# we plot both the UMAP projections and the mapped clusters
#------------------------------------------------------------------------------#
seu <- RunPCA(brainCounts)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.7, verbose = FALSE)
seuData <- FetchData(seu, c("UMAP_1","UMAP_2","seurat_clusters")) %>%
           group_by(seurat_clusters) %>%
           mutate(xs = mean(UMAP_1), ys = mean(UMAP_2))

#------------------------------------------------------------------------------#
# This section just rebuilds a data frame to make it easier to plot with ggplot
#------------------------------------------------------------------------------#
coordSeu <- GetTissueCoordinates(seu)
seuData <- cbind(seuData,coordSeu[,c("x","y")])
seu_col <- length(unique(seuData$seurat_clusters))
seu_pal <- colorRampPalette(brewer.pal(8, "Accent"))


coord_seu <- ggplot(seuData, aes(x,y,col = as.factor(seurat_clusters))) +
            geom_point(size = 0.25, alpha = 0.65) +
            theme_void() +
            scale_color_manual(values = seu_pal(seu_col)) +
            theme(legend.text = element_text(size = 12),
                  legend.position = "left",
                  plot.title = element_text(size=15)) +
            guides(colour = guide_legend(override.aes = list(size=5)))+
            labs(colour = "Cluster nr.", title = "Seurat - Clusters")






library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tvR)
library(sp)
library(grid)





input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)
time <-vector("list",length(input))
count <- 1
n <-c(7,7,7,7,5,5,5,5,7,7,7,7)
for(i in input){
  s <- Sys.time()
  output <- paste0(i,"/Seurat")
  if(!dir.exists(output)){
      dir.create(output)
  }

  sec <- Read10X(i)
  img <- Read10X_Image(i)
  sec <- CreateSeuratObject(counts = sec, assay = "Spatial")
  img <- img[Cells(x = sec)]
  DefaultAssay(object = img) <- "Spatial"
  sec[["slice1"]] <- img

  sec <- NormalizeData(sec)
  sec <- FindVariableFeatures(sec, selection.method = "vst", nfeatures = 2000)
  sec <- ScaleData(sec)

  sec <- RunPCA(sec, npcs = 30)
  sec <- FindNeighbors(sec, dims = 1:30,verbose =F)
  for(k in seq(0.5,1,by = 0.01)){
      sec <- FindClusters(sec, resolution = k, verbose = FALSE)
      cl <- FetchData(sec, c("seurat_clusters"))$seurat_clusters
      if(length(levels(cl))==n[count]){
        break
      }
  }

  tag <- strsplit(i, "/")[[1]]
  tag <- tag[length(tag)]
  file <- paste0(i,"/Seurat/",tag,"_seurat.rda")
  save(sec,file = file)
  e <- Sys.time()
  time[[count]] <- e - s
  count <- count +1

}
save(time,file = "Seurat_time.Rda")
