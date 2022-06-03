#-----------------------------/Running Seurat/---------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to Seurat runs for the purpose of
# performance comaprison. Here you will find the code related to Slide-seq V2
# and the DLPFC Visium data sets
#------------------------------------------------------------------------------#

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(grid)

#------------------------------------------------------------------------------#
# Slide-seq V2
#------------------------------------------------------------------------------#

input <- "~/group/slide_seqV2/"

slideBeads <-paste0(input,list.files(input, pattern ="locations.csv"))
slideCounts <- paste0(input,list.files(input, pattern ="digital_expression.txt.gz"))

slideTag <- sapply(strsplit(slideBeads,"V2/"),"[[",2)
slideTag <- gsub("_bead_locations.csv","",slideTag)


output <- paste0(input,"/SeuratBenchMarking/")
if(!dir.exists(output)){
    dir.create(output)
}

time <-vector("list",length(slideBeads))

for(i in seq_along(slideBeads)){
#for(i in c(10,11)){
  s <- Sys.time()
  #brainCoord <- utils::read.csv(slideBeads[i], header=F)
  #brainCoord <- brainCoord[-1,]
  #colnames(brainCoord)<- c("barcodes","xcoord","ycoord")
  brainCounts <- read.table(slideCounts[i], header = TRUE )
  rownames(brainCounts) <- brainCounts[,1]
  brainCounts <- brainCounts[,-1]
  seu <- CreateSeuratObject(brainCounts, assay ="Spatial")

  beadBrain <- ReadSlideSeq(slideBeads[i])
  if(dim(beadBrain@coordinates)[2] == 3){
      colnames(beadBrain@coordinates) <- c("x","y","cells")
      rownames(beadBrain@coordinates) <- beadBrain@coordinates$cells
      DefaultAssay(object = beadBrain) <- "Spatial"
      seu[["slice1"]] <- beadBrain
  } else {
    beadBrain <- beadBrain[Cells(x = brainCounts)]
    DefaultAssay(object = beadBrain) <- "Spatial"
    seu[["slice1"]] <- beadBrain
  }


  seu <- SCTransform(seu,variable.features.n = 2000, assay = "Spatial")


  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:30)
  seu <- FindNeighbors(seu, dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.7, verbose = FALSE)
  seuData <- FetchData(seu, c("UMAP_1","UMAP_2","seurat_clusters")) %>%
             group_by(seurat_clusters) %>%
             mutate(xs = mean(UMAP_1), ys = mean(UMAP_2))


  coordSeu <- GetTissueCoordinates(seu)
  seuData <- cbind(seuData,coordSeu[,c("x","y")])
  fileOut <- paste0(output,slideTag[i],"_SSV2_BM.Rda")
  save(seuData,file= fileOut)
  e <- Sys.time()
  time[[i]] <- e - s

}
save(time, file = paste0(output,"Seurat_SSV2_Time.Rda"))


#seuData$seurat_clusters <- as.numeric(as.character(seuData$seurat_clusters)) + 1
#seu_col <- length(unique(seuData$seurat_clusters))
#seu_pal <- colorRampPalette(brewer.pal(8, "Accent"))


#coord_seu <- ggplot(seuData, aes(x,y,col = as.factor(seurat_clusters))) +
#            geom_point(size = 0.25, alpha = 0.65) +
#            theme_void() +
#            scale_color_manual(values = seu_pal(seu_col)) +
#            theme(legend.text = element_text(size = 12),
#                  legend.position = "left",
#                  plot.title = element_text(size=15)) +
#            guides(colour = guide_legend(override.aes = list(size=5)))+
#            labs(colour = "Cluster nr.", title = "Seurat - Clusters")




#------------------------------------------------------------------------------#
# Visium DLPFC
#------------------------------------------------------------------------------#

### Input directories containing Visium data

# input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)
# time <-vector("list",length(input))
# count <- 1
# n <-c(7,7,7,7,5,5,5,5,7,7,7,7)
# for(i in input){
#   s <- Sys.time()
#   output <- paste0(i,"/Seurat")
#   if(!dir.exists(output)){
#       dir.create(output)
#   }
#
#   sec <- Read10X(i)
#   img <- Read10X_Image(i)
#   sec <- CreateSeuratObject(counts = sec, assay = "Spatial")
#   img <- img[Cells(x = sec)]
#   DefaultAssay(object = img) <- "Spatial"
#   sec[["slice1"]] <- img
#
#   #sec <- NormalizeData(sec)
#   #sec <- FindVariableFeatures(sec, selection.method = "vst", nfeatures = 2000)
#   #sec <- ScaleData(sec)
#   sec <- SCTransform(sec,variable.features.n = 2000, assay = "Spatial")
#   sec <- RunPCA(sec, npcs = 30)
#   sec <- FindNeighbors(sec, dims = 1:30,verbose =F)
#   for(k in seq(0.5,1,by = 0.01)){
#       sec <- FindClusters(sec, resolution = k, verbose = FALSE)
#       cl <- FetchData(sec, c("seurat_clusters"))$seurat_clusters
#       if(length(levels(cl))==n[count]){
#         break
#       }
#   }
#
#   tag <- strsplit(i, "/")[[1]]
#   tag <- tag[length(tag)]
#   file <- paste0(i,"/Seurat/",tag,"_seurat.rda")
#   save(sec,file = file)
#   e <- Sys.time()
#   time[[count]] <- e - s
#   count <- count +1
#
# }
# ## Saving run time
# save(time,file = "~/group/visium/DLPFC_globus/Seurat_Visium_Time.Rda")
