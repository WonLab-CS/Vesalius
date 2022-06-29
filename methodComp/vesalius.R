#-----------------------------/Running Vesalius/---------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to Vesalius runs for the purpose of
# performance comaprison. Here you will find the code related to Slide-seq V2
# the DLPFC Visium data sets and the seqScope data sets.
#------------------------------------------------------------------------------#

library(vesalius)
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


output <- paste0(input,"/VesaliusBenchMarking/")
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


  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 2000)
  seu <- ScaleData(seu)


  seu <- rgbUMAP(seu,pcs=30, conserveSparse = FALSE)

  seu <- buildImageArray(seu,
                         filterGrid =0.01,
                         filterThreshold = 0.9975,
                         resolution = 40, cores = 1)





  seu <- regulariseImage(seu,lambda = 5,niter = 50, normalise=T)

  seu <- iterativeSegmentation.array(seu, colDepth = 11,
                                          smoothIter = 20,
                                          method = c("iso","box"),
                                          sigma=1,box = 15,
                                          useCenter = T)



  seu <- isolateTerritories.array(seu,captureRadius = 0.008,minBar = 50)

  fileOut <- paste0(output,slideTag[i],"_SSV2_BM.Rda")
  save(seu,file = fileOut)
  e <- Sys.time()
  time[[i]] <- e - s

}
save(time, file = paste0(output,"Vesalius_SSV2_Time.Rda"))


## Running Vesalius in lower resolution
#Here, we demonstrate running Vesalius over 12 10X Visium data sets taken from
#[here](https://www.nature.com/articles/s41593-020-00787-0). All files were
#downloaded and distributed into seperate directories. We assess run time for
#each data set and compared performance using an ARI score. Please note that
#the comprative analysis code can be found
#[here](https://github.com/patrickCNMartin/Vesalius/blob/main/methodComp/methodComp.R).



#------------------------------------------------------------------------------#
# Visium DLPFC
#------------------------------------------------------------------------------#

### Input directories containing Visium data

# input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)
# time <-vector("list",length(input))
# count <- 1
# n <-c(7,7,7,7,5,5,5,5,7,7,7,7)
# for(i in input){
#
#   s <- Sys.time()
#   output <- paste0(i,"/Vesalius")
#   if(!dir.exists(output)){
#       dir.create(output)
#   }
#   sec <- Read10X(i)
#   img <- Read10X_Image(i)
#   sec <- CreateSeuratObject(counts = sec, assay = "Spatial")
#   img <- img[Cells(x = sec)]
#   DefaultAssay(object = img) <- "Spatial"
#   sec[["slice1"]] <- img
#
#   sec <- NormalizeData(sec)
#   sec <- FindVariableFeatures(sec, nfeatures = 2000)
#   sec <- ScaleData(sec)
#
#
#
#   sec <- rgbUMAP(sec, pcs = 30, conserveSparse =F)
#
#   image <- buildImageArray(sec,filterGrid=1, resolution = 100, filterThreshold = 1, invert =F)
#   image <- equalizeHistogram(image, sleft = 15,sright =15,invert =F)
#
#   image <- iterativeSegmentation.array(image,
#       colDepth = seq(15,n[count], by =-2),
#       smoothIter = 1,
#       method = c("median","iso"),
#       sigma=1,
#       box = 13,
#       useCenter = T,
#       invert =F)
#   image <- isolateTerritories.array(image, captureRadius = 0.035,minBar =5)
#
#   tag <- strsplit(i, "/")[[1]]
#   tag <- tag[length(tag)]
#
#   file <- paste0(i,"/Vesalius/",tag,"_vesalius.rda")
#   save(image,file = file)
#   e <- Sys.time()
#   time[[count]] <- e - s
#   count <- count +1
#
# }
#
# ## Saving run time
# save(time,file = "~/group/visium/DLPFC_globus/Vesalius_Visium_Time.Rda")
