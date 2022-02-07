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

  seu <- buildImageArray(seu,invert=F,
                         filterThreshold = 0.9975,
                         resolution = 40, cores = 1)


  seu <- equalizeHistogram(seu,sleft = 2.5, sright=2.5,invert =F)


  seu <- regulariseImage(seu,lambda = 10,niter = 200, normalise=T,invert=F)

  seu <- iterativeSegmentation.array(seu, colDepth = 6,
                                          smoothIter = 20,
                                          method = c("iso","median"),
                                          sigma=1.5,box = 10,
                                          useCenter = T,
                                          invert =F)



  seu <- isolateTerritories.array(seu,captureRadius = 0.008,minBar = 40)

  fileOut <- paste0(output,slideTag[i],"_SSV2_BM.Rda")
  save(seu,file = fileOut)
  e <- Sys.time()
  time[[i]] <- e - s

}
save(time, file = paste0(output,"Vesalius_SSV2_Time.Rda"))





#------------------------------------------------------------------------------#
# Visium DLPFC
#------------------------------------------------------------------------------#

### Input directories containing Visium data

input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)
time <-vector("list",length(input))
count <- 1

for(i in input){

  s <- Sys.time()
  output <- paste0(i,"/Vesalius")
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



  sec <- rgbUMAP(sec, pcs = 30, conserveSparse =F)

  image <- buildImageArray(sec, resolution = 100, filterThreshold = 1, invert =F)
  image <- equalizeHistogram(image, sleft = 15,sright =15,invert =F)

  image <- iterativeSegmentation.array(image,
      colDepth = seq(11,7, by =-2),
      smoothIter = 1,
      method = c("median","iso"),
      sigma=1,
      box = 13,
      useCenter = T,
      invert =F)
  image <- isolateTerritories.array(image, captureRadius = 0.035,minBar =5)

  tag <- strsplit(i, "/")[[1]]
  tag <- tag[length(tag)]

  file <- paste0(i,"/Vesalius/",tag,"_vesalius.rda")
  save(image,file = file)
  e <- Sys.time()
  time[[count]] <- e - s
  count <- count +1

}

## Saving run time
save(time,file = "~/group/visium/DLPFC_globus/Vesalius_Visium_Time.Rda")

#------------------------------------------------------------------------------#
# Seq Scope
#------------------------------------------------------------------------------#


input <- "~/group/seqScope/"

seqScope <-paste0(input,list.files(input, pattern =".rds"))
seqScopeTag <- sapply(strsplit(seqScope,"Scope/"),"[[",2)
seqScopeTag <- gsub(".rds","",seqScopeTag)


output <- paste0(input,"/VesaliusBenchMarking/")
if(!dir.exists(output)){
    dir.create(output)
}

time <-vector("list",length(seqScope))

for(i in seq_along(seqScope)){
  s <- Sys.time()
  seu <- readRDS(seqScope[i])


  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 2000)
  seu <- ScaleData(seu)


  seu <- rgbUMAP(seu,pcs=30, conserveSparse = FALSE)
  # Adding a slight jitter to avoid overlapping coordinates
  seu$x <-jitter(seu$x/20, factor=0.1)
  seu$y <-jitter(seu$y/20, factor=0.1)

  seu <- buildImageArray(seu,invert=F,
                         filterThreshold = 0.99,
                         resolution = 40, cores = 1)


  seu <- equalizeHistogram(seu,sleft = 2.5, sright=2.5,invert =F)


  seu <- regulariseImage(seu,lambda = 10,niter = 200, normalise=T,invert=F)

  seu <- iterativeSegmentation.array(seu, colDepth = 6,
                                          smoothIter = 20,
                                          method = c("iso","median"),
                                          sigma=1.5,box = 10,
                                          useCenter = T,
                                          invert =F)



  seu <- isolateTerritories.array(seu,captureRadius = 0.008,minBar = 40)

  fileOut <- paste0(output,seqScopeTag[i],"_SeqScope_BM.Rda")
  save(seu,file = fileOut)
  e <- Sys.time()
  time[[i]] <- e - s
  rm(seu)
  gc()
}
save(time, file = paste0(output,"Vesalius_SeqScope_Time.Rda"))
