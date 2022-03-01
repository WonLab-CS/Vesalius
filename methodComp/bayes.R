#---------------------------/Running BayesSpace/-------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to BayesSpace runs for the purpose of
# performance comaprison. Here you will find the code related to Slide-seq V2
# and the DLPFC Visium data sets
#------------------------------------------------------------------------------#

### NOTE To run BayesSpace on slide_seqV2 please install BayesSpace via:
library(devtools)
install_github("patrickCNMartin/BayesSpace")
## This contain a modified version of BayesSpace that accomodates Slide-seq data

library(BayesSpace)
library(SingleCellExperiment)
library(Seurat)

#------------------------------------------------------------------------------#
# Slide-seq V2 - Hippocampus
#------------------------------------------------------------------------------#

set.seed(1)

slideTag <- c("Puck_200115_08","Puck_190926_03")

slideBeads <-c("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv",
               "~/group/slide_seqV2/Puck_190926_03_bead_locations.csv")

slideCounts <- c("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz",
                 "~/group/slide_seqV2/Puck_190926_03.digital_expression.txt.gz")


input <- "~/group/slide_seqV2/"

#slideBeads <-paste0(input,list.files(input, pattern ="locations.csv"))
#slideCounts <- paste0(input,list.files(input, pattern ="digital_expression.txt.gz"))

#slideTag <- sapply(strsplit(slideBeads,"V2/"),"[[",2)
#slideTag <- gsub("_bead_locations.csv","",slideTag)


output <- paste0(input,"BayesSpaceBenchMarking/")
if(!dir.exists(output)){
    dir.create(output)
}

## removing all the intermeidate files
## set to out put dir


time <-vector("list",length(slideBeads))
#for(i in seq_along(slideBeads)){
for(i in 2){
  s <- Sys.time()
  coord <- utils::read.csv(slideBeads[i], header=F)
  coord <- coord[-1,]
  colnames(coord)<- c("barcodes","xcoord","ycoord")
  rownames(coord) <- coord$barcodes

  raw <- read.table(slideCounts[i], header = TRUE )
  rownames(raw) <- raw[,1]
  raw <- raw[,-1]
  raw <- as.matrix(raw)

  sce <- SingleCellExperiment(list(counts=as(raw, "dgCMatrix")))


  coord <- coord[colnames(raw),c("ycoord","xcoord")]
  colnames(coord)<- c("row","col")
  colData(sce) <- DataFrame(coord)

  genes <- data.frame(rownames(raw), rownames(raw))
  colnames(genes) <- c("gene_id","gene_name")
  rownames(genes) <- rownames(raw)
  rowData(sce) <- genes
  meta <- list(sample ="Puck_200115_08", dataset = "SSv2",
            BayesSpace.data = list(platform = "SS",is.enhanced=FALSE))
  metadata(sce) <- meta

  sce <- spatialPreprocess(sce, platform="SS",skip.PCA=FALSE,
                               n.PCs=30, n.HVGs=2000, log.normalize=TRUE)
  #sce <- qTune(sce, qs=seq(10,30), platform="SS", d=30)

  sce <- spatialCluster(sce,q=15,d=30,nrep=50000, gamma=3, save.chain=TRUE,platform = "SS")
  time[[i]] <- Sys.time() -s
  # BayesSpace Run time
  bayesData <- as.data.frame(colData(sce))
  save(bayesData, file= paste0(output,slideTag[i],"_SSV2_BM.Rda"))
}
save(time, file = paste0(output,"Bayes_SSV2_time.Rda"))





#------------------------------------------------------------------------------#
# Visium DLPFC
#------------------------------------------------------------------------------#
# Directories containing Visium data
# input <- list.dirs("~/group/visium/DLPFC_globus/",recursive =F)
#
# time <-vector("list",length(input))
# count <- 1
# n <-c(7,7,7,7,5,5,5,5,7,7,7,7)
# for(i in input){
#
#   s <- Sys.time()
#   output <- paste0(i,"/BayesSpace")
#   if(!dir.exists(output)){
#       dir.create(output)
#   }
#   tag <- strsplit(i, "/")[[1]]
#   tag <- tag[length(tag)]
#   ## The file name requirements between Seurat and BayesSpace are different
#   ## I'm not going to start renaming everything just for that
#   ## SO loading with Seurat and then parsing to BayesSpace
#
#   sec <- Read10X(i)
#   img <- Read10X_Image(i)
#
#
#   sce <- SingleCellExperiment(list(counts=sec))
#
#   coord <- img@coordinates[,c("row","col")]
#   coord <- coord[colnames(sec),]
#   colData(sce) <- DataFrame(coord)
#
#   genes <- data.frame(rownames(sec), rownames(sec))
#   colnames(genes) <- c("gene_id","gene_name")
#   rownames(genes) <- rownames(sec)
#   rowData(sce) <- genes
#   meta <- list(sample =tag, dataset = "DPLFC",
#           BayesSpace.data = list(platform = "Visium",is.enhanced=FALSE))
#   metadata(sce) <- meta
#
#   #set.seed(149)
#   sce <- spatialPreprocess(sce, platform="Visium",skip.PCA=FALSE,
#                              n.PCs=15, n.HVGs=2000, log.normalize=TRUE)
#
#   sce <- spatialCluster(sce,q=n[count],d=15,nrep=50000, gamma=3, save.chain=TRUE,platform = "Visium")
#
#   file <- paste0(i,"/BayesSpace/",tag,"_BayesSpace.rda")
#   save(sce,file = file)
#   e<- Sys.time()
#   time[[count]] <- e - s
#   count <- count +1
#
# }
# #Save Bayes Space Full run time
# save(time,file = "~/group/visium/DLPFC_globus/BayesSpace_time.Rda")
