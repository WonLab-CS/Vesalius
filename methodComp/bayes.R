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
slideTag <- c("Puck_200115_08")

slideBeads <-c("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv")

slideCounts <- c("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz")


coord <- utils::read.csv(slideBeads, header=F)
coord <- coord[-1,]
colnames(coord)<- c("barcodes","xcoord","ycoord")
rownames(coord) <- coord$barcodes

raw <- read.table(slideCounts, header = TRUE )
rownames(raw) <- raw[,1]
raw <- raw[,-1]
raw <- as.matrix(raw)
s <- Sys.time()
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
sce <- qTune(sce, qs=seq(10,30), platform="SS", d=30)

sce <- spatialCluster(sce,q=30,d=30,nrep=50000, gamma=3, save.chain=TRUE,platform = "SS")
time <- Sys.time() -s
# BayesSpace Run time
save(time, file = "BayesSSV2_time.Rda")


bayesData <- as.data.frame(colData(bayes))
bayes_col <- length(unique(bayesData$spatial.cluster))
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


#------------------------------------------------------------------------------#
# Slide-seq V2 - Embryo
#------------------------------------------------------------------------------#

set.seed(1)
slideTag <- c("Puck_200115_08")

slideBeads <-c("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv")

slideCounts <- c("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz")


coord <- utils::read.csv(slideBeads, header=F)
coord <- coord[-1,]
colnames(coord)<- c("barcodes","xcoord","ycoord")
rownames(coord) <- coord$barcodes

raw <- read.table(slideCounts, header = TRUE )
rownames(raw) <- raw[,1]
raw <- raw[,-1]

s <- Sys.time()
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
sce <- qTune(sce, qs=seq(10,30), platform="SS", d=30)

sce <- spatialCluster(sce,q=30,d=30,nrep=50000, gamma=3, save.chain=TRUE,platform = "SS")
time <- Sys.time() -s
# BayesSpace Run time
save(time, file = "BayesSSV2Embryo_time.Rda")


bayesData <- as.data.frame(colData(bayes))
bayes_col <- length(unique(bayesData$spatial.cluster))
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














#------------------------------------------------------------------------------#
# Visium DLPFC
#------------------------------------------------------------------------------#
# Directories containing Visium data
input <- list.dirs("visium/DLPFC",recursive =F)

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
#Save Bayes Space Full run time
save(time,file = "BayesSpace_time.Rda")
