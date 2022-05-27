#-----------------------------/Running Seurat/---------------------------------#
#------------------------------------------------------------------------------#
# Running Seurat on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(grid)

#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#

# Out put directory
input <- "/home/pcnmartin/Vesalius/Simulation"
output <- paste0(input,"/seuratSim")
if(!dir.exists(output)){
    dir.create(output)
}

# Time benchmarking file
time <- paste0(input,"/seuratSim/time.txt")
if(!file.exists(time)){
    file.create(time)
}

# Get files
simFiles <- list.files(pattern = ".csv",full.names = TRUE)
fileTag <- list.files(pattern = ".csv",full.names = FALSE)
counts <- read.table("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz", header = TRUE )
rownames(counts) <- counts[,1]
counts <- counts[,-1]

# run benchmarking

for(i in seq_along(files)){
    sim <- read.table(simFiles[i], sep = ",", header=T)
    subCounts <- counts[,sim$barcodes]
    #----------------------------------------------------------------------------#
    # Rename barcodes to avoid potential duplicated names
    #----------------------------------------------------------------------------#
    colnames(subCounts) <- paste0("bar_",seq_len(ncol(subCounts)))


    rownames(sim) <- sim$simBarcode
    ss <- new(Class = 'SlideSeq',
              assay = "Spatial",
              coordinates = sim[,c("x","y")])

    rownames(ss@coordinates) <- sim$simBarcode

    st <- CreateSeuratObject(counts, assay ="Spatial")
    ss <- ss[Cells(x = st)]
    DefaultAssay(object = ss) <- "Spatial"
    st[["slice1"]] <- ss
    st <- AddMetaData(st,metadata = sim$ter,col.name = "Territory")
    s <- Sys.time()
    st <- SCTransform(st,assay ="Spatial") %>%
           RunPCA(dims =1:30) %>%
           RunUMAP(dims =1:30) %>%
           FindNeighbors() %>%
           FindClusters(st,verbose = FALSE)




    seuData <- FetchData(st, c("UMAP_1","UMAP_2","seurat_clusters")) %>%
               group_by(seurat_clusters) %>%
               mutate(xs = mean(UMAP_1), ys = mean(UMAP_2))
    coordSeu <- GetTissueCoordinates(st)
    seuData <- cbind(seuData,coordSeu[,c("x","y")])
    t <- Sys.time() - s
    ftime <- paste0(fileTag[i],",",as.numeric(t),",",units(t),"\n")
    cat(ftime,file=time,append=TRUE)

    ari <- adjustedRandIndex(tmpPred$seurat_clusters,tmpSim$territory)
    vi <- vi.dist(tmpPred$seurat_clusters,tmpSim$territory)

    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)

    fileOut <- paste0(output,"/Seurat_",fileTags[i])
    write.table(seuData,file =fileOut,sep =",",quote=F)
    rm(seuData); gc()
}
