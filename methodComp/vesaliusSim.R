#-----------------------------/Running Vesalius/---------------------------------#
#------------------------------------------------------------------------------#
# Running Vesalius on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#

library(vesalius)
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(grid)


library(mclust)
library(mcclust)


#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#

# Out put directory
input <- "/home/pcnmartin/Vesalius/Simulation"
output <- paste0(input,"/vesaliusSim")
if(!dir.exists(output)){
    dir.create(output)
}

# Time benchmarking file
time <- paste0(input,"/vesaliusSim/time.txt")
if(!file.exists(time)){
    file.create(time)
}

# performance benchmarking file
perf <- paste0(input,"/vesaliusSim/performance.txt")
if(!file.exists(perf)){
    file.create(perf)
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
    if(grepl(pattern = "dot", x = simFiles[i])){
        colDepth <- 6
        iter <- 5
    } else {
        colDepth <- c(81,27,9,3)
        iter <- 15
    }
    #----------------------------------------------------------------------------#
    # Rename barcodes to avoid potential duplicated names
    #----------------------------------------------------------------------------#
    colnames(subCounts) <- paste0("bar_",seq_len(ncol(subCounts)))


    rownames(sim) <- sim$simBarcode
    ss <- new(Class = 'SlideSeq',
              assay = "Spatial",
              coordinates = sim[,c("x","y")])

    rownames(ss@coordinates) <- sim$simBarcode

    st <- CreateSeuratObject(subCounts, assay ="Spatial")
    ss <- ss[Cells(x = st)]
    DefaultAssay(object = ss) <- "Spatial"
    st[["slice1"]] <- ss
    st <- AddMetaData(st,metadata = sim$ter,col.name = "Territory")

    s <- Sys.time()

    ves <- NormalizeData(st) %>%
           FindVariableFeatures(nfeatures=2000) %>%
           ScaleData()%>%
           rgbUMAP()%>%
           buildImageArray(resolution=50,filterThreshold=1,cores =1)%>%
           #equalizeHistogram(sleft =2.5,sright=2.5)%>%
           regulariseImage(lambda = 5,niter=200)%>%
           iterativeSegmentation.array(colDepth=colDepth,
                                       smoothIter = iter,
                                       method = c("iso","box"),
                                       sigma=6,box = 15,
                                       useCenter = T) %>%
          isolateTerritories.array(captureRadius=0.1,minBar=1) %>%
          filter(tile==1) %>%
          distinct(barcodes, .keep_all =TRUE)

    t <- Sys.time() - s
    ftime <- paste0(fileTag[i],",",as.numeric(t),",",units(t),"\n")
    cat(ftime,file=time,append = TRUE)

    aligned <- match(tmpPred$barcodes, tmpSim$barcodes)
    simTer <- as.character(tmpSim$territory)[aligned]
    ari <- adjustedRandIndex(tmpPred$territory,simTer)
    vi <- vi.dist(tmpPred$territory,simTer)

    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)

    fileOut <- paste0(output,"/Vesalius_",fileTags[i])
    write.table(ves,file =fileOut,sep =",",quote=F)
    rm(ves); gc()
}


pdf("test.pdf", width = 12, height=5)
g <- ggplot(ves, aes(x,y, col = as.factor(territory))) + geom_point()
g1 <- ggplot(sim, aes(x,y, col = as.factor(ter))) + geom_point()
g + g1
dev.off()
