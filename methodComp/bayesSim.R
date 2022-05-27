#---------------------------/Running BayesSpace/-------------------------------#
#------------------------------------------------------------------------------#
# Running BayesSpace on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#

### NOTE To run BayesSpace on slide_seqV2 please install BayesSpace via:
library(devtools)
install_github("patrickCNMartin/BayesSpace")
## This contain a modified version of BayesSpace that accomodates Slide-seq data

library(BayesSpace)
library(SingleCellExperiment)
library(Seurat)

# Out put directory
input <- "/home/pcnmartin/Vesalius/Simulation"
output <- paste0(input,"/bayesSim")
if(!dir.exists(output)){
    dir.create(output)
}

# Time benchmarking file
time <- paste0(input,"/bayesSim/time.txt")
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
    if(grepl(pattern = "dot", x = simFiles[i])){
        # Dotted regime should contain 10 total territories
        # one background and 9 "spots"
        q <- 10
    } else {
        q <- 3
    }

    subCounts <- counts[,sim$barcodes]
    #----------------------------------------------------------------------------#
    # Rename barcodes to avoid potential duplicated names
    #----------------------------------------------------------------------------#
    sce <- SingleCellExperiment(list(counts=as(as.matrix(subCounts), "dgCMatrix")))
    coord <- sim[colnames(sce),c("y","x")]
    colnames(coord)<- c("row","col")
    colData(sce) <- DataFrame(coord)

    genes <- data.frame(rownames(subCounts), rownames(subCounts))
    colnames(genes) <- c("gene_id","gene_name")
    rownames(genes) <- rownames(subCounts)
    rowData(sce) <- genes
    meta <- list(sample ="BayesSpaceSim", dataset = "Simulated",
                 BayesSpace.data = list(platform = "SS",is.enhanced=FALSE))
    metadata(sce) <- meta
    s <- Sys.time()
    sce <- spatialPreprocess(sce,platform="SS",skip.PCA=FALSE,
                               n.PCs=30, n.HVGs=2000, log.normalize=TRUE)
      #sce <- qTune(sce, qs=seq(10,30), platform="SS", d=30)

    sce <- spatialCluster(sce,q=q,platform = "SS")
    bayes <- as.data.frame(colData(sce))

    t <- Sys.time() - s
    ftime <- paste0(fileTag[i],",",as.numeric(t),",",units(t),"\n")
    cat(ftime,file=time,append=TRUE)

    aligned <- match(rownames(tmpPred), tmpSim$barcodes)
    simTer <- as.character(tmpSim$territory)[aligned]
    ari <- adjustedRandIndex(tmpPred$spatial.cluster,simTer)
    vi <- vi.dist(tmpPred$spatial.cluster,simTer)

    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)

    fileOut <- paste0(output,"/BayesSpace_",fileTags[i])
    write.table(bayes,file =fileOut,sep =",",quote=F)
    rm(bayes); gc()
}
