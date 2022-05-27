#-----------------------------/Running Giotto/---------------------------------#
#------------------------------------------------------------------------------#
# Running Giotto on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#
library(Giotto)
library(Seurat)
library(ggplot2)
library(patchwork)
library(mclust)
library(mcclust)


#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#

# Out put directory
input <- "/home/pcnmartin/Vesalius/Simulation"
output <- paste0(input,"/giottoSim")
if(!dir.exists(output)){
    dir.create(output)
}

# Time benchmarking file
time <- paste0(input,"/giottoSim/time.txt")
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
    if(grepl(pattern = "dot", x = simFiles[i])){
        k <- 10
    } else {
        k <- 3
    }
    #----------------------------------------------------------------------------#
    # Rename barcodes to avoid potential duplicated names
    #----------------------------------------------------------------------------#
    colnames(subCounts) <- paste0("bar_",seq_len(ncol(subCounts)))


    s <- Sys.time()
    instruc <- createGiottoInstructions(save_plot=F, show_plot=F, save_dir=output)
    giotto <- createGiottoObject(raw_exprs = subCounts,
                                       spatial_locs = sim[,c("x","y")],
                                       instructions = instruc,
                                       cell_metadata = sim[,c("x","y","ter")])

    giotto <- filterGiotto(gobject = giotto,
                           expression_threshold = 1,
                           expression_values = c('raw'))

    giotto <- normalizeGiotto(gobject = giotto, scalefactor = 6000)

    giotto <- addStatistics(gobject = giotto)



    giotto <- calculateHVG(gobject = giotto)

    gene_metadata <- fDataDT(giotto)
    featgenes <- gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

    giotto <- runPCA(gobject = giotto,
                     genes_to_use = featgenes,
                     scale_unit = F,
                     center=T,
                     method="irlba")


    giotto <- runUMAP(giotto, dimensions_to_use = 1:30)

    giotto <- createSpatialNetwork(gobject=giotto,
                                   method='kNN',
                                   k=k,
                                   maximum_distance_knn=400,
                                   name='spatial_network')
    spatial_genes <- silhouetteRank(giotto)
    ext_spatial_genes <-spatial_genes[1:2000,]$gene
    HMRF_spatial_genes <- doHMRF(gobject=giotto,
                                 expression_values='scaled',
                                 spatial_genes=ext_spatial_genes,
                                 k=k,
                                 spatial_network_name="spatial_network",
                                 betas=c(0, 10, 5),
                                 output_folder=output)

    giotto <- addHMRF(gobject=giotto,
                      HMRFoutput=HMRF_spatial_genes,
                      k=3,
                      betas_to_add=c(0, 10, 20, 30, 40),
                      hmrf_name='HMRF')
    giotto <- giotto@cell_metadata


    t <- Sys.time() - s
    ftime <- paste0(fileTag[i],",",as.numeric(t),",",units(t),"\n")
    cat(ftime,file=time,append=TRUE)

    ari <- adjustedRandIndex(tmpPred$territory,tmpPred$HMRF_k3_b.40)
    vi <- vi.dist(tmpPred$territory,tmpPred$HMRF_k3_b.40)

    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)

    fileOut <- paste0(output,"/Giotto_",fileTags[i])
    write.table(ves,file =fileOut,sep =",",quote=F)
    rm(giotto); gc()

    frem <- list.files(pattern =".txt")
    for(f in frem){file.remove(f)}
    unlink("result.spatial.zscore", recursive =TRUE)
}
