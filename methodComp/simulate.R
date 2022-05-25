#-----------------------------/Simulated Data/---------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to data simumlation in high res data sets
# We used slide-seq data as basis for data simulation.
#------------------------------------------------------------------------------#
library(vesalius)
library(BayesSpace)
library(Seurat)
library(Giotto)


library(imagerExtra)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tvR)
library(sp)
library(grid)
library(SingleCellExperiment)

#library(devtools)
#install_github("patrickCNMartin/BayesSpace")
## This contain a modified version of BayesSpace that accomodates Slide-seq data

# Standard libraries - just in case you are a command line fiend
library(utils)
library(stats)
library(graphics)
library(grDevices)
### Set seed
set.seed(1)

#------------------------------/ Functions /-----------------------------------#

globaliseTerritories <- function(sim){

    ter <- paste0(sim$ter,"_", sim$cells)
    allTer <- unique(ter)
    ter <- seq_along(allTer)[match(ter,allTer)]
    sim$cells <- ter
    return(sim)

}

pureDist <- function(n_celltypes,n_territories, xlim = c(1,3000), ylim=c(1,3000),n_cells=6000){
    #--------------------------------------------------------------------------#
    # First lets set up the the coordinates dor each tissue.
    # I think I will only do the split along the x axis
    #--------------------------------------------------------------------------#
    if(n_celltypes!=n_territories){
        stop("Number of cells don't match number territories in pure Dist")

    }
    xlims <- seq(1,max(xlim),length.out = n_territories + 1)

    simulated <- data.frame("x" = runif(n = n_cells,min = min(xlim),max = max(xlim)),
                            "y" = runif(n = n_cells,min = min(ylim), max = max(ylim)),
                            "ter" = rep(0,n_cells),
                            "cells" = rep(0,n_cells))

    #--------------------------------------------------------------------------#
    # Set territories
    # There is only one cell per territory
    #--------------------------------------------------------------------------#
    for(i in seq(1,length(xlims)-1)){
        simulated$ter[simulated$x >= xlims[i] & simulated$x <= xlims[i+1]] <- i
        simulated$cells[simulated$x >= xlims[i] & simulated$x <= xlims[i+1]] <- 1
    }
    zeroes <- which(simulated$cells == 0)
    if(length(zeroes)>0){
        for(i in zeroes){
            simulated$cells[i] <- sample(seq(1,n_celltypes),1)
        }
    }
    simulated <- globaliseTerritories(simulated)
    return(simulated)

}


uniformDist <- function(n_celltypes,n_territories,xlim = c(1,3000), ylim=c(1,3000),n_cells=6000){
    #--------------------------------------------------------------------------#
    # get lims and set up df
    #--------------------------------------------------------------------------#
    xlims <- seq(1,max(xlim),length.out = n_territories + 1)

    simulated <- data.frame("x" = runif(n = n_cells,min = min(xlim),max = max(xlim)),
                            "y" = runif(n = n_cells,min = min(ylim), max = max(ylim)),
                            "ter" = rep(0,n_cells),
                            "cells" = rep(0,n_cells))
    #--------------------------------------------------------------------------#
    # Now we assign territories and cells by sampling idices
    #--------------------------------------------------------------------------#
    for(i in seq(1, length(xlims)-1)){
        simulated$ter[simulated$x >= xlims[i] & simulated$x <= xlims[i+1]] <- i
        tmp <- which(simulated$ter == i)

        for(j in seq(1,n_celltypes)){
            ter <- sample(x=tmp,size = floor(sum(simulated$ter == i)/n_celltypes), replace=FALSE)
            simulated$cells[ter] <- j
            tmp <- tmp[!tmp %in% ter]
        }
    }
    #--------------------------------------------------------------------------#
    # Due to rounding there might still be some "0" valued territories
    # so we will randomly assign a cell and territory to those ones
    #--------------------------------------------------------------------------#
    zeroes <- which(simulated$cells == 0)
    if(length(zeroes)>0){
        for(i in zeroes){
            simulated$cells[i] <- sample(seq(1,n_celltypes),1)
        }
    }
    simulated <- globaliseTerritories(simulated)
    return(simulated)

}


expDist <- function(n_celltypes,n_territories,xlim = c(1,3000), ylim=c(1,3000),n_cells=6000){
    divs <- exp(1)^seq(1,n_celltypes)
    divs <- round(divs/sum(divs), digits = 2)
    mat <- matrix(0,ncol =n_territories, nrow = n_celltypes)
    for(i in seq_len(n_celltypes)){
        mat[i,] <- seq(divs[i], divs[1+n_celltypes-i], length.out = n_territories)

    }
    xlims <- seq(1,max(xlim),length.out = n_territories + 1)

    simulated <- data.frame("x" = runif(n = n_cells,min = min(xlim),max = max(xlim)),
                            "y" = runif(n = n_cells,min = min(ylim), max = max(ylim)),
                            "ter" = rep(0,n_cells),
                            "cells" = rep(0,n_cells))
    for(i in seq(1, length(xlims)-1)){
        simulated$ter[simulated$x >= xlims[i] & simulated$x <= xlims[i+1]] <- i
        tmp <- which(simulated$ter == i)
        for(j in seq(1,n_celltypes)){

            ter <- sample(x=tmp,size = floor(sum(simulated$ter == i)* mat[j,i]), replace=FALSE)
            simulated$cells[ter] <- j
            tmp <- tmp[!tmp %in% ter]
          }
    }
    zeroes <- which(simulated$cells == 0)
    if(length(zeroes)>0){
        for(i in zeroes){
            simulated$cells[i] <- sample(seq(1,n_celltypes),1)
        }
    }
    simulated <- globaliseTerritories(simulated)
    return(simulated)
}

dottedDist <- function(n_celltypes,
                       n_territories,
                       randomise = TRUE,
                       xlim = c(1,3000),
                       ylim=c(1,3000),
                       n_cells=6000,
                       radius = c(50, 500)){
    #--------------------------------------------------------------------------#
    # we first create a back ground that will serve as a canvas
    # this back ground can contain more than one cell type
    #--------------------------------------------------------------------------#
    simulated <- data.frame("x" = runif(n = n_cells,min = min(xlim),max = max(xlim)),
                            "y" = runif(n = n_cells,min = min(ylim), max = max(ylim)),
                            "ter" = rep(1,n_cells),
                            "cells" = rep(0,n_cells))
    seeds <- simulated[sample(seq(1,n_cells), size = n_territories -1 , replace= FALSE),]
    if(randomise){
        rad <- sample(seq(radius[1L],radius[2L],l = n_territories -1))
        for(i in seq(1, length(rad))){
            simulated$ter[sqrt(((simulated$x - seeds$x[i])^2 + (simulated$y - seeds$y[i])^2)) <= rad[i]] <- i+1
            tmp <- which(simulated$ter == i+1)
            types <- sample(seq(1,n_celltypes),size = 1)
            for(j in seq(1,types)){

                ter <- sample(x=tmp,size = floor(sum(simulated$ter == i+1)/types), replace=FALSE)
                simulated$cells[ter] <- j
                tmp <- tmp[!tmp %in% ter]
              }
        }
        #----------------------------------------------------------------------#
        # Assigning random number of cell types to background
        # Okay this is a bit messy - I agree
        #----------------------------------------------------------------------#
        types <- sample(seq(1,n_celltypes),size = 1)
        tmp <- which(simulated$ter == 1)
        for(j in seq(1,types)){

            ter <- sample(x=tmp,size = floor(sum(simulated$ter == 1)/types), replace=FALSE)
            simulated$cells[ter] <- j
            tmp <- tmp[!tmp %in% ter]
        }

    } else {
        rad <- seq(radius[1L],radius[2L], l = n_territories -1)
        for(i in seq(1, length(rad))){
            simulated$ter[sqrt(((simulated$x - seeds$x[i])^2 + (simulated$y - seeds$y[i])^2)) <= rad[i]] <- i+1
            tmp <- which(simulated$ter == i+1)

            for(j in seq(1,n_celltypes)){
                ter <- sample(x=tmp,size = round(sum(simulated$ter == i+1)/n_celltypes), replace=FALSE)
                simulated$cells[ter] <- j
                tmp <- tmp[!tmp %in% ter]
              }
        }

        tmp <- which(simulated$ter == 1)
        for(j in seq(1,types)){

            ter <- sample(x=tmp,size = round(sum(simulated$ter == 1)/n_celltypes), replace=FALSE)
            simulated$cells[ter] <- j
            tmp <- tmp[!tmp %in% ter]
        }
    }
    zeroes <- which(simulated$cells == 0)
    if(length(zeroes)>0){
        for(i in zeroes){
            simulated$cells[i] <- sample(seq(1,n_celltypes),1)
        }
    }
    simulated <- globaliseTerritories(simulated)
    return(simulated)
}

simulateCells <- function(cells,
                          method = c("pure","uni","exp","dot")){
    return(NULL)

}

reSample <- function(counts,n_c,d_t,celltypes,cells){
  #----------------------------------------------------------------------------#
  # for some reason some cell compositions just don't work and I have no idea
  # why so I will reSample cells and try again
  # Hopefully it won't lead to anything weird
  #----------------------------------------------------------------------------#
  tmpCells <- celltypes[sample(seq_along(celltypes),n_c)]
  sim <- simulateCells(n_t=3,
                         n_c=n_c,
                         d_t=d_t,
                         n_cells=6000,
                         xlims=c(1,3000),
                         ylims=c(1,3000),
                         celltypes = tmpCells,
                         cells = cells)
  counts <- counts[,sim$barcodes]
  #----------------------------------------------------------------------------#
  # Rename barcodes to avoid potential duplicated names
  #----------------------------------------------------------------------------#
  colnames(counts) <- paste0("bar_",seq_len(ncol(counts)))
  rownames(counts)<- rownames(brainCounts)
  sim$barcodes <- paste0("bar_",seq_len(nrow(sim)))
  rownames(sim) <- sim$barcodes
  #----------------------------------------------------------------------------#
  # Running sims with various tools
  # Just coercing seurat to use slide seq data
  #----------------------------------------------------------------------------#
  ss <- new(Class = 'SlideSeq',
            assay = "Spatial",
            coordinates = sim[,c("xcoord","ycoord")])

  rownames(ss@coordinates) <- sim$barcodes

  st <- CreateSeuratObject(counts, assay ="Spatial")
  ss <- ss[Cells(x = st)]
  DefaultAssay(object = ss) <- "Spatial"
  st[["slice1"]] <- ss
  st <- AddMetaData(st,metadata = sim$territory,col.name = "Territory")

  seu <- tryCatch(suppressWarnings(SCTransform(st,assay = "Spatial")),
                  error = function(cond){
                      return(NULL)
                  })

  if(is.null(seu)){
      seu <- reSample(counts,n_c,d_t,celltypes,cells)
      return(seu)
  } else{
      return(list("cells" = tmpCells,"sim"=sim,"st" =st,"counts" =counts))
  }

}

#------------------------------------------------------------------------------#

# Loading RCTD cell type Annotations and getting homotypic beads
ref <- readRDS("/isilonsund/NextGenSeqData/project-data/hkim/ST_project/Slide-seq/Slide-seq_hippocampus.rds")
cells <- unique(unlist(strsplit(ref@meta.data$celltype, "\\+")))
## NOTE I will only run 3 cells types at the moment.
## Could run more for sure but for the sake of simplicity
method <- "homotypic"

## If we only want to homotypic cells
celltypes <- paste0(cells,"+",cells)

bar_cells <- ref@meta.data[ref@meta.data$celltype %in% celltypes,]
bar_cells <- data.frame("barcodes"=rownames(bar_cells),"celltype"=bar_cells$celltype)

brainCounts <- read.table("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz", header = TRUE )
rownames(brainCounts) <- brainCounts[,1]
brainCounts <- brainCounts[,-1]



## Let's run this bad boy or girl or which ever gender makes you a happy human
input <- "/home/pcnmartin/Vesalius"
output <- paste0(input,"/Simulation")
if(!dir.exists(output)){
    dir.create(output)
}
#------------------------------------------------------------------#
# First simulate territories
# n_t = number of territories
# n_c = number of cell types
# d_t = distribution type
# n_cells = number of cells to simulate !!! Will likely be less due to rounding
# xlims and ylims = spatial coodindate boundaries
# celltypes = cell types used for simulation
# cells = cell info
# tmpCells = which cells were sampled
#------------------------------------------------------------------#
time <- vector("list", 10)

for(k in seq(1,10)){

  n_c <-c(9,12,15,3,3,4,5)
  d_t <-c(rep("uniform",3),"pure",rep("exp",3))
  tmpCells <- list()
  fileTags <- c()
  simList <- list()

  for(j in seq_along(n_c)){
      fileTags <- c(fileTags,paste0("nt3_nc",n_c[j],"_",d_t[j],"_",method,".csv"))
      tmpCells[[j]] <- celltypes[sample(seq_along(celltypes),n_c[j])]

      simList[[j]] <- simulateCells(n_t=3,
                           n_c=n_c[j],
                           d_t=d_t[j],
                           n_cells=6000,
                           xlims=c(1,3000),
                           ylims=c(1,3000),
                           celltypes = tmpCells[[j]],
                           cells = bar_cells)
  }


  time[[k]] <- vector("list",4)
  names(time[[k]]) <- c("Seurat","Vesalius","BayesSpace","Giotto")
  for(i in seq_along(simList)){

  #for(i in seq(6,length(simList))){
    #----------------------------------------------------------------------------#
    # Lets sim
    #----------------------------------------------------------------------------#
    sim <- simList[[i]]
    #----------------------------------------------------------------------------#
    # Next get counts from count df
    #----------------------------------------------------------------------------#
    counts <- brainCounts[,sim$barcodes]

    #----------------------------------------------------------------------------#
    # Rename barcodes to avoid potential duplicated names
    #----------------------------------------------------------------------------#
    colnames(counts) <- paste0("bar_",seq_len(ncol(counts)))
    rownames(counts)<- rownames(brainCounts)
    sim$barcodes <- paste0("bar_",seq_len(nrow(sim)))
    rownames(sim) <- sim$barcodes
    #----------------------------------------------------------------------------#
    # Running sims with various tools
    # Just coercing seurat to use slide seq data
    #----------------------------------------------------------------------------#
    ss <- new(Class = 'SlideSeq',
              assay = "Spatial",
              coordinates = sim[,c("xcoord","ycoord")])

    rownames(ss@coordinates) <- sim$barcodes

    st <- CreateSeuratObject(counts, assay ="Spatial")
    ss <- ss[Cells(x = st)]
    DefaultAssay(object = ss) <- "Spatial"
    st[["slice1"]] <- ss
    st <- AddMetaData(st,metadata = sim$territory,col.name = "Territory")
    #----------------------------------------------------------------------------#
    # Running Seurat
    # We first test if this sim works - some of them dont work
    # and I'll be honnest I have no idea why...
    # Seurat and Giotto seem to crash on these ones. I will resample until
    # a Seurat run goes to completion and update the cells and sim list
    # I don't see any issue in the data sim output it looks all the same as
    # the other ones...And it is not related to the number of unique barcodes
    # used for sampling.
    #----------------------------------------------------------------------------#
    tmp <- tryCatch(suppressWarnings(SCTransform(st,assay = "Spatial")),
                  error = function(cond){
                      return(NULL)
                  })
    if(is.null(tmp)){
        tmp <- reSample(brainCounts,n_c[i],d_t[i],celltypes,bar_cells)
        sim <- tmp$sim
        simList[[i]] <- sim
        counts <- tmp$counts
        tmpCells[[i]] <- tmp$cells
        st <- tmp$st

    }
    s <- Sys.time()
    seu <- SCTransform(st,assay ="Spatial") %>%
           RunPCA(dims =1:30) %>%
           RunUMAP(dims =1:30) %>%
           FindNeighbors()
    closestRes <- c()
    for(res in seq(0.1,1.4,by = 0.05)){
        seu <- FindClusters(seu, resolution = res, verbose = FALSE)
        cl <- length(levels(FetchData(seu, c("seurat_clusters"))$seurat_clusters))
        closestRes <- c(closestRes,abs(n_c[i] - cl))

    }
    res <- seq(0.1,1.4,by = 0.05)[closestRes == min(closestRes)]
    seu <- FindClusters(seu, resolution =res[1L] ,
                        verbose = FALSE)

    seuData <- FetchData(seu, c("UMAP_1","UMAP_2","seurat_clusters")) %>%
               group_by(seurat_clusters) %>%
               mutate(xs = mean(UMAP_1), ys = mean(UMAP_2))
    coordSeu <- GetTissueCoordinates(seu)
    seuData <- cbind(seuData,coordSeu[,c("x","y")])
    time[[k]]$Seurat <- Sys.time() - s
    fileOut <- paste0(output,"/Seurat_Sim_rep",k,"_",fileTags[i])
    write.table(seuData,file =fileOut,sep =",",quote=F)
    rm(seuData,seu); gc()
    #----------------------------------------------------------------------------#
    # Running Vesalius
    #----------------------------------------------------------------------------#
    #colDepth <- seq(n_c[i],3)
    s <- Sys.time()
    colDepth <- c(81,27,9,3)
    iter <- 15
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
          isolateTerritories.array(captureRadius=0.1,minBar=10) %>%
          filter(tile==1) %>%
          distinct(barcodes, .keep_all =TRUE)
    time[[k]]$Vesalius <- Sys.time() - s
    fileOut <- paste0(output,"/Vesalius_Sim_rep",k,"_",fileTags[i])
    write.table(ves,file =fileOut,sep =",",quote=F)
    rm(ves); gc()
    #----------------------------------------------------------------------------#
    # Running BayesSpace
    #----------------------------------------------------------------------------#

    sce <- SingleCellExperiment(list(counts=as(as.matrix(counts), "dgCMatrix")))
    coord <- sim[colnames(sce),c("ycoord","xcoord")]
    colnames(coord)<- c("row","col")
    colData(sce) <- DataFrame(coord)

    genes <- data.frame(rownames(counts), rownames(counts))
    colnames(genes) <- c("gene_id","gene_name")
    rownames(genes) <- rownames(counts)
    rowData(sce) <- genes
    meta <- list(sample ="BayesSpaceSim", dataset = "Simulated",
                 BayesSpace.data = list(platform = "SS",is.enhanced=FALSE))
    metadata(sce) <- meta
    s <- Sys.time()
    sce <- spatialPreprocess(sce,platform="SS",skip.PCA=FALSE,
                             n.PCs=30, n.HVGs=2000, log.normalize=TRUE)
    #sce <- qTune(sce, qs=seq(10,30), platform="SS", d=30)

    sce <- spatialCluster(sce,q=3,platform = "SS")
    bayes <- as.data.frame(colData(sce))
    time[[k]]$BayesSpace <- Sys.time() - s
    fileOut <- paste0(output,"/BayesSpace_Sim_rep",k,"_",fileTags[i])
    write.table(bayes,file =fileOut,sep =",",quote=F)
    rm(bayes,sce);gc()
    #----------------------------------------------------------------------------#
    # Running Giotto
    #----------------------------------------------------------------------------#
    s <- Sys.time()
    instruc <- createGiottoInstructions(save_plot=F, show_plot=F, save_dir=output)
    giotto <- createGiottoObject(raw_exprs = counts,
                                       spatial_locs = sim[,c("xcoord","ycoord")],
                                       instructions = instruc,
                                       cell_metadata = sim[,c("xcoord","ycoord","territory")])

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
                                   k=3,
                                   maximum_distance_knn=400,
                                   name='spatial_network')
    spatial_genes <- silhouetteRank(giotto)
    ext_spatial_genes <-spatial_genes[1:2000,]$gene
    HMRF_spatial_genes <- doHMRF(gobject=giotto,
                                 expression_values='scaled',
                                 spatial_genes=ext_spatial_genes,
                                 k=3,
                                 spatial_network_name="spatial_network",
                                 betas=c(0, 10, 5),
                                 output_folder=output)

    giotto <- addHMRF(gobject=giotto,
                      HMRFoutput=HMRF_spatial_genes,
                      k=3,
                      betas_to_add=c(0, 10, 20, 30, 40),
                      hmrf_name='HMRF')
    giotto <- giotto@cell_metadata
    time[[k]]$Giotto <- Sys.time() - s
    fileOut <- paste0(output,"/Giotto_Sim_rep",k,"_",fileTags[i])
    write.table(giotto,file =fileOut,sep =",",quote=F)
    # Removing all the stuff giotto outputs...
    rm(giotto);gc()
    frem <- list.files(pattern =".txt")
    for(f in frem){file.remove(f)}
    unlink("result.spatial.zscore", recursive =TRUE)
  }
  ## save at the very end incase there are some re-runs
  ## I don't know why some samples don't work...
  save(tmpCells,simList, file = paste0("Simulated_cells_rep",k,"_",method,".Rda"))
}
save(time,file= paste0("SimulationTimes_",method,".Rda"))
