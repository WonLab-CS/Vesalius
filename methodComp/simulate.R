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
pureDist <- function(n_celltypes,n_territories){

    if(n_celltypes!=n_territories){
        stop("Number of cells don't match number territories in pure Dist")
        return(NULL)
    } else {
        dist <- diag(n_celltypes)
        return(dist)

    }
}

##A bit odd as this need to be able to split cells between territories
uniformDist <- function(n_celltypes,n_territories){
      n_groups <- n_celltypes %/% n_territories

      if(n_groups ==0){
          stop("Not enough different cell types to populate each territory")
      } else {
          mat <- matrix(0,ncol =n_territories,nrow=n_territories*n_groups)
          locs <- seq(1,n_groups)
          for(i in seq_len(n_territories)){
              mat[locs,i] <- 1/ n_groups
              locs <- locs + n_groups
          }
      }
      return(mat)

}


expDist <- function(n_celltypes,n_territories){
    divs <- exp(1)^seq(1,n_celltypes)
    divs <- divs/sum(divs)
    mat <- matrix(0,ncol =n_territories, nrow = n_celltypes)
    for(i in seq_len(n_celltypes)){
        mat[i,] <- seq(divs[i], divs[1+n_celltypes-i], length.out = n_territories)

    }
    return(mat)
}

simulateCells <- function(n_t,n_c,d_t,n_cells,xlims,ylims,celltypes,cells){

  #------------------------------------------------------------------#
  # we create the territory splits
  # we only going to split along the x axis for now
  # It shouldnt change to much tbh
  # creating coordinates for each territory
  # Here im just going to randomly distribute all cells
  # Not sure if it matters how many cells there are in each territory
  # There should be roughly the same amount anyway
  #------------------------------------------------------------------#
  bounds <- round(seq(xlims[1],xlims[2], length.out = n_t + 1))
  coord <- data.frame("xcoord" = jitter(sample(seq(xlims[1],xlims[2]),size = n_cells, replace =T)),
                      "ycoord" = jitter(sample(seq(ylims[1],ylims[2]),size = n_cells, replace =T)))
  #------------------------------------------------------------------#
  # Next we create cell sampling templates
  # No exp for now... I have no idea how to do that yet
  #------------------------------------------------------------------#

  dist <- switch(d_t,
                "pure" = pureDist(n_c,n_t),
                "uniform" = uniformDist(n_c,n_t),
                "exp" = expDist(n_c,n_t))
  if(is.null(dist))return(NULL)
  #------------------------------------------------------------------#
  # Using these proportion templates to sample barcodes
  # and generate coordinates
  # for now we will used imposed cell types.. you can just do the sampling
  # out of function
  #------------------------------------------------------------------#
  #ct <- sample(celltypes,n_c, replace = FALSE)

  territory_barcodes <- vector("list", n_t)
  for(l in seq_len(ncol(dist))){
      tmp <- c()
      for(m in seq_len(nrow(dist))){
          r <- nrow(filter(coord,xcoord >= bounds[l] & xcoord < bounds[l+1]))

          n_points <- r * dist[m,l]
          if(m == nrow(dist)){
              n_points <- r - length(tmp)
          }
          tmp <- c(tmp,sample(cells$barcodes[cells$celltype==celltypes[m]],n_points, replace =TRUE))
      }

      territory_barcodes[[l]] <- data.frame("barcodes" = tmp,
                                            filter(coord,xcoord >= bounds[l] & xcoord < bounds[l+1]),
                                            "territory" = rep(l,length(tmp)))
  }

  territory_barcodes <- do.call("rbind",territory_barcodes)

  return(territory_barcodes)
}



#------------------------------------------------------------------------------#

# Loading RCTD cell type Annotations and getting homotypic beads
ref <- readRDS("/isilonsund/NextGenSeqData/project-data/hkim/ST_project/Slide-seq/Slide-seq_hippocampus.rds")
cells <- unique(unlist(strsplit(ref@meta.data$celltype, "\\+")))
## NOTE I will only run 3 cells types at the moment.
## Could run more for sure but for the sake of simplicity
homo <- paste0(cells,"+",cells)

homo_cells <- ref@meta.data[ref@meta.data$celltype %in% homo,]
homo_cells <- data.frame("barcodes"=rownames(homo_cells),"celltype"=homo_cells$celltype)

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


n_c <-c(9,12,15,3,3,4,5)
d_t <-c(rep("uniform",3),"pure",rep("exp",3))
tmpCells <- list()
fileTags <- c()
simList <- list()

for(i in seq_along(n_c)){
    fileTags <- c(fileTags,paste0("nt3_nc",n_c[i],"_",d_t[i],".csv"))
    tmpCells[[i]] <- homo[sample(seq_along(homo),n_c[i])]

    simList[[i]] <- simulateCells(n_t=3,
                         n_c=n_c[i],
                         d_t=d_t[i],
                         n_cells=6000,
                         xlims=c(1,3000),
                         ylims=c(1,3000),
                         celltypes = tmpCells[[i]],
                         cells = homo_cells)
}

save(tmpCells, file = "Simulated_cells.Rda")


for(i in seq_along(simList)){
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
  #----------------------------------------------------------------------------#
  seu <- SCTransform(st,assay = "Spatial")%>%
         RunPCA(dims =1:30) %>%
         RunUMAP(dims =1:30) %>%
         FindNeighbors()
  for(k in seq(0.1,1,by = 0.05)){
      seu <- FindClusters(seu, resolution = k, verbose = FALSE)
      cl <- length(levels(FetchData(seu, c("seurat_clusters"))$seurat_clusters))
      closestRes <- abs(n_c[i] - cl)

  }
  seu <- FindClusters(seu, resolution = seq(0.1,1,by = 0.05)[closestRes == min(closestRes)],
                      verbose = FALSE)

  seuData <- FetchData(seu, c("UMAP_1","UMAP_2","seurat_clusters")) %>%
             group_by(seurat_clusters) %>%
             mutate(xs = mean(UMAP_1), ys = mean(UMAP_2))
  coordSeu <- GetTissueCoordinates(seu)
  seuData <- cbind(seuData,coordSeu[,c("x","y")])
  fileOut <- paste0(output,"/Seurat_Sim_",fileTags[i])
  write.table(seuData,file =fileOut,sep =",",quote=F)
  rm(seuData,seu); gc()
  #----------------------------------------------------------------------------#
  # Running Vesalius
  #----------------------------------------------------------------------------#
  ves <- NormalizeData(st) %>%
         FindVariableFeatures(nfeatures=2000) %>%
         ScaleData()%>%
         rgbPCA()%>%
         buildImageArray(resolution=50,filterThreshold=1,keep_edge=T,cores =5)%>%
         equalizeHistogram(sleft =2.5,sright=2.5)%>%
         regulariseImage(lambda = 10,niter=200)%>%
         iterativeSegmentation.array(colDepth=3,
                                     smoothIter = 30,
                                     method = c("iso","box"),
                                     sigma=1.5,box = 30,
                                     useCenter = T) %>%
        isolateTerritories.array(captureRadius=0.1,minBar=10) %>%
        filter(tile==1) %>%
        distinct(barcodes, .keep_all =TRUE)
  fileOut <- paste0(output,"/Vesalius_Sim_",fileTags[i])
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

  sce <- spatialPreprocess(sce,platform="SS",skip.PCA=FALSE,
                           n.PCs=30, n.HVGs=2000, log.normalize=TRUE)
  #sce <- qTune(sce, qs=seq(10,30), platform="SS", d=30)

  sce <- spatialCluster(sce,q=3,platform = "SS")
  bayes <- as.data.frame(colData(sce))
  fileOut <- paste0(output,"/BayesSpace_Sim_",fileTags[i])
  write.table(bayes,file =fileOut,sep =",",quote=F)
  rm(bayes,sce);gc()
  #----------------------------------------------------------------------------#
  # Running Giotto
  #----------------------------------------------------------------------------#
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
  fileOut <- paste0(output,"/Giotto_Sim_",fileTags[i])
  write.table(giotto,file =fileOut,sep =",",quote=F)
  # Removing all the stuff giotto outputs...
  rm(giotto);gc()
  frem <- list.files(pattern =".txt")
  for(i in frem){file.remove(i)}
  unlink("result.spatial.zscore", recursive =TRUE)
}
