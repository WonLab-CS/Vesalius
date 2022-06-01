#-----------------------------/Simulated Data/---------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to data simumlation in high res data sets
# We used slide-seq data as basis for data simulation.
#------------------------------------------------------------------------------#

library(Seurat)




library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(dplyr)


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
        simulated$cells[simulated$x >= xlims[i] & simulated$x <= xlims[i+1]] <- i
    }
    zeroes <- which(simulated$cells == 0)
    if(length(zeroes)>0){
        for(i in zeroes){
            simulated$cells[i] <- sample(seq(1,n_celltypes),1)
        }
    }
    #simulated <- globaliseTerritories(simulated)
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
    #simulated <- globaliseTerritories(simulated)
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

simulateCells <- function(c,s,co){
    cells <- c
    sim <- s
    counts <- co
    #--------------------------------------------------------------------------#
    # Assuming that we receive a data frame with barcodes and the cell types
    # associated to each cell type
    #--------------------------------------------------------------------------#
    if(length(unique(sim$cells)) > length(unique(cells$celltype))){
      tmp <- sample(unique(cells$celltype), size = length(unique(cells$celltype)),replace = FALSE)
      tmp <- c(tmp,sample(unique(cells$celltype), size = length(unique(sim$cells)) - length(unique(cells$celltype)),replace = FALSE))

      sim$celltype <- tmp[sim$cells]
      sim$barcodes <- "none"
    } else {
      sim$celltype <- sample(unique(cells$celltype), size = length(unique(sim$cells)),replace = FALSE)[sim$cells]
      sim$barcodes <- "none"
    }

    for(i in seq_along(unique(sim$celltype))){
        barcodes <- cells$barcodes[cells$celltype == unique(sim$celltype)[i]]
        simBarcodes <- sim$celltype == unique(sim$celltype)[i]
        sim$barcodes[simBarcodes] <- sample(barcodes,size = sum(simBarcodes), replace = TRUE)
    }
    sim$simBarcode <- paste0("bar_",seq_len(nrow(sim)))
    return(sim)

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
bar_cells <- bar_cells[bar_cells$celltype %in% names(table(bar_cells$celltype))[table(bar_cells$celltype)>50],]

brainCounts <- read.table("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz", header = TRUE )
rownames(brainCounts) <- brainCounts[,1]
brainCounts <- brainCounts[,-1]

# PARAMS
n_runs <- 10
n_territories <- c(3,3,3,3,3,3,3,10)
n_cells <- c(3,3,3,4,4,5,5,3)
method <- c("pure","uni","exp","uni","exp","uni","exp","dot")

## Let's run this bad boy or girl or which ever gender makes you a happy human
## we are just going to create data sets that will be used for each tool
input <- "/home/pcnmartin/Vesalius"
output <- paste0(input,"/Simulation")
if(!dir.exists(output)){
    dir.create(output)
}

for(i in seq_len(n_runs)){
    for(j in seq_along(n_territories)){
        sim <- switch(method[j],
                      "pure" = pureDist(n_cells[j],n_territories[j]),
                      "uni" = uniformDist(n_cells[j],n_territories[j]),
                      "exp" = expDist(n_cells[j],n_territories[j]),
                      "dot" = dottedDist(n_cells[j],n_territories[j]))
        sim <- simulateCells(bar_cells,sim, brainCounts)
        write.table(sim,
                    file = paste0(output,"/Simulation_",method[j],"_",n_territories[j],"_",n_cells[j],"_",i,".csv"),
                    sep =",",
                    quote = FALSE,
                    row.names = FALSE)
        rm(sim);gc()
    }
}
