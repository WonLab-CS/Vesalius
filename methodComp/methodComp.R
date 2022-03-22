#-----------------------------/Method comparison/------------------------------#
#------------------------------------------------------------------------------#
# This file contains all the code used to compare methods - Table of content:
# * Time Plots For Visum
# * ARI for Visium
# * Example data set in visium for viz - using the common example 151673
# * Plotting slide-seqV2 output for BayesSpace and Seurat
# * Simulation scores for all regimes
# * Simulation Run times
# * Concat cell types use in each sim round
# * Simulation Plots for all regimes
# * Simulation plots for the manuscript
# * Concat all diff gene expression files into a single file
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Time Plots for Visium Figure S18
#------------------------------------------------------------------------------#
library(ggplot2)
files <- list.files(pattern = "time.Rda")

convertToSecs <- function(x,u){
    m <- which(u == "mins")
    h <- which(u == "hours")

    x[m] <- x[m] * 60
    x[h] <- x[h] * 60 * 60
    return(x)
}

df <- data.frame("time" = numeric(), "method" = character())
for(i in seq_along(files)){
    tmp <- get(load(files[i]))
    units <- sapply(tmp, units)
    time <- unlist(tmp)
    name <- gsub("_time.Rda","", files[i])
    df <- rbind(df, data.frame(convertToSecs(time,units), rep(name,12)))

}
colnames(df) <- c("time" , "method" )
df$method <- as.factor(df$method)
cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,4,7,8)]

g <- ggplot(df, aes(x = method, y= time,fill = method)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size =15)) +
  #guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(x = "", y = "Time in Seconds",fill = " ")

pdf("MethodComp_time.pdf",width = 7, height = 3.5)
print(g)
dev.off()


#------------------------------------------------------------------------------#
# Adjusted Rand Index for Visium - Figure S18
#------------------------------------------------------------------------------#
library(mclust)
library(spatialLIBD)
library(dplyr)
library(Giotto)
library(Seurat)
library(BayesSpace)

labels <- fetch_data(type = 'sce')


labels <- colData(labels)
labels <- as.data.frame(labels)
labels <- split(labels, labels$sample_name)

sample <- names(labels)

input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)


performance <- data.frame("method" = character(), "ARI" = numeric())
for(i in seq_along(input)){
    # get labels
    tmp_labels <- labels[names(labels) == sample[i]]
    tmp_labels <- tmp_labels[[1]]
    ## First Get Vesalius
    vesalius <- paste0(input[i],"/Vesalius/",sample[i],"_vesalius.rda")
    vesalius <- get(load(vesalius))
    vesalius <- filter(vesalius,tile == 1) %>% distinct(barcodes,.keep_all = TRUE)
    aligned <- match(vesalius$barcodes, tmp_labels$barcode)
    tmp <- as.character(tmp_labels$layer_guess)[aligned]
    vesalius <- vesalius$territory
    vesalius <- adjustedRandIndex(tmp,vesalius)
    ### Seurat
    seurat <- paste0(input[i],"/Seurat/",sample[i],"_seurat.rda")
    seurat <- get(load(seurat))
    seurat <- FetchData(seurat,c("seurat_clusters"))
    aligned <- match(rownames(seurat), tmp_labels$barcode)
    tmp <- as.character(tmp_labels$layer_guess)[aligned]
    seurat <- seurat$seurat_clusters
    seurat <- adjustedRandIndex(tmp,seurat)
    ## Bayes space
    bayes <- paste0(input[i],"/BayesSpace/",sample[i],"_BayesSpace.rda")
    bayes <- as.data.frame(colData(get(load(bayes))))
    aligned <- match(rownames(bayes), tmp_labels$barcode)
    tmp <- as.character(tmp_labels$layer_guess)[aligned]
    bayes <- bayes$spatial.cluster
    bayes <- adjustedRandIndex(tmp,bayes)

    ## giotto
    giotto <- paste0(input[i],"/Giotto/",sample[i],"_giotto.rda")
    giotto <- get(load(giotto))
    giotto <- giotto@cell_metadata
    aligned <- match(giotto$cell_ID, tmp_labels$barcode)
    tmp <- as.character(tmp_labels$layer_guess)[aligned]
    giotto <- unlist(giotto[,12])
    giotto <- adjustedRandIndex(tmp,giotto)

    performance <- rbind(performance, data.frame(
        "method" = c("Vesalius","Seurat","BayesSpace","Giotto"),
        "ARI" = c(vesalius,seurat,bayes,giotto)
    ))

}

performance$method <- as.factor(performance$method)
cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,4,7,8)]

g <- ggplot(performance, aes(x = method, y= ARI,fill = method)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size =15)) +
  #guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(x = "", y = "ARI",fill = " ")

pdf("MethodComp_ARI.pdf",width = 7, height = 3.5)
print(g)
dev.off()



#------------------------------------------------------------------------------#
# Example data set in visium for viz - using the common example 151673
# Figure S18
#------------------------------------------------------------------------------#
library(ggplot2)
library(patchwork)
library(dplyr)
library(spatialLIBD)
library(RColorBrewer)
library(Seurat)
library(BayesSpace)
library(Giotto)


labels <- fetch_data(type = 'sce')
labels <- colData(labels)
labels <- as.data.frame(labels)
labels <- filter(labels, sample_name == "151673") %>% select(row,col,layer_guess,layer_guess_reordered)

input <- list.dirs("~/group/visium/DLPFC_globus",recursive =F)
input <- input[9]

cols <- length(levels(labels$layer_guess_reordered))
pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
cols <- pal(cols)#[sample(1:cols)]


ground <- ggplot(labels, aes(col,row,col = layer_guess_reordered)) +
     geom_point(size = 2.3,alpha = 1)+
     theme_void()+
     scale_color_manual(values = cols)+
     theme(legend.text = element_text(size = 12),
           legend.title = element_text(size=12),
           plot.title = element_text(size =15),
           legend.position = "right",
           plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"))+
     labs(title = "Ground Truth", col = " ")

vesalius <- paste0(input,"/Vesalius/151673_vesalius.rda")
vesalius <- get(load(vesalius))
vesalius <- filter(vesalius,tile == 1) %>% distinct(barcodes,.keep_all = TRUE)
sorted_labels <- order(levels(as.factor(vesalius$territory)))
sorted_labels[length(sorted_labels)] <- "isolated"
vesalius$territory <- factor(vesalius$territory, levels = sorted_labels)

cols <- length(levels(vesalius$territory))
pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
cols <- pal(cols)[sample(1:cols)]


ves <- ggplot(vesalius, aes(x,y,col = territory)) +
     geom_point(size = 1.5,alpha = 1)+
     theme_void()+
     scale_color_manual(values = cols)+
     theme(legend.text = element_text(size = 12),
           legend.title = element_text(size=12),
           plot.title = element_text(size =15),
           legend.position = "right")+
     labs(title = "Vesalius", col = "Territory")

seurat <- paste0(input,"/Seurat/151673_seurat.rda")
seurat <- get(load(seurat))
coord <- GetTissueCoordinates(seurat)
cluster <- FetchData(seurat,c("seurat_clusters"))
seurat <- cbind(coord,cluster)
seurat$seurat_clusters <- as.factor(as.numeric(as.character(seurat$seurat_clusters))+1)

cols <- length(levels(seurat$seurat_clusters))
pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
cols <- pal(cols)#[sample(1:cols)]


seu <- ggplot(seurat, aes(imagecol,imagerow,col = as.factor(seurat_clusters))) +
     geom_point(size = 1.5,alpha = 1)+
     theme_void()+
     scale_color_manual(values = cols)+
     theme(legend.text = element_text(size = 12),
           legend.title = element_text(size=12),
           plot.title = element_text(size =15),
           legend.position = "right")+
     labs(title = "Seurat", col = "Cluster")

giotto <- paste0(input,"/Giotto/151673_giotto.rda")
giotto <- get(load(giotto))
giotto <- giotto@cell_metadata
cols <- length(unique(giotto$HMRF_k7_b.40))
pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
cols <- pal(cols)#[sample(1:cols)]


gio <- ggplot(giotto, aes(array_col,array_row,col = as.factor(HMRF_k7_b.40))) +
     geom_point(size = 1.5,alpha = 1)+
     theme_void()+
     scale_color_manual(values = cols)+
     theme(legend.text = element_text(size = 12),
           legend.title = element_text(size=12),
           plot.title = element_text(size =15),
           legend.position = "right")+
     labs(title = "Giotto", col = "Spatial Domain")

bayes <- paste0(input,"/BayesSpace/151673_BayesSpace.rda")
bayes <- as.data.frame(colData(get(load(bayes))))

cols <- length(unique(bayes$spatial.cluster))
pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
cols <- pal(cols)#[sample(1:cols)]


bay <- ggplot(bayes, aes(col,row,col = as.factor(spatial.cluster))) +
     geom_point(size = 1.5,alpha = 1)+
     theme_void()+
     scale_color_manual(values = cols)+
     theme(legend.text = element_text(size = 12),
           legend.title = element_text(size=12),
           plot.title = element_text(size =15),
           legend.position = "right")+
     labs(title = "BayesSpace", col = "Spatial Cluster")

layout <- "
ABC
ADE"
pdf("MethodComp_Sample.pdf",width = 16, height = 8)
ground + ves + seu + gio + bay + plot_layout(design = layout,width = c(2,1,1))
dev.off()

#------------------------------------------------------------------------------#
# Plotting slide seq output for BayesSpace and Seurat
# Figure 2b and Figure 2d
# Figure S3 to Figure S6
#------------------------------------------------------------------------------#
bayes <- get(load("~/group/slide_seqV2/BayesSpaceBenchMarking/Puck_200115_08_SSV2_BM.Rda"))
bayes$row <- as.numeric(bayes$row)
bayes$col <- as.numeric(bayes$col)

ter_col <- length(unique(bayes$spatial.cluster))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
b1 <- ggplot(data = bayes,aes(x = col,y=row, col = as.factor(spatial.cluster)))+
      geom_point(size = 0.2, alpha = 0.65)+
      theme_void() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "BayesSpace", col = "Spatial Cluster")
pdf("~/Vesalius/BayesSpaceBrain.pdf", width = 6, height=7)
b1
dev.off()
b1s <- ggplot(data = bayes,aes(x = col,y=row, col = as.factor(spatial.cluster)))+
      geom_point(size = 0.5, alpha = 0.65)+
      facet_wrap(~spatial.cluster)+
      theme_light() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "BayesSpace", col = "Spatial Cluster")
pdf("~/Vesalius/BayesSpaceBrainSplit.pdf", width = 16, height=18)
b1s
dev.off()



bayes <- get(load("~/group/slide_seqV2/BayesSpaceBenchMarking/Puck_190926_03_SSV2_BM.Rda"))
bayes$row <- as.numeric(bayes$row)
bayes$col <- as.numeric(bayes$col)

ter_col <- length(unique(bayes$spatial.cluster))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
b2 <- ggplot(data = bayes,aes(x = col,y=row, col = as.factor(spatial.cluster)))+
      geom_point(size = 0.2, alpha = 0.65)+
      theme_void() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "BayesSpace", col = "Spatial Cluster")
pdf("~/Vesalius/BayesSpaceEmbryo.pdf", width = 6, height=7)
b2
dev.off()
b2s <- ggplot(data = bayes,aes(x = col,y=row, col = as.factor(spatial.cluster)))+
      geom_point(size = 0.5, alpha = 0.65)+
      theme_light() +
      facet_wrap(~spatial.cluster)+
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "BayesSpace", col = "Spatial Cluster")
pdf("~/Vesalius/BayesSpaceEmbryoSplit.pdf", width = 16, height=18)
b2s
dev.off()

seurat <- get(load("~/group/slide_seqV2/SeuratBenchMarking/Puck_200115_08_SSV2_BM.Rda"))
ter_col <- length(unique(seurat$seurat_clusters))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
s1 <- ggplot(data = seurat,aes(x = x,y=y, col = as.factor(seurat_clusters)))+
      geom_point(size = 0.2, alpha = 0.65)+
      theme_void() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "Seurat", col = "Seurat Cluster")
pdf("~/Vesalius/SeuratBrain.pdf", width = 6, height=7)
s1
dev.off()
s1s <- ggplot(data = seurat,aes(x = x,y=y, col = as.factor(seurat_clusters)))+
      geom_point(size = 0.5, alpha = 0.65)+
      theme_light() +
      facet_wrap(~seurat_clusters)+
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "Seurat", col = "Seurat Cluster")
pdf("~/Vesalius/SeuratBrainSplit.pdf", width = 16, height=18)
s1s
dev.off()

seurat <- get(load("~/group/slide_seqV2/SeuratBenchMarking/Puck_190926_03_SSV2_BM.Rda"))
ter_col <- length(unique(seurat$seurat_clusters))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
s2 <- ggplot(data = seurat,aes(x = x,y=y, col = as.factor(seurat_clusters)))+
      geom_point(size = 0.2, alpha = 0.65)+
      theme_void() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "Seurat", col = "Seurat Cluster")
pdf("~/Vesalius/SeuratEmbryo.pdf", width = 6, height=7)
s2
dev.off()
s2s <- ggplot(data = seurat,aes(x = x,y=y, col = as.factor(seurat_clusters)))+
      geom_point(size = 0.5, alpha = 0.65)+
      theme_light() +
      facet_wrap(~seurat_clusters)+
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom")+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "Seurat", col = "Seurat Cluster")
pdf("~/Vesalius/SeuratEmbryoSplit.pdf", width = 16, height=18)
s2s
dev.off()

#------------------------------------------------------------------------------#
# Simulation scores for all regimes
# Figure 2g
#------------------------------------------------------------------------------#
files <- list.files(pattern = ".csv")
gt <- list.files(pattern = "homotypic.Rda")
gt <- gt[gt!= "SimulationTimes_homotypic.Rda"]

## reorder because 1 and 10
gt <- gt[c(1,3:10,2)]

library(mclust)
library(mcclust)
library(dplyr)
library(ggplot2)

performance <- data.frame("Method" = character(),
                          "Regime" = character(),
                          "Rep" = numeric(),
                          "ARI" = numeric(),
                          "viDist" = numeric())

for(i in seq_along(gt)){
    ## We can assign names to the simList as we know what they are
    load(gt[i])
    n_c <-c(9,12,15,3,3,4,5)
    d_t <-c(rep("uniform",3),"pure",rep("exp",3))
    names(simList) <- paste0("nt3_nc",n_c,"_",d_t)
    for(j in seq_along(simList)){
        tmpSim <- simList[[j]]
        tmpSim$barcodes <- paste0("bar_",seq_len(nrow(tmpSim)))
        fileTag <- sapply(strsplit(files,"_"),"[[",3)
        repsByRegime <- files[fileTag == paste0("rep",i) &
                              grepl(names(simList)[j],files)]
        for(k in seq_along(repsByRegime)){
            if(grepl("BayesSpace",repsByRegime[k])){
                tmpPred <- read.csv(repsByRegime[k],header=T)
                aligned <- match(rownames(tmpPred), tmpSim$barcodes)
                simTer <- as.character(tmpSim$territory)[aligned]
                ari <- adjustedRandIndex(tmpPred$spatial.cluster,simTer)
                vi <- vi.dist(tmpPred$spatial.cluster,simTer)
                tmpPred <- data.frame("BayesSpace",names(simList)[j],i,ari,vi)
                colnames(tmpPred)<- c("Method","Regime","Rep","ARI","viDist")
                performance <- rbind(performance,tmpPred)

            }else if(grepl("Giotto",repsByRegime[k])){
                tmpPred <- read.csv(repsByRegime[k],header=T)
                ari <- adjustedRandIndex(tmpPred$territory,tmpPred$HMRF_k3_b.40)
                vi <- vi.dist(tmpPred$territory,tmpPred$HMRF_k3_b.40)
                tmpPred <- data.frame("Giotto",names(simList)[j],i,ari,vi)
                colnames(tmpPred)<- c("Method","Regime","Rep","ARI","viDist")
                performance <- rbind(performance,tmpPred)

            }else if(grepl("Seurat",repsByRegime[k])){
                tmpPred <- read.csv(repsByRegime[k],header=T)
                #aligned <- match(tmpPred$barcodes, tmpSim$barcodes)
                #simTer <- as.character(tmpSim$territory)[aligned]
                ari <- adjustedRandIndex(tmpPred$seurat_clusters,tmpSim$territory)
                vi <- vi.dist(tmpPred$seurat_clusters,tmpSim$territory)
                tmpPred <- data.frame("Seurat",names(simList)[j],i,ari,vi)
                colnames(tmpPred)<- c("Method","Regime","Rep","ARI","viDist")
                performance <- rbind(performance,tmpPred)

            }else if(grepl("Vesalius",repsByRegime[k])){
                tmpPred <- read.csv(repsByRegime[k],header=T)
                aligned <- match(tmpPred$barcodes, tmpSim$barcodes)
                simTer <- as.character(tmpSim$territory)[aligned]
                ari <- adjustedRandIndex(tmpPred$territory,simTer)
                vi <- vi.dist(tmpPred$territory,simTer)
                tmpPred <- data.frame("Vesalius",names(simList)[j],i,ari,vi)
                colnames(tmpPred)<- c("Method","Regime","Rep","ARI","viDist")
                performance <- rbind(performance,tmpPred)

            } else {
                message("Yeah nah don't know what that is...")
            }
        }
    }
}



performance$Method <- as.factor(performance$Method)
cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,4,7,8)]

g <- ggplot(performance, aes(x = Method, y= ARI,fill = Method)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size =15)) +
  #guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(x = "", y = "ARI",fill = " ")

pdf("MethodCompSim_ARI.pdf",width = 7, height = 3.5)
print(g)
dev.off()

g1 <- ggplot(performance, aes(x = Method, y= viDist,fill = Method)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size =15)) +
  #guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(x = "", y = "Variation of information",fill = " ")

pdf("MethodCompSim_viDist.pdf",width = 7, height = 3.5)
print(g1)
dev.off()

#------------------------------------------------------------------------------#
# Simulation run times Figure 2h
#------------------------------------------------------------------------------#
convertToSecs <- function(x){
    u <- units(x)
    if(u == "mins"){
        return(as.numeric(x) *60)
    }else if(u == "hours"){
        return(as.numeric(x) * 60 * 60)
    } else {
        stop("what this?")
    }

}

time <- get(load("SimulationTimes_homotypic.Rda"))
df <- data.frame("time" = numeric(), "method" = character(),"Replicate" = numeric())


for(i in seq_along(time)){
    for(j in seq_along(time[[i]])){
        tmp <- convertToSecs(time[[i]][[j]])
        tmp <- data.frame("time" = tmp, "method" = names(time[[i]])[j],"Replicate" = i)
        df <- rbind(df, tmp)
    }
}

#### changing to minutes for now
df$time <- df$time / 60

df$method <- as.factor(df$method)
cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,4,7,8)]

g <- ggplot(df, aes(x = method, y= time,fill = method)) +
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size =15)) +
  #guides(colour = guide_legend(override.aes = list(size=4)))+
  labs(x = "", y = "Time in Minutes",fill = " ")

pdf("~/Vesalius/MethodCompSim_time.pdf",width = 7, height = 3.5)
print(g)
dev.off()

#------------------------------------------------------------------------------#
# concat cell types use in each sim round
# Suplemntary material
#------------------------------------------------------------------------------#

gt <- list.files(pattern = "homotypic.Rda")
gt <- gt[gt!= "SimulationTimes_homotypic.Rda"]

## reorder because 1 and 10
gt <- gt[c(1,3:10,2)]
cells <- data.frame("Replicate" = numeric(),
                    "Regime" = character(),
                    "n_cells" = numeric(),
                    "territory" = numeric(),
                    "cellType" = character())
for(i in seq_along(gt)){
  load(gt[i])
  n_c <-c(9,12,15,3,3,4,5)
  d_t <-c(rep("uniform",3),"pure",rep("exp",3))
  regime <- paste0("nt3_nc",n_c,"_",d_t)
  for(j in seq_along(tmpCells)){
    if(d_t[j] == "exp"){
      df <- data.frame("Replicate" = rep(i,n_c[j]*3),
                       "Regime" = rep(regime[j],n_c[j]*3),
                       "n_cells" = rep(n_c[j],n_c[j]*3),
                       "territory" = rep(1:3,each=n_c[j]),
                       "cellType" = rep(tmpCells[[j]], times = 3))
    } else {
      df <- data.frame("Replicate" = rep(i,n_c[j]),
                       "Regime" = rep(regime[j],n_c[j]),
                       "n_cells" = rep(n_c[j],n_c[j]),
                       "territory" = rep(1:3,each=n_c[j]/3),
                       "cellType" = tmpCells[[j]])
    }

    cells <- rbind(cells,df)
  }

}
write.csv(cells, file = "~/Vesalius/Simulation_CellTypes.csv")

#------------------------------------------------------------------------------#
# Simulation Plots for all regimes
# Figure S7 to Figure S17
# Figure 2i
#------------------------------------------------------------------------------#
library(ggplot2)
library(patchwork)
library(dplyr)
library(grid)
library(RColorBrewer)

files <- list.files(pattern = ".csv")


plots <- vector("list", length(files))
plotNames <- rep("Simulation", length(files))
for(i in seq_along(files)){
    tmp <- read.csv(files[i], header=T)
    tag <- sapply(strsplit(files[i],"_"),"[[",1)
    title <- gsub("_"," ",files[i])
    title <- gsub(".csv","",title)
    plotNames[i] <- title
    if(tag == "Seurat"){
        seurat <- tmp[,c("x","y","seurat_clusters")]
        seurat$seurat_clusters <- as.factor(as.numeric(as.character(seurat$seurat_clusters))+1)

        cols <- length(levels(seurat$seurat_clusters))
        pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
        cols <- pal(cols)#[sample(1:cols)]


        seu <- ggplot(seurat, aes(x,y,col = as.factor(seurat_clusters))) +
             geom_point(size = 1.5,alpha = 1)+
             theme_void()+
             scale_color_manual(values = cols)+
             theme(legend.text = element_text(size = 10),
                   legend.title = element_text(size=10),
                   plot.title = element_text(size =12),
                   legend.position = "right",
                   plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
             labs(title = title, col = "Cluster")+
             guides(colour = guide_legend(override.aes = list(size=3)))
        plots[[i]] <- seu
    } else if(tag == "Vesalius"){

        vesalius <- tmp %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)


        cols <- length(unique(vesalius$territory))
        pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
        cols <- pal(cols)[sample(1:cols)]


        ves <- ggplot(vesalius, aes(x,y,col = as.factor(territory))) +
             geom_point(size = 1.5,alpha = 1)+
             theme_void()+
             scale_color_manual(values = cols)+
             theme(legend.text = element_text(size = 10),
                   legend.title = element_text(size=10),
                   plot.title = element_text(size =12),
                   legend.position = "right",
                   plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
             labs(title = title, col = "Territory")+
             guides(colour = guide_legend(override.aes = list(size=3)))
        plots[[i]] <- ves
    }else if(tag == "BayesSpace"){
        bayes <- tmp[,c("col","row","spatial.cluster")]
        cols <- length(unique(bayes$spatial.cluster))
        pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
        cols <- pal(cols)#[sample(1:cols)]


        bay <- ggplot(bayes, aes(col,row,col = as.factor(spatial.cluster))) +
             geom_point(size = 1.5,alpha = 1)+
             theme_void()+
             scale_color_manual(values = cols)+
             theme(legend.text = element_text(size = 10),
                   legend.title = element_text(size=10),
                   plot.title = element_text(size =12),
                   legend.position = "right",
                   plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
             labs(title = title, col = "Spatial Cluster")+
             guides(colour = guide_legend(override.aes = list(size=3)))
        plots[[i]] <- bay
    }else if(tag == "Giotto"){
      # If everything is NA assign a single territory.
      # Giotto does not work in some case and I have no idea why.
      if(all(is.na(tmp$HMRF_k3_b.40))){
          tmp$HMRF_k3_b.40 <- 1
      }
      cols <- length(unique(tmp$HMRF_k3_b.40))
      pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
      cols <- pal(cols)#[sample(1:cols)]


      gio <- ggplot(tmp, aes(xcoord,ycoord,col = as.factor(HMRF_k3_b.40))) +
           geom_point(size = 1.5,alpha = 1)+
           theme_void()+
           scale_color_manual(values = cols)+
           theme(legend.text = element_text(size = 10),
                 legend.title = element_text(size=10),
                 plot.title = element_text(size =12),
                 legend.position = "right",
                 plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
           labs(title = title, col = "Spatial Domain")+
           guides(colour = guide_legend(override.aes = list(size=3)))
        plots[[i]] <- gio
    } else{
      message("No idea what this is...")
    }
}
names(plots) <- plotNames


# Get main manuscript plots
purePlots <- plots[grepl("rep6",plotNames) & grepl("pure", plotNames)]
uniformPlots <- plots[grepl("rep6",plotNames) & grepl("uniform", plotNames) & grepl("nc9",plotNames)]
expPlots <- plots[grepl("rep6",plotNames) & grepl("exp", plotNames) & grepl("nc4",plotNames)]

tmpPlots <- c(purePlots,uniformPlots,expPlots)

## Figure 2i
pdf("Simulation.pdf",width = 19,height=15)
plotCols <- 4
plotRows <- 3

grid.newpage()
pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y)}

for (i in seq(1,plotCols*plotRows)) {
  curRow <- ceiling(i/plotCols)
  curCol <- (i-1) %% plotCols + 1
  print(tmpPlots[[i]], vp = vplayout(curRow, curCol ))

}
dev.off()

#------------------------------------------------------------------------------#
# Simulation plots for manuscript
# NOTE this uses the for loop shown above - we are just selecting a subset of
# of those plots to be shown in the main manuscript
# Figure S7 to S17
#------------------------------------------------------------------------------#

reorder <- c(4,7,1,2,3,5,6)
reorder <- c(reorder,reorder+7,reorder+14,reorder+21)

#for(k in seq(1,10)){
for(k in 10){
  locs <- sapply(strsplit(plotNames," "),"[[",3)
  tmpPlots <- plots[locs == paste0("rep",k)]
  tmpPlots <- tmpPlots[reorder]
  png(paste0("~/Vesalius/Simulation_rep",k,".png"),width =1800,height=2400)
  plotCols <- 4
  plotRows <- 7

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y)}

  for (i in seq_along(tmpPlots)) {
    curCol <- ceiling(i/plotRows)
    curRow <- (i-1) %% plotRows + 1
    print(tmpPlots[[i]], vp = vplayout(curRow,curCol ))

  }
  dev.off()
}




#------------------------------------------------------------------------------#
# Concat all diff gene expression files into a single file
# Not so much method comp but hey it's there
# make sure you can differentiate each folder
# Suplemntary material
#------------------------------------------------------------------------------#
files <- list.files(pattern=".csv")
ves <- files[grepl("ves",files)]
vesDEG <- list()
for(i in seq_along(ves)){
    if(any(grepl(",",readLines(ves[i],5)))){
        sep <- ","
        tmp <- read.csv(ves[i],sep=sep, header=T)
        tmp$celltype <- NA
        tmp$group_1 <- NA
        tmp$group_2 <- NA
        tmp <- tmp[,-1]
        tag <- sapply(strsplit(ves[i],".csv"),"[[",1)
        tmp$tissue <- tag

    } else {
        sep <- " "
        tmp <- read.table(ves[i],sep=sep, header=T)
        tag <- sapply(strsplit(ves[i],".csv"),"[[",1)
        tmp$tissue <- tag
    }

    vesDEG[[i]]<- tmp
}
vesDEG <- do.call("rbind",vesDEG)
write.csv(vesDEG,file = "~/vesalius_DEG_concat.csv")

seu <- files[!grepl("ves",files)]
seuDEG <- list()
for(i in seq_along(seu)){
    tmp <- read.csv(seu[i], header=T)
    if(ncol(tmp)==8){
        tmp <- tmp[,-1]
        tag <- sapply(strsplit(seu[i],".csv"),"[[",1)
        tmp$tissue <- tag
    } else {
      tmp$cluster <- NA
      tmp$gene <- tmp[,1]
      tmp <- tmp[,-1]
      tag <- sapply(strsplit(seu[i],".csv"),"[[",1)
      tmp$tissue <- tag
    }

    seuDEG[[i]] <- tmp
}
seuDEG <- do.call("rbind", seuDEG)
write.csv(seuDEG,file = "~/seurat_DEG_concat.csv")
