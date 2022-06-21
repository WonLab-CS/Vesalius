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
# scoring of STAGATE
#------------------------------------------------------------------------------#
input <- "~/group/slide_seqV2"
files <- list.files(paste0(input,"/stagateSim"), pattern =".csv",full.names=T)
fileTag <- list.files(paste0(input,"/stagateSim"), pattern =".csv",full.names=F)
perf <- paste0(input,"/stagateSim/performance.txt")
if(!file.exists(perf)){
    file.create(perf)
}
for(i in seq_along(files)){
    tmp <- read.csv(files[i], header=T)
    ari <- adjustedRandIndex(tmp$ter,tmp$louvain)
    vi <- vi.dist(tmp$ter,tmp$louvain)
    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)
}

perf <- read.table(paste0(input,"/stagateSim/performance.txt"), sep =",")
tag <- sapply(strsplit(perf$V1,"_"),"[",4:6)
tag <- apply(tag,2, paste0, sep =" ", collapse ="")
perf$V1 <- tag

#------------------------------------------------------------------------------#
# scoring of SpaGCN
#------------------------------------------------------------------------------#
library(mclust)
library(mcclust)
input <- "~/group/slide_seqV2"
files <- list.files(paste0(input,"/spagcnSim"), pattern =".csv",full.names=T)
fileTag <- list.files(paste0(input,"/spagcnSim"), pattern =".csv",full.names=F)
perf <- paste0(input,"/spagcnSim/performance.txt")
if(!file.exists(perf)){
    file.create(perf)
}
for(i in seq_along(files)){
    tmp <- read.csv(files[i], header=T)
    ari <- adjustedRandIndex(tmp$ter,tmp$refined_pred)
    vi <- vi.dist(tmp$ter,tmp$refined_pred)
    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)
}

#perf <- read.table(paste0(input,"/spagcnSim/performance.txt"), sep =",")
#tag <- sapply(strsplit(perf$V1,"_"),"[",4:6)
#tag <- apply(tag,2, paste0, sep =" ", collapse ="")
#perf$V1 <- tag

#------------------------------------------------------------------------------#
# scoring of SEDR
#------------------------------------------------------------------------------#
library(mclust)
library(mcclust)
input <- "~/group/slide_seqV2"
files <- list.files(paste0(input,"/sedrSim"), pattern =".csv",full.names=T)
fileTag <- list.files(paste0(input,"/sedrSim"), pattern =".csv",full.names=F)
perf <- paste0(input,"/sedrSim/performance.txt")
if(!file.exists(perf)){
    file.create(perf)
}
for(i in seq_along(files)){
    tmp <- read.csv(files[i], header=T)
    ari <- adjustedRandIndex(tmp$ter,tmp$SEDR_leiden)
    vi <- vi.dist(tmp$ter,tmp$SEDR_leiden)
    fperf <- paste0(fileTag[i],",",ari,",",vi,"\n")
    cat(fperf, file = perf, append = TRUE)
}


#------------------------------------------------------------------------------#
# Combine performances
#------------------------------------------------------------------------------#
library(stringr)
input <- "~/group/slide_seqV2"
methods <- c("vesaliusSim","bayesSim","giottoSim","seuratSim","stagateSim","spagcnSim","sedrSim")

perf <- paste0(input,"/",methods,"/performance.txt")
performance <- data.frame("method" = character(),
                          "regime" = character(),
                          "ARI" = numeric(),
                          "VI" = numeric())
for(i in seq_along(perf)){
    tmp <- read.table(perf[i], sep =",")
    if(methods[i] == "stagateSim" | methods[i] == "spagcnSim" | methods[i] == "sedrSim"){
      tag <- gsub(".csv","",tmp$V1)
      tag <- sapply(strsplit(tmp$V1,"_"),"[",3:5)
      tag <- apply(tag,2, paste0, sep =" ", collapse ="")
    } else {
      tag <- gsub(".csv","",tmp$V1)
      tag <- sapply(strsplit(tmp$V1,"_"),"[",2:4)
      tag <- apply(tag,2, paste0, sep =" ", collapse ="")
    }
    meth <- str_to_title(rep(gsub("Sim","",methods[i]),nrow(tmp)))
    df <- data.frame("method" = meth,
                     "regime" = str_to_title(tag),
                     "ARI" = tmp$V2,
                     "VI" = tmp$V3)
    performance <- rbind(performance,df)

}
performance$method <- as.factor(performance$method)
performance$regime <- as.factor(performance$regime)


#------------------------------------------------------------------------------#
# Plotting performance
#------------------------------------------------------------------------------#
library(ggplot2)
library(ggpubr)

cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,3,4,5,6,7,8)]

ari <- ggplot(performance, aes(method,ARI,fill = method)) +
       geom_boxplot()+
       ylim(-0.2,1.2)+
       scale_fill_manual(values=cols)+
       theme_bw()+
       facet_wrap(~regime)+
       stat_compare_means(method = "anova", label.y = 1)+      # Add global p-value
       stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Vesalius",hide.ns = TRUE)+
       theme(legend.text = element_text(size = 12),
                           axis.text = element_text(size = 15),
                           axis.title = element_text(size = 15),
                           plot.tag = element_text(size=15),
                           plot.title = element_text(size=15),
                           legend.title = element_text(size=15),
                           strip.text.x = element_text(size = 15),
                          axis.text.x=element_text(angle=45,hjust=1)) +
       guides(colour = guide_legend(override.aes = list(size=5)))+
       labs(fill = " ",y = "Adjusted Rand Index",x="")

vi <- ggplot(performance, aes(method,VI,fill = method)) +
      geom_boxplot()+
      ylim(-0.2,7.5)+
      scale_fill_manual(values=cols)+
      theme_bw()+
      facet_wrap(~regime)+
      stat_compare_means(method = "anova", label.y = 6)+      # Add global p-value
      stat_compare_means(label = "p.signif", method = "t.test",
                    ref.group = "Vesalius",hide.ns = TRUE)+
      theme(legend.text = element_text(size = 15),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            plot.tag = element_text(size=15),
            plot.title = element_text(size=15),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15),
          axis.text.x=element_text(angle=45,hjust=1)) +
      guides(colour = guide_legend(override.aes = list(size=5)))+
      labs(fill = " ",y = "Variation of Information",x="")

pdf("ari.pdf", width = 15, height = 9)
print(ari)
dev.off()

pdf("vi.pdf", width = 15, height = 9)
print(vi)
dev.off()

#------------------------------------------------------------------------------#
# Simulation scores for all regimes
# Figure 2g
#------------------------------------------------------------------------------#



performance$Method <- as.factor(performance$Method)


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
