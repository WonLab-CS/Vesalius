#-----------------------------/Method comparison/------------------------------#
#------------------------------------------------------------------------------#
# This file contains all the code used to compare methods - Table of content:
# * BayesSpace and Seurat in slide_seqV2
# * scoring of STAGATE, SPAGCN and SEDR
# * Plotting performance (ARI and VI)
# * Plotting run time
# * Plotting Simulation runs
# * Concet DEG from Vesalius and Seurat
# * Concat of ssimulation cell types
#------------------------------------------------------------------------------#





#------------------------------------------------------------------------------#
# Plotting slide seq output for BayesSpace and Seurat
# Figure 2b and Figure 2d
# Figure S3 to Figure S6
#------------------------------------------------------------------------------#
library(spdep)
library(RColorBrewer)
library(ggplot2)

bayes <- get(load("~/group/slide_seqV2/BayesSpaceBenchMarking/Puck_200115_08_SSV2_BM.Rda"))
bayes$row <- as.numeric(bayes$row)
bayes$col <- as.numeric(bayes$col)

# coord <- Rotation(bayes[,c("col","row")],angle = 190)
# colnames(coord) <- c("x","y")
# bayes <- cbind(bayes,coord)

ter_col <- length(unique(bayes$spatial.cluster))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
ter_col <- ter_pal(ter_col)
b1 <- ggplot(data = bayes,aes(x = col,y=row, col = as.factor(spatial.cluster)))+
      geom_point(size = 0.2, alpha = 0.65)+
      coord_polar()+
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
      theme_bw() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
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
      theme_bw() +
      facet_wrap(~spatial.cluster)+
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = "BayesSpace", col = "Spatial Cluster")
pdf("~/Vesalius/BayesSpaceEmbryoSplit.pdf", width = 16, height=18)
b2s
dev.off()

seurat <- get(load("~/group/slide_seqV2/SeuratBenchMarking/Puck_200115_08_SSV2_BM.Rda"))
ter_col <- length(unique(seurat$seurat_clusters))
ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
#ter_pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
#          "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
ter_col <- ter_pal(ter_col)[sample(seq_len(ter_col),ter_col)]
#ter_col <- viridis(ter_col)
s1 <- ggplot(data = seurat,aes(x = x,y=y, col = as.factor(seurat_clusters)))+
      geom_point(size = 0.2, alpha = 1)+
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
      theme_bw() +
      facet_wrap(~seurat_clusters)+
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
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
      facet_wrap(~seurat_clusters)+
      theme_bw() +
      scale_color_manual(values = ter_col)+
      theme(legend.text = element_text(size = 12),
            legend.title = element_text(size=12),
            plot.title = element_text(size =15),
            legend.position = "bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
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
       stat_compare_means(method = "anova", label.y = 1.15,label.x ="Seurat")+      # Add global p-value
       stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Vesalius",hide.ns = TRUE,label.y=1)+
       theme(legend.text = element_text(size = 12),
                           axis.text = element_text(size = 10),
                           axis.title = element_text(size = 10),
                           plot.tag = element_text(size=15),
                           plot.title = element_text(size=15),
                           legend.title = element_text(size=15),
                           strip.text.x = element_text(size = 15),
                          axis.text.x=element_text(angle=45,hjust=1)) +
       guides(colour = guide_legend(override.aes = list(size=5)))+
       labs(fill = " ",y = "Adjusted Rand Index",x="")

vi <- ggplot(performance, aes(method,VI,fill = method)) +
      geom_boxplot()+
      ylim(-0.2,8)+
      scale_fill_manual(values=cols)+
      theme_bw()+
      facet_wrap(~regime)+
      stat_compare_means(method = "anova", label.y = 7.5,label.x ="Seurat")+      # Add global p-value
      stat_compare_means(label = "p.signif", method = "t.test",
                    ref.group = "Vesalius",hide.ns = TRUE,label.y = 6.7)+
      theme(legend.text = element_text(size = 15),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            plot.tag = element_text(size=15),
            plot.title = element_text(size=15),
            legend.title = element_text(size=15),
            strip.text.x = element_text(size = 15),
          axis.text.x=element_text(angle=45,hjust=1)) +
      guides(colour = guide_legend(override.aes = list(size=5)))+
      labs(fill = " ",y = "Variation of Information",x="")

pdf("ari.pdf", width = 9, height = 9)
print(ari)
dev.off()

pdf("vi.pdf", width = 9, height = 9)
print(vi)
dev.off()



#------------------------------------------------------------------------------#
# Simulation run times Figure 2h
#------------------------------------------------------------------------------#

library(stringr)
input <- "~/group/slide_seqV2"
methods <- c("vesaliusSim","bayesSim","giottoSim","seuratSim","stagateSim","spagcnSim","sedrSim")

time <- paste0(input,"/",methods,"/time.txt")
runTime <- data.frame("method" = character(),
                          "regime" = character(),
                          "time" = numeric(),
                          "unit" = character())
for(i in seq_along(time)){
    tmp <- read.table(time[i], sep =",")
    colnames(tmp) <- c("regime","time","unit")
    tag <- gsub(".csv","",tmp$regime)
    tag <- sapply(strsplit(tmp$regime,"_"),"[",2:4)
    tag <- apply(tag,2, paste0, sep =" ", collapse ="")
    ## Convert to munites
    if(all(tmp$unit == "hours")){
        tmp$time <- tmp$time * 60
        tmp$unit <- "minutes"
    } else if(all(tmp$unit == "sec")){
        tmp$time <- tmp$time/60
        tmp$unit <- "minutes"
    }
    meth <- str_to_title(rep(gsub("Sim","",methods[i]),nrow(tmp)))
    df <- data.frame("method" = meth,
                     tmp)
    runTime <- rbind(runTime,df)

}
runTime$method <- as.factor(runTime$method)
runTime$regime <- as.factor(runTime$regime)

library(ggplot2)
library(ggpubr)

cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,3,4,5,6,7,8)]
runTime$time <- log(runTime$time)
ti <- ggplot(runTime, aes(method,time,fill = method)) +
       geom_boxplot()+
       scale_fill_manual(values=cols)+
       theme_bw()+
       #facet_wrap(~regime)+
       #stat_compare_means(method = "anova", label.y = 1.15)+      # Add global p-value
       #stat_compare_means(label = "p.signif", method = "t.test",
        #             ref.group = "Vesalius",hide.ns = TRUE,label.y=1)+
       theme(legend.text = element_text(size = 12),
              legend.position = "bottom",
                           axis.text = element_text(size = 15),
                           axis.title = element_text(size = 15),
                           plot.tag = element_text(size=15),
                           plot.title = element_text(size=15),
                           legend.title = element_text(size=15),
                           strip.text.x = element_text(size = 15),
                          axis.text.x=element_text(angle=45,hjust=1)) +
       guides(colour = guide_legend(override.aes = list(size=5)))+
       labs(fill = " ",y = "Log(Time in minutes)",x="")

pdf("~/Vesalius/MethodCompSim_time.pdf",width = 7, height = 4.5)
print(ti)
dev.off()


#------------------------------------------------------------------------------#
# Simulation Plots for all regimes
# Figure S7 to Figure S17
# Figure 2i
#------------------------------------------------------------------------------#
library(ggplot2)
library(stringr)
library(patchwork)
library(dplyr)
library(grid)
library(RColorBrewer)

input <- "~/group/slide_seqV2"
methods <- c("vesaliusSim","bayesSim","giottoSim","seuratSim","stagateSim","spagcnSim","sedrSim")

sim <- list.files("~/Vesalius/Simulation", pattern = ".csv",full.names=T)

predictions <- paste0(input,"/",methods)
predictions <- lapply(predictions,function(x){return(list.files(x,pattern=".csv",full.names=T))})
names(predictions) <- str_to_title(gsub("Sim","",methods))



generateSimPlots <- function(sim,predictions){
    n_plots <- length(sim) + sum(sapply(predictions, length))
    plotList <- vector("list", n_plots)
    count <- 1
    for(i in seq_along(sim)){
        simTmp <- read.csv(sim[i], header=T)
        simTag <- strsplit(sim[i],"Simulation/")[[1]][2]
        tag <- gsub(".csv","",simTag)
        tag <- sapply(strsplit(tag,"_"),"[",2:5)
        tag <- apply(tag,2, paste0, sep =" ", collapse ="")
        ter_col <- length(unique(simTmp$ter))
        if(ter_col <= 8){
          ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
        } else {
          ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
        }

        ter_col <- ter_pal(ter_col)
        g <- ggplot(simTmp, aes(x,y, col = as.factor(ter)))+
                             geom_point(size =0.25)+
                             scale_color_manual(values = ter_col)+
                             theme_void()+
                             theme(legend.text = element_text(size = 8),
                                   legend.title = element_text(size=8),
                                   plot.title = element_text(size =10),
                                   legend.position = "left")+
                             guides(colour = guide_legend(override.aes = list(size=3)))+
                             labs(title = paste("Ground Truth",str_to_title(tag)), col = "Ter.")
        plotList[[count]] <-g
        count <- count +1
        for(j in seq_along(predictions)){
            predTmp <- grep(pattern = simTag,x=predictions[[j]], value = T)
            predTmp <- read.csv(predTmp,header=T)
            if(names(predictions)[j] == "Vesalius"){
              ter_col <- length(unique(predTmp$territory))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              g <- ggplot(predTmp, aes(x,y, col = as.factor(territory)))+
                                   geom_point(size =0.25)+
                                   scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste(names(predictions)[j],str_to_title(tag)), col = "Ter.")
              plotList[[count]] <-g
              count <- count +1

            } else if(names(predictions)[j] == "Seurat"){
              ter_col <- length(unique(predTmp$seurat_clusters))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              g <- ggplot(predTmp, aes(x,y, col = as.factor(seurat_clusters)))+
                                   geom_point(size =0.25)+
                                   scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste(names(predictions)[j],str_to_title(tag)), col = "Ter.")
              plotList[[count]] <-g
              count <- count +1

            } else if(names(predictions)[j] == "Bayes"){
              ter_col <- length(unique(predTmp$spatial.cluster))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              g <- ggplot(predTmp, aes(col,row, col = as.factor(spatial.cluster)))+
                                   geom_point(size =0.25)+
                                   scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste0(names(predictions)[j],"Space ",str_to_title(tag)), col = "Ter.")
              plotList[[count]] <-g
              count <- count +1
            } else if(names(predictions)[j] == "Giotto"){
              ter_col <- length(unique(predTmp[,ncol(predTmp)]))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              if(grepl("dot",tag)){
                g <- ggplot(predTmp, aes(x,y, col = as.factor(HMRF_k6_b.40)))

              }else{
                g <- ggplot(predTmp, aes(x,y, col = as.factor(HMRF_k3_b.40)))

              }

                g <- g + geom_point(size =0.25)+
                         scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste(names(predictions)[j],str_to_title(tag)), col = "Ter.")
              plotList[[count]] <-g
              count <- count +1

            }else if(names(predictions)[j] == "Stagate"){
              ter_col <- length(unique(predTmp$louvain))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              g <- ggplot(predTmp, aes(x,y, col = as.factor(louvain)))+
                                   geom_point(size =0.25)+
                                   scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste(names(predictions)[j],str_to_title(tag)), col = "Ter.")
                plotList[[count]] <-g
                count <- count +1
            } else if(names(predictions)[j] == "Spagcn"){
              ter_col <- length(unique(predTmp$refined_pred))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              g <- ggplot(predTmp, aes(x,y, col = as.factor(refined_pred)))+
                                   geom_point(size =0.25)+
                                   scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste(names(predictions)[j],str_to_title(tag)), col = "Ter.")
                plotList[[count]] <-g
                count <- count +1

            } else if (names(predictions)[j] == "Sedr"){
              ter_col <- length(unique(predTmp$SEDR_leiden))
              if(ter_col <= 8){
                ter_pal <- colorRampPalette(brewer.pal(ter_col, "Accent"))
              } else {
                ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
              }

              ter_col <- ter_pal(ter_col)
              g <- ggplot(predTmp, aes(x,y, col = as.factor(SEDR_leiden)))+
                                   geom_point(size =0.25)+
                                   scale_color_manual(values = ter_col)+
                                   theme_void()+
                                   theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size=8),
                                         plot.title = element_text(size =10),
                                         legend.position = "none")+
                                   guides(colour = guide_legend(override.aes = list(size=3)))+
                                   labs(title = paste(names(predictions)[j],str_to_title(tag)), col = "Ter.")
              plotList[[count]] <-g
              count <- count +1
            } else {
              stop("don't know what this is")
            }
        }
    }
    return(plotList)
}


plots <- generateSimPlots(sim, predictions)
count <- 1

for(k in seq(1,8)){


  png(paste0("~/Vesalius/Simulation_rep",k,".png"),width =2400,height=2400)
  plotCols <- 10
  plotRows <- 8

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y)}

  for (i in seq_len(plotCols * plotRows)) {
    curCol <- ceiling(i/plotRows)
    curRow <- (i-1) %% plotRows + 1
    print(plots[[count]], vp = vplayout(curRow,curCol ))
    count <- count +1
  }
  dev.off()
}


simSub <- sim[c(10,12,27,39,41,55,62,74)]
simSub <-simSub[c(1,4,5,6)]
plotSub <- generateSimPlots(simSub, predictions)
count <- 1
pdf(paste0("~/Vesalius/Simulation_summary.pdf"),width = 20, height=12)
plotCols <- 8
plotRows <- 4

grid.newpage()
pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y)}

for (i in seq_len(plotCols * plotRows)) {
    curCol <- ceiling(i/plotCols)
    curRow <- (i-1) %% plotCols + 1
    print(plotSub[[count]], vp = vplayout(curCol,curRow ))
    count <- count +1
}
dev.off()




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

#------------------------------------------------------------------------------#
# Concat cell types from sims
#------------------------------------------------------------------------------#
library(tidyverse)
sim <- list.files("~/Vesalius/Simulation", pattern = ".csv",full.names=T)
df <- vector("list", length(sim))
for(i in seq_along(sim)){
    tmp <- read.csv(sim[i],header=T)
    tmp <- tmp %>% group_by(ter) %>% distinct(celltype,.keep_all=T)
    tag <- gsub(".csv","",sim[i])
    tag <- gsub("/home/pcnmartin/Vesalius/Simulation/","",tag)
    tmp$sim <- tag
    df[[i]] <- select(tmp,c(sim,ter,cells, celltype)) %>% arrange(ter, group_by=sim)
}
df <- do.call("rbind",df)
