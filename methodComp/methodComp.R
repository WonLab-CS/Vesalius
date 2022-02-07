#-----------------------------/Method comparison/------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code use to generate comparison plots
# This includes:
# computational time
# ARI scores
# Example sample 151673
#------------------------------------------------------------------------------#

##### Getting time plots
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


### ARI of data sets
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



#### example plots
library(ggplot2)
library(patchwork)
library(dplyr)
library(spatialLIBD)
library(RColorBrewer)
library(Seurat)
library(BayesSpace)
library(Giotto)

## using the common example 151673
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




### Simulation plots
files <- list.files(pattern = ".csv")
plots <- list()

for(i in seq_along(files)){
    tmp <- read.csv(files[i], header=T)
    tag <- sapply(strsplit(files[i],"_"),"[[",1)
    title <- gsub("_"," ",files[i])
    title <- gsub(".csv","",title)
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
             theme(legend.text = element_text(size = 12),
                   legend.title = element_text(size=12),
                   plot.title = element_text(size =15),
                   legend.position = "right")+
             labs(title = title, col = "Cluster")
        plots[[i]] <- seu
    } else if(tag == "Vesalius"){

        vesalius <- tmp %>% filter(tile == 1) %>% distinct(barcodes,.keep_all = TRUE)
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
             labs(title = title, col = "Territory")
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
             theme(legend.text = element_text(size = 12),
                   legend.title = element_text(size=12),
                   plot.title = element_text(size =15),
                   legend.position = "right")+
             labs(title = title, col = "Spatial Cluster")
        plots[[i]] <- bay
    }else if(tag == "Giotto"){
      cols <- length(unique(tmp$HMRF_k3_b.40))
      pal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
      cols <- pal(cols)#[sample(1:cols)]


      gio <- ggplot(tmp, aes(xcoord,ycoord,col = as.factor(HMRF_k3_b.40))) +
           geom_point(size = 1.5,alpha = 1)+
           theme_void()+
           scale_color_manual(values = cols)+
           theme(legend.text = element_text(size = 12),
                 legend.title = element_text(size=12),
                 plot.title = element_text(size =15),
                 legend.position = "right")+
           labs(title = title, col = "Spatial Domain")
        plots[[i]] <- gio
    } else{
      message("No idea what this is...")
    }
}

pdf("Simulation.pdf",width = 20,height=35)


plotCols <- 4
plotRows <- 7

grid.newpage()
pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y)}

for (i in seq(1,plotCols*plotRows)) {
  curRow <- ceiling(i/plotCols)
  curCol <- (i-1) %% plotCols + 1
  print(plots[[i]], vp = vplayout(curRow, curCol ))

}
dev.off()
