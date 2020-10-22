################################################################################
############################   ST package        ###############################
################################################################################

#/---------------------/ Plotting territories Functions /----------------------/

    ### Plotting

plotClusters <- function(img, invert = TRUE){
      clusters <- unique(img$cluster)
      par(mfrow=c(3,3),mar = c(3,3,3,3))

      plot(0, type="n", xlim=c(min(img$xcoord),max(img$xcoord)),
                        ylim=c(min(img$ycoord),max(img$ycoord)),axes=F )
      for(i in seq_along(clusters)){
        tmp <- img[img$cluster == clusters[i],]
        cols <- c(median(tmp$R),median(tmp$G),median(tmp$B))
        points(tmp$xcoord,tmp$ycoord,col=rgb(1-cols[1],1-cols[2],1-cols[3],alpha=0.75), pch=19,cex=0.9)
                              #plot(tmp$xcoord,tmp$ycoord,col=rgb(1-cols[1],1-cols[2],1-cols[3]), pch=19,cex=0.4)
      }
      for(i in seq_along(clusters)){
          tmp <- img[img$cluster == clusters[i],]
          cols <- c(median(tmp$R),median(tmp$G),median(tmp$B))
          #points(tmp$xcoord,tmp$ycoord,col=rgb(1-cols[1],1-cols[2],1-cols[3]), pch=19,cex=0.4)
          plot(tmp$xcoord,tmp$ycoord,xlim=c(min(img$xcoord),max(img$xcoord)),axes=F,
                            ylim=c(min(img$ycoord),max(img$ycoord)),
                            col=rgb(1-cols[1],1-cols[2],1-cols[3], alpha = 0.75), pch=19,cex=0.9)
      }
}



plotTerritories <- function(img,cluster,cex=1){
      par(mar=c(3,3,3,3))
      xlims <- c(min(img$xcoord),max(img$xcoord))
      ylims <- c(min(img$ycoord),max(img$ycoord))

      plot(img$xcoord,img$ycoord,xlim=xlims,ylim=ylims,
          col = rgb(1-img$R,1-img$G,1-img$B), pch=19, cex=cex,axes=F)

      tmp <- img[img$cluster == cluster,]
      cols <- c(median(tmp$R),median(tmp$G),median(tmp$B))

      plot(tmp$xcoord,tmp$ycoord,,xlim=xlims,ylim=ylims,
           col=rgb(1-cols[1],1-cols[2],1-cols[3], alpha = 0.75),
           pch=19,cex=cex,axes=F, main = cluster)

      tercol <- rainbow(length(unique(tmp$territories)))
      centers <- split(tmp,tmp$territories)
      centers <- lapply(centers, function(tmp){
                        x <- median(tmp$xcoord)
                        y <- median(tmp$ycoord)
                        return(list("x" = x,"y" = y))
                  })
      plot(tmp$xcoord,tmp$ycoord,,xlim=xlims,ylim=ylims,
           col = tercol[tmp$territories], pch = 19,cex=cex,axes=F)

      for(i in seq_along(names(centers))){
          text(names(centers)[i],x=centers[[i]]$x,y=centers[[i]]$y,cex=2)
      }
}





## Plotting RGB colour code to barcode coordinates
## rgb = output of assignRGBtoPixel - data frame with barcode name,
# x/y coordinates  and RGB colour code
## This is the quick and dirty method for plotting
## The aesthetics depend of r base plot aesthetics.
## invert = if colour code should be inverted (from rgb to 1-rgb)

plotRGBPCA <- function(rgb,invert = FALSE){
  ## Get plot limits
  xmin <- min(as.numeric(as.character(rgb$xcoord)))
  xmax <- max(as.numeric(as.character(rgb$xcoord)))
  ymin <- min(as.numeric(as.character(rgb$ycoord)))
  ymax <- max(as.numeric(as.character(rgb$ycoord)))

  ## Remove plot inner margins
  par(mar = c(0,0,0,0))
  if(invert){
    ## define plotting area
    plot(0,type ="n", axes = F, xlab = "X coordinates",ylab ="Y coordinates",
         xlim = c(xmin,xmax),
         ylim = c(ymin,ymax))


    ## Adding point to plot
    ## point type is defined by pch
    ## point transparency defined by alpha
    points(x = as.numeric(as.character(rgb$xcoord)),
           y = as.numeric(as.character(rgb$ycoord)),pch =20 ,
           col = rgb(1-as.numeric(as.character(rgb$R)),
                     1-as.numeric(as.character(rgb$G)),
                     1-as.numeric(as.character(rgb$B)),alpha =0.75))

  } else {
    ## define plotting area
    ## define plotting area
    plot(0,type ="n", axes = F, xlab = "X coordinates",ylab ="Y coordinates",
         xlim = c(xmin,xmax),
         ylim = c(ymin,ymax))
    ## Adding point to plot
    ## point type is defined by pch
    ## point transparency defined by alpha
    points(x = as.numeric(as.character(rgb$xcoord)),
           y = as.numeric(as.character(rgb$ycoord)),pch =19 ,
           col = rgb(as.numeric(as.character(rgb$R)),
                     as.numeric(as.character(rgb$G)),
                     as.numeric(as.character(rgb$B)),alpha =0.75))
    }


}

## Plotting RGB colour code to matrix coordinates
## rgb = output of assignRGBtoPixelBlockQuick - data frame with barcode names,
# column and row coordinates in pixel matrix and RGB colour code
## This is the quick and dirty method for plotting
## The aesthetics depend of r base plot aesthetics.
## invert = if colour code should be inverted (from rgb to 1-rgb)
plotRGBBlock <- function(rgb,invert = TRUE){
    ## removing inner margins and defining plotting area
    par(mar =c(0,0,0,0))
    plot(0,type ="n", axes = F, xlab = "",ylab ="",
      xlim = c(min(as.numeric(as.character(rgb$idx))),max(as.numeric(as.character(rgb$idx)))),
      ylim = c(min(as.numeric(as.character(rgb$idy))),max(as.numeric(as.character(rgb$idy)))))




  if(invert){
    # plotting rectangles for each pixel
    # pixel shift towards the left and upwards
    rect(as.numeric(as.character(rgb$idx)),
               as.numeric(as.character(rgb$idy)),
               as.numeric(as.character(rgb$idx+1)),
               as.numeric(as.character(rgb$idy+1)),
               col = rgb(1-as.numeric(as.character(rgb[,"R"])),
                         1-as.numeric(as.character(rgb[,"G"])),
                         1-as.numeric(as.character(rgb[,"B"])),alpha =1),
               border =NA)
  } else {
    rect(as.numeric(as.character(rgb$idx)),
               as.numeric(as.character(rgb$idy)),
               as.numeric(as.character(rgb$idx+1)),
               as.numeric(as.character(rgb$idy+1)),
               col = rgb(as.numeric(as.character(rgb[,"R"])),
                         as.numeric(as.character(rgb[,"G"])),
                         as.numeric(as.character(rgb[,"B"])),alpha =1),
               border =NA)
  }




}


plotSubCluster <- function(subCluster){
    ## Extract Data
    ## We will consider that you are running the for loop out side for now

    sub <- subCluster$subCluster
    FetchData(sub,
            c("UMAP_1", "UMAP_2","seurat_clusters"
              )) %>%
    group_by(seurat_clusters) %>%
    mutate(xu = median(UMAP_1), yu = median(UMAP_2))  -> all_p


  #colset <- wes_palette("Darjeeling1",5)


  p1 <- ggplot(all_p, aes(UMAP_1, UMAP_2, colour = seurat_clusters)) +
    geom_point(alpha = 0.5, size = 2.2) +
    theme_bw() +
    theme(legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=5)))+
    geom_text(mapping = aes(xu, yu, label = seurat_clusters),
              data = unique(all_p[, c("seurat_clusters", "xu", "yu")]),
              inherit.aes = F) +
    labs(colour = "Cluster nr.")

  p2 <- SpatialDimPlot(sub, stroke = 0.1,pt.size.factor = 2.2)

  return(list("p1" = p1,"p2" = p2))
}
