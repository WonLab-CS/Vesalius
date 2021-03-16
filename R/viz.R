################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------------/Data Viz /------------------------------------#


imagePlot <- function(image, as.cimg = TRUE){
    if(as.cimg){
        plot(as.cimg(image[,c("x","y","cc","value")]))
    } else {
        fgcol <- select(image, c("cc","value"))
        fgcol <- data.frame(fgcol$value[fgcol1$cc ==1],
                             fgcol$value[fgcol1$cc ==2],
                             fgcol$value[fgcol1$cc ==3])
        fgcol <- rgb(fgcol[,1],fgcol[,2],fgcol[,3])
        image <- image %>% filter(cc==1)
        image <- ggplot(image,aes(x,y))+
                     geom_raster(aes(fill=fgcol))+
                     scale_fill_identity()+
                     theme_minimal()+
                     labs(title = "Vesalius Image")+
                     theme(plot.title = element_text(hjust = 0.5))
         return(image)
    }
}





territoryPlot <- function(territories, split = FALSE,randomise = TRUE){
    #--------------------------------------------------------------------------#
    # Dirty ggplot - this is just a quick and dirty plot to show what it look like
    # At the moment - I thinking the user should make their own
    # Not a prority for custom plotting functions
    #--------------------------------------------------------------------------#
    ter <- territories %>% filter(tile==1) %>% distinct(barcodes, .keep_all =TRUE)

    #--------------------------------------------------------------------------#
    # My pure hatred for the standard ggplot rainbow colours has forced me
    # to use this palette instead - Sorry Hadely
    #--------------------------------------------------------------------------#
    ter_col <- length(unique(ter$territory))
    ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
    if(randomise){
        ter_col <- sample(ter_pal(ter_col),ter_col)
    } else {
        ter_col <- ter_pal(ter_col)
    }
    if(split){
        ter <- ggplot(ter, aes(x,y,col = as.factor(territory))) +
               geom_point(size= 0.25, alpha = 0.65)+
               facet_wrap(~territory)+
               theme_minimal() +
               scale_color_manual(values = ter_col)+
               theme(legend.text = element_text(size = 12),
                     axis.text = element_text(size = 12),
                     axis.title = element_text(size = 12),
                     plot.tag = element_text(size=20)) +
               guides(colour = guide_legend(override.aes = list(size=3)))+
               labs(colour = "Territory nr.", title = "Vesalius - Territories",
                     x = "X coordinates", y = "Y coordinates")
    } else {
      ter <- ggplot(ter, aes(x,y,col = as.factor(territory))) +
                 geom_point(size= 0.25, alpha = 0.65)+
                 theme_minimal() +
                 scale_color_manual(values = ter_col)+
                 theme(legend.text = element_text(size = 12),
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 12),
                       plot.tag = element_text(size=20)) +
                 guides(colour = guide_legend(override.aes = list(size=3)))+
                 labs(colour = "Territory nr.", title = "Vesalius - Territories",
                                        x = "X coordinates", y = "Y coordinates")
    }

    return(ter)
}
