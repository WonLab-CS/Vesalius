################################################################################
############################   ST package        ###############################
################################################################################

#' imagePlot - Vesalius image plotting
#' @param image a Vesalius data frame containing at least
#' barcodes, x, y, cc, value
#' @param as.cimg logical - If TRUE, will plot directly using cimg plotting.
#' If FALSE, will produce a ggplot object.
#' @param cex numeric - font and point size resizing factor.
#' @details Vesalius provides basic plotting functions. This function
#' only shows Vesalius images (i.e. using true embedded colours or segmented
#' colours).
#'
#' Currently, we advise to use as.cimg as false for custom plots.
#'
#' As Vesalius nearly always returns data frames, custom plots are left to
#' the discretion of the user.
#' @return cimg plot or ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' imagePlot(image)
#' # as ggplot
#' g <- imagePlot(image,as.cimg = F)
#' }

imagePlot <- function(image,
                      as.cimg = TRUE,
                      cex=1){
    if(as.cimg){
        plot(as.cimg(image[,c("x","y","cc","value")]))
    } else {
        fgcol <- select(image, c("cc","value"))
        fgcol <- data.frame(fgcol$value[fgcol$cc ==1],
                            fgcol$value[fgcol$cc ==2],
                            fgcol$value[fgcol$cc ==3])
        fgcol <- rgb(fgcol[,1],fgcol[,2],fgcol[,3])
        image <- image %>% filter(cc==1)
        image <- ggplot(image,aes(x,y))+
                     geom_raster(aes(fill=fgcol))+
                     scale_fill_identity()+
                     theme_classic()+
                     theme(axis.text = element_text(size = cex *1.2),
                           axis.title = element_text(size = cex *1.2),
                           plot.tag = element_text(size=cex * 2),
                           plot.title = element_text(size=cex * 1.5)) +
                    labs(title = "Vesalius - Image",
                           x = "X coordinates", y = "Y coordinates")

         return(image)
    }
}




#' territoryPlot plotting Vesalius territories
#' @param territories a Vesalius data frame conatining territory information
#' (i.e. containing the terriotry column)
#' @param split logical - If TRUE, territories will be plotted in
#' separate panels
#' @param randomise logical - If TRUE, colour palette will be randomised.
#' @param cex numeric describing font size multiplier.
#' @param cex.pt numeric describing point size multiplier.
#' @details Territory plots show all territories in false colour after they
#' have been isolated from a Vesalius image.
#' @return a ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' g <- territoryPlot(image,cex = 12, cex.pt = 1)
#' }

territoryPlot <- function(territories,
                          split = FALSE,
                          randomise = TRUE,
                          cex=1,
                          cex.pt=0.25){
    #--------------------------------------------------------------------------#
    # Dirty ggplot - this is just a quick and dirty plot to show what it look
    # like
    # At the moment - I thinking the user should make their own
    # Not a prority for custom plotting functions
    #--------------------------------------------------------------------------#
    ter <- territories %>% filter(tile==1) 
    #--------------------------------------------------------------------------#
    # Changing label order because factor can suck ass sometimes
    #--------------------------------------------------------------------------#

    sorted_labels <- order(levels(as.factor(ter$territory)))
    if(any(grepl("isolated",ter$territory))){
      sorted_labels[length(sorted_labels)] <- "isolated"
    }

    ter$territory <- factor(ter$territory, levels = sorted_labels)


    #--------------------------------------------------------------------------#
    # My pure hatred for the standard ggplot rainbow colours has forced me
    # to use this palette instead - Sorry Hadely
    #--------------------------------------------------------------------------#
    ter_col <- length(levels(ter$territory))
    ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))

    if(randomise){
        ter_col <- sample(ter_pal(ter_col),ter_col)
    } else {
        ter_col <- ter_pal(ter_col)
    }
    if(split){
        terPlot <- ggplot(ter, aes(x,y,col = territory)) +
               geom_point(size= cex.pt, alpha = 0.65)+
               facet_wrap(~territory)+
               theme_classic() +
               scale_color_manual(values = ter_col)+
               theme(legend.text = element_text(size = cex *1.2),
                     axis.text = element_text(size = cex *1.2),
                     axis.title = element_text(size = cex *1.2),
                     plot.tag = element_text(size=cex * 2),
                     plot.title = element_text(size=cex * 1.5),
                     legend.title = element_text(size=cex * 1.2)) +
               guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
               labs(colour = "Territory nr.", title = "Vesalius - Territories",
                     x = "X coordinates", y = "Y coordinates")
    } else {
      terPlot <- ggplot(ter, aes(x,y,col = territory)) +
                 geom_point(size= cex.pt, alpha = 0.65)+
                 theme_classic() +
                 scale_color_manual(values = ter_col)+
                 theme(legend.text = element_text(size = cex *1.2),
                       axis.text = element_text(size = cex *1.2),
                       axis.title = element_text(size = cex *1.2),
                       plot.tag = element_text(size=cex * 2),
                       plot.title = element_text(size=cex * 1.5),
                       legend.title = element_text(size=cex * 1.2)) +
                 guides(colour = guide_legend(override.aes = list(size=cex * 0.3)))+
                 labs(colour = "Territory nr.", title = "Vesalius - Territories",
                                        x = "X coordinates", y = "Y coordinates")
    }

    return(terPlot)
}



# Not in use at the moment
# Might add in futur iterations of the packages
#.cellProportion <- function(image){
#    prop <- FetchData(image,"seurat_clusters") %>% table %>% as.data.frame
#    colnames(prop) <- c("Cell","Cell Proportion")

#    ter_col <- nrow(prop)
#    ter_pal <- colorRampPalette(brewer.pal(8, "Accent"))
#    ter_col <- ter_pal(ter_col)

#    treemap(prop,
#            index="Cell",
#            vSize="Cell Proportion",
#            type="index",
#            palette = ter_col
#            )
#}




#' viewGeneExpression - plot gene expression.
#' @param image a Vesalius data frame containing barcodes, x, y, cc, value,
#' cluster, and territory.
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param ter integer - Territory ID in which gene expression will be viewed. If
#' NULL, all territories will be selected.
#' @param gene character - gene of interest (only on gene at a time)
#' @param normalise logical - If TRUE, gene expression values will be min/max
#' normalised.
#' @param cex numeric - font size modulator
#' @details Visualization of gene expression in all or in selected territories.
#' Gene expression is shown "as is". This means that if no transformation
#' is applied to the data then normalized raw count will be shown.
#'
#' If normalization and scaling have been applied, normalized counts will be
#' shown.
#'
#' This also applies to any data transformation applied on the count matrix
#' or the Seurat object.
#'
#' We recommend using Seurat Objects.
#'
#' NOTE : You might be promoted to use geom_tiles instead of geom_raster.
#' However - this will increase file size to unreasonable heights.
#' @return a ggplot object (geom_raster)
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' # In all points
#' g <- viewGeneExpression(image,vesalius, genes = "Cst3")
#' # In a specific territory
#' g1 <- viewGeneExpression(image, vesalius, ter = 1, genes = "Cst3")
#' }


viewGeneExpression <- function(image,
                               counts,
                               ter = NULL,
                               gene = NULL,
                               normalise = TRUE,
                               cex =10){
    #--------------------------------------------------------------------------#
    # We will assume that you parse the image with all pixel for now
    # if no territory is specified gene expression on all
    # First let's get the territory
    #--------------------------------------------------------------------------#
    if(!is.null(ter)){
      image <- filter(image, territory %in% as.character(ter)) %>%
               select(c("barcodes","x","y"))

    }
    #--------------------------------------------------------------------------#
    # Next lets get the unique barcodes associated to the selected territories
    # Not great but it gets the job done I guess
    #--------------------------------------------------------------------------#
    barcodes <- image %>% distinct(barcodes, .keep_all = FALSE)
    barcodes <- barcodes$barcodes



    #--------------------------------------------------------------------------#
    # Extract counts from the seurat object
    # This will need to be cleaned up later
    #--------------------------------------------------------------------------#

    if(is(counts) == "Seurat"){
        counts <- subset(counts,cells = barcodes)
        counts <- GetAssayData(counts, slot = "data")
    } else {
        counts <- counts[,barcodes]
    }
    type <- "Expression"

    #--------------------------------------------------------------------------#
    # Getting genes - For now it will stay as null
    # I could add multiple gene viz and what not
    # And yes I know it would be better to be consitent between tidyr and base
    #--------------------------------------------------------------------------#
    if(is.null(gene)){
        stop("Please specifiy which gene you would like to visualize")
    } else {
        #----------------------------------------------------------------------#
        # will need some regex here probs
        #----------------------------------------------------------------------#
        inCount <- rownames(counts) == gene
        #----------------------------------------------------------------------#
        # Just in case the gene is not present
        # this should probably be cleaned up
        #----------------------------------------------------------------------#
        if(sum(inCount) == 0){
            warning(paste(gene," is not present in count matrix.
                          Returning NULL"),immediate.=T)
            return(NULL)
        }

        counts <- counts[inCount,]

        counts <- data.frame(names(counts),counts)
        colnames(counts) <-c("barcodes","counts")
        if(normalise){
            counts$counts <- (counts$counts - min(counts$counts)) /
                             (max(counts$counts) - min(counts$counts))
            type <- "Norm. Expression"
        }

    }
    #--------------------------------------------------------------------------#
    # Adding count values to iamge
    #--------------------------------------------------------------------------#
    image <- right_join(image,counts,by = "barcodes")

    ge <- ggplot(image,aes(x,y))+
          geom_raster(aes(fill = counts))+
          scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")))+
          theme_classic()+
          theme(axis.text = element_text(size = cex ),
                axis.title = element_text(size = cex),
                legend.title = element_text(size = cex),
                legend.text = element_text(size = cex),
                plot.title = element_text(size=cex)) +
          labs(title = gene,fill = type,
               x = "X coordinates", y = "Y coordinates")
    return(ge)
}


#' viewLayerExpression - plot gene expression in layered Vesalius territories
#' @param image a Vesalius data frame containing barcodes, x, y, cc, value,
#' cluster, and territory.
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param gene character - gene of interest (only on gene at a time)
#' @param normalise logical - If TRUE, gene expression values will be min/max
#' normalised.
#' @param cex numeric - font size modulator
#' @details Visualization of gene expression a layered territory.
#' After isolation or a territory, Vesalius can apply territory layering. The
#' expression shown here describes the mean expression in each layer.
#'
#' Gene expression is shown "as is". This means that if no transformation
#' is applied to the data then normalized raw count will be used for each layer.
#'
#' If normalization and scaling have been applied, normalized counts will be
#' used for each layer
#'
#' This also applies to any data transformation applied on the count matrix
#' or the Seurat object.
#'
#' We recommend using Seurat Objects.
#'
#' NOTE : You might be promoted to use geom_tiles instead of geom_raster.
#' However - this will increase file size to unreasonable heights.
#'
#' @return a ggplot object (geom_raster)
#' @examples
#' \dontrun{
#' data("vesalius")
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' layer <- layerTerritory.edge(image, seedTerritory = 1)
#' g <- viewLayerExpression(layer, vesalius, gene = "Cst3")
#' }


viewLayerExpression <- function(image,
                                counts,
                                gene = NULL,
                                normalise =TRUE,
                                cex =10){

    #------------------------------------------------------------------------#
    # Next lets get the unique barcodes associated to the selected territories
    # Not great but it gets the job done I guess
    #------------------------------------------------------------------------#
    barcodes <- image %>% distinct(barcodes, .keep_all = FALSE)
    barcodes <- barcodes$barcodes

    #------------------------------------------------------------------------#
    # Extract counts from the seurat object
    # This will need to be cleaned up later
    #------------------------------------------------------------------------#
    counts <- subset(counts,cells = barcodes)
    if(is(counts) == "Seurat"){
        counts <- GetAssayData(counts, slot = "data")
    }


    #--------------------------------------------------------------------------#
    # Getting genes - For now it will stay as null
    # I could add multiple gene viz and what not
    # And yes I know it would be better to be consitent between tidyr and base
    #--------------------------------------------------------------------------#
    if(is.null(gene)){
        stop("Please specifiy which gene you would like to visualize")
    } else {
        #----------------------------------------------------------------------#
        # will need some regex here probs
        #----------------------------------------------------------------------#
        inCount <- rownames(counts) == gene
        #----------------------------------------------------------------------#
        # Just in case the gene is not present
        # this should probably be cleaned up
        #----------------------------------------------------------------------#
        if(sum(inCount) == 0){
            warning(paste(gene," is not present in count matrix.
                          Returning NULL"),immediate.=T)
            return(NULL)
        }

        counts <- counts[inCount,]

        counts <- data.frame(names(counts),counts)
        colnames(counts) <-c("barcodes","counts")


    }
    #--------------------------------------------------------------------------#
    # Adding count values to iamge
    #--------------------------------------------------------------------------#
    image <- right_join(image,counts,by = "barcodes")

    image <- image %>% group_by(layer) %>% mutate(counts = mean(counts))
    #--------------------------------------------------------------------------#
    # Normalising counts in image
    # There might be a cleaner way of doing This
    # Should I norm over all counts?
    #--------------------------------------------------------------------------#

    if(normalise){
        image$counts <- (image$counts - min(image$counts)) /
                         (max(image$counts) - min(image$counts))
    }
    ge <- ggplot(image,aes(x,y))+
          geom_raster(aes(fill = counts))+
          scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")))+
          theme_classic()+
          theme(axis.text = element_text(size = cex ),
                axis.title = element_text(size = cex),
                legend.title = element_text(size = cex),
                legend.text = element_text(size = cex),
                plot.title = element_text(size=cex)) +
          labs(title = gene,fill = "Mean Expression",
               x = "X coordinates", y = "Y coordinates")
    return(ge)
}




#' viewCellExpression - plot gene expression in cell of choice
#' @param image a Vesalius data frame containing barcodes, x, y, cc, value,
#' cluster, and territory.
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param cells - character vector containing cells of interest.
#' @param gene character - gene of interest (only on gene at a time)
#' @param ter1  vector/integer - Territory ID in which gene expression will be viewed.
#' See details.
#' @param ter2 vector/integer - Territory ID in which gene expression will be viewed.
#' See details.
#' @param normalise logical - If TRUE, gene expression values will be min/max
#' normalised.
#' @param cex numeric - font size modulator
#' @details Visualization of gene expression in selected territories.
#' The purpose of this function is to contrast the expression of genes between
#' different territories.
#' The cells argument should be supplied as a vector of barcodes containing
#' all barcodes of the cell type of interest. Theoretically, you can supply
#' more than one cell type in this vector but Vesalius will treat them as single
#' cell type.
#' If you do not provide any territory and leave the \code{ter1} and \code{ter2}
#' as NULL, this function will plot the expression of your gene of interest
#' in your cell type of interest accross the whole ST assay.
#' Supplying territories to only one of \code{ter1} or \code{ter2} will plot
#' cell type of interest in one set of territories.
#' Supplying territories to both will contrast expression between territories
#'
#' @return a ggplot object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' # For the example we will take a random selection of barcodes
#' cells <- GetAssayData(vesalius, slot = "data")
#' cells <- sample(colnames(cells), 200, replace =F)
#' g <- viewCellExpression(image, vesalius,cells,
#'                         gene = "Cst3",ter1 = 1,ter2=2)
#' }

viewCellExpression <- function(image,
                               counts,
                               cells,
                               gene = NULL,
                               ter1=NULL,
                               ter2=NULL,
                               normalise=FALSE,
                               cex =10){

    #--------------------------------------------------------------------------#
    # First lets do some checks
    # make sure you have at least 1 gene and 1 territory
    #--------------------------------------------------------------------------#
    if(is.null(gene)){
        stop("Please supply the gene you would like to visualise!")
    }
    if(is.null(ter1) && is.null(ter2)){
      CellCoord <- image %>%
                   filter(barcodes %in% cells & tile ==1)%>%
                   droplevels()
      CellCoord$territory <- "Territory 1"
    } else if(!is.null(ter1) && is.null(ter2)){
      CellCoord <- image %>%
                   filter(barcodes %in% cells & tile ==1 &
                          territory %in% ter1) %>%
                   droplevels()
      CellCoord[CellCoord$territory %in% ter1,"territory"] <- "Territory 1"

    } else if(is.null(ter1) && !is.null(ter2)){
      CellCoord <- image %>%
                   filter(barcodes %in% cells & tile ==1 &
                          territory %in% ter2) %>%
                   droplevels()
      CellCoord[CellCoord$territory %in% ter2,"territory"] <- "Territory 2"
    } else {
      CellCoord <- image %>%
                   filter(barcodes %in% cells & tile ==1 &
                          territory %in% c(ter1,ter2)) %>%
                   droplevels()
      CellCoord[CellCoord$territory %in% ter1,"territory"] <- "Territory 1"
      CellCoord[CellCoord$territory %in% ter2,"territory"] <- "Territory 2"
    }
    #--------------------------------------------------------------------------#
    # All non cell barcodes - we still want to plot them so you can at least
    # your cells are with respect to the other beads
    #--------------------------------------------------------------------------#
    NonCellCoord <- image %>% filter(!barcodes %in% cells & tile ==1)%>%
      droplevels()

    #--------------------------------------------------------------------------#
    # Get your count data - using normalized data
    #--------------------------------------------------------------------------#
    if(is(counts)=="Seurat"){
        counts <- GetAssayData(counts, slot = "data")
    }

    #--------------------------------------------------------------------------#
    # Lets prep the data - counts and coords
    #--------------------------------------------------------------------------#
    counts <- counts[gene,cells]

    counts <- counts[!is.na(match(names(counts),CellCoord$barcodes))]
    counts <- as.data.frame(cbind(CellCoord,counts))
    if(normalise){
        counts$counts <- (counts$counts - min(counts$counts)) /
                         (max(counts$counts) - min(counts$counts))
    }
    #--------------------------------------------------------------------------#
    # Lets make the ggplot - dark background for now...
    # If compplaints I will change this
    #--------------------------------------------------------------------------#
    ge <- ggplot()+
          geom_point(data = NonCellCoord ,
                     aes(x,y),
                     alpha =0.75,
                     fill="#2e2e2e",size =cex* 0.01,show.legend=F)+
          geom_point(data = counts,
                     aes(x=x,y=y,shape = territory,col =counts),
                     size = cex *0.125)+
          scale_color_gradientn(colors = rev(brewer.pal(11,"Spectral")))+
          theme_void()+
          theme(panel.background = element_rect(fill = "#474747"),
                legend.title = element_text(size = cex),
                legend.text = element_text(size = cex),
                plot.margin = margin(1, 1, 1, 1, "cm")) +
          labs(col = "Expression")+
          guides(shape = guide_legend(override.aes = list(size=cex * 0.5)))
    return(ge)
}
