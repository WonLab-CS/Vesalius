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

imagePlot <- function(vesalius,
                      dims = seq(1,3),
                      embedding = "last",
                      as.cimg = TRUE,
                      cex=1){
    #--------------------------------------------------------------------------#
    # Can we somehow vizualise more than 3 dims? how would we go about doing that
    #--------------------------------------------------------------------------#
    if(length(dims)!=3 & length(dims)!=1){
        stop("Only grey scale or RGB images are supported!")
    }
    if(as.cimg){
        image <- .vesToC(vesalius,dims,embed = embedding, verbose=FALSE)

        image <- imappend(image,"cc")
        plot(image)
    } else {
        image <- .vesToC(vesalius,dims,embed = embedding)

        image <- as.data.frame(imappend(image,"cc"))
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

territoryPlot <- function(vesalius,
                          trial = "last",
                          split = FALSE,
                          randomise = TRUE,
                          cex=10,
                          cex.pt=0.25){
    #--------------------------------------------------------------------------#
    # Dirty ggplot - this is just a quick and dirty plot to show what it look
    # like
    # At the moment - I thinking the user should make their own
    # Not a prority for custom plotting functions
    #--------------------------------------------------------------------------#

    if(!is.null(vesalius@territories) & trial == "last"){
      trial <- colnames(vesalius@territories)[ncol(vesalius@territories)]
      ter <- vesalius@territories[,c("x","y",trial)]
      colnames(ter) <- c("x","y","territory")
      ter$territory <- as.factor(ter$territory)
      legend <- sapply(strsplit(trial,"_"),"[[",1)
    } else if(!is.null(vesalius@territories) & trial != "last") {
      if(length(grep(x = colnames(vesalius@territories),pattern = trial))==0){
          stop(paste(deparse(substitute(trial)),"is not in territory data frame"))
      }
      ter <- vesalius@territories[,c("x","y",trial)]
      colnames(ter) <- c("x","y","territory")
      ter$territory <- as.factor(ter$territory)
      legend <- sapply(strsplit(trial,"_"),"[[",1)

    } else {
      stop("No territories have been computed!")
    }

    #--------------------------------------------------------------------------#
    # Changing label order because factor can suck ass sometimes
    #--------------------------------------------------------------------------#

    #sorted_labels <- order(levels(ter$territory))
    #ter$territory <- factor(ter$territory) %>% fct_reorder(sorted_labels)


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
               labs(colour = legend, title = paste("Vesalius",trial),
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
                 labs(colour = legend, title = paste("Vesalius",trial),
                                        x = "X coordinates", y = "Y coordinates")
    }

    return(terPlot)
}




#' viewGeneExpression - plot gene expression.
#' @param image a Vesalius data frame containing barcodes, x, y, cc, value,
#' cluster, and territory.
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param ter integer - Territory ID in which gene expression will be viewed. If
#' NULL, all territories will be selected.
#' @param genes character - gene of interest (only on gene at a time)
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


viewGeneExpression <- function(vesalius,
                               genes = NULL,
                               normMethod = "last",
                               trial = "last",
                               territory1 = NULL,
                               territory2 = NULL,
                               cells = NULL,
                               norm = TRUE,
                               as.layer = FALSE,
                               cex =10){
    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    if(normMethod == "last"){
        counts <- as.matrix(vesalius@counts[[length(vesalius@counts)]])

    }else{
        if(length(grep(x=names(vesalius@counts), pattern = normMethod))==0){
            stop(paste0(deparse(substitute(normMethod)),"is not in count list!"))
        }
        counts <- as.matrix(vesalius@counts[[normMethod]])
    }
    #--------------------------------------------------------------------------#
    # Next let's get the territory data
    # it occurs to me that If done propely I could remove the view layer
    # expression function
    #--------------------------------------------------------------------------#
    if(!is.null(vesalius@territories) & trial == "last"){
      trial <- colnames(vesalius@territories)[ncol(vesalius@territories)]
      ter <- vesalius@territories[,c("barcodes","x","y",trial)]
      colnames(ter) <- c("barcodes","x","y","trial")
      #------------------------------------------------------------------------#
      # Getting and setting territory categories
      #------------------------------------------------------------------------#
      if(is.null(territory1) & is.null(territory2)){
        ter <- select(ter, c("barcodes","x","y","trial"))
      }else if(!is.null(territory1) & is.null(territory2)){
        ter$trial[!ter$trial %in% territory1 ] <- "other"
        if(!is.null(cells)){
            ter$trial[ter$trial != "other" & ter$barcodes %in% cells] <-
            paste0(territory1, collapse=" ")
        }

      }else if(is.null(territory1) & !is.null(territory2)){
        ter$trial[!ter$trial %in% territory2 ] <- "other"
        if(!is.null(cells)){
            ter$trial[ter$trial != "other" & ter$barcodes %in% cells] <-
            paste0(territory2, collapse=" ")
        }
      }else {
        ter$trial[!ter$trial %in% c(territory1,territory2) ] <- "other"
        if(!is.null(cells)){
            ter$trial[ter$trial %in% territory1 & ter$barcodes %in% cells] <-
            paste0(territory1, collapse=" ")
            ter$trial[ter$trial %in% territory2 & ter$barcodes %in% cells] <-
            paste0(territory2, collapse=" ")
        }
      }



    } else if(!is.null(vesalius@territories) & trial != "last"){
      if(!any(grepl(x = colnames(vesalius@territories),pattern = trial))){
          stop(paste(deparse(substitute(trial)),"is not in territory data frame"))
      }
      ter <- vesalius@territories[,c("barcodes","x","y",trial)]
      colnames(ter) <- c("barcodes","x","y","trial")
      #------------------------------------------------------------------------#
      # Getting and setting territory categories
      #------------------------------------------------------------------------#
      if(is.null(territory1) & is.null(territory2)){
        ter <- select(ter, c("barcodes","x","y","trial"))
      }else if(!is.null(territory1) & is.null(territory2)){
        ter$trial[!ter$trial %in% territory1 ] <- "other"

      }else if(is.null(territory1) & !is.null(territory2)){
        ter$trial[!ter$trial %in% territory1 ] <- "other"
      }else {
        ter$trial[!ter$trial %in% c(territory1,territory2) ] <- "other"
      }

    }else{
      ter <- vesalius@tiles %>% filter(origin == 1)
    }

    #--------------------------------------------------------------------------#
    # Getting genes
    # We will facet wrap if there is more than one.
    #--------------------------------------------------------------------------#
    if(is.null(genes)){
        stop("Please specifiy which gene you would like to visualize")
    } else {
        #----------------------------------------------------------------------#
        # First lets make a list to store our counts and then plots
        #----------------------------------------------------------------------#
        geneList <- vector("list",length(genes))
        names(geneList) <- genes
        #----------------------------------------------------------------------#
        # Next we loop,extract, combine
        #----------------------------------------------------------------------#
        for(i in seq_along(genes)){
            gene <- rownames(counts) == genes[i]
            if(sum(gene) == 0){
                warning(paste(genes[i]," is not present in count matrix.
                              Returning NULL"),immediate.=T)
                geneList[[i]] <- NULL
                next()
            }

            gene <- counts[gene,]
            gene <- data.frame(names(gene),gene)
            colnames(gene) <- c("barcodes","gene")

            gene <- left_join(gene,ter, by = c("barcodes")) %>% na.exclude()
            if(norm){
                gene$gene <- (gene$gene - min(gene$gene)) /
                                 (max(gene$gene) - min(gene$gene))
                type <- "Norm. Expression"
            }else{
                type <- "Expression"
            }
            geneList[[i]] <- gene
        }
    }


    #--------------------------------------------------------------------------#
    # Now we can loop of gene list
    #--------------------------------------------------------------------------#
    for(i in seq_along(geneList)){
        #----------------------------------------------------------------------#
        # Okay let's check if there are territories to prioritise
        #----------------------------------------------------------------------#

        if(any(grepl("trial",colnames(ter)))){
          #--------------------------------------------------------------------#
          # Extract data and convert count to mean count if as.layer=T
          #--------------------------------------------------------------------#
          other <- filter(geneList[[i]],trial == "other") %>% droplevels()
          territory <- filter(geneList[[i]],trial != "other") %>% droplevels()

          if(as.layer){
            territory <- territory %>% group_by(trial) %>% mutate(gene = mean(gene))
          }
          #--------------------------------------------------------------------#
          # Create background
          #--------------------------------------------------------------------#
          if(nrow(other)>0){
            ge <- ggplot() + geom_point(data = other ,
                       aes(x,y),
                       alpha =0.75,
                       fill="#2e2e2e",size =cex* 0.01,show.legend=T)
          } else{
            ge <- ggplot()
          }
          #--------------------------------------------------------------------#
          # Create ggplot with foreground
          #--------------------------------------------------------------------#
          if(is.null(cells)){
            ge <- ge + geom_point(data = territory,
                             aes(x=x,y=y,col =gene),
                             size = cex *0.1)
          } else {
            ge <- ge + geom_point(data = territory,
                             aes(x=x,y=y,shape = trial,col =gene),
                             size = cex *0.1)
          }

          ge <- ge +
                scale_color_gradientn(colors = rev(brewer.pal(11,"Spectral")))+
                theme_classic()+
                theme(#panel.background = element_rect(fill = "#474747"),
                      legend.title = element_text(size = cex),
                      legend.text = element_text(size = cex),
                      plot.margin = margin(1, 1, 1, 1, "cm")) +
                labs(col = type,title = names(geneList)[i])+
                guides(shape = guide_legend(override.aes = list(size=cex * 0.5)))

            geneList[[i]] <- ge
        }else{
          #--------------------------------------------------------------------#
          # Create simple ggplot as there is not background here
          #--------------------------------------------------------------------#
          ge <- ggplot()+
                geom_point(data = geneList[[i]],
                           aes(x=x,y=y,col =gene),
                           size = cex *0.1)+
                scale_color_gradientn(colors = rev(brewer.pal(11,"Spectral")))+
                theme_classic()+
                theme(#panel.background = element_rect(fill = "#474747"),
                      legend.title = element_text(size = cex),
                      legend.text = element_text(size = cex),
                      plot.margin = margin(1, 1, 1, 1, "cm")) +
                labs(col = type,title = names(geneList)[i])+
                guides(shape = guide_legend(override.aes = list(size=cex * 0.5)))
            geneList[[i]] <- ge
        }
    }
    #--------------------------------------------------------------------------#
    # returning ggplot or ggplot list with ggarrange
    #--------------------------------------------------------------------------#
    if(length(geneList)==1){
        geneList <- geneList[[1L]]
    }else{
        geneList <- ggarrange(plotlist = geneList,ncol = floor(sqrt(length(geneList))))
    }
    return(geneList)
}
