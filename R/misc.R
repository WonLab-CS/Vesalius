################################################################################
############################   ST package        ###############################
################################################################################

# Misc functions function that may or may not be dropped in the future.
# Who knows? I certainly don't. But you designed and wrote this packages.
# That means diddly squat!! Since when do programmers know what they are doing?

#' subSetTerritories gets all barcodes associated with a territory from a
#' Seurat object.
#' @param territories vector of barcode values (generally character strings)
#' @param seurat Seurat object containing all barcode values from ST assay
#' @details Essentially a wrapper function to the Seurat \code{subset} function.
#' The function serves mainly the purpose of a place holder function for
#' future iterations of Vesalius.
#' @return Seurat object only containing the desired barcodes.
#' @examples
#' \dontrun{

#' }

subSetTerritories <- function(territories,seurat){
    #--------------------------------------------------------------------------#
    # Simplified version for now
    # It might be worth while getting away from seurat later
    # essentially this is a template function
    #--------------------------------------------------------------------------#

    seurat <- subset(seurat, cells = territories)
    return(seurat)
}


#' getSeuratCoordinates get barcode coordinates from a seurat object
#' @param seurat Seurat object containing all barcode values from ST assay
#' @details Essentially a wrapper function to the Seurat
#' \code{GetTissueCoordinates} function.
#' Mainly serve as a way if generalising output format.
#' @return Seurat object only containing the desired barcodes.
#' @examples
#' \dontrun{
#' data(Vesalius)
#' }
getSeuratCoordinates <- function(seurat){
  ret <-GetTissueCoordinates(seurat)
  if(sum(colnames(ret) %in% c("imagerow","imagecol"))==2){
      colnames(ret) <- c("y","x")
  }
  return(ret)
}



# not in use
.bindCounts <- function(territories,counts){
    #--------------------------------------------------------------------------#
    # Binding all count matrices together and filling all the "gaps"
    #--------------------------------------------------------------------------#

    cells <- unlist(lapply(territories,colnames))
    counts <- counts[,cells]
    return(counts)
}


# Extracting count values for each cell in a cluster/territory
# not in use
.getGeneCounts <- function(object,by){
    #--------------------------------------------------------------------------#
    # This is just some cleaning and data extraction
    # Note that this code relies on Seurat code
    # This will need to be changed when refactoring
    #--------------------------------------------------------------------------#
    if(by == "cluster"){
      clusters <- FetchData(object,"seurat_clusters")
      #------------------------------------------------------------------------#
      # Get all barcodes associated with each cluster
      #------------------------------------------------------------------------#
      barcodes <- lapply(unique(clusters$seurat_clusters),function(idx,obj){
                return(WhichCells(obj,idx))
      }, object)
      #------------------------------------------------------------------------#
      # Get counts associated to each clusters
      #------------------------------------------------------------------------#
      counts <- lapply(barcodes, function(bar,obj){
                return(subset(obj, cells = bar))
      })
      #------------------------------------------------------------------------#
      # Rebuild subsetted count matrix
      ## Will need to change this for more felxibility
      # Set to default assay or reduction based
      #------------------------------------------------------------------------#
      newCount <- lapply(object, function(x){
                  return(GetAssayData(x,slot ="data"))
      })
    } else if(by == "territory"){
      #------------------------------------------------------------------------#
      # Just return all cells for that territory
      # Normalised counts ! This is important and will need to change the code
      # accordingly
      # This just considers normalised data
      #------------------------------------------------------------------------#
      newCount <- GetAssayData(object, slot ="data")
    } else {
      #------------------------------------------------------------------------#
      # Placeholder for now
      #------------------------------------------------------------------------#
      newCount <- NULL
    }


    return(newCount)

}


### might be dropped in final version
### Not in use
.addTerritories <- function(dat,coordinates = NULL, global = TRUE){

    ters <- names(dat)

    dat <- lapply(seq_along(ters), function(idx,ters,dat){
                  dat[[idx]]$territory <- ters[idx]
                  return(dat[[idx]])
    },ters,dat)

    if(!is.null(coordinates)){
        dat <- mapply(function(dat,coordinates){
                      tmp <- cbind(dat,coordinates[,c("x","y")])
                      return(tmp)
        },dat,coordinates, SIMPLIFY = FALSE)
    }


    dat <- do.call("rbind",dat)
    if(global){
      dat <- .globaliseTerritories(dat, seurat = TRUE)
    }

    return(dat)
}





# Used to convert territories per cluster to territories across the whole
# ST array
.globaliseTerritories <- function(img,seurat=FALSE){
    if(!seurat){
      imgTmp <- img %>% filter(territory != "isolated")
      ter <- paste0(imgTmp$cluster,"_", imgTmp$territory)
      allTer <- unique(ter)
      ter <- seq_along(allTer)[match(ter,allTer)]
      img$territory[img$territory != "isolated"] <- ter
      return(img)
    } else {
      imgTmp <- img %>% filter(territory != "isolated")
      ter <- paste0(img$seurat_clusters,"_", img$territory)
      allTer <- unique(ter)
      ter <- seq_along(allTer)[match(ter,allTer)]
      img$seurat_clusters[img$territory != "isolated"] <- ter
      return(img)
    }


}


.selectSimilar <- function(img,cols,segment,threshold = "auto"){
  img <- img %>% { . - imfill(dim=dim(img),val=cols[segment,2:4]) } %>%
         imsplit("c") %>%
         enorm()
  img <- !threshold(img,threshold)
  img <- as.cimg(img)
  return(img)
}
.detectEdges <- function(img){
  img <- img %>% imgradient("xy") %>%
         enorm() %>%
         add() %>%
         sqrt()
  return(img)
}

.pmap <- function(edges,method = c("inverse","none")){
    pmap <- switch(method[1],
                   "inverse" = 1 /(1+ edges),
                   "none" = edges)
    return(pmap)
}

.watershed <- function(img,pmap){
  #----------------------------------------------------------------------------#
  # First create empty image to fill
  #----------------------------------------------------------------------------#
  seed <- imfill(dim=dim(pmap))
  #----------------------------------------------------------------------------#
  # Next we can get back ground and foreground pixels
  # basing this on the pixSet object describing similar "colours"
  # Fill in seed
  #----------------------------------------------------------------------------#
  back <- which(img == 0, arr.ind = TRUE)
  back <- back[sample(seq_len(nrow(back)),1),]
  seed[back[1],back[2],back[3],back[4]]<- 0
  fore <- which(img == 1, arr.ind = TRUE)
  fore <- fore[sample(seq_len(nrow(fore)),1),]
  seed[fore[1],fore[2],fore[3],fore[4]]<- 1
  #----------------------------------------------------------------------------#
  # watershed using seed image and priority map built using the pixset
  # describing similar colours
  #----------------------------------------------------------------------------#
  wt <- as.cimg(watershed(seed, pmap))
  return(wt)
}

.removeWaterPool <- function(water, image){
  #----------------------------------------------------------------------------#
  # Now that we have one territory, we want to remove from the original pool
  # and start again with the next water pool
  #----------------------------------------------------------------------------#
  w <- which(water == 2,arr.ind = TRUE)
  n <- which(water == 1,arr.ind = TRUE)
  image[w[,1],w[,2],w[,3],w[,4]] <- 2
  image[n[,1],n[,2],n[,3],n[,4]] <- 1
  return(image)
}



.vesToC <- function(object, dims,correctBackground =TRUE,verbose =TRUE){
    #--------------------------------------------------------------------------#
    # Get relevant slots and prepare for joining
    #--------------------------------------------------------------------------#
    tiles <- object@tiles
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    #--------------------------------------------------------------------------#
    imageList <- list()
    .vtc(verbose)
    for(i in seq_along(dims)){
      embeds <- object@embeddings[,dims[i]]
      embeds <- data.frame(names(embeds),embeds)
      colnames(embeds) <- c("barcodes",as.character(dims[i]))
      cimg <- right_join(tiles,embeds, by= "barcodes")
      colnames(cimg) <- c("barcodes","x","y","value")
      cimg <- na.exclude(cimg)
      if(correctBackground){
          cimgTmp <- cimg %>%
             select(c("x","y","value")) %>%
             suppressWarnings %>%
             as.cimg %>%
             as.data.frame
          nonImg <- paste0(cimg$x,"_",cimg$y)
          inImg <- paste0(cimgTmp$x,"_",cimgTmp$y)
          cimgTmp[!inImg %in% nonImg,"value"] <- median(cimg$value)
          imageList[[i]] <- suppressWarnings(as.cimg(cimgTmp[,c("x","y","value")]))
      } else {
         imageList[[i]] <- suppressWarnings(as.cimg(cimg[,c("x","y","value")]))
      }


    }
    return(imageList)
}

.vesToDF <- function(object, dims,verbose =TRUE){
  #--------------------------------------------------------------------------#
  # Get relevant slots and prepare for joining
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  #--------------------------------------------------------------------------#
  # generate a list of images based on the number of dims
  # Remember that any time you do anything to an image
  # it is always applied to each "channel" separately - it will be so
  # much easier to just consider everything as gray scale for this
  #--------------------------------------------------------------------------#
  imageDF <- list()
  .vtdf(verbose)
  for(i in seq_along(dims)){
    embeds <- object@embeddings[,dims[i]]
    embeds <- data.frame(names(embeds),embeds)
    colnames(embeds) <- c("barcodes",as.character(dims[i]))
    cimg <- right_join(tiles,embeds, by= "barcodes")
    cimg$cc <- i
    colnames(cimg) <- c("barcodes","x","y","value","cc")
    cimg <- na.exclude(cimg)
    cimg <- cimg[,c("barcodes","x","y","cc","value")]
    imageDF[[i]] <- cimg
   }
  imageDF <- do.call("rbind", imageDF)
  return(imageDF)
}

.cToVes <- function(cimg,object,dims){
    #--------------------------------------------------------------------------#
    # Get stuff out
    #--------------------------------------------------------------------------#
    tiles <- object@tiles
    embeds <- object@embeddings
    #--------------------------------------------------------------------------#
    # Always going to be a gray scale image.
    # Colour are only used when during viz
    # This allows any arbitrary number of embeds
    #--------------------------------------------------------------------------#
    for(i in seq_along(dims)){
      img <- as.data.frame(cimg[[i]])
      barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                  distinct(barcodes,value) %>% na.exclude()

      locs <- match(rownames(embeds),barcodes$barcodes)
      embeds[locs,dims[i]] <- barcodes$value

    }
    object@embeddings <- embeds

    return(object)
}

.dfToVes <- function(df, object,dims){
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  embeds <- object@embeddings
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  # this is probably a shit way to do it
  # re think for later 
  #--------------------------------------------------------------------------#
  for(i in seq_along(dims)){
    img <- filter(df, cc == i)
    barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                distinct(barcodes,value) %>% na.exclude()

    locs <- match(rownames(embeds),barcodes$barcodes)
    embeds[locs,dims[i]] <- barcodes$value

  }
  object@embeddings <- embeds

  return(object)
}

.vesToSIS <- function(object,dims){
    tiles <- object@tiles
    #--------------------------------------------------------------------------#
    # generate a list of images based on the number of dims
    # Remember that any time you do anything to an image
    # it is always applied to each "channel" separately - it will be so
    # much easier to just consider everything as gray scale for this
    # also this is a different image format so we will adpt the function
    # we had above.
    #--------------------------------------------------------------------------#
    imageList <- list()
    for(i in seq_along(dims)){
      embeds <- object@embeddings[,dims[i]]
      embeds <- data.frame(names(embeds),embeds)
      colnames(embeds) <- c("barcodes",as.character(dims[i]))
      sis <- right_join(tiles,embeds, by= "barcodes") %>% na.exclude()
      colnames(sis) <- c("barcodes","x","y","value")
      #------------------------------------------------------------------------#
      # Using the cimg format just because they have optimised format conversion
      #------------------------------------------------------------------------#
      sis <- suppressWarnings(as.cimg(sis[,c("x","y","value")]))
      sis <- as.matrix(sis)
      #------------------------------------------------------------------------#
      # hopefully this will do the trick - OpenImage require 3d Arrays
      # so we will create 3 instances of the same gray scale array
      #------------------------------------------------------------------------#
      img <- c(sis,sis,sis)
      dim(img) <- c(nrow(sis),ncol(sis),3)
      imageList[[i]] <- img
    }
    return(imageList)
}

.SISToVes <- function(image,object,dims){
  #--------------------------------------------------------------------------#
  # Get stuff out
  #--------------------------------------------------------------------------#
  tiles <- object@tiles
  embeds <- object@embeddings
  #--------------------------------------------------------------------------#
  # Always going to be a gray scale image.
  # Colour are only used when during viz
  # This allows any arbitrary number of embeds
  #--------------------------------------------------------------------------#
  for(i in seq_along(dims)){
    img <- .SISToDF(image[[i]])
    barcodes <- left_join(tiles,img, by = c("x","y")) %>%
                distinct(barcodes,value)
    locs <- match(barcodes$barcodes,rownames(embeds))
    embeds[locs,dims[i]] <- barcodes$value

  }
  object@embeddings <- embeds

  return(object)
}

.SISToDF <- function(image, is.cimg = TRUE){
    image <- image$AP_image_data
    y <- rep(seq(1,ncol(image)), each = nrow(image))
    x <- rep(seq(1,nrow(image)), times = ncol(image))
    value <- as.vector(image[seq_len(nrow(image)),seq_len(ncol(image)),1])

    df <- data.frame("x" = x,"y" = y, "value" = value)
    return(df)
}
