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
#' data(Vesalius)
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
#' The function serves mainly the purpose of a place holder function for 
#' future iterations of Vesalius. 
#' @return Seurat object only containing the desired barcodes. 
#' @examples 
#' \dontrun{
#' data(Vesalius)
#' }
getSeuratCoordinates <- function(seurat){
  return(GetTissueCoordinates(seurat))
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
