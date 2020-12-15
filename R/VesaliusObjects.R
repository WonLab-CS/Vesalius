################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/ Vesalius Objects /---------------------------------#
## Class union for count matrices
setClassUnion(name = 'countMatrix', members = c("matrix", "dgCMatrix"))


#' Vesalius Object
#'
#' Vesalius objects are place holders for data required to run Vesalius.
#'
#' @slot countMatrix matrix containing gene (rows) counts and barcodes (columns)
#' @slot cooridinates data.frame containing x/y cooridinates for each barcode (barcode , xcoord,ycoord)
#' @slot SO seurat object containing gene count and barcode locations (DefaultAssay = "Spatial")
#' @slot rgb List containing the rgb code for each bead and each slice.
#' @slot nn List of nearest neighbors in 2D space.
#' @slot territory  data.frame containing barcodes, cooridinates, RGB code
#' as well as territories and sub territories associated with each barcode.
#' @slot subTerritory data.frame containing a subset of full data set.
#' @slot subCount matrix containing a subset of full data set.
#' @slot subSO seurat object containing a subset of full data set.
#' @slot combine List containing territories and subterritories to extract and combine for downsteam analysis.
#'
#' @name Vesalius-class
#' @rdname Vesalius-class
#' @exportClass Vesalius
#'

setClass("Vesalius",
    slots = c(countMatrix = "countMatrix",
              cooridinates = "data.frame",
              SO = "SeuratObject",
              rgb = "list",
              nn = "list",
              territory = "data.frame",
              subTerritory  = "data.frame",
              subCount = "countMatrix",
              subSO = "SeuratObject",
              combine = "list"
              )
)


Vesalius <- function(countMatrix = NULL, cooridinates = NULL, SO = NULL){
    ## It's check time!
    if(!is.null(countMatrix)){
      if(is(countMatrix) != "matrix" & is(countMatrix) != "dgCMatrix"){
          stop("Unrecognised Object type for countMatrix")
      }
    }
    if(!is.null(cooridinates)){
      if(is(cooridinates) != "data.frame"){
        stop("Unrecognised Object type for cooridinates")
      }
    }
    if(!is.null(SO)){
      if(is(SO) != "SeuratObject"){
        stop("Unrecognised Object type for SO")
      }
    }

    ##
    if(!all(colnames(countMatrix) %in% cooridinates$barcodes)){
        warning("Barcodes in countMatrix and cooridinates do not match.")
    }
    if(!all(colnames(countMatrix) %in% SO@barcodes)){
        warning("Barcodes in countMatrix and cooridinates do not match.")
    }

    return(new("Vesalius",
              countMatrix = countMatrix,
              cooridinates = cooridinates,
              SO = SO))


}
