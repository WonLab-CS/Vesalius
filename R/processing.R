################################################################################
###############################   Vesalius      ################################
################################################################################

#/----------------------/ Pre-procssing Functions /----------------------------/

#' filter count Matrix - removing low information beads
#'@param slide count matrix with beads/barcodes as columns and genes as rows
#'@param min.genes minimum number of different genes a barcodes should have
#'@param min.count minimum number of counts per barcode
#'@return A reduced version of the count matrix with barcodes as columns and
#'genes as rows
filterCountMatrix <- function(slide, min.genes = 50, min.count = 100 ){
    ## Extracting values
    count <- apply(slide,2,function(x){
                  return(sum(x))
    })
    genes <- apply(slide,2,function(x){
                  return(sum(x!=0))
    })

    cells_to_keep <- which(count >= min.count & genes >= min.genes)

    # Subsetting slide
    slide <- slide[,cells_to_keep]
    ## checking for zero gene count
    zero <- as.vector(as.matrix(apply(slide,1,sum)))
    slide <- slide[zero != 0,]

    return(slide)

}


