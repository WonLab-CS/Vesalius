################################################################################
############################   ST package        ###############################
################################################################################

#/----------------------/ Pre-procssing Functions /----------------------------/

## Removing low genes and low counts
## Assuming that this a count matrix
# slide = gene count matrix (digital expression)
# min.genes = number of genes to have for a barcode to be kept
# min.count = number of gene counts to have for a barcode to be kept
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

## normalise data
# slide = gene count matrix (digital expression)
normalizeSlide <- function(slide, method = c("minmax","quantile")){
    norm <- switch(method[1],
                   minmax = .minmax(slide),
                   quantile = normalize.quantiles(slide))
    return(norm)
}

### min max norm
# slide = gene count matrix (digital expression)
# Function used in normalizeSlide (line 45 - 51)
.minmax <- function(slide){
    mi <- min(slide)
    ma <- max(slide)
    slide <- (slide - mi)/ (ma - mi)
    return(slide)
}


## Naive histogram trimmer (Not used)
## channel = colour channel value
## trim = quantile of values to remove of each side
trimHistogram <- function(channel, trim = 0.0125){
    trim <- as.numeric(as.character(channel)) > quantile(as.numeric(as.character(channel)),0 + trim) &
             as.numeric(as.character(channel)) < quantile(as.numeric(as.character(channel)),1 - trim)
    return(trim)
}




trimHistogram <- function(channel, trim = 0.0125){
    trim <- as.numeric(as.character(channel)) > quantile(as.numeric(as.character(channel)),0 + trim) &
             as.numeric(as.character(channel)) < quantile(as.numeric(as.character(channel)),1 - trim)
    return(trim)
}
