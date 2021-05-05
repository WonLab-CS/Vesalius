################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Spatial Clustering/---------------------------------#

spatialCluster <- function(image, count,seedTerritory = NULL,
    expressionThreshold = 0.75,verbose =T){
    #--------------------------------------------------------------------------#
    # First let's get the beads associated with territory
    #--------------------------------------------------------------------------#

}

regulariseGeneArray <- function(geneArray,gaussian = 1,inflationFactor = 1000,lambda =1){
    #--------------------------------------------------------------------------#
    # First we assume that files are loaded outside of this function
    # Let's inflate count values
    #--------------------------------------------------------------------------#
    geneArray$counts <- geneArray$counts * inflationFactor
    #geneArray$counts <- (geneArray$counts - min(geneArray$counts)) /
    #(max(geneArray$counts) - min(geneArray$counts))
    #--------------------------------------------------------------------------#
    # now that we have this we can then just convert to a cimg format
    #--------------------------------------------------------------------------#
    img <- geneArray[,c("x","y","cc","counts")]
    colnames(img)<- c("x","y","cc","value")
    img <- as.cimg(img)

    #--------------------------------------------------------------------------#
    # if we want to add a little blur just in case
    #--------------------------------------------------------------------------#
    if(gaussian>0){
      img <- isoblur(img, sigma = gaussian)
    }
    #--------------------------------------------------------------------------#
    # Let's use the tvR package
    #--------------------------------------------------------------------------#
    img <- denoise2(as.matrix(img), lambda = 1,
    method ="TVL2.FiniteDifference",normalize =TRUE)
    return(img)

}

generateGeneImageArray <- function(image, count,minCount = 500,allPixel = FALSE,cores = 1){
    #--------------------------------------------------------------------------#
    # Going to build it as list for now
    #--------------------------------------------------------------------------#
    numCount <- Matrix::rowSums(count) >= minCount
    count <- count[numCount,]
    #--------------------------------------------------------------------------#
    # next we create a pseudo array for each gene
    # everything is save in a seperate file otherwise RAM usage is going to
    #go of the charts
    #--------------------------------------------------------------------------#
    parallel::mclapply(seq_len(nrow(count)),function(idx,img,count, allPixel){
                  print(idx)

                  locs <- colnames(count)[count[idx,] >0]
                  if(allPixel){
                    img <- filter(img, barcodes %in% locs & cc == 1)

                  } else {
                    img <- filter(img, barcodes %in% locs & tile == 1 & cc == 1)

                  }
                  locs <- locs[locs %in% img$barcodes]
                  localCounts <- count[idx,locs]
                  img$counts <- localCounts
                  img <- regulariseGeneArray(img)
                  write.csv(img, file = paste0(rownames(count)[idx],"_ia.csv"),
                  row.names=F, col.names=F, quote =F)
                  #rm(list= ls());gc()
    },img = image , count = count, allPixel = allPixel ,mc.cores =cores)



}

compareGeneArrays <- function(dir,method = "RMSE", fuzz = 1){
    allGenes <- list.files(dir,pattern = ".csv")

    tags <- gsub("_ia.csv","",allGenes)

    simMatrix <- matrix(0,ncol = length(allGenes), nrow = length(allGenes))
    colnames(simMatrix) <- tags
    rownames(simMatrix) <- tags

    for(i in seq_along(allGenes)){
        g1 <- as.matrix(read.csv(allGenes[i]))
        dim(g1) <- c(dim(g1),1)
        g1 <- image_read(g1) %>% image_convert(colorspace = "gray")
        for(j in seq_along(allGenes)){
            g2 <- as.matrix(read.csv(allGenes[j]))
            dim(g2) <- c(dim(g1),1)
            g2 <- image_read(g2) %>% image_convert(colorspace = "gray")

            comp <- image_compare(g1,g2, metric = method,fuzz = 1)
            simMatrix[i,j] <- attributes(comp)$distortion
        }
    }
    return(simMatrix)
}

#generateGeneImageArray(imageBrain, counts,cores =5)
#simMat <- compareGeneArrays(dir())
#save(simMat)
