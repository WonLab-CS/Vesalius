################################################################################
###############################   Vesalius      ################################
################################################################################


runTransform <- function(vesalius,
                         method = "FFT",
                         chunkSize = 1,
                         cores =1){
  .simpleBar(verbose)

    if(is(vesalius)[1L] == "vesaliusObject"){
      images <- .vesToC(object = vesalius,embed = embedding, dims = dims)
    } else {
      stop("Unsupported format to runTransform function")
    }
    #--------------------------------------------------------------------------#
    # Chunking image into smaller images for each embedding
    # the idea is that we can do this to look at the image at a smaller scale
    # I think this will need to be cleaned up to make sure that we look at Cells
    # and not just sub sections of the images
    #--------------------------------------------------------------------------#
    if(chunkSize <1){
        images <- parallel::mclapply(images,.chunkImage,chunkSize= chunkSize,
                                    mc.cores = cores)
    } else if(chuckSize >1){
        stop("chunkSize to large! Must be between 0 and 1")
    }
    #--------------------------------------------------------------------------#
    
}


.chunkImage <- function(image,chunkSize =1){
    #--------------------------------------------------------------------------#
    # Get dims and convert to cut index
    #--------------------------------------------------------------------------#
    rows <- round(seq(1,dim(image)[1L],l=(1/chunkSize)+1))
    cols <- round(seq(1,dim(image)[2L],l=(1/chunkSize)+1))
    #--------------------------------------------------------------------------#
    # sub image an create list of subimage
    # Images are represented as y = columns and x = rows for some reason
    #--------------------------------------------------------------------------#
    subImages <- list()
    count <- 1
    for(i in seq(1, length(rows)-1)){
        for(j in seq(1,length(cols)-1)){
            subImages[[count]] <- imsub(image,
                                        x >= rows[j] & x <= rows[j+1],
                                        y >= cols[i] & y <= cols[i+1])
            print(paste0(cols[j],"_",cols[j+1]))
            print(paste0(rows[i],"_",rows[i+1]))
            print(count)
            count <- count +1
        }
    }
    return(subImages)

}
