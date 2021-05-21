################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Territory identity/---------------------------------#

#' extractIdentity compute differentially expressed genes for each territory
#' @param image a dataframe conatining territory information
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param method character describing the statistical test to use in order to extract
#' Differentially expressed genes
#' @param logFC numeric describing minimum log fold change value for diff gene expression
#' @param pval numeric bound between 0 and 1 defining minimum p value for Differentially
#' expressed genes
#' @param minPct numeric defining the minimum percentage of cells that should contain
#' any given gene.
#' @param dilationFactor positive Integer describing the extent of territory growth
#' @param cores numeric describing the number of cores used for the analyis.
#' @details To add
#' @return A data.frame (tibble) containing differentially expressed genes as well
#' p.value, logFC, seedPct (percentage of cells containing gene in first group), queryPct
#' (percentage of cells containing gene in second group), seedTerritory (territory used as group 1)
#' queryTerritory ( territory value that is compared to seed territory)

extractAllMarkers <- function(image,counts,method = "wilcox",
  logFC = 0.25, pval = 0.05,minPct = 0.05,minCell = 10,dilationFactor=0,cores = 1){
    .simpleBar(verbose)
    .checkCounts(verbose)
    #--------------------------------------------------------------------------#
    # Check class of counts and determine approach based on that
    #--------------------------------------------------------------------------#
    if(is(counts) == "Seurat"){
        if(DefaultAssay(counts) == "Spatial"){
            counts <- counts@assays$Spatial@counts
        } else if(DefaultAssay(counts) == "SCT"){
            counts <- counts@assays$SCT@counts
        }
    }
    #--------------------------------------------------------------------------#
    # Let's split the image into seed territories
    #--------------------------------------------------------------------------#
    territories <- split(image,territory)
    seedTerritory <- names(territories)
    #--------------------------------------------------------------------------#
    # For each territory, dilate and subset
    #--------------------------------------------------------------------------#
    for(i in seq_along(territories)){

        #----------------------------------------------------------------------#
        # Dilation factor can only be postive - now erosion for now
        #----------------------------------------------------------------------#
        if(dilationFactor==0){
            seed <- territories[[i]] %>%
              filter(tile ==1 & territory %in% seedTerritory[i])
        } else {
            seed <- territories[[i]] %>% filter(territory %in% seedTerritory[i])
            seed <- .territoryDilation(seed,dilationFactor,image,verbose)
            seed <- filter(seed, tile ==1)
        }
        #----------------------------------------------------------------------#
        # Now that we have the seed barcodes , we can take all the other
        # barcodes as the query image. This ensure that we don't have
        # any overlap between territories even affter dilation
        #----------------------------------------------------------------------#
        query <-image%>%filter(tile ==1 & !barcodes %in% seed$barcodes & cc ==1)

        #----------------------------------------------------------------------#
        # aight let's get dem juicy counts out for both seed and query
        #----------------------------------------------------------------------#
        seed <- counts[,colnames(counts) %in% unique(seed$barcodes)]
        query <- counts[,colnames(counts) %in% unique(query$barcodes)]
        #----------------------------------------------------------------------#
        # Let's do some diff expression analysis
        #----------------------------------------------------------------------#
        deg <- .VesaliusDEG(seed,query,seedTerritory[i],"Remaining",method,
            logFC,pval,minPct,minCell,verbose)
        territories[[i]] <- deg
      }
      #------------------------------------------------------------------------#
      # find nulls
      #------------------------------------------------------------------------#
      nulls <- sapply(territories,is.null)
      territories <- territories[!nulls]
      territories <- do.call("rbind",territories)
      return(territories)

  }


### might not use this function
### I honestly think it is a pain to just do one at a time
### But if needed I can implement it - if required and asked for
extractMarkers <- function(image,counts,seed = NULL,query = NULL, method = "wilcox",
     logFC = 0.25, pval = 0.05,minPct = 0.05,minCell = 10,
     dilationFactorSeed = 0,dilationFactorQuery = 0,verbose=TRUE){
    .simpleBar(verbose)
    .checkCounts(verbose)
    #--------------------------------------------------------------------------#
    # Check class of counts and determine approach based on that
    #--------------------------------------------------------------------------#
    if(is(counts) == "Seurat"){
        if(DefaultAssay(counts) == "Spatial"){
            counts <- counts@assays$Spatial@counts
        } else if(DefaultAssay(counts) == "SCT"){
            counts <- counts@assays$SCT@counts
        }
    }

    #--------------------------------------------------------------------------#
    # setting up seed territoy counts
    #--------------------------------------------------------------------------#
    .seedSelect(verbose)
    if(is.null(seed)){
        stop("Please supply a seed as numeric values describing territory identity")
    }
    #--------------------------------------------------------------------------#
    # getting tmp to use as template when getting query territorires
    #--------------------------------------------------------------------------#
    tmp <- seed
    #--------------------------------------------------------------------------#
    # Dilating if needed
    #--------------------------------------------------------------------------#
    if(dilationFactorSeed == 0){
        seedID <- paste0(seed, sep =" ",collapse="")
        seed <- image %>% filter(tile == 1) %>% filter(territory %in% seed)
        seed <- counts[,colnames(counts) %in% unique(seed$barcodes)]
    } else {
        seedID <- paste0(seed, sep =" ",collapse="")
        seed <- image %>% filter(territory %in% seed)
        seed <- .territoryDilation(seed,dilationFactorSeed,image, verbose)
        seed <- counts[,colnames(counts) %in% unique(seed$barcodes)]
    }

    #--------------------------------------------------------------------------#
    # setting up query territories - if it is set to null then
    # take all the other territories otherwise territories supplied by user
    #--------------------------------------------------------------------------#
    .querySelect(verbose)

    if(is.null(query)){
        queryID <- "Remaining"
        query <- image %>% filter(!territory %in% tmp)
    } else {
        queryID <- paste0(query, sep =" ",collapse="")
        query <- territories %>% filter(territory %in% query)
    }
    if(dilationFactorQuery == 0){
        queryID <- paste0(query, sep =" ",collapse="")
        query <- image %>% filter(tile == 1)
        query <- counts[,colnames(counts) %in% unique(query$barcodes)]
        query <- counts[,!colnames(query) %in% colnames(seed)]
    } else {
        queryID <- paste0(query, sep =" ",collapse="")
        query <- image %>% filter(territory %in% query)
        query <- .territoryDilation(query,dilationFactorQuery,image, verbose)
        query <- filter(query, tile == 1)
        query <- counts[,colnames(counts) %in% unique(query$barcodes)]
        query <- counts[,!colnames(query) %in% colnames(seed)]
    }

    #--------------------------------------------------------------------------#
    # running grouped analysis for DEG
    #--------------------------------------------------------------------------#
    .degProg(verbose)
    deg <- .VesaliusDEG(seed,query,seedID,queryID,method,
        logFC,pval,minPct,minCell,verbose)
    cat("\n")
   .simpleBar(verbose)
   return(markers)

}





.territoryDilation <- function(territories,dilationFactor,image,verbose){

  if(dilationFactor<0){
      stop("Territory Erosion not supported - dilationFactor < 0")
  } else if(dilationFactor >0){
      .dilate(verbose)
      #------------------------------------------------------------------#
      # First we define territory limits and add a little on each
      # side - this ensures that we won't be clipping any parts of the
      # territory
      #------------------------------------------------------------------#
      seed <- territories %>% mutate(value=1)

      dilationFactor <- ifelse(dilationFactor <= 0,1,dilationFactor)
      ymin <- ifelse((min(seed$y) - dilationFactor *2) <=0,1,
         min(seed$y) - dilationFactor *2)
      xmin <- ifelse((min(seed$x) - dilationFactor *2) <=0,1,
         min(seed$x) - dilationFactor *2)
      ymax <- max(seed$y) + dilationFactor * 2
      xmax <- max(seed$x) + dilationFactor * 2
      #------------------------------------------------------------------#
      # Now we convert add buffer boundaeries, covert to grey scale and
      # dilate the territory
      #------------------------------------------------------------------#
      seed <- seed %>% select(c("x","y","cc","value")) %>%
              rbind(.,c(xmin,ymin,1,0),c(xmax,ymax,1,0)) %>%
              as.cimg() %>% grayscale() %>% grow(dilationFactor)
      #------------------------------------------------------------------#
      # Next we rebuild the image data frame with dilated
      #------------------------------------------------------------------#
      seed %>% as.data.frame()
      seed <- inner_join(seed,image,by = c("x","y"))%>%
             select(c("barcodes","x","y","cc.y","z","tile",
                      "cluster","territory"))%>%
             filter(cc.y ==1)
      colnames(seed) <- c("barcodes","x","y","cc","value","tile",
                         "cluster","territory")

     return(seed)
  }
}




.VesaliusDEG <- function(seed,query,seedID,queryID,method,logFC,
    pval,minPct,minCell,verbose =T){
    .degEachProg(seedID,queryID,verbose)
    #--------------------------------------------------------------------------#
    # We assume here that we are parsing cleaned up version of each object
    # This is strictly just computing DEG between two groups
    #
    # Just in case there are not enough cells
    #--------------------------------------------------------------------------#
    dimSeed <- dim(seed)
    dimQuery <- dim(query)
    if(is.null(dimSeed)){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    }else if(is.null(dimQuery)){
        warning(paste0("Territory ",queryID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    } else if(dimSeed[2L] < minCell){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    } else if(dimQuery[2L] < minCell){
        warning(paste0("Territory ",queryID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    }
    #--------------------------------------------------------------------------#
    # Computing DEG metrics
    # Pseudo count = 1
    #--------------------------------------------------------------------------#
    seedPct <- Matrix::rowSums(seed >0) / ncol(seed)
    queryPct <- Matrix::rowSums(query >0) / ncol(query)
    FC <- log(Matrix::rowMeans(seed) +1) - log(Matrix::rowMeans(query)+1)
    #--------------------------------------------------------------------------#
    # Dropping genes that don't fit the logFC and pct criteria
    # At least I won't be computing anything that doesn't need to be
    # I say that while creating a whole bunch of variables
    #--------------------------------------------------------------------------#
    keep <- (seedPct >= minPct | queryPct >= minPct) & abs(FC) >= logFC
    genes <- rownames(seed)[keep]
    seed <- seed[keep,]
    query <- query[keep,]
    seedPct <- seedPct[keep]
    queryPct <- queryPct[keep]
    FC <- FC[keep]

    #--------------------------------------------------------------------------#
    # testing diff gene expression
    # Rebuilding a matrix just in case you only have one gene to test
    # Not going to rewrite the sapply for that
    #--------------------------------------------------------------------------#
    if(is.null(dim(seed))){
      idx <- 1
      seed <- matrix(seed,nrow=1)
      query <- matrix(query, nrow=1)
    } else {
      idx <- seq_len(nrow(seed))
    }

    deg <- sapply(idx,function(idx,seed,query,method){
                  res <- switch(method,
                         "wilcox" = wilcox.test(seed[idx,],query[idx,])$p.value,
                         "t.test" = t.test(seed[idx,],query[idx,])$p.value)
                  return(res)
    },seed,query,method)


    #--------------------------------------------------------------------------#
    # rebuilding data.frame and filtering p values
    # filtering on pval or pvaladjusted?
    #--------------------------------------------------------------------------#
    deg <- tibble("genes" = genes,"p.value" = deg,"p.value.adj" = p.adjust(deg),
                  "seedPct" = seedPct,
                  "queryPct" = queryPct,"logFC" = FC,
                  "seedTerritory" = rep(seedID,length(deg)),
                  "queryTerritory" = rep(queryID,length(deg)))

    deg <- deg %>% filter(p.value <= pval)
    return(deg)
}



# Extracting count values for each cell in a cluster/territory
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
                  return(x@assays$Spatial@counts)
      })
    } else if(by == "territory"){
      #------------------------------------------------------------------------#
      # Just return all cells for that territory
      # Normalised counts ! This is important and will need to change the code accordingly
      # This just considers normalised data
      #------------------------------------------------------------------------#
      newCount <- object@assays$Spatial@counts
    } else {
      #------------------------------------------------------------------------#
      # Placeholder for now
      #------------------------------------------------------------------------#
      newCount <- NULL
    }


    return(newCount)

}


.bindCounts <- function(territories,counts){
    #--------------------------------------------------------------------------#
    # Binding all count matrices together and filling all the "gaps"
    #--------------------------------------------------------------------------#

    cells <- unlist(lapply(territories,colnames))
    counts <- counts[,cells]
    return(counts)
}


### Cant remember why I took this approach
### Lets keep it ultra simple for now
### Dont create this intermediate object only if required
extractTerritories <- function(image,seurat,seedID = NULL,dilationFactor = 0,
                               minCell = 10, verbose = TRUE,cores = 1){
    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # if combine is null we consider that we want to have all
    # territories prepared for clustering analysis
    #--------------------------------------------------------------------------#
    if(is.null(seedID)){
        .extractTerProg("all",verbose)
        territories <- image %>% filter(tile == 1) %>% distinct(barcodes, .keep_all = FALSE)
        territories <- split(territories, territories$territory)
        #----------------------------------------------------------------------#
        # Revert back to fine grain if inage array was reduced
        #----------------------------------------------------------------------#
        if(dilationFactor != 0){
           cat("No seedID specified - Territory Dilation will not be applied.\n")
        }
        cells <- sapply(territories,nrow) > minCell
        territories <- territorries[cells]
        territories <- lapply(territories,"$",territory)
        territories <- parallel::mclapply(territories,subSetTerritories,seurat,mc.cores= cores)
        .simpleBar(verbose)
        return(territories)
    } else {
        #----------------------------------------------------------------------#
        # Filtering only the barcodes that are in the territories of interest
        # To ADD check for combine - numeric values
        #----------------------------------------------------------------------#
        .extractTerProg(seedID,verbose)
        territories <- image %>% filter(territory %in% seedID)
        territories <- .territoryDilation(territories,dilationFactor,image,verbose)
        territories <- filter(territories,tile == 1)  %>%
            distinct(barcodes,.keep_all=F)
        territories <- territories$barcodes

        #----------------------------------------------------------------------#
        # Revert back to fine graind and return null if territory does not
        # Contain enough cells - in this case throw in a warning
        #----------------------------------------------------------------------#
        #territories <- .fineGrain(territories,minCell)

        if(length(territories) < minCell){

            warning("Territory selection does not contain enough cells - NULL returned",
                  immediate. = TRUE)
            .simpleBar(verbose)
            return(NULL)
        }else{

            territories <- subSetTerritories(territories,seurat)
            .simpleBar(verbose)
            return(territories)
        }
    }


}


subSetTerritories <- function(territories,seurat){
    #--------------------------------------------------------------------------#
    # Simplified version for now
    # It might be worth while getting away from seurat later
    # essentially this is a template function
    #--------------------------------------------------------------------------#

    seurat <- subset(seurat, cells = territories)
    return(seurat)
}

compareClusters <- function(ref,seedCluster,queryCluster,seed = NULL,query = NULL,method = "wilcox",
  logFC = 0.25, pval = 0.05,minPct = 0.05,minCell = 10,verbose=TRUE){
      .simpleBar(verbose)
      #------------------------------------------------------------------------#
      # First let's set up seed cells
      # For now we will consider that if seed is null then you want to compare
      # to everything. That is the same as comparing territories
      #------------------------------------------------------------------------#
      .seedSelect(verbose)
      if(!is.null(seed)){
        seedCluster <- FetchData(seedCluster,c("seurat_clusters")) %>% filter(seurat_clusters %in% seed)
      } else {
        seedCluster <- FetchData(seedCluster,c("seurat_clusters"))
      }

      seedCluster <- data.frame(rownames(seedCluster),seedCluster)
      colnames(seedCluster) <- c("barcodes","territory")
      #------------------------------------------------------------------------#
      # Same thing with query cells
      # For now we will consider that if seed is null then you want to compare
      # to everything. That is the same as comparing territories
      #------------------------------------------------------------------------#
      .querySelect(verbose)
      if(!is.null(query)){
          queryCluster <- FetchData(queryCluster,c("seurat_clusters")) %>% filter(seurat_clusters %in% query)
      } else {
          queryCluster <- FetchData(queryCluster,c("seurat_clusters"))
      }

      queryCluster <- data.frame(rownames(queryCluster),queryCluster)
      colnames(queryCluster) <- c("barcodes","territory")
      #------------------------------------------------------------------------#
      # Now we can subset the refence
      # First we check if it is seurat or count matrix
      #------------------------------------------------------------------------#
      .checkCounts(verbose)
      if(is(ref) == "Seurat"){
          if(DefaultAssay(ref) == "Spatial"){
              ref <- ref@assays$Spatial@counts
          } else if(DefaultAssay(ref) == "SCT"){
              ref <- ref@assays$SCT@counts
          }
      }
      #------------------------------------------------------------------------#
      # get counts for seed and query
      #------------------------------------------------------------------------#
      .degProg(verbose)
      deg<- .VesaliusDEG.grouped(counts =ref,seedTer =seedCluster,
                                 queryTer = queryCluster,
                                 logFC= logFC,pval = pval,
                                 minPct = minPct,minCell = minCell,method=method,
                                 verbose = verbose)
      cat("\n")
      .simpleBar(verbose)
      return(deg)
  }


  compareLayers <- function(layers,counts,l1 = NULL, l2 = NULL, method = "wilcox",
    logFC = 0.25, pval = 0.05,minPct = 0.05,minCell = 10,verbose=TRUE){
      #------------------------------------------------------------------------#
      # Ayo now we comapre layers
      # First thing is to get layers so lets start with layer 1
      #------------------------------------------------------------------------#
      .simpleBar(verbose)
      if(!is.null(l1)){

          l1 <- list(filter(layers, layer %in% l1))

      } else {

          l1 <- split(layers,layers$layer)
      }
      #------------------------------------------------------------------------#
      # Get layer 2
      #------------------------------------------------------------------------#
      if(!is.null(l2)){

          l2 <- list(filter(layers, layer %in% l2))
      } else {

          l2 <- split(layers,layers$layer)
      }
      #------------------------------------------------------------------------#
      # Prep counts for deg
      #------------------------------------------------------------------------#
      .checkCounts(verbose)
      if(is(counts) == "Seurat"){
          if(DefaultAssay(counts) == "Spatial"){
              counts <- counts@assays$Spatial@counts
          } else if(DefaultAssay(counts) == "SCT"){
              counts <- counts@assays$SCT@counts
          }
      }
      #------------------------------------------------------------------------#
      # Now we can loop over layers
      #------------------------------------------------------------------------#
      deg <- vector("list", length(l1) * length(l2))
      counter <- 1
      for(i in seq_along(l1)){
          for(j in seq_along(l2)){
              tmp1 <- counts[,colnames(counts %in% unique(l1[[i]]$barcodes))]
              tmp2 <- counts[,colnames(counts %in% unique(l2[[j]]$barcodes))]
              tmpID1 <- paste0(l1[[i]]$layer,collapse = "",sep = " ")
              tmpID2 <- paste0(l2[[j]]$layer,collapse = "",sep = " ")
              deg[[counter]] <- .VesaliusDEG(tmp1,tmp2,tmpID1,tmpID2,
                                                method = method,
                                                logFC=logFC,pval = pval,
                                                minPct= minPct,minCell=minCell,
                                                verbose=verbose)
          }
      }
      deg <- do.call("rbind",deg)
      .simpleBar(verbose)
      return(deg)
  }
