################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Territory identity/---------------------------------#

#' extractIdentity compute differentially expressed genes for each territory
#' @param territories dataframe conatining territory information
#' @param counts count matrix - either matrix, sparse matrix or seurat object
#' This matrix should contain genes as rownames and cells/barcodes as colnames
#' @param method character describing the statistical test to use in order to extract
#' Differentially expressed genes
#' @param logFC numeric describing minimum log fold change value for diff gene expression
#' @param pval numeric bound between 0 and 1 defining minimum p value for Differentially
#' expressed genes
#' @param minPct numeric defining the minimum percentage of cells that should contain
#' any given gene.
#' @param extractOver character describing how Differentially expressed genes
#' should be computed. Three option are available: "territory","cluster" or "all".See details
#' @param between logical describing if territories should be compared between each other
#' or only against all other cells.
#' @param cores numeric describing the number of cores used for the analyis.
#' @details To add
#' @return A data.frame (tibble) containing differentially expressed genes as well
#' p.value, logFC, seedPct (percentage of cells containing gene in first group), queryPct
#' (percentage of cells containing gene in second group), seedTerritory (territory used as group 1)
#' queryTerritory ( territory value that is compared to seed territory)
extractAllMarkers <- function(territories,counts, method = "wilcox",
  logFC = 0.25, pval = 0.05,minPct = 0.05,minCell = 10, extractOver = "territory",between = TRUE,cores = 1){
    #--------------------------------------------------------------------------#
    # Lets assume that you have dataframe or tibble
    # We can do the split in side it will make thing cleaner for
    # code conversion later down the line
    #--------------------------------------------------------------------------#
    countType <- class(counts)

    #--------------------------------------------------------------------------#
    # Switch depending on how we want to split
    # all checks can be done internally
    #--------------------------------------------------------------------------#
    markers <- switch(extractOver,
                      "territory" = .territoryMarkers(territories,counts,method,between,
                                     logFC,pval,minPct,cores),
                      "cluster" = .clusterMarkers(territories,counts,method,between,
                                   logFC,pval,minPctcores),
                      "all" = .allMarkers(territories,counts,method,between,
                              logFC,pval,minPct,cores))


   return(markers)

}

### might not use this function
### I honestly think it is a pain to just do one at a time
### But if needed I can implement it - if required and asked for
extractMarkers <- function(territories,counts,seedTer = NULL,queryTer = NULL, method = "wilcox",
  logFC = 0.25, pval = 0.05,minPct = 0.05,minCell = 10,verbose=TRUE){
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
    if(is.null(seedTer)){
        stop("Please supply seedTer as numeric values describing territory identity")
    }
    seed <- territories %>% filter(tile == 1) %>% filter(territory %in% seedTer)

    #--------------------------------------------------------------------------#
    # setting up query territories - if it is set to null then
    # take all the other territories otherwise territories supplied by user
    #--------------------------------------------------------------------------#
    .querySelect(verbose)

    if(is.null(queryTer)){
        query <- territories %>% filter(tile == 1) %>% filter(!territory %in% seedTer)
    } else {
        query <- territories %>% filter(tile == 1) %>% filter(territory %in% queryTer)
    }
    #--------------------------------------------------------------------------#
    # running grouped analysis for DEG
    #--------------------------------------------------------------------------#
    .degProg(verbose)
    markers <- .VesaliusDEG.grouped(counts,seed,query,
                logFC,pval,minPct,minCell,method,verbose)
    cat("\n")
   .simpleBar(verbose)
   return(markers)

}



.territoryMarkers <- function(territories, counts, method,between = TRUE,
                             logFC,pval,minPct,minCell,verbose=TRUE,cores=1){
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
    # Next if between == TRUE we want to compare territories between each other
    #--------------------------------------------------------------------------#

    if(between){
        territories <- .VesaliusDEG.each(territories,counts,logFC,
                                         pval,minPct,minCell,method,verbose,cores)

    } else {
        #----------------------------------------------------------------------#
        # splitting only so i can use mclapply
        #----------------------------------------------------------------------#
        territories <- territories %>% distinct(barcodes,.keep_all = TRUE)
        territories <- split(territories, territories$territory)
        territories <- parallel::mclapply(territories,.VesaliusDEG.all,
                                          counts,logFC,pval,minPct,minCell,
                                          method,verbose,
                                          mc.cores = cores)

        territories <- do.call("rbind",territories)
    }
    .simpleBar(verbose)
    return(tibble(territories))


}

.VesaliusDEG.grouped <- function(counts,seedTer,queryTer,logFC,pval,
                                 minPct,minCell,method,verbose){
    #--------------------------------------------------------------------------#
    # generating messages
    #--------------------------------------------------------------------------#

    seedID <- paste(unique(seedTer$territory),collapse = " ")
    ter <- paste(unique(queryTer$territory),collapse = " ")
    #--------------------------------------------------------------------------#
    # Don't need to subset seed as it has already been done
    # subset the query variable
    #--------------------------------------------------------------------------#
    .degEachProg(seedID,ter,verbose)
    seed <- counts[,colnames(counts) %in% seedTer$barcodes]
    query <- counts[,colnames(counts) %in% queryTer$barcodes]

    #--------------------------------------------------------------------------#
    # Just in case there are not enough cells
    #--------------------------------------------------------------------------#
    dimSeed <- dim(seed)
    dimQuery <- dim(query)
    if(is.null(dimSeed)| is.null(dimQuery)){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    } else if(dimSeed[2L] < minCell | dimQuery[2L] < minCell){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
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
    #--------------------------------------------------------------------------#
    deg <- tibble("genes" = genes,"p.value" = deg,"p.value.adj" = p.adjust(deg),"seedPct" = seedPct,
                  "queryPct" = queryPct,"logFC" = FC,
                  "seedTerritory" = rep(seedID,length(deg)),
                  "queryTerritory" = rep(ter,length(deg)))

    deg <- deg %>% filter(p.value <= pval)
    return(deg)

}


.VesaliusDEG.each <- function(territory,counts,logFC,pval,minPct,minCell,method,verbose,cores){
    #--------------------------------------------------------------------------#
    # Selecting territories from counts - each comparing to each individual
    # territory - I don't think I'm going to add "all" in this section
    # That can be added to the upper level function
    #--------------------------------------------------------------------------#
    territory <- territory %>% distinct(barcodes,.keep_all = TRUE)
    splitTerritory <- split(territory, territory$territory)
    deg <- vector("list",length(splitTerritory))

    for(i in seq_along(splitTerritory)){
      #------------------------------------------------------------------------#
      # lapply over all the other territories - this is were cores come into play
      #------------------------------------------------------------------------#
      query <- splitTerritory[!names(splitTerritory) %in% names(splitTerritory)[i]]
      query <- parallel::mclapply(names(query),.VesaliusDEG.int,
                                  query = query,seed = splitTerritory,
                                  seedID = i,method = method ,counts = counts,
                                  logFC = logFC,pval = pval,minPct = minPct,
                                  minCell = minCell,verbose = verbose,
                                  mc.cores = cores)
      cat("\n")
      deg[[i]] <- do.call("rbind",query)

    }
    deg <- do.call("rbind",deg)
    return(deg)
}

.VesaliusDEG.int <- function(ter,query,seed,seedID,method,counts,logFC,pval,
                             minPct,minCell,verbose=TRUE){
    #--------------------------------------------------------------------------#
    # Don't need to subset seed as it has already been done
    # subset the query variable
    # TO SANATISE! REFACTOR THE SHIT OUT OF THIS
    #--------------------------------------------------------------------------#
    .degEachProg(seedID,ter,verbose)
    seed <- counts[,colnames(counts) %in% seed[[seedID]]$barcodes]
    query <- counts[,colnames(counts) %in% query[[ter]]$barcodes]

    #--------------------------------------------------------------------------#
    # Just in case there are not enough cells
    #--------------------------------------------------------------------------#
    dimSeed <- dim(seed)
    dimQuery <- dim(query)
    if(is.null(dimSeed)| is.null(dimQuery)){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    } else if(dimSeed[2L] < minCell | dimQuery[2L] < minCell){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
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
    #--------------------------------------------------------------------------#
    deg <- tibble("genes" = genes,"p.value" = deg,"p.value.adj" = p.adjust(deg),
                  "seedPct" = seedPct,
                  "queryPct" = queryPct,"logFC" = FC,
                  "seedTerritory" = rep(seedID,length(deg)),
                  "queryTerritory" = rep(ter,length(deg)))

    deg <- deg %>% filter(p.value <= pval)
    return(deg)
}

.VesaliusDEG.all <- function(territory,counts,logFC,pval,minPct,minCell,
                             method = "wilcox",verbose=TRUE){
    #--------------------------------------------------------------------------#
    # Selecting territories from counts - All = Against all
    #--------------------------------------------------------------------------#
    .degEachProg(unique(territory$territory),"All",verbose)
    seed <- counts[,colnames(counts) %in% territory$barcodes]
    query <- counts[,!colnames(counts) %in% territory$barcodes]

    #--------------------------------------------------------------------------#
    # Just in case there are not enough cells
    #--------------------------------------------------------------------------#
    dims <- dim(seed)
    if(is.null(dims)){
        warning(paste0("Territory ",unique(territory$territory)," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    } else if(dims[2]<minCell){
        warning(paste0("Territory ",unique(territory$territory)," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    }

    #--------------------------------------------------------------------------#
    # Computing DEG metrics
    # Pseudo count = 1
    #--------------------------------------------------------------------------#
    seedPct <- (Matrix::rowSums(seed >0) / ncol(seed))
    queryPct <- (Matrix::rowSums(query >0) / ncol(query))

    FC <- log(Matrix::rowMeans(seed) +1) - log(Matrix::rowMeans(query)+1)

    #--------------------------------------------------------------------------#
    # Dropping genes that don't fit the logFC and pct criteria
    # At least I won't be computing anything that doesn't need to be
    # I say that while creating a whole bunch of variables
    # also I want to return all diff expressed genes even if that means
    # that the genes are down regulated
    # You can do the filtering afterwards and that is just as interesting
    #--------------------------------------------------------------------------#

    keep <- (seedPct >= minPct | queryPct >= minPct) & abs(FC) >= logFC
    genes <- rownames(seed)[keep]
    seed <- seed[keep,]
    query <- query[keep,]
    seedPct <- seedPct[keep]
    queryPct <- queryPct[keep]
    FC <- FC[keep]
    territory <- unique(territory$territory)

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
    # rebuilding data.frame
    #--------------------------------------------------------------------------#
    deg <- tibble("genes" = genes,"p.value" = deg,
                  "p.value.adj" = p.adjust(deg),"seedPct" = seedPct,
                  "queryPct" = queryPct,"logFC" = FC,
                  "seedTerritory" = territory,
                  "queryTerritory" = rep("all",length(deg)))

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
extractTerritories <- function(img,seurat,terIdent = NULL,combine = FALSE,
                               minCell = 10, verbose = TRUE,cores = 1){
    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # if combine is null we consider that we want to have all
    # territories prepared for clustering analysis
    #--------------------------------------------------------------------------#
    if(is.null(terIdent)){
        .extractTerProg("all",verbose)
        territories <- img %>% filter(tile == 1) %>% distinct(barcodes, .keep_all = FALSE)
        territories <- split(territories, territories$territory)
        #----------------------------------------------------------------------#
        # Revert back to fine grain if inage array was reduced
        #----------------------------------------------------------------------#
        #territories <- parallel::mclapply(territories,.fineGrain,minCell,mc.cores = cores)
        #territories <- territories[!sapply(territories,is.null)]
        cells <- sapply(territories,nrow) > minCell
        territories <- territorries[cells]
        territories <- lapply(territories,"$",territory)
        territories <- parallel::mclapply(territories,.subSetTerritories,seurat,mc.cores= cores)
        .simpleBar(verbose)
        return(territories)
    } else {
        #----------------------------------------------------------------------#
        # Filtering only the barcodes that are in the territories of interest
        # To ADD check for combine - numeric values
        #----------------------------------------------------------------------#
        territories <- img %>% filter(tile==1 & territory %in% terIdent)%>%
                       distinct(barcodes, .keep_all = FALSE)
        .extractTerProg(terIdent,verbose)
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

            territories <- .subSetTerritories(territories,seurat)
            .simpleBar(verbose)
            return(territories)
        }
    }


}


.subSetTerritories <- function(territories,seurat){
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
          } else if(DefaultAssay(ref) == "SCT"){
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
              #TODO change when DEG functions are refactored
              ## This sucks - it is so messy
              
              deg[[counter]] <- .VesaliusDEG.int(j,l2,l1,i,
                                                method = method,counts =counts,
                                                logFC=logFC,pval = pval,
                                                minPct= minPct,minCell=minCell,
                                                verbose=verbose)
          }
      }
      deg <- do.call("rbind",deg)
      return(deg)
  }
