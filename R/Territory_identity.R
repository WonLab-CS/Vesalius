################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Territory identity/---------------------------------#

#' extractIdentity compute differentially expressed genes for each territory

extractIdentity <- function(territories,counts, method = "wilcox",
  logFC = 0.5, pvalThreshold = 0.05,min.pct = 0.1, extractOver = "territory",between = TRUE){
    #--------------------------------------------------------------------------#
    ## First run I will only consider the list object
    ## Get the data out then we can think about refactoring
    ## Consider you have a full list after computing clusters
    #--------------------------------------------------------------------------#

    markers <- switch(extractOver,
                      "territory" = .territoryMarkers(territories,counts,method,logFC,pvalThreshold,min.pct,between),
                      "cluster" = .clusterMarkers(territories,counts,method,logFC,pvalThreshold,min.pct,between),
                      "all" = .allMarkers(territories,counts,method,logFC,pvalThreshold,min.pct,between))


  return(markers)

}



.territoryMarkers <- function(territories,counts,method,logFC,pvalThreshold,min.pct,between = TRUE){
    #--------------------------------------------------------------------------#
    # Get counts for territories
    #--------------------------------------------------------------------------#
    all.counts <- lapply(territories, .getGeneCounts,"territory")

    #--------------------------------------------------------------------------#
    # Comparing all territories between each other
    #--------------------------------------------------------------------------#
    all.deg <- vector("list", length(all.counts))
    for(i in seq_along(all.counts)){
        message(paste0("Territoty ",i))
        #----------------------------------------------------------------------#
        # Comparing to all other territories individually
        #----------------------------------------------------------------------#
        seedCount <- all.counts[[i]]
        otherTer <- seq_along(all.counts)[seq_along(all.counts) != i]
        #----------------------------------------------------------------------#
        # Comparing seed territory to all other territories after combining
        # Might need to think about a better way of describing
        #----------------------------------------------------------------------#
        queryCount <- all.counts[otherTer]
        queryCount <- .bindCounts(queryCount,counts)
        deg <- .VesaliusDEGGlobal(seedCount,queryCount,method,i,"all")

        if(between){
            for(j in otherTer){
                #--------------------------------------------------------------#
                # Getting all genes with at least one count in either group
                #--------------------------------------------------------------#
                queryCount <- all.counts[[j]]
                degAll <- .VesaliusDEGGlobal(seedCount,queryCount,method,i,j)
            }
            deg <- rbind(deg,degAll)
        }
    }
    #--------------------------------------------------------------------------#
    # Filtering based on provided thresholds
    #--------------------------------------------------------------------------#
    deg <- deg %>% filter(seedPct >= min.pct | queryPct >= min.pct &
                          logFC >= logFC &
                          p.value <= pvalThreshold)

    return(deg)
}


.VesaliusDEGGlobal <- function(seedCount,queryCount,method,i,j){
    geneUnion <- union(rownames(seedCount), rownames(queryCount))

    #------------------------------------------------------------------#
    # Overlap between genes - commonGenes will stay as is
    # diff.genes we be converted to a zero matrix
    #------------------------------------------------------------------#

    seedDiff <- setdiff(rownames(seedCount),geneUnion)
    if(length(seedDiff)>0){
      seedDiffMat <-  Matrix(0,nrow = length(seedDiff), ncol = ncol(seedCount))
      rownames(seedDiffMat) <- seedDiff
      seedCount <- rbind(seedCount,seedDiffMat)
    }
    seedCount <- seedCount[order(rownames(seedCount)),]

    queryDiff <- setdiff(rownames(queryCount), geneUnion)
    if(length(queryDiff)>0){
      queryDiffMat <-  Matrix(0,nrow = length(queryDiff), ncol = ncol(queryCount))
      rownames(queryDiffMat) <- queryDiff
      queryCount <- rbind(queryCount,queryDiffMat)
    }

    queryCount <- queryCount[order(rownames(queryCount)),]
    #------------------------------------------------------------------#
    # computing percantge of cells containing each gene
    # comouting log Fold change
    #------------------------------------------------------------------#

    seedPct <- Matrix::rowSums(seedCount >0) / ncol(seedCount)
    queryPct <- Matrix::rowSums(queryCount >0) / ncol(queryCount)
    logFC <- log(Matrix::rowMeans(seedCount) - Matrix::rowMeans(queryCount))


    #------------------------------------------------------------------#
    # Diff genes expression in common gene pool
    #------------------------------------------------------------------#
    deg <- pblapply(rownames(seedCount),.VesaliusDEG,seedCount,queryCount,method)
    deg <- do.call("rbind",deg)

    #------------------------------------------------------------------#
    # Adding percentage and fold change
    # Adding territory info
    #------------------------------------------------------------------#

    deg$seedPct <- seedPct
    deg$queryPct <- queryPct
    deg$logFC <- logFC
    deg$seedTerritory <- i
    deg$compTerritory <- j

    return(deg)
}

.VesaliusDEG <- function(genes,group1,group2,method){
    #--------------------------------------------------------------------------#
    # Selecting genes in common for analysis
    #--------------------------------------------------------------------------#
    group1 <- group1[rownames(group1) == genes,]
    group2 <- group2[rownames(group2) == genes,]



    #--------------------------------------------------------------------------#
    # testing diff gene expression
    #--------------------------------------------------------------------------#
    deg <- switch(method,
                  "wilcox" = wilcox.test(group1,group2)$p.value,
                  "t.test" = t.test(group1,group2)$p.value)

    #--------------------------------------------------------------------------#
    # rebuilding data.frame
    #--------------------------------------------------------------------------#
    deg <- data.frame("genes" = genes,"p.value" = deg)


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
