################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Territory identity/---------------------------------#

#' extractAllMarkers compute differential gene expression for each territory
#' @param image a dataframe containing territory information
#' @param counts seurat object containing counts. Alternatively, matrix or
#' sparse matrix. Colnames as barcodes and rownames as genes.
#' @param method character describing the statistical test to use in order to
#' extract differential gene expression (currently only wilcox and t.test)
#' @param logFC numeric describing minimum log fold change value for
#' differential gene expression. Default set at 0.25.
#' @param pval numeric for pval threshold. Default set at 0.05
#' @param minPct numeric defining the minimum percentage of cells that should
#' contain any given gene.
#' @param minBar integer defining minimum number of barcodes in a territory.
#' @param morphologyFactor Integer or vector of integers describing growth or
#' shrink rate of territories.
#' @param cores numeric describing the number of cores used for the analysis.
#' @param verbose logical - progress message output
#' @details extractAllMarkers function compares a territory after any
#' morphological operator to all remaining barcodes.
#' If one of territory sets does not contain enough
#' barcodes, \code{extractMarkers} will return NULL. For example, if you dilate
#' a territory (let's say you want to grow each territory by 5 pixels -
#' morphologyFactor = 5),the seed territory will be dilated and compared to all
#' remaining barcodes that not part of this new dilated territory.
#'
#' For further details on territory Morphology - refer to
#' \code{territoryMorphology}
#' @return A data.frame (tibble) containing differential gene expression as well
#' p.value,
#' logFC,
#' seedPct (percentage of cells containing gene in first group),
#' queryPct (percentage of cells containing gene in second group),
#' seedTerritory (territory used as group 1)
#' queryTerritory (territory used as group 2)
#' @examples
#'\dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' mark <- extractAllMarkers(image, vealius )
#'
#' }

extractAllMarkers <- function(image,counts,method = "wilcox",
  logFC = 0.25, pval = 0.05,minPct = 0.05,minBar = 10,morphologyFactor=0,
  cores = 1,verbose = TRUE){
    .simpleBar(verbose)
    .checkCounts(verbose)
    #--------------------------------------------------------------------------#
    # Get counts
    #--------------------------------------------------------------------------#
    if(is(counts)=="Seurat"){
      counts <- GetAssayData(counts, slot = "data")
    }

    #--------------------------------------------------------------------------#
    # Let's split the image into seed territories
    #--------------------------------------------------------------------------#
    territories <- split(image,image$territory)
    tmp <- vector("list", length(territories))
    seedTerritory <- names(territories)
    #--------------------------------------------------------------------------#
    # For each territory, dilate and subset
    #--------------------------------------------------------------------------#
    for(i in seq_along(territories)){

        #----------------------------------------------------------------------#
        # Dilation factor can only be postive - now erosion for now
        #----------------------------------------------------------------------#
        if(morphologyFactor==0){
            seed <- territories[[i]] %>%
              filter(tile ==1 & territory %in% seedTerritory[i])
        } else {
            seed <- territories[[i]] %>% filter(territory %in% seedTerritory[i])
            seed <- territoryMorphing(seed,morphologyFactor,image,verbose)
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
            logFC,pval,minPct,minBar,verbose)
        #if(is.null(deg)) deg <- "empty"
        tmp[[i]] <- deg
      }
      #------------------------------------------------------------------------#
      # find nulls
      #------------------------------------------------------------------------#
      nulls <- sapply(tmp,is.null)

      tmp <- tmp[!nulls]
      tmp <- do.call("rbind",tmp)
      return(tmp)

  }

#' extractMarkers compute differential gene expression between selected
#' territories.
#' @param image a dataframe containing territory information
#' @param counts seurat object containing counts. Alternatively, matrix or
#' sparse matrix. Colnames as barcodes and rownames as genes.
#' @param seed Integer or vector of integers describing territories to be
#' included in group 1 for differential gene expression analysis.
#' @param query Integer or vector of integers describing territories to be
#' included in group 2 for differential gene expression analysis. Default = NULL
#' @param method character describing the statistical test to use in order to
#' extract differantial gene expression (currently only wilcox and t.test)
#' @param logFC numeric describing minimum log fold change value for
#' differential gene expression. Default set at 0.25.
#' @param pval numeric for pval threshold. Default set at 0.05
#' @param minPct numeric defining the minimum percentage of cells that should
#' contain any given gene.
#' @param minBar integer defining minimum number of barcodes in a territory.
#' @param morphologyFactorSeed Integer or vector of integers describing growth
#' or shrink rate of seed territories
#' @param morphologyFactorQuery Integer or vector of integers describing growth
#' or shrink rate of query territories
#' @param verbose logical - progress message output
#' @details extractMarkers compares a set of selected territories to another set
#' of selected territories. If one of territory sets does not contain enough
#' barcodes, \code{extractMarkers} will return NULL.
#'
#' A seed territory (or territories) should always be
#' provided. Query territories can be left as the default NULL. In this case,
#' the seed territory will be compared to all remaining barcodes that are NOT
#' present in the seed. Otherwise, seed territories will be compared to query
#' territories after any morphological operations
#' (see \code{territoryMorphology}.
#'
#' extractMarkers provides a way to manipulate each set of territories
#' independtly as described in \code{morphologyFactorSeed} and
#' \code{morphologyFactorQuery.}
#'
#' @return A data.frame (tibble) containing differential gene expression as well
#' p.value,
#' logFC,
#' seedPct (percentage of cells containing gene in first group),
#' queryPct (percentage of cells containing gene in second group),
#' seedTerritory (territory used as group 1)
#' queryTerritory (territory used as group 2)
#' @examples
#'\dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' markers <- extractMarkers(image, vesalius, seed = 1, query = 2)
#' }


extractMarkers <- function(image,counts,seed = NULL,query = NULL,
     method = "wilcox",logFC = 0.25, pval = 0.05,minPct = 0.05,minBar = 10,
     morphologyFactorSeed = 0,morphologyFactorQuery = 0,verbose=TRUE){
    .simpleBar(verbose)
    .checkCounts(verbose)
    #--------------------------------------------------------------------------#
    # Get counts
    #--------------------------------------------------------------------------#
    if(is(counts)=="Seurat"){
      counts <- GetAssayData(counts, slot = "data")
    }

    #--------------------------------------------------------------------------#
    # setting up seed territory counts
    #--------------------------------------------------------------------------#
    .seedSelect(verbose)
    if(is.null(seed)){
        stop("Please supply a seed as numeric values describing territory
             identity")
    }
    #--------------------------------------------------------------------------#
    # getting tmp to use as template when getting query territories
    #--------------------------------------------------------------------------#
    tmp <- seed
    #--------------------------------------------------------------------------#
    # Dilating if needed
    #--------------------------------------------------------------------------#
    if(morphologyFactorSeed == 0){
        seedID <- paste0(seed, sep =" ",collapse="")
        seed <- image %>% filter(tile == 1) %>% filter(territory %in% seed)
        seed <- counts[,colnames(counts) %in% unique(seed$barcodes)]
    } else {
        seedID <- paste0(seed, sep =" ",collapse="")
        seed <- image %>% filter(territory %in% seed)
        seed <- territoryMorphing(seed,morphologyFactorSeed,image, verbose)
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
        query <- image %>% filter(territory %in% query)
    }
    if(morphologyFactorQuery == 0){
        query <- image %>% filter(tile == 1)
        query <- counts[,colnames(counts) %in% unique(query$barcodes)]
        query <- counts[,!colnames(query) %in% colnames(seed)]
    } else {
        query <- image %>% filter(territory %in% query)
        query <- territoryMorphing(query,morphologyFactorQuery,image, verbose)
        query <- filter(query, tile == 1)
        query <- counts[,colnames(counts) %in% unique(query$barcodes)]
        query <- counts[,!colnames(query) %in% colnames(seed)]
    }

    #--------------------------------------------------------------------------#
    # running grouped analysis for DEG
    #--------------------------------------------------------------------------#
    .degProg(verbose)
    deg <- .VesaliusDEG(seed,query,seedID,queryID,method,
        logFC,pval,minPct,minBar,verbose)
    cat("\n")
   .simpleBar(verbose)
   return(deg)

}


#' extractClusterMarkers compute differential gene expression between clusters
#' and remaning barcodes
#' @param cluster a Seurat object containing clusters of an isolated territory
#' @param counts seurat object containing counts. Alternatively, matrix or
#' sparse matrix. Colnames as barcodes and rownames as genes.
#' @param method character describing the statistical test to use in order to
#' extract differential gene expression (currently only wilcox and t.test)
#' @param logFC numeric describing minimum log fold change value for
#' differential gene expression. Default set at 0.25.
#' @param pval numeric for pval threshold. Default set at 0.05
#' @param minPct numeric defining the minimum percentage of cells that should
#' contain any given gene.
#' @param minBar integer defining minimum number of barcodes in a territory.
#' @param verbose logical - progress message output
#' @details extractClusterMarkers compares clusters to all remaining barcodes.
#' If a territory is isolated (see \code{extractTerritories}) for further
#' analysis, Seurat can provide a cluster analysis of the isolated territory.
#' Seurat will compared clusters between each other within the isolated
#' territory. However, in some cases, it can be useful to compared the barcodes
#' present in a Seurat cluster to all remaining barcodes. This provides overall
#' differences in expression within the cluster as opposed to differential
#' expression between isolated clusters.
#'
#' @return A data.frame (tibble) containing differential gene expression as well
#' p.value,
#' logFC,
#' seedPct (percentage of cells containing gene in first group),
#' queryPct (percentage of cells containing gene in second group),
#' seedTerritory (territory used as group 1)
#' queryTerritory (territory used as group 2)
#' @examples
#'\dontrun{
#' data("vesalius")
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' ter <- extractTerritories(image, vesalius, seedID = 1)
#' ter <- ter %>% SCTransform(assay = "Spatial") %>%
#'              RunPCA(dims = 1:30) %>% RunUMAP(dims = 1:30) %>%
#'              FindNeighbors()%>% FindClusters(resolution = 0.3)
#' clusterMarkers <- extractClusterMarkers(ter,vesalius)
#' }
extractClusterMarkers <- function(cluster,counts,
     method = "wilcox",logFC = 0.25, pval = 0.05,minPct = 0.05,minBar = 10,
     verbose=TRUE){
    .simpleBar(verbose)
    .checkCounts(verbose)
    #--------------------------------------------------------------------------#
    # Get counts
    #--------------------------------------------------------------------------#
    if(is(counts)=="Seurat"){
      counts <- GetAssayData(counts, slot = "data")
    }

    #--------------------------------------------------------------------------#
    # Get data relating the seurat clusters
    # And coordinates - Not in use at the moment but would be useful if
    # I want to add some morphology for each cluster later on
    #--------------------------------------------------------------------------#
    #coord <-.getSeuratCoordinates(CACluster)
    cluster <- FetchData(cluster, c("seurat_clusters"))

    #--------------------------------------------------------------------------#
    # Now we can loop over each cluster and compare each cluster to everything
    # sounds stupid but other wise you won't always get the proper markers
    # You get subtle difference when isolating territories and that is great
    # but sometimes you have to compare to everything else for markers
    #--------------------------------------------------------------------------#
    clustNum <- unique(cluster$seurat_clusters)
    deg <- vector("list", length(clustNum))
    for(i in seq_along(clustNum)){
        tmp <- cluster %>% filter(seurat_clusters == clustNum[i]) %>% rownames()
        seed <- counts[,colnames(counts) %in% tmp]
        query <- counts[,!colnames(counts) %in% tmp]
        .degProg(verbose)
        deg[[i]] <- .VesaliusDEG(seed,query,clustNum[i],"Remaining",method,
            logFC,pval,minPct,minBar,verbose)
    }
    deg <- do.call("rbind", deg)
    cat("\n")
    .simpleBar(verbose)
    return(deg)



}



# Internal differantial gene expression function. Essentially all other DEG
# Functions will call this one function to do diff analysis.
# seed = group1 count data
# query = group 2 count data
# seedID = territory ID's for group 1
# queryID = territory ID's for group 2
# method = DEG stat method
# logFC = fold change threshold
# pval = p value threshold
# minPct = minimum percentage of barcodes that should contain a given gene
# minBar = minimum number of barcodes present in a territory
# verbose  = progress message output

.VesaliusDEG <- function(seed,query,seedID,queryID,method,logFC,
    pval,minPct,minBar,verbose =T){
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
    } else if(dimSeed[2L] < minBar){
        warning(paste0("Territory ",seedID," does not contain enough cells.\n
        Territory will be skipped"),call.=FALSE)
        return(NULL)
    } else if(dimQuery[2L] < minBar){
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





#' extractTerritories provides a Seurat object from a Vesalius Territory
#' @param image a data.frame containing territory information
#' @param seurat Seurat object. Contains all barcodes present on the spatial
#' transcrtiptomic assay.
#' @param seedID integer or vector of integers describing territories that
#' should be extracted.
#' @param morphologyFactor Integer or vector of integers describing growth or
#' shrink rate of territories. Only applied if seedID is provided.
#' @param minBar integer describing minimum number of barcodes in a territory
#' set.
#' @param verbose logical - progress message output
#' @param cores numeric describing the number of cores used for the analysis.
#' @details Once territories have been isolated from Vesalius images, the
#' \code{extractTerritories} provides a convenient way to extract a set of
#' territories, apply morphological operators, and return a Seurat object for
#' further analysis. This Seurat Object can be analysed using Seurat's
#' recommended pipeline.
#' @return A Seurat oject containing barcodes taken from territory set.
#' @examples
#' \dontrun{
#' data(vesalius)
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' # No dilation
#' ter <- extractTerritories(image, vesalius, seedID = 1)
#' # With dilation
#' ter <- extractTerritories(image, vesalius, seedID = 1, morphologyFactor = 5)
#' }
#'

extractTerritories <- function(image,seurat,seedID = NULL,morphologyFactor = 0,
                               minBar = 10, verbose = TRUE,cores = 1){
    .simpleBar(verbose)
    #--------------------------------------------------------------------------#
    # if combine is null we consider that we want to have all
    # territories prepared for clustering analysis
    #--------------------------------------------------------------------------#
    if(is.null(seedID)){
        .extractTerProg("all",verbose)
        territories <- image %>% filter(tile == 1) %>%
            distinct(barcodes, .keep_all = FALSE)
        territories <- split(territories, territories$territory)
        #----------------------------------------------------------------------#
        # Revert back to fine grain if inage array was reduced
        #----------------------------------------------------------------------#
        if(morphologyFactor != 0){
           cat("No seedID specified -
               Territory Dilation will not be applied.\n")
        }
        cells <- sapply(territories,nrow) > minBar
        territories <- territorries[cells]
        territories <- lapply(territories,"$",territory)
        territories <- parallel::mclapply(territories,subSetTerritories,seurat,
                                          mc.cores= cores)
        .simpleBar(verbose)
        return(territories)
    } else {
        #----------------------------------------------------------------------#
        # Filtering only the barcodes that are in the territories of interest
        # To ADD check for combine - numeric values
        #----------------------------------------------------------------------#
        .extractTerProg(seedID,verbose)
        territories <- image %>% filter(territory %in% seedID)
        territories <- territoryMorphing(territories,morphologyFactor,image,
                                         verbose)
        territories <- filter(territories,tile == 1)  %>%
            distinct(barcodes,.keep_all=F)
        territories <- territories$barcodes

        #----------------------------------------------------------------------#
        # Revert back to fine graind and return null if territory does not
        # Contain enough cells - in this case throw in a warning
        #----------------------------------------------------------------------#
        #territories <- .fineGrain(territories,minBar)

        if(length(territories) < minBar){

            warning("Territory selection does not contain enough cells -
                    NULL returned",
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



#' compareClusters computes differential gene expression between seurat clusters
#' @param counts seurat object containing ALL counts. Alternatively, matrix or
#' sparse matrix of all counts. Colnames as barcodes and rownames as genes.
#' @param seedCluster Seurat Object of seed clusters
#' @param queryCluster Seurat Object of query clusters
#' @param seed integer describing the Seurat cluster (in seed - group 1) to be
#' used for Differential expression analysis. Default NULL
#' @param query integer describing the Seurat cluster (in query - group 2) to be
#' used for Differential expression analysis. Default NULL
#' @param method character describing the statistical test to use in order to
#' extract differential gene expression (currently only wilcox and t.test)
#' @param logFC numeric describing minimum log fold change value for
#' differential gene expression. Default set at 0.25.
#' @param pval numeric for pval threshold. Default set at 0.05
#' @param minPct numeric defining the minimum percentage of cells that should
#' contain any given gene.
#' @param minBar integer defining minimum number of barcodes in a territory.
#' @param verbose logical - progress message output.
#' @details Territory isolation provides a convenient way to investigate
#' the finer detail of anatomy and cellular spatial distribution. In some case,
#' two territories separated in 2D space may contain barcodes that cluster
#' together and are described by similar cell type labels. The
#' \code{compareClusters} functions compares clusters in different territories
#' between each other. The count values are the overall count values that have
#' been normalized, scaled and centered all together. This ensure that the
#' barcodes are comparable. The territory isolation and clustering serves as a
#' way to isolate a specific set of barcodes that will be used for differential
#' gene expression analysis.
#' @return A data.frame (tibble) containing differential gene expression as well
#' p.value,
#' logFC,
#' seedPct (percentage of cells containing gene in first group),
#' queryPct (percentage of cells containing gene in second group),
#' seedTerritory (territory used as group 1)
#' queryTerritory (territory used as group 2)
#' @examples
#'\dontrun{
#' data("vesalius")
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' ter1 <- extractTerritories(image, vesalius, seedID = 1)
#' ter1 <- ter1 %>% SCTransform(assay = "Spatial") %>%
#'              RunPCA(dims = 1:30) %>% RunUMAP(dims = 1:30) %>%
#'              FindNeighbors()%>% FindClusters(resolution = 0.3)

#' ter2 <- extractTerritories(image, vesalius, seedID = 2)
#' ter2 <- ter2 %>% SCTransform(assay = "Spatial") %>%
#'              RunPCA(dims = 1:30) %>% RunUMAP(dims = 1:30) %>%
#'              FindNeighbors()%>% FindClusters(resolution = 0.3)
#' clusterComp <- compareClusters(vesalius, seedCluster  = ter1,
#' queryCluster = ter2)
#' # Comparing Specific clusters
#' clusterComp <- compareClusters(vesalius, seedCluster  = ter1,
#' queryCluster = ter2, seed = 1, query = 2)
#' }

compareClusters <- function(counts,seedCluster,queryCluster,seed = NULL,
                            query = NULL,method = "wilcox",
  logFC = 0.25, pval = 0.05,minPct = 0.05,minBar = 10,verbose=TRUE){
      .simpleBar(verbose)
      #------------------------------------------------------------------------#
      # Get counts
      #------------------------------------------------------------------------#
      .checkCounts(verbose)
      if(is(counts)=="Seurat"){
        counts <- GetAssayData(counts, slot = "data")
      }
      #------------------------------------------------------------------------#
      # First let's set up seed cells
      # For now we will consider that if seed is null then you want to compare
      # to everything. That is the same as comparing territories
      #------------------------------------------------------------------------#
      .seedSelect(verbose)
      if(!is.null(seed)){
        seedCluster <- FetchData(seedCluster,c("seurat_clusters")) %>%
        filter(seurat_clusters %in% seed)
        seedCluster <- counts[,rownames(seedCluster)]
        seedID <- paste0(seed, collapse ="",sep=" ")
      } else {
        seedCluster <- FetchData(seedCluster,c("seurat_clusters"))
        seedCluster <- counts[,rownames(seedCluster)]
        seedID <- "All"
      }

      #------------------------------------------------------------------------#
      # Same thing with query cells
      # For now we will consider that if seed is null then you want to compare
      # to everything. That is the same as comparing territories
      #------------------------------------------------------------------------#
      .querySelect(verbose)
      if(!is.null(query)){
          queryCluster <- FetchData(queryCluster,c("seurat_clusters")) %>%
          filter(seurat_clusters %in% query)
          queryCluster <- counts[,rownames(queryCluster)]
          queryID <- paste0(query,sep = " ", collapse ="")
      } else {
          queryCluster <- FetchData(queryCluster,c("seurat_clusters"))
          queryCluster <- counts[,rownames(queryCluster)]
          queryID <- "All"
      }

      #------------------------------------------------------------------------#
      # Setting up dat in right format to be parsed to VesDeg
      #------------------------------------------------------------------------#

      .degProg(verbose)
      deg<- .VesaliusDEG(seed = seedCluster, query= queryCluster,
                         seedID =seedID,
                         queryID = queryID,
                         logFC= logFC,pval = pval,
                         minPct = minPct,minBar = minBar,
                         method=method,
                         verbose = verbose)
      cat("\n")
      .simpleBar(verbose)
      return(deg)
  }


#' compareLayers computes differential gene expression between territory layers
#' @param layers data.frame containing layered territory information
#' (See \code{layerTerritory.edge} and \code{layerTerritory.concave})
#' @param counts seurat object containing counts. Alternatively, matrix or
#' sparse matrix. Colnames as barcodes and rownames as genes.
#' @param l1 integer or vector of integers indicating layers to be contained in
#' group 1. Default NULL - If NULL will compare all layers independently
#' @param l2 integer or vector of integers indicating layers to be contained in
#' group 2. If NULL will compare all layers independently
#' @param method character describing the statistical test to use in order to
#' extract differential gene expression (currently only wilcox and t.test)
#' @param logFC numeric describing minimum log fold change value for
#' differential gene expression. Default set at 0.25.
#' @param pval numeric for pval threshold. Default set at 0.05
#' @param minPct numeric defining the minimum percentage of cells that should
#' contain any given gene.
#' @param minBar integer defining minimum number of barcodes in a territory.
#' @param verbose logical - progress message output.
#' @details Territory isolation provides a convenient way to investigate
#' the finer detail of anatomy and cellular spatial distribution. Image
#' representation of territories provides a convenient way to manipulate
#' territories such as dividing a territory into layers. Vesalius provides two
#' method to layer a territory : \code{layerTerritory.edge} and
#' \code{layerTerritory.concave}. The output of both function returns a data
#' frame with layer information. The \code{compareLayers} compares a layer or
#' layer sets between each other. If both l1 and l2 layers are left to NULL,
#' each layer will be compared to each other layer on at a time.
#' To compare one layer to all other layers, specific layer ID's should be
#' specified.
#' @return A data.frame (tibble) containing differential gene expression as well
#' p.value,
#' logFC,
#' seedPct (percentage of cells containing gene in first group),
#' queryPct (percentage of cells containing gene in second group),
#' seedTerritory (territory used as group 1)
#' queryTerritory (territory used as group 2)
#' @examples
#'\dontrun{
#' data("vesalius")
#' # Seurat pre-processing
#' image <- NormalizeData(vesalius)
#' image <- FindVariableFeatures(image, nfeatures = 2000)
#' image <- ScaleData(image)
#' # converting to rgb
#' image <- rgbPCA(image,slices = 1)
#' image <- buildImageArray(image, sliceID=1)
#' # One segmentation round
#' image <- iterativeSegmentation.array(image)
#' image <- isolateTerritories.array(image, minBar = 5)
#' layer <- layerTerritory.edge(image, seedTerritory = 1)
#' layerComp <- compareLayers(layer, vesalius, l1 = 1, l2 = 2)
#' }

compareLayers <- function(layers,counts,l1 = NULL, l2 = NULL, method = "wilcox",
    logFC = 0.25, pval = 0.05,minPct = 0.05,minBar = 10,verbose=TRUE){
      #------------------------------------------------------------------------#
      # Ayo now we comapre layers
      # First thing is to get layers so lets start with layer 1
      # Will need to update so there is no overlap.
      # relevant if you set one of the l1/l2 layers as NULL
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
      #------------------------------------------------------------------------#
      # Get counts
      #------------------------------------------------------------------------#
      if(is(counts)=="Seurat"){
        counts <- GetAssayData(counts, slot = "data")
      }

      #------------------------------------------------------------------------#
      # Now we can loop over layers
      #------------------------------------------------------------------------#
      deg <- vector("list", length(l1) * length(l2))
      counter <- 1
      for(i in seq_along(l1)){
          for(j in seq_along(l2)){

              tmp1 <- counts[,colnames(counts) %in% unique(l1[[i]]$barcodes)]
              tmp2 <- counts[,colnames(counts) %in% unique(l2[[j]]$barcodes)]
              tmpID1 <- paste0(unique(l1[[i]]$layer),collapse = "",sep = " ")
              tmpID2 <- paste0(unique(l2[[j]]$layer),collapse = "",sep = " ")
              deg[[counter]] <- .VesaliusDEG(tmp1,tmp2,tmpID1,tmpID2,
                                                method = method,
                                                logFC=logFC,pval = pval,
                                                minPct= minPct,minBar=minBar,
                                                verbose=verbose)
          }
      }
      deg <- do.call("rbind",deg)
      cat("\n")
      .simpleBar(verbose)
      return(deg)
  }
