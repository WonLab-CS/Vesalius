################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Territory identity/---------------------------------#


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


extractMarkers <- function(vesalius,
                           trial = "last",
                           normMethod = "last",
                           seed = NULL,
                           query = NULL,
                           cells = NULL,
                           method = c("wilcox", "t.test", "chisq"),
                           logFC = 0.25,
                           pval = 0.05,
                           minPct = 0.05,
                           minBar = 10,
                           verbose = TRUE) {
    .simpleBar(verbose)
    .consSparse(verbose)

    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    if (normMethod == "last") {
        counts <- as.matrix(vesalius@counts[[length(vesalius@counts)]])

    } else {
        if (length(grep(x = names(vesalius@counts),
                pattern = normMethod)) == 0) {
            stop(paste0(deparse(substitute(normMethod)),
                "is not in count list!"))
        }
        counts <- as.matrix(vesalius@counts[[normMethod]])
    }
    #--------------------------------------------------------------------------#
    # Next let's get the territory data
    #--------------------------------------------------------------------------#

    if (!is.null(vesalius@territories) && trial == "last") {
      trial <- colnames(vesalius@territories)[ncol(vesalius@territories)]
      DEG_trial <- paste0("DEG_", trial)
      ter <- vesalius@territories[, c("barcodes", "x", "y", trial)]
      colnames(ter) <- c("barcodes", "x", "y", "trial")
    } else if (!is.null(vesalius@territories) && trial != "last") {
      if (!any(grepl(x = colnames(vesalius@territories), pattern = trial))) {
          stop(paste(deparse(substitute(trial)),
            "is not in territory data frame"))
      }
      DEG_trial <- paste0("DEG_",trial)
      ter <- vesalius@territories[, c("barcodes", "x", "y", trial)]
      colnames(ter) <- c("barcodes", "x", "y", "trial")
    } else {
      stop("No territories have been computed! Cannot perform DEG.")
    }

    #------------------------------------------------------------------------#
    # Getting and setting territory categories
    #------------------------------------------------------------------------#
    if (is.null(seed) && is.null(query)) {
        #----------------------------------------------------------------------#
        # if no seed nor query is specified - assume that we want everything
        #----------------------------------------------------------------------#
        seed <- split(ter$barcodes, ter$trial)
        query <- lapply(seed, function(bar, ter) {
            return(ter$barcodes[!ter$barcodes %in% bar])
        }, ter = ter)
        names(query) <- rep("remaining", length(query))
        if (!is.null(cells)) {
            seed <- lapply(seed,function(s, c) {
              return(s[s %in% c])
              }, c = cells)
            if (length(seed) == 0) {
              stop("No cells of interest are present in seed territory") 
            }
            query <- lapply(query,function(s,c){
                return(s[s %in% c])
                },c = cells)
            if (length(query) == 0) {
              stop("No cells of interest are present in query territory")
            }
        }

    } else if (!is.null(seed) && is.null(query)) {
        #----------------------------------------------------------------------#
        # if only seed we compare seed to everything else
        # Get initial seed territory
        #----------------------------------------------------------------------#
        seedID <- paste0(seed, collapse=" ")
        seed <- ter[ter$trial %in% seed, "barcodes"]
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        queryID <- "remaining"
        query <- ter[!ter$barcodes %in% seed, "barcodes"]
        #----------------------------------------------------------------------#
        # If cell barcodes are present we filter based on cells as well
        #----------------------------------------------------------------------#
        if (!is.null(cells)) {
            seed <- seed[seed %in% cells]
            if (length(seed) == 0) {
              stop("No cells of interest are present in seed territory")
            }
            query <- query[query %in% cells]
            if (length(query) == 0) {
              stop("No cells of interest are present in query territory")
            }
        }


    } else if (is.null(seed) && !is.null(query)) {
        #----------------------------------------------------------------------#
        # if only query we compare query to everything else
        # Get initial query territory
        #----------------------------------------------------------------------#
        queryID <- paste0(query, collapse = " ")
        query <- ter[ter$trial %in% query, "barcodes"]
        #----------------------------------------------------------------------#
        # Filter seed based on query
        #----------------------------------------------------------------------#
        seed <- "remaning"
        seed <- ter[!ter$barcodes %in% query, "barcodes"]
        #----------------------------------------------------------------------#
        # If cell barcodes are present we filter based on cells as well
        #----------------------------------------------------------------------#
        if (!is.null(cells)) {
            seed <- seed[seed %in% cells]
            if (length(seed) == 0) {
              stop("No cells of interest are present in seed territory")
            }
            query <- query[query %in% cells]
            if (length(query) == 0) {
              stop("No cells of interest are present in query territory")
            }
        }

    } else {
        #----------------------------------------------------------------------#
        # if get both filter based on both
        #----------------------------------------------------------------------#
        seedID <- paste0(seed, collapse = " ")
        seed <- ter[ter$trial %in% seed, "barcodes"]
        #----------------------------------------------------------------------#
        # Filter query based on seed
        #----------------------------------------------------------------------#
        queryID <- paste0(query, collapse = " ")
        query <- ter[ter$trial %in% query, "barcodes"]
        #----------------------------------------------------------------------#
        # If cell barcodes are present we filter based on cells as well
        #----------------------------------------------------------------------#
        if (!is.null(cells)) {
            seed <- seed[seed %in% cells]
            if (length(seed) == 0) {
              stop("No cells of interest are present in seed territory")
            }
            query <- query[query %in% cells]
            if (length(query) == 0) {
              stop("No cells of interest are present in query territory")
            }
        }
    }
    #--------------------------------------------------------------------------#
    # Now we have some barcodes let's run DEG
    # We need to check what case we have
    # if list then loop between all territories - this includes layers
    # Else we can just parse the counts
    #--------------------------------------------------------------------------#
    if (any(is(seed) == "list") && any(is(query) == "list")) {
        #----------------------------------------------------------------------#
        # Now we can compared each territory to the rest
        # we won't will skip same territory comparison
        #----------------------------------------------------------------------#
        deg <- vector("list", length(seed))
        .degProg(verbose)
        for (i in seq_along(seed)) {

            seedTmp <- counts[, colnames(counts) %in% seed[[i]]]
            seedID <- names(seed)[i]

            queryTmp <- counts[, colnames(counts) %in% query[[i]]]
            queryID <- names(query)[i]


            deg[[i]] <- .VesaliusDEG(seedTmp,
                                queryTmp,
                                seedID,
                                queryID,
                                method,
                                logFC,
                                pval,
                                minPct,
                                minBar,
                                verbose)
        }
        deg <- do.call("rbind", deg)
    } else {
        #----------------------------------------------------------------------#
        # Getting count out of matrix and parsing to Internal
        #----------------------------------------------------------------------#
        seed <- counts[, colnames(counts) %in% seed]
        query <- counts[, colnames(counts) %in% query]
        .degProg(verbose)
        deg <- .VesaliusDEG(seed,
                            query,
                            seedID,
                            queryID,
                            method,
                            logFC,
                            pval,
                            minPct,
                            minBar,
                            verbose)
    }
    deg <- list(deg)
    names(deg) <- DEG_trial
    vesalius <- .updateVesalius(vesalius = vesalius,
                                data = deg,
                                slot = "DEG",
                                commit = as.list(match.call()),
                                defaults = as.list(args(extractMarkers)),
                                append = TRUE)


    #--------------------------------------------------------------------------#
    # running grouped analysis for DEG
    #--------------------------------------------------------------------------#

    cat("\n")
   .simpleBar(verbose)
   return(vesalius)

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
# Don't know what to do with these ones yet.
# Not sure I want to keep it for the next version.
# I am think to simplify - we are not trying to make a full pipeline
# Let that be done externaly.
.extractClusterMarkers <- function(cluster,counts,
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
    # This sections need to be rethought - logFC
    # Depending on what the data is you will need a different fold change
    # method. This section Need to be thought through also this will need to
    # change if not using seurat
    #--------------------------------------------------------------------------#
    #seedPct <- Matrix::rowSums(seed >0) / ncol(seed)
    seedPct <- rowSums(seed >0) / ncol(seed)
    #queryPct <- Matrix::rowSums(query >0) / ncol(query)
    queryPct <- rowSums(query >0) / ncol(query)
    #FC <- log(Matrix::rowMeans(seed) +1) - log(Matrix::rowMeans(query)+1)
    #FC <- Matrix::rowMeans(seed)- Matrix::rowMeans(query)
    FC <- rowMeans(seed)- rowMeans(query)
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

                  res <- switch(method[1L],
                         "wilcox" = wilcox.test(seed[idx,],query[idx,])$p.value,
                         # maybe add eq of variance for internal function
                         "t.test" = t.test(seed[idx,],query[idx,])$p.value,
                         "chisq" = chisq.test(cbind(table(seed[idx,]),
                                                    table(query[,idx])))$p.value)
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

    deg <- deg %>% filter(p.value.adj <= pval)
    return(deg)
}
