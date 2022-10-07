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


extract_markers <- function(vesalius,
  trial = "last",
  norm_method = "last",
  seed = NULL,
  query = NULL,
  cells = NULL,
  method = c("wilcox",
    "t.test",
    "chisq",
    "fisher.exact",
    "DEseq2",
    "edgeR",
    "ArchR"),
  log_fc = 0.25,
  pval = 0.05,
  min_pct = 0.05,
  min_spatial_index = 10,
  verbose = TRUE) {
    simple_bar(verbose)

    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    counts <- check_norm(vesalius, norm_method, method, verbose)
    #--------------------------------------------------------------------------#
    # Next let's get the territory data
    #--------------------------------------------------------------------------#
    ter <- check_territories(vesalius, trial)
    deg_trial <- create_trial_tag(names(vesalius@DEG), "DEG")
    #--------------------------------------------------------------------------#
    # Getting and setting territory categories
    #--------------------------------------------------------------------------#
    buffer <- deg_group_dispatch(ter, seed, query, cells, verbose)
    #--------------------------------------------------------------------------#
    # buffer will contain seed group(s) and query group(s) as well ids for each
    # group. We can loop over each to get Differentially expressed genes
    #--------------------------------------------------------------------------#
    message_switch("deg_prog", verbose)
    deg <- vector("list", length(seed))
    for (i in seq_along(buffer$seed)) {
      seed <- counts[, colnames(counts) %in% buffer$seed[[i]]]
      seed_id <- buffer$seed_id[i]

      query <- counts[, colnames(counts) %in% buffer$query[[i]]]
      query_id <- buffer$query_id[i]

      deg[[i]] <- vesalius_deg(seed,
        query,
        seed_id,
        query_id,
        method,
        log_fc,
        pval,
        min_pct,
        min_spatial_index,
        verbose)
    }
    deg <- do.call("rbind", deg)
    deg <- list(deg)
    names(deg) <- deg_trial
    vesalius <- update_vesalius(vesalius = vesalius,
      data = deg,
      slot = "DEG",
      commit = as.list(match.call()),
      defaults = as.list(args(extract_markers)),
      append = TRUE)


    #--------------------------------------------------------------------------#
    # running grouped analysis for DEG
    #--------------------------------------------------------------------------#

  cat("\n")
   simple_bar(verbose)
   return(vesalius)
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

vesalius_deg <- function(seed,
  query,
  seed_id,
  query_id,
  method,
  log_fc,
  pval,
  min_pct,
  min_spatial_index,
  verbose = TRUE) {
    message_switch("deg_prog_each", verbose, seed = seed_id, query = query_id)
    #--------------------------------------------------------------------------#
    # We assume here that we are parsing cleaned up version of each object
    # This is strictly just computing DEG between two groups
    #
    # Just in case there are not enough cells
    #--------------------------------------------------------------------------#
    seed <- check_min_spatial_index(seed, min_spatial_index, seed_id)
    query <- check_min_spatial_index(query, min_spatial_index, query_id)
    
    #--------------------------------------------------------------------------#
    # testing diff gene expression
    # Rebuilding a matrix just in case you only have one gene to test
    # Not going to rewrite the sapply for that
    #--------------------------------------------------------------------------#
    params <- list("log_fc" = log_fc, "pval" = pval, "min_pct" = min_pct)
    deg <- switch(EXPR = method,
      "wilcox" = vesalius_deg_wilcox(seed, query, params),
      "t.test" = vesalius_deg_ttest(seed, query, params),
      "chisq" = vesalius_deg_chisq(seed, query, params),
      "fisher.exact" = vesalius_deg_fisher(seed, query, params),
      "DESeq2" = vesalius_deg_deseq2(seed, query, params),
      "edgeR" = vesalius_deg_edger(seed,query, params),
      "ArchR" = vesalius_deg_archr(seed, query, params))


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

vesalius_deg_wilcox <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  pvals <- sapply(seq_len(nrow(buffer$seed)), function(idx, seed, query) {
    return(wilcox.test(seed[idx,],query[idx,])$p.value)
  }, seed = buffer$seed, query = buffer$query)
  deg <- tibble("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

vesalius_deg_ttest <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  pvals <- sapply(seq_len(nrow(buffer$seed)), function(idx, seed, query) {
    return(t.test(seed[idx,],query[idx,])$p.value)
  }, seed = buffer$seed, query = buffer$query)
  deg <- tibble("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

vesalius_deg_chisq <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  pvals <- sapply(seq_len(nrow(buffer$seed)), function(idx, seed, query) {
    dat <- cbind(table(seed[idx, ]), table(query[, idx]))
    if (!identical(dim(dat),c(2,2)) {
      warning("Count data is not binary! 
        Contingency table for Chi-squared test not in 2 x 2 format")
    }
    return(chisq.test(dat)$p.value)
  }, seed = buffer$seed, query = buffer$query)
  deg <- tibble("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

vesalius_deg_fisher <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  pvals <- sapply(seq_len(nrow(buffer$seed)), function(idx, seed, query) {
    dat <- cbind(table(seed[idx, ]), table(query[, idx]))
    if (!identical(dim(dat),c(2,2)) {
      stop("Count data is not binary! 
        Cannot create 2 x 2 contingency table for fisher's exact test")
    }
    return(fisher.test(dat)$p.value)
  }, seed = buffer$seed, query = buffer$query)
  deg <- tibble("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

get_deg_metrics <- function(seed, query, params) {
  #--------------------------------------------------------------------------#
  # this asumes that we receive normalised counts
  # NEED to check and update for TFID and SCTRANSFORM
  #--------------------------------------------------------------------------#
  seed_pct <- rowSums(seed > 0) / ncol(seed)
  query_pct <- rowSums(query > 0) / ncol(query)
  fc <- rowMeans(seed) - rowMeans(query)
  #--------------------------------------------------------------------------#
  # Dropping genes that don't fit the logFC and pct criteria
  #--------------------------------------------------------------------------#
  keep <- (seed_pct >= params$min_pct ||
    query_pct >= params$min_pct) &&
    abs(fc) >= params$log_fc
    return(list("genes" = rownames(seed)[keep],
      "seed" = seed[keep, ],
      "query" = query[keep, ],
      "seed_pct" = seed_pct[keep],
      "query_pct" = query_pct[keep],
      "fc" = fc))

}
