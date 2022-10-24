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
#' @export


extract_markers <- function(vesalius_assay,
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
    "logit"),
  log_fc = 0.25,
  pval = 0.05,
  min_pct = 0.05,
  min_spatial_index = 10,
  verbose = TRUE,
  cores = 1,
  ...) {
    simple_bar(verbose)
    args <- list(...)
    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    counts <- check_norm(vesalius_assay, norm_method, method, verbose)
    #--------------------------------------------------------------------------#
    # Next let's get the territory data
    #--------------------------------------------------------------------------#
    ter <- check_territories(vesalius_assay, trial)
    seed <- check_group_value(ter, seed)
    query <- check_group_value(ter, query)
    deg_trial <- create_trial_tag(names(vesalius_assay@DEG), "DEG")
    #--------------------------------------------------------------------------#
    # Getting and setting territory categories
    #--------------------------------------------------------------------------#
    buffer <- dispatch_deg_group(ter, seed, query, cells, verbose)
    if (is.null(buffer$seed) || is.null(buffer$query)){
      return(NULL)
    }
    #--------------------------------------------------------------------------#
    # buffer will contain seed group(s) and query group(s) as well ids for each
    # group. We can loop over each to get Differentially expressed genes
    #--------------------------------------------------------------------------#
    message_switch("in_assay", verbose,
      comp_type = "Computing DEGs",
      assay = get_assay_names(vesalius_assay))
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
        verbose,
        args)
    }
    deg <- do.call("rbind", deg)
    deg <- list(deg)
    names(deg) <- deg_trial
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = deg,
    slot = "DEG",
    append = TRUE)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(extract_markers))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      get_assay_names(vesalius_assay))
    simple_bar(verbose)
    return(vesalius_assay)
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
  verbose = TRUE,
  args) {
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
    #--------------------------------------------------------------------------#
    params <- list("log_fc" = log_fc, "pval" = pval, "min_pct" = min_pct)
    deg <- switch(EXPR = method[1L],
      "wilcox" = vesalius_deg_wilcox(seed, seed_id, query, query_id, params),
      "t.test" = vesalius_deg_ttest(seed, seed_id, query, query_id, params),
      "chisq" = vesalius_deg_chisq(seed, seed_id, query, query_id, params),
      "fisher.exact" = vesalius_deg_fisher(seed, seed_id, query, query_id,
        params),
      "DESeq2" = vesalius_deg_deseq2(seed, seed_id, query, query_id, params),
      "edgeR" = vesalius_deg_edger(seed, seed_id, query, query_id, params),
      "logit" = vesalius_deg_logit(seed, seed_id, query, query_id, params))
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
    if (!identical(dim(dat),c(2,2))){
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
    if (!identical(dim(dat),c(2,2))) {
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

vesalius_deg_deseq2 <- function(seed, seed_id, query, query_id, params) {
  if (!require("DESeq2")) {
    stop("DESeq2 not installed!")
  }
  buffer <- get_deg_metrics(seed, query, params)
  seed <- buffer$seed
  query <- buffer$query
  #--------------------------------------------------------------------------#
  # we might need to add some more flexibility here by parsing DEseq specific 
  # paramters... 
  # maybe using the glmGamPoi option to increase speed ? 
  #--------------------------------------------------------------------------#
  deg <- format_counts_for_deseq2(seed, query)
  #--------------------------------------------------------------------------#
  # run deseq with recommended params for single cell data
  #--------------------------------------------------------------------------#
  deg <- DESeq2::DEseq(deg,
    test = "LRT",
    useT = TRUE,
    minmu = 1e-6,
    minReplicatesForReplace = Inf)
  deg <- DESeq2::results(
    object = deg,
    contrast = c("group", "seed", "query"),
    alpha = params$pval
  )
  deg <- tibble("genes" = rownames(deg),
    "p_value" = deg$pvalue,
    "p_value_adj" = deg$padj,
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = deg$log2FoldChange,
    "seed" = rep(seed_id, nrow(deg)),
    "query" = rep(query_id, nrow(deg)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
  
}

vesalius_deg_edger <- function(seed, seed_id, query, query_id, params) {
  if (!require("edgeR")) {
    stop("edgeR not installed!")
  }
  buffer <- get_deg_metrics(seed, query, params)
  seed <- buffer$seed
  query <- buffer$query
  #--------------------------------------------------------------------------#
  # Formatting to edgeR specifications 
  #--------------------------------------------------------------------------#
  deg <- format_counts_for_edger(seed, query)
  #--------------------------------------------------------------------------#
  # run edgeR - not sure if I could optimise these paramters
  # recommaned default - not much is said about single cell
  # TODO look into best paramters for edgeR for single cell data 
  #--------------------------------------------------------------------------#
  deg <- calcNormFactors(deg)
  design <- model.matrix(~group)
  deg <- estimateDisp(deg, design)
  fit <- glmFit(deg, design)
  lrt <- glmLRT(fit, coef = 2)
  deg <- topTags(lrt)
  deg <- tibble("genes" = deg$Symbol,
    "p_value" = deg$PValue,
    "p_value_adj" = deg$FDR,
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = deg$logFC,
    "seed" = rep(seed_id, nrow(deg)),
    "query" = rep(query_id, nrow(deg)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

vesalius_deg_logit <- function(seed, seed_id, query, query_id, params) {
  #--------------------------------------------------------------------------#
  # Setting up data for loigt regression as suggested by the signac package 
  # I don't want to rebuild a seurat object as they do in their vignette
  #--------------------------------------------------------------------------#
  buffer <- get_deg_metrics(seed, query, params)
  merged <- format_counts_for_logit(buffer$seed, buffer$query) 
  pvals <- sapply(seq_len(nrow(merged)), function(idx, merged) {
    model_data <- cbind("gene" = merged$merged[idx, ], merged$seed_query_info)
    gene_model <- as.formula("group ~ gene")
    gene_model <- glm(formula = gene_model,
      data = model_data,
      family = "binomial")
    null_model <- as.formula("group ~ 1")
    null_model <- glm(formula = null_model,
      data = model_data,
      family = "binomial")
    return(lrtest(gene_model, null_model)$Pr[2])
  }, merged = merged)
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
  # this can be handle by edgeR as well but let's stay consistent here
  # and keep this approach.
  # Maybe it would be a good idea to look into what is the best method
  # to select genes for DEG - with spatial componnent
  #--------------------------------------------------------------------------#
  keep <- (seed_pct >= params$min_pct |
    query_pct >= params$min_pct) &
    abs(fc) >= params$log_fc
    return(list("genes" = rownames(seed)[keep],
      "seed" = seed[keep, ],
      "query" = query[keep, ],
      "seed_pct" = seed_pct[keep],
      "query_pct" = query_pct[keep],
      "fc" = fc[keep]))

}
