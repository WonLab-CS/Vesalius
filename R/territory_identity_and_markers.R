################################################################################
############################   ST package        ###############################
################################################################################

#-------------------------/Territory identity/---------------------------------#


#' identify_markers computes differential observation expression 
#' between selected territories.
#' @param vesalius_assay a vesalius_assay
#' @param trial character string - which territory trial that 
#' should be used to select
#' territorires. Default is last one computed
#' @param norm_method charcater string - which normalisation method should 
#' be used. 
#' @param seed Integer or vector of integers describing territories to be
#' included in group 1 for differential gene expression analysis.
#' @param query Integer or vector of integers describing territories to be
#' included in group 2 for differential gene expression analysis. Default = NULL
#' @param cells character vector containing barcodes of cells of interest. 
#' @param method character describing the statistical test to use in order to
#' extract differantial gene expression.
#' Select from:
#' "wilcox", "t.test", "chisq", "fisher.exact", "DEseq2", "QLF", "LRT","logit"
#' @param log_fc numeric describing minimum log fold change value for
#' differential gene expression. Default set at 0.25.
#' @param pval numeric for pval threshold. Default set at 0.05
#' @param min_pct numeric defining the minimum percentage of cells that should
#' contain any given gene. Deault set at 0.05
#' @param min_spatial_index integer defining minimum number of 
#' barcodes in a territory.
#' @param verbose logical - progress message output
#' @param ... other parameters parsed to DESeq2 or edgeR (not functional)
#' @details Identifying markers is a key aspect of spatial data analysis. 
#' This functions let's you select which territory trial you which to use.
#' Note that this can be any territory trial that you have run, including 
#' color segments, isolated territories, dilated or eroded territories and
#' layered territories. By default, \code{identify_markers} takes the last 
#' one that has been computed.
#' 
#' The normalisation method refers to the normalisation method applied to
#' the count matrices. If you use, DESeq2, QLF (edgeR) or LRT (edgeR), 
#' raw counts will be selected and will ignore any cother command. 
#' This is a requirement for both of these packages.
#' 
#' If you have some cells you are interested in comparing between territories,
#' you can simply parse a character vector containing the barcodes of your 
#' cells of interest. \code{identify_markers} will automatically retrieve cells
#' in each territory and only use these cell for comparison.
#' 
#' Once the Differentially expressed genes/oberservations have been computed 
#' they are stored in the vesalius_assay object. This allows you to run multiple
#' trial and have all these trials sorted within your object.
#' 
#' To retrieve them from the object, you can use \code{\link{get_markers}}
#'
#' @return a vesalius_assay object
#' @examples
#' \dontrun{
#' data(vesalius)
#' # First we build a simple object
#' ves <- build_vesalius_object(coordinates, counts)
#' # We can do a simple run
#' ves <- build_vesalius_embeddings(ves)
#'
#' # simple smoothing
#' ves <- smooth_image(ves, dimensions = seq(1, 30))
#' 
#' # quick segmentation
#' ves <- segment_image(ves, dimensions = seq(1, 30))
#' 
#' # isolate territories
#' ves <- isolate_territories(ves)
#' 
#' # identify markers
#' ves <- identify_markers(ves, seed = c(3,5), query = 8)
#' deg <- get_markers(ves)
#'}
#' @export


identify_markers <- function(vesalius_assay,
  trial = "last",
  norm_method = "last",
  seed = NULL,
  query = NULL,
  cells = NULL,
  method = "wilcox",
  log_fc = 0.25,
  pval = 0.05,
  min_pct = 0.05,
  min_spatial_index = 10,
  verbose = TRUE,
  ...) {
    simple_bar(verbose)
    args <- list(...)
    if (!is(vesalius_assay, "vesalius_assay")) {
      stop("Unsupported format to identify_markers function \n
      Please use a vesalius_assay object")
    }
    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    # qand check if DEG method is present 
    #--------------------------------------------------------------------------#
    method <- check_deg_method(method)
    counts <- check_norm(vesalius_assay, norm_method, method, verbose)
    #--------------------------------------------------------------------------#
    # Next let's get the territory data
    #--------------------------------------------------------------------------#
    ter <- check_territory_trial(vesalius_assay, trial)
    seed <- check_group_value(ter, seed)
    query <- check_group_value(ter, query)
    #--------------------------------------------------------------------------#
    # Getting and setting territory categories
    #--------------------------------------------------------------------------#
    buffer <- dispatch_deg_group(ter, seed, query, cells, verbose)
    if (is.null(buffer$seed[[1L]]) || is.null(buffer$query[[1L]])) {
      return(vesalius_assay)
    }
    #--------------------------------------------------------------------------#
    # buffer will contain seed group(s) and query group(s) as well ids for each
    # group. We can loop over each to get Differentially expressed genes
    # The loop is here because we assume that if both seed and query are null
    # the user wants to compare everytyhing to everything.
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
    names(deg) <- tail(create_trial_tag(names(vesalius_assay@DEG), "DEG"), 1)
    vesalius_assay <- update_vesalius_assay(vesalius_assay = vesalius_assay,
    data = deg,
    slot = "DEG",
    append = TRUE)
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(identify_markers))
    vesalius_assay <- commit_log(vesalius_assay,
      commit,
      get_assay_names(vesalius_assay))
    simple_bar(verbose)
    return(vesalius_assay)
}

#' Internal differantial gene expression function.
#' 
#' This function dispatches the groups to the various methods 
#' That are avialbale 
#' @param seed = group1 count data
#' @param query = group 2 count data
#' @param seed_id = territory ID's for group 1
#' @param query_id = territory ID's for group 2
#' @param method = DEG stat method
#' @param log_fc = fold change threshold
#' @param pval = p value threshold
#' @param min_pct = minimum percentage of barcodes that should contain a
#'  given gene
#' @param min_spatial_index = minimum number of barcodes present in a territory
#' @param verbose  = progress message output
#' @param args arguments parse to (...) in upper level function (not functional)

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
    # return NULL if either one is NULL
    #--------------------------------------------------------------------------#
    if (is.null(seed) || is.null(query)) {
      return(NULL)
    }
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
      "QLF" = vesalius_deg_edger(seed, seed_id, query, query_id, params, "QLF"),
      "LRT" = vesalius_deg_edger(seed, seed_id, query, query_id, params, "LRT"),
      "logit" = vesalius_deg_logit(seed, seed_id, query, query_id, params))
    return(deg)
}

#' wilcox rank sum test for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @importFrom stats wilcox.test
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
vesalius_deg_wilcox <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
  pvals <- sapply(seq_len(nrow(buffer$seed)), function(idx, seed, query) {
    return(suppressWarnings(wilcox.test(seed[idx, ], query[idx, ])$p.value))
  }, seed = buffer$seed, query = buffer$query)
  effect_size <- sapply(pvals, compute_effect_size,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

#' t.test for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @importFrom stats t.test
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
vesalius_deg_ttest <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
  pvals <- sapply(seq_len(nrow(buffer$seed)), function(idx, seed, query) {
    return(suppressWarnings(t.test(seed[idx, ], query[idx, ])$p.value))
  }, seed = buffer$seed, query = buffer$query)
  effect_size <- sapply(pvals, compute_effect_size,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

#' chisq test for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @importFrom stats chisq.test
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
vesalius_deg_chisq <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
  #buffer <- check_binary_nature(buffer)
  seed <- apply(buffer$seed, 1, sum)
  query <- apply(buffer$query,1, sum)
  pvals <- suppressWarnings(chisq.test(x = seed, y = query))$p.value
  effect_size <- compute_effect_size(pvals,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("cells" = paste(buffer$genes,collapse = " "),
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

#' Fisher's excat test for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @importFrom stats fisher.test
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
vesalius_deg_fisher <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
  #buffer <- check_binary_nature(buffer)
  seed <- apply(buffer$seed, 1, sum)
  query <- apply(buffer$query,1, sum)
  pvals <- suppressWarnings(fisher.test(x = seed, y = query))$p.value
  effect_size <- compute_effect_size(pvals,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("cells" = paste(buffer$genes,collapse = " "),
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

#' DESeq negatove binomail for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @importFrom DESeq2 DESeq results
#' @importFrom dplyr filter
vesalius_deg_deseq2 <- function(seed, seed_id, query, query_id, params) {
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
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
  deg <- DESeq2::DESeq(deg,
    test = "LRT",
    useT = TRUE,
    minmu = 1e-6,
    minReplicatesForReplace = Inf,
    fitType = "local",
    reduced = ~1,
    quiet = TRUE)
  deg <- DESeq2::results(
    object = deg,
    contrast = c("group", "seed", "query"),
    alpha = params$pval)
  effect_size <- sapply(deg$pvalue, compute_effect_size,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("genes" = rownames(deg),
    "p_value" = deg$pvalue,
    "p_value_adj" = deg$padj,
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = deg$log2FoldChange,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, nrow(deg)),
    "query" = rep(query_id, nrow(deg)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

#' edgeR functions for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @param type either "QLF" for quasi-likelihood F-test or 
#' "LRT" fpr likelihood ratio test
#' @importFrom edgeR calcNormFactors estimateDisp
#' @importFrom edgeR glmQLFit glmQLFTest glmFit glmLRT topTags
#' @importFrom stats model.matrix
#' @importFrom dplyr filter
vesalius_deg_edger <- function(seed,
  seed_id,
  query,
  query_id,
  params,
  type) {
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
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
  # at the moment running try - some subset dont work not sure why
  # error is not informative. will skip territories if that occurs
  #--------------------------------------------------------------------------#
  deg <- edgeR::calcNormFactors(deg)
  design <- stats::model.matrix(~deg@.Data[[2]]$group)
  deg <- edgeR::estimateDisp(deg, design)
  if (type == "QLF") {
    fit <- edgeR::glmQLFit(deg, design)
    mod_fit <- edgeR::glmQLFTest(fit)
  } else {
    fit <- edgeR::glmFit(deg, design)
    mod_fit <- edgeR::glmLRT(fit)
  }
  deg <- edgeR::topTags(mod_fit, n = nrow(seed))@.Data[[1]]
  ori_ord <- match(buffer$genes, rownames(deg))
  deg <- deg[ori_ord, ]
  effect_size <- sapply(deg$PValue, compute_effect_size,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("genes" = rownames(deg),
    "p_value" = deg$PValue,
    "p_value_adj" = deg$FDR,
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = deg$logFC,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, nrow(deg)),
    "query" = rep(query_id, nrow(deg)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}

#' Logistic Regression for DEG
#' @param seed count matrix for group 1
#' @param seed_id territory used in group 1
#' @param query count matrix for group 2
#' @param query_id territory used in group 2
#' @param params parameter value list (pval, log_fc, min_pct)
#' @importFrom stats glm as.formula
#' @importFrom lmtest lrtest
#' @importFrom dplyr filter
#' @importFrom stats p.adjust
vesalius_deg_logit <- function(seed, seed_id, query, query_id, params) {
  #--------------------------------------------------------------------------#
  # Setting up data for loigt regression as suggested by the signac package 
  # I don't want to rebuild a seurat object as they do in their vignette
  # NOTE that there will be a lot of work to do here!!!!
  # It does not work completly for ATAC data. 
  #--------------------------------------------------------------------------#
  buffer <- get_deg_metrics(seed, query, params)
  if (is.null(buffer)){
    return(NULL)
  }
  merged <- format_counts_for_logit(buffer$seed, buffer$query)
  pvals <- sapply(seq_len(nrow(merged[[1L]])), function(idx, merged) {
    model_data <- cbind("gene" = merged$merged[idx, ], merged$seed_query_info)
    gene_model <- as.formula("group ~ gene")
    gene_model <- suppressWarnings(glm(formula = gene_model,
      data = model_data,
      family = "binomial"))
    null_model <- as.formula("group ~ 1")
    null_model <- suppressWarnings(glm(formula = null_model,
      data = model_data,
      family = "binomial"))
    return(lrtest(gene_model, null_model)$Pr[2])
  }, merged = merged)
  effect_size <- sapply(pvals, compute_effect_size,
    seed = ncol(buffer$seed),
    query = ncol(buffer$query))
  deg <- data.frame("genes" = buffer$genes,
    "p_value" = pvals,
    "p_value_adj" = p.adjust(pvals),
    "seed_pct" = buffer$seed_pct,
    "query_pct" = buffer$query_pct,
    "fold_change" = buffer$fc,
    "effect_size" = effect_size,
    "seed" = rep(seed_id, length(pvals)),
    "query" = rep(query_id, length(pvals)))
  deg <- filter(deg, p_value_adj <= params$pval)
  return(deg)
}


#' get DEG metrics from groups
#' @param seed count matrix for group 1
#' @param query count matrix for group 2
#' @param params parameter list (pval,log_fc, min_pct)
#' @details Computes basic metrics such as fold change 
#' and percent of cells containing each gene. This functions 
#' is used to filter out genes so we don't run differential 
#' expression on all genes. We also compute an effct size estimate 
#' by running a power analysis with 2 uneven groups 
#' importFrom pwr pwr.2p2n.test
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
    keep <- which(((seed_pct >= params$min_pct |
        query_pct >= params$min_pct) &
        abs(fc) >= params$log_fc))
    if(length(keep) == 0) {
        return(NULL)
    } else {
        seed_format <- matrix(seed[keep, ],nrow = length(keep))
        colnames(seed_format) <- colnames(seed)
        rownames(seed_format) <- rownames(seed)[keep]
        query_format <- matrix(query[keep, ],nrow = length(keep))
        colnames(query_format) <- colnames(query)
        rownames(query_format) <- rownames(query)[keep]
        return(list("genes" = rownames(seed)[keep],
        "seed" = seed_format,
        "query" = query_format,
        "seed_pct" = seed_pct[keep],
        "query_pct" = query_pct[keep],
        "fc" = fc[keep]))
    }
}

#' compute_effect_size
#' @param pval pvalue for a given gene
#' @param seed number of cells in seed
#' @param query number of cells in query
#' @details If pval is 0 we convert that to a very small number
#' pwr does not take 0 as input values.
#' DESeq might return NA pvals - just skip the effect size and 
#' return NA
#' @return efect size estimate for an unbalenced 
#' power analysis.
#' importFrom pwr pwr.2p2n.test
compute_effect_size <- function(pval, seed, query) {
  if (is.na(pval)) {
      return(NA)
  }
  pval <- max(c(pval, 1e-100))
  pval <- min(c(pval, 0.8))
  effect_size <- pwr::pwr.2p2n.test(h = NULL,
    n1 = seed,
    n2 = query,
    sig.level = pval,
    power = 0.8)$h
  return(effect_size)
}