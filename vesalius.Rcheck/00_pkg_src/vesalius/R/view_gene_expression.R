
#' view_gene_expression
#' 
#' View gene expression in spatial omics data, in specific territories or
#' the expression of genes in a subset of cells.
#' @param vesalius_assay a vesalius_assay object
#' @param genes character vector containing the genes you wish to 
#' visualise.
#' @param norm_method character string - which count matrix should be used.
#' @param trial character string describing which segmentation trial
#' to use. Default is "last" which is the last trial used.
#' @param territory_1 integer or vector of integers descrbing territories in
#' group 1 (see details)
#' @param territory_2 integer or vector of integers descrbing territories in
#' group 2 (see details)
#' @param cells charactor vector containing barcodes/spatial_indices
#' associated with cell types of interest (see details)
#' @param norm logical indicating if expression should be min/max normalised
#' @param as_layer logical indicating if expression should represented as
#' a territory/ layer.
#' @param cex numeric - font size modulator
#' @param cex_pt numeric point size
#' @param alpha point transparency
#' @param return_as_list logical - should plot be returned as simple list
#' or as a ggplot object (single gene)/ patchwork object (multiple genes)
#' @details Vesalius offers a plotting function that allows you to 
#' visualise gene expression. 
#' 
#' This function offers multiple options depending on what you provide.
#' 
#' 1. Overall expression
#' You can visualise the overall expression pattern of a set of genes
#' by providing a vesalius_assay object containing counts.
#' If \code{as_layer} is set to FALSE this will show the expression at
#' each sptial index indivdually. If set to TRUE, this will show the
#' average expression on a gene in all territories present.
#' 
#' 2. Expression in a territory
#' You can visualise the expression of a gene in an isolated territory.
#' 
#' 3. Expression of cells in one or more territory
#' If you want to visualise the expression of specific cells, you can
#' parse a character vector containing your cells of interest. This
#' function will automatically subset the relevant territory data and
#' show the expression only in the spatial indeces hat are associated
#' with your cell type of interest. You can use this option to contrast
#' the expression of cells between territories by also providing which
#' territories you wish to contrast (`territory_1` and `territory_2`).
#' If only a single territory is provided, vesalius will only shows cell
#' in that territory.
#' 
#' If you provide more than one gene, `view_gene_expression` will return 
#' a ggarrange list containing all your genes as individual plots. 
#' @return a ggplot object
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
#' # view over all expression
#' p <- view_gene_expression(ves, genes = "Malat1")
#' p1 <- view_gene_expression(ves, genes = "Malat1", as_layer = TRUE)
#' 
#' # view expression in isolated territory 
#' p2 <- view_gene_expression(ves, genes = "Malat1", territory_1 = 5)
#' 
#' # view expression of cells
#' cells <- sample(colnames(get_counts(ves)),300)
#' p3 <- view_gene_expression(ves,
#'  genes = "Malat",
#'  cells = cells,
#'  territory_1 = 5,
#'  terriotry_2 = 8)
#'}
#' @export
#' @importFrom infix %||%
#' @importFrom dplyr left_join filter group_by mutate
#' @importFrom stats na.exclude
#' @importFrom ggplot2 ggplot geom_point aes facet_wrap theme_classic
#' @importFrom ggplot2 scale_color_gradientn theme element_text
#' @importFrom ggplot2 guides guide_legend labs margin element_rect
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
#' @importFrom patchwork wrap_plots
view_gene_expression <- function(vesalius_assay,
  genes = NULL,
  norm_method = "last",
  trial = "last",
  territory_1 = NULL,
  territory_2 = NULL,
  cells = NULL,
  norm = TRUE,
  as_layer = FALSE,
  with_background = FALSE,
  cex = 10,
  cex_pt = 1,
  alpha = 0.75,
  max_size = 5,
  return_as_list = FALSE) {
    #--------------------------------------------------------------------------#
    # First lets get the norm method out and the associated counts
    #--------------------------------------------------------------------------#
    counts <- check_norm(vesalius_assay, norm_method)
    territories <- check_territory_trial(vesalius_assay, trial) %>%
      dispatch_territory(territory_1, territory_2, cells)
    #--------------------------------------------------------------------------#
    # Getting genes
    # We will facet wrap if there is more than one.
    #--------------------------------------------------------------------------#
    genes <- genes %||%
      stop("Please specifiy which gene(s) you would like to visualize")
    #----------------------------------------------------------------------#
    # First lets make a list to store our counts and then plots
    #----------------------------------------------------------------------#
    gene_list <- vector("list", length(genes))
    names(gene_list) <- genes
    #----------------------------------------------------------------------#
    # Next we loop,extract, combine
    #----------------------------------------------------------------------#
    for (i in seq_along(genes)) {
      gene <- rownames(counts) == genes[i]
      if (sum(gene) == 0) {
        warning(paste(genes[i], "is not present in count matrix.
          Returning NULL"), immediate. = TRUE)
        gene_list[[i]] <- NULL
        next()
      }
      gene <- counts[gene, ]
      gene <- data.frame(names(gene), gene)
      colnames(gene) <- c("barcodes", "gene")
      gene <- left_join(gene, territories, by = c("barcodes")) %>%
        na.exclude()
      gene_list[[i]] <- gene
    }

    #----------------------------------------------------------------------#
    # if nothing to plot
    #----------------------------------------------------------------------#
    nulls <- sapply(gene_list, is.null)
    if (all(nulls)) {
        return(NULL)
    }
    
    #--------------------------------------------------------------------------#
    # Now we can loop of gene list
    #--------------------------------------------------------------------------#
    for (i in seq_along(gene_list)) {
      #------------------------------------------------------------------------#
      # Extract data and convert count to mean count if as_layer=T
      #------------------------------------------------------------------------#
      if (is.null(gene_list[[i]])) next()

      other <- filter(gene_list[[i]], trial == "other") %>% droplevels()
      territory <- filter(gene_list[[i]], trial != "other") %>% droplevels()
      if (as_layer) {
          territory <- territory %>%
            group_by(trial) %>%
            mutate(gene = mean(gene)) %>%
            ungroup()
      }
      #------------------------------------------------------------------------#
      # checking if layer and normalising post layering if layering required 
      #------------------------------------------------------------------------#
      if (norm) {
        territory$gene <- min_max(territory$gene)
        type <- "Norm. Expression"
      } else {
        type <- "Expression"
      }
      #------------------------------------------------------------------------#
      # Create background
      #------------------------------------------------------------------------#
      if (nrow(other) > 0) {
        ge <- ggplot() +
          geom_point(data = other, aes(x, y),
            alpha = alpha,
            fill = "#2e2e2e",
            size = cex_pt,
            show.legend = TRUE)
      } else {
        ge <- ggplot()
      }
      #------------------------------------------------------------------------#
      # Create ggplot with foreground
      #------------------------------------------------------------------------#
      if (is.null(cells)) {
        ge <- ge + geom_point(data = territory,
          aes(x = x, y = y, col = gene), size = cex_pt * 2, alpha = alpha)
      } else {
        ge <- ge + geom_point(data = territory,
          aes(x = x, y = y, shape = trial, col = gene),
          size = cex_pt * 2,
          alpha = alpha)
      }

      ge <- ge +
        scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
        theme_classic() +
        theme(legend.title = element_text(size = cex),
          legend.text = element_text(size = cex),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
        labs(col = type, title = names(gene_list)[i]) +
        guides(shape = guide_legend(override.aes = list(size = cex * 0.5)))
      
      if (with_background){
        ge <- ge +
          theme(panel.background = element_rect(
            fill = rev(brewer.pal(11, "Spectral"))[1]))
      }
      gene_list[[i]] <- ge
    }
    #--------------------------------------------------------------------------#
    # returning ggplot or ggplot list with ggarrange
    #--------------------------------------------------------------------------#
    if (is.null(gene_list)) {
      return(NULL)
    } else {
      gene_list <- gene_list[!sapply(gene_list, is.null)]
    }

    if (length(gene_list) == 1) {
      gene_list <- gene_list[[1L]]
    } else if (!return_as_list) {
      gene_list <- ggarrange(plotlist = gene_list,
        ncol = min(c(floor(sqrt(length(gene_list))), max_size)),
        nrow = min(c(floor(sqrt(length(gene_list))), max_size)))
    }
    return(gene_list)
}
