# view_gene_expression

View gene expression in spatial omics data, in specific territories or
the expression of genes in a subset of cells.

## Usage

``` r
view_gene_expression(
  vesalius_assay,
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
  return_as_list = FALSE
)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- genes:

  character vector containing the genes you wish to visualise.

- norm_method:

  character string - which count matrix should be used.

- trial:

  character string describing which segmentation trial to use. Default
  is "last" which is the last trial used.

- territory_1:

  integer or vector of integers descrbing territories in group 1 (see
  details)

- territory_2:

  integer or vector of integers descrbing territories in group 2 (see
  details)

- cells:

  charactor vector containing barcodes/spatial_indices associated with
  cell types of interest (see details)

- norm:

  logical indicating if expression should be min/max normalised

- as_layer:

  logical indicating if expression should represented as a territory/
  layer.

- with_background:

  logical - include background in plot

- cex:

  numeric - font size modulator

- cex_pt:

  numeric point size

- alpha:

  point transparency

- max_size:

  numeric - maximum point size

- return_as_list:

  logical - should plot be returned as simple list or as a ggplot object
  (single gene)/ patchwork object (multiple genes)

## Value

a ggplot object

## Details

Vesalius offers a plotting function that allows you to visualise gene
expression.

This function offers multiple options depending on what you provide.

1\. Overall expression You can visualise the overall expression pattern
of a set of genes by providing a vesalius_assay object containing
counts. If `as_layer` is set to FALSE this will show the expression at
each sptial index indivdually. If set to TRUE, this will show the
average expression on a gene in all territories present.

2\. Expression in a territory You can visualise the expression of a gene
in an isolated territory.

3\. Expression of cells in one or more territory If you want to
visualise the expression of specific cells, you can parse a character
vector containing your cells of interest. This function will
automatically subset the relevant territory data and show the expression
only in the spatial indeces hat are associated with your cell type of
interest. You can use this option to contrast the expression of cells
between territories by also providing which territories you wish to
contrast (\`territory_1\` and \`territory_2\`). If only a single
territory is provided, vesalius will only shows cell in that territory.

If you provide more than one gene, \`view_gene_expression\` will return
a ggarrange list containing all your genes as individual plots.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vesalius)
# First we build a simple object
ves <- build_vesalius_object(coordinates, counts)
# We can do a simple run
ves <- build_vesalius_embeddings(ves)

# simple smoothing
ves <- smooth_image(ves, dimensions = seq(1, 30))

# quick segmentation
ves <- segment_image(ves, dimensions = seq(1, 30))

# isolate territories
ves <- isolate_territories(ves)

# view over all expression
p <- view_gene_expression(ves, genes = "Malat1")
p1 <- view_gene_expression(ves, genes = "Malat1", as_layer = TRUE)

# view expression in isolated territory 
p2 <- view_gene_expression(ves, genes = "Malat1", territory_1 = 5)

# view expression of cells
cells <- sample(colnames(get_counts(ves)),300)
p3 <- view_gene_expression(ves,
 genes = "Malat",
 cells = cells,
 territory_1 = 5,
 terriotry_2 = 8)
} # }
```
