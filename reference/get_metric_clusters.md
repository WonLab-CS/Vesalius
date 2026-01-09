# Cluster query cells based on which reference cells they tend to mapped to

Cluster query cells based on which reference cells they tend to mapped
to

## Usage

``` r
get_metric_clusters(
  vesalius_assay,
  use_cost = "feature",
  cluster_method = "hclust",
  trial = NULL,
  group_identity = NULL,
  ref_cells = NULL,
  query_cells = NULL,
  top_nn = 30,
  h = 0.75,
  k = NULL,
  nn = 30,
  resolution = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- vesalius_assay:

  vesalius_assay object after mapping a query onto a reference.

- use_cost:

  character vector describing which cost matrices should be used to
  compare cells

- cluster_method:

  character string - which method should be used for clustering (hclust,
  louvain, leiden)

- trial:

  character string defining which trial should be used for clustering if
  any. If NULL, will search for "Cells".

- group_identity:

  character vector - which specific substes of trial should be used for
  clustering By default will use all labels present.

- ref_cells:

  character vector with reference cell barcodes (by default will use all
  barcodes)

- query_cells:

  character with query cell barcodes (by default will use all barcodes)

- top_nn:

  int - how many cells should be used to define clustering similarity
  (see details)

- h:

  numeric - normalized height to use as hclust cutoff \[0,1\]

- k:

  int - number of cluster to obtain from hclust

- nn:

  int - number of nearest neighbors to use when creating graph for
  community clustering algorithms

- resolution:

  numeric - clustering resolution to be parsed to community clustering
  algorithms

- verbose:

  logical - print output message

- ...:

  additional arguments

## Value

vesalius_assay with clustering results

## Details

Once we have mapped cells between sample, we can identify which cells
tend to map to the same group of cells. To achieve this, we first create
a cost matrix that will serve as a basis to find similar-mapping
instances. The cost matrix can be constructed from any cost matrix that
was used during the mapping phase. Next, for each query cell we extract
the top_nn cells in the reference with lowest cost. Using the ordered
index as a character label, we compute a jaccard index between
overlapping labels. Query cells with a high jaccard index tend to map to
the same reference cells. We then use the reciprocal to define a
distance between cells and cluster cells based on this distance. The
same approach is used for every clustering method provided. This
function will add a new column with the metric clustering results.
