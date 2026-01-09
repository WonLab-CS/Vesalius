# dispatch barcodes to subset cost for match clustering

dispatch barcodes to subset cost for match clustering

## Usage

``` r
dispatch_cost_groups(
  vesalius_assay,
  cost,
  trial = NULL,
  group_identity = NULL,
  ref_cells = NULL,
  query_cells = NULL
)
```

## Arguments

- vesalius_assay:

  vesalius_assay object post cell mapping

- cost:

  cost matrix

- trial:

  character - trial name

- group_identity:

  character - name of column containing group to be used for cluster
  (i.e. territories, segments, layers etc)

- ref_cells:

  character vector - reference cell barcodes

- query_cells:

  character vector - query cell barcodes

## Value

character string of barcodes
