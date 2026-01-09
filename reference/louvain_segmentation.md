# louvain segmentation

using leiden clustering to cluster colors

## Usage

``` r
louvain_segmentation(
  vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 0.01,
  embedding = "last",
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- dimensions:

  embedding dimensions used for clustering

- col_resolution:

  clustering resolution used for leiden

- embedding:

  embedding type used for clustering

- verbose:

  logical if progress message should outputed

## Value

list with updated segmented embedding values and segment territories.
