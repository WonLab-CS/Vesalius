# kmeans segmentation function

kmeans segmentation function

## Usage

``` r
kmeans_segmentation(
  vesalius_assay,
  dimensions = seq(1, 3),
  col_resolution = 10,
  embedding = "last",
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  a vesalius_assay

- dimensions:

  numeric vector of latent space dimensions to use.

- col_resolution:

  integer or vector of positive integers. Colour depth used for
  segmentation.

- embedding:

  character string describing which embedding should be used.

- verbose:

  logical - progress message output.

## Value

list containing segmented image as an active embedding and territory
cluster for all barcodes.

## Details

Run an interaive kmeans segmentation with the possibility to run
multiple rounds of smoothing.
