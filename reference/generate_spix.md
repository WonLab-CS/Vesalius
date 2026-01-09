# Generate super pixels from ST

Generate super pixels from ST

## Usage

``` r
generate_spix(
  vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "kmeans",
  col_resolution = 10,
  compactness = 1,
  scaling = 0.5,
  threshold = 0.9,
  index_selection = "bubble",
  max_iter = 10000,
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- dimensions:

  numeric vector of latent space dimensions to use.

- embedding:

  character string describing which embedding should be used.

- method:

  character string for which method should be used for segmentation.
  Select from "slic", "leiden_slic","louvain_slic"

- col_resolution:

  numeric colour resolution used for segmentation. (see details)

- compactness:

  numeric - factor defining super pixel compaction.

- scaling:

  numeric - scaling image ration during super pixel segmentation.

- threshold:

  numeric \[0,1\] - correlation threshold between nearest neighbors when
  generating segments from super pixels.

- index_selection:

  character - method for selecting initial indices

- max_iter:

  int - max number of kmean iterations for slic segments

- verbose:

  logical - progress message output.

## Value

a vesalius_assay object
