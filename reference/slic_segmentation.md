# SLIC segmentation

SLIC segmentation

## Usage

``` r
slic_segmentation(
  vesalius_assay,
  dimensions,
  col_resolution,
  embedding,
  index_selection = "bubble",
  compactness = 1,
  scaling = 0.3,
  threshold = 0.9,
  max_iter = 1000,
  verbose
)
```

## Arguments

- vesalius_assay:

  vesalius_assay object

- dimensions:

  dimensions to use

- col_resolution:

  color resolution for segmentation

- embedding:

  embedding to use

- index_selection:

  method for selecting initial indices

- compactness:

  compactness factor

- scaling:

  scaling factor

- threshold:

  correlation threshold

- max_iter:

  int - max number of kmeans iterations

- verbose:

  logical - progress messages

## Value

segmented territories
