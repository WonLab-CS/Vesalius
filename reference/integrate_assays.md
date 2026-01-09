# integrate 2 vesalius assays

integrate 2 vesalius assays

## Usage

``` r
integrate_assays(
  mapped,
  reference,
  method = "CCAIntegration",
  nfeatures = 2000,
  signal = "variable_features",
  dimensions = 30,
  infer = TRUE,
  use_counts = "raw",
  labels_mapped = NULL,
  labels_reference = NULL,
  regenerate_tiles = TRUE,
  tensor_resolution = 1,
  filter_grid = 1,
  filter_threshold = 1,
  verbose = TRUE
)
```

## Arguments

- mapped:

  vesalius_assay object - mapped veslius assay (map_assays)

- reference:

  vesalius_assay object - reference vesalius_assay (seed)

- method:

  character - count integration method (methods provided by Seurat v5)

- nfeatures:

  integer - number of variable features to sue during integration.

- signal:

  character - defining which signal should be returned:
  variable_features, all_features or custom gene list.

- dimensions:

  interger - number of dimensions integrated latent space dimensions.

- infer:

  logical - back infer original counts by reversing reduced dimensional
  space roations.

- use_counts:

  character - which count matrix to use during integration

- labels_mapped:

  character - which columns in the mapped data assay should be merged
  with the reference data (see details)

- labels_reference:

  character - which columns in the reference data assay should be merged
  with the mapped data (see details)

- regenerate_tiles:

  logical - should tiles be regenrated from integrated coordinates

- tensor_resolution:

  numeric - tensor resolution for tile generation

- filter_grid:

  numeric - filter grid parameter

- filter_threshold:

  numeric - filter threshold parameter

- verbose:

  logical - should progressed message be printed

## Value

a vesalius_assay object

## Details

After mapping coordinates from a query onto a reference, vesalius
provides a way to then integrate the assays together. This function
will:

\* Integrate Counts using Seurat \* Merge coordinates (adding a jitter
to avoid overlapping coordinayes) \* Merge territories (pair-wise
merging using labels_mapped and labels_reference - everything else will
have a separate column).

We also infer log nomalized counts from CCA latent space.

The final output of this function is a vesalius_assay object containing
coordinates from the mapped and reference, integrated latent space
(e.g.CCA) integrated counts, merged territories, an additional meta
data.

It should be noted that this object does not contain tiles. To use this
vesalius_assay as any other, add tiles by using
[`generate_tiles`](wonlab-cs.github.io/Vesalius/reference/generate_tiles.md)
