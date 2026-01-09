# embed latent space

Embed latent space into grey color scale.

## Usage

``` r
embed_latent_space(
  counts,
  assay,
  dim_reduction,
  dimensions,
  features = NULL,
  remove_lsi_1,
  verbose
)
```

## Arguments

- counts:

  Seurat object containing counts (generally normalised)

- assay:

  charcter string of the assay being used

- dim_reduction:

  dimensionality reduction method that will be used Select from PCA,
  PCA_L, UMAP, LSI, LSI_UMAP

- dimensions:

  numeric for number of dimeniosn top retain after dimensionality
  reduction

- features:

  custom features used for dim reduction

- remove_lsi_1:

  logical if first dimension of LSI embedding should be removed (will
  soon be depreciated)

- verbose:

  logical if progress messages should be outputed or not

## Value

data frame of normalised embedding values.

## Details

General method dispatch function for dim reduction methods
