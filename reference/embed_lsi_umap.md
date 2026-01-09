# embed lsi

embed in grey scale using latent semantic indexing followed by UMAP

## Usage

``` r
embed_lsi_umap(
  counts,
  dimensions,
  features = NULL,
  remove_lsi_1,
  verbose = TRUE
)
```

## Arguments

- counts:

  Seurat object containing normalised counts

- dimensions:

  numeric for number of latent space dimensions to use

- features:

  custom vector of features

- remove_lsi_1:

  logical if first LSI dimenions should be removed

- verbose:

  logical if progress messages should be outputed

## Value

normalised LSI embedding matrix
