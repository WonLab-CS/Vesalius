# embed PCA

embed in grey scale using PCA embeddings

## Usage

``` r
embed_pca(counts, dimensions, features = NULL, verbose = TRUE)
```

## Arguments

- counts:

  Seurat object containing normalised counts

- dimensions:

  number dimension to retain from PCA

- features:

  custom vector of features

- verbose:

  logical if progress messages should be outputed

## Value

normalised PCA embedding matrix
