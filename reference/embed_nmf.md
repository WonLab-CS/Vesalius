# embed nmf

embed in grey scale using NMF embeddings

## Usage

``` r
embed_nmf(counts, dimensions, verbose = TRUE)
```

## Arguments

- counts:

  Seurat object containing normalised counts

- dimensions:

  number dimension to retain from NMF

- verbose:

  logical if progress messages should be outputed

## Value

normalised NMF embedding matrix
