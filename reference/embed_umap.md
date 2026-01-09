# embed umap

embed in gray scale using UMAP projections

## Usage

``` r
embed_umap(counts, dimensions, features = NULL, verbose)
```

## Arguments

- counts:

  Seurat object containing normalised counts

- dimensions:

  number of PCs to use for the UMAP projections

- features:

  custom vector of features

- verbose:

  logical if progress messages should be outputed

## Value

normalised UMAP projection matrix

## Details

Note that while you can select any number of dimensions the number of
UMAP dimensions will always be 3.
