# SCTransform

SCTransform pre-processing from Seurat

## Usage

``` r
int_sctransform(counts, nfeatures)
```

## Arguments

- counts:

  seurat object containing counts

- nfeatures:

  number of top variable features to select

## Value

list with seurat object used later and normalised counts to be stored in
a vesalius object
