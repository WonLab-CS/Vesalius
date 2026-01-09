# embed PCA loading values

embed in grey scale using PCA Loading value

## Usage

``` r
embed_pcal(counts, dimensions, features = NULL, verbose = TRUE)
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

normalised PCA loading matrix

## Details

This approach is a slightly different as it takes the loading value
associted to each gene in a barcode and sums the absolute value of each
of those values. Once all genes in all barcodes have been summed, we
normalise the latent space and return the matrix.
