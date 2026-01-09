# convert ves embedding to image

convert ves embedding to image

## Usage

``` r
future_ves_to_cimg(i, embeddings, dims, tiles, full_image = TRUE)
```

## Arguments

- i:

  index of embedding to use

- embeddings:

  matrix - embedding matrix

- dims:

  dimensions to use

- tiles:

  tile data frame used to reconstruct images

- full_image:

  logical - should the background be returned as well

## Value

cimg object of embedding

## Details

using this as a way to run this section in parallel way to slow
otherwise. The back ground represents all pixels that are not part of
the Spatial data but constiture the "rest" of the pixels in the image.
This tends to happen when you have a non rectangular assay that needs to
be fitted into a n \* p or n \* p \* d array.
