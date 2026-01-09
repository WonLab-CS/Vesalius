# reduce tensor resoltuon

reduce the size of the image tensor by merging barcodes together after
compressing coordinates values

## Usage

``` r
reduce_tensor_resolution(coordinates, tensor_resolution = 1)
```

## Arguments

- coordinates:

  data frame with barcode coordinates

- tensor_resolution:

  numeric (range 0 - 1) describing the compression ratio to be applied
  to the final image. Default = 1

## Value

a data frame with barcodes, x and coordinates

## Details

While we compress the coordinates, we retain all barcodes. Each barcode
that overlap with each other are marged together. Their respective
counts will also be merged together. This allows us to retain all
barcodes for downstream analysis.
