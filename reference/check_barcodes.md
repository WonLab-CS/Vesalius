# checking overlap between barcodes in counts and coordinates

checking overlap between barcodes in counts and coordinates

## Usage

``` r
check_barcodes(mat_barcodes, coordinates, throw = TRUE)
```

## Arguments

- mat_barcodes:

  character vector containing barcode names in matrix (count matrix or
  embedding matrix)

- coordinates:

  character vector containing barcode names in tile data frame

- throw:

  logical if warning should be thrown or not

## Value

overlapping location between counts and coordinates

## Details

Will throw in a warning if the overlap is not perfect.
