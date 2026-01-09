# filter grid

filtered stray barcodes/spatial indices.

## Usage

``` r
filter_grid(coordinates, filter_grid)
```

## Arguments

- coordinates:

  data frame containing barcodes, / y coordinates of each barcode.

- filter_grid:

  numeric describing size of the grid to use as proporiton of array
  size.

## Value

a data.frame containing barcodes, x and y coordinates.

## Details

we create a grid over the data an check which grid section contain a
number of barcodes lower than a quantile threshold. Those barcodes will
be removed. As it stands, I am not satisfied with this function. I think
it is too restrictive and will most likely not be applicable to highly
standardised assay.
