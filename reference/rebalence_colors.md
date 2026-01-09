# rebalence colors

rebalence colors

## Usage

``` r
rebalence_colors(coordinates, dimensions, method = "minmax")
```

## Arguments

- coordinates:

  data frame containing barcodes, x/y coord, origin, and color value.

- dimensions:

  number dimensions select for plotting

- method:

  character string: min max or truncate

## Value

data frame with \[0,1\] bound values

## Details

This function is use to re-bound values between 0 and 1. Some image
processing steps may lead to negative values being introduced. Imager
handles this without an issue but for plotting we want to make sure that
all values are indeed bound between 0 and 1. We will try to different
methods either we truncate the value or we min max norm.
