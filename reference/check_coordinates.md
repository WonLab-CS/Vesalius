# check if coordinates are of the correct type and format and adjust coordinate value to remove white edge space

check if coordinates are of the correct type and format and adjust
coordinate value to remove white edge space

## Usage

``` r
check_coordinates(coordinates, assay, verbose)
```

## Arguments

- coordinates:

  coordinate data

- assay:

  string - assay name

- verbose:

  logical if progress message should be printed

## Value

coordinates data.frame or error

## Details

adjusts coordinates by either snapping coordinates to origin or min max
normalisation of coordinates. Might add polar for future tests.
