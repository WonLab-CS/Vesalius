# retrieve the points contained in the edge of each territory

retrieve the points contained in the edge of each territory

## Usage

``` r
unpack_territory_path(trial, tiles, method = "none")
```

## Arguments

- trial:

  name of territory trial that should be selected

- tiles:

  vesalius tiles

- method:

  string - method for how territories edges should be selected

## Value

a data frame containing edge of each territory.

## Details

Here we are using convex as start point. Essentially, we order the
coordinates based on their polar coordinates using the median coordinate
as the center point.
