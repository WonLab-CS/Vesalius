# filter tiles

filter tiles based on tile area and if they share an edge with the
tesselation box

## Usage

``` r
filter_tiles(tesselation, coordinates, filter_threshold)
```

## Arguments

- tesselation:

  data.frame output from the deldir function

- coordinates:

  data frame with original coordinates

- filter_threshold:

  numeric describing the quantile threshold value to use for area
  filtering

## Value

a list with 2 data frame. 1 with filtered tesselation results 2 filtered
coordinate file.

## Details

Here we want to filter based on the size of the tile under the
assumption that very large tiles are probably due to unpexted space. The
issue is that if you don't apply this threshold, subtle pixel patterns
are lost wuhtin the sea of pixels. If the ueser really dont want to
filter anything out then you don't We will set a filter threshold pretty
high though to make sure that it does get filter out as default.
