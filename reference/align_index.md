# assign coordinates to matched indices

assign coordinates to matched indices

## Usage

``` r
align_index(matched_index, coord, jitter = 0)
```

## Arguments

- matched_index:

  data.frame containing matching pairs of coordinates

- coord:

  data.frame containing coordinate data

- jitter:

  numeric - amount of jitter to add to duplicated coordinates

## Value

adjusted coordinate data.frame where each point receives the coordinates
of its best match.
