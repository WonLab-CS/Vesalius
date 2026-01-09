# adjust count

adjust counts after reducing the resolution of the image tensor or after
filtering stray beads

## Usage

``` r
adjust_counts(coordinates, counts, throw = TRUE, verbose = TRUE)
```

## Arguments

- coordinates:

  data frame containing coordinates after reducing resolution and
  compressing cooridnates

- counts:

  count matrix

- throw:

  logical - throwing warning message for unshared barcodes

- verbose:

  logical if progress messaged should be outputed.

## Value

a count matrix with adjusted count values

## Details

This function will check the coordinate file to see if any barcodes have
been merged together. If so, the counts will be adjusted by taking the
average count value accross all barcodes that have been merged together.
