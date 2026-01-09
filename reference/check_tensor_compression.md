# check tensor compression

check tensor compression

## Usage

``` r
check_tensor_compression(locs)
```

## Arguments

- locs:

  rle output of coordinate compression

## Value

rle out of coordinate compression

## Details

Here we just want to check if the output is reasonable or not
Essentially if you end up with very little barcodes, this could mean
that the user compressed to much. Set warning if compression is set to
less 10
