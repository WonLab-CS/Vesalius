# Merge coordinates from 2 vesalius assays

Merge coordinates from 2 vesalius assays

## Usage

``` r
merge_coordinates(matched, reference, barcodes)
```

## Arguments

- matched:

  matched vesalius_assay object (query object)

- reference:

  reference vesalius_assay object

- barcodes:

  barcodes to use for matching and merging

## Value

merged coordinate data frame (barcodes, x, y)

## Details

Simple merge followed a duplication check. Duplicated coordinates
produce errors when running tesselation. Just adding some noise
