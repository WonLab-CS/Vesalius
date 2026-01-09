# rasterise tiles

fill tiles with pixel - rasterisation

## Usage

``` r
rasterise(filtered)
```

## Arguments

- filtered:

  data.frame with voronoi tile coordinates

## Value

a data frame barcodes and their associated pixels.

## Details

Here, we take our tile cooridnates and fill them with pixels.
Essentially, each voronoi tile can be discretised into a series of
pixels and we achieve this by reconstructing a polygon from the
tesselation coordinates and finding all discrete value that fall within
this polygon.

Note that the polygon coordinates need to be "convexified". Essentially,
the order of the coordinates maters here so we order the coordinates
before recontructing the polygon (see
[`convexify`](wonlab-cs.github.io/Vesalius/reference/convexify.md))
