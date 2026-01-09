# distance pooling beads of colour segment into seperate territories

distance pooling beads of colour segment into seperate territories

## Usage

``` r
distance_pooling(img, capture_radius, min_spatial_index)
```

## Arguments

- img:

  data frame contain all barcodes of a single sgement

- capture_radius:

  numeric proportion of max distance between beads to use as distance
  threshold between beads

- min_spatial_index:

  numeric minimum number of beads that should be contained in a
  territory.

## Details

Beads that are too far away or bead cluster that are below the minimum
number of spatial indices will all be pooled under the isolated label.
Note that this label is used across color segments.
