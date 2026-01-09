# Build a new vesalius_assay object using mapping information

Build a new vesalius_assay object using mapping information

## Usage

``` r
build_mapped_assay(mapped, seed_assay, query_assay, meta_labels, jitter = 0)
```

## Arguments

- mapped:

  mapping results (see output of point_mapping)

- seed_assay:

  vesalius_assay object used as reference

- query_assay:

  vesalius_assay that was mapped onto the reference

- meta_labels:

  character - name of column to be transfered to new object

- jitter:

  numeric - how much coordiate jitter should be added to the coordinates
  to avoid duplication. If 0, no jitter will be added.

## Value

vesalius_assay object with new coordinates.

## Details

This function will only return a vesalius_assay with new coordinates and
mapping inforation such as cost. The rest will remain the same. We are
only transfering coordinates here. No integration.
