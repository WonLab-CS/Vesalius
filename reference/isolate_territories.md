# isolating territories from vesalius image segments

isolating territories from vesalius image segments

## Usage

``` r
isolate_territories(
  vesalius_assay,
  method = "distance",
  trial = "last",
  capture_radius = 0.05,
  global = TRUE,
  min_spatial_index = 10,
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  vesalius_Assay object

- method:

  character describing barcode pooling method. Currently, only
  "distance" availble

- trial:

  character string describing which segmentation trial to use. Default
  is "last" which is the last segmentation trial used.

- capture_radius:

  numeric - proportion of maximum distance between barcodes that will be
  used to pool barcodes together (range 0 - 1).

- global:

  logical - If TRUE, territories will be numbered across all colour
  segments. If FALSE, territories will be numbered within each colour
  segment.

- min_spatial_index:

  integer - minimum number of barcodes/spots/beads required in each
  territory

- verbose:

  logical - progress message output.

## Value

a vesalius_assay object

## Details

Image segments can be further subdivided into 2D seperated territorires.
This is accomplished by pooling barcodes that are associated with a
colour cluster into territories based on the distance between each
barcode.

First, `isolate_territories` considers the maximum distance between all
beads. The `capture_radius` will define which proportion of this
distance should be considered.

Second, a seed barcode will be selected and all barcodes that are within
the capture distance of the seed barcode with be pooled together. This
process is then applied on barcodes that have been selected in this
manner. The process is repeated until all barcodes have been pooled into
a territory. If there are still barcodes remaining, a new seed barcode
is selected and the whole process is repeated. NOTE : Territory
isolation is applied to each colour segment independently.

If a territory does not contain enough barcodes, it will be pooled into
the isolated territory. This territory contains all isolated territories
regardless of colour cluster of origin.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vesalius)
# First we build a simple object
ves <- build_vesalius_object(coordinates, counts)
# We can do a simple run
ves <- build_vesalius_embeddings(ves)

# simple smoothing
ves <- smooth_image(ves, dimensions = seq(1, 30))

# quick segmentation
ves <- segment_image(ves, dimensions = seq(1, 30))

# isolate territories
ves <- isolate_territories(ves)
} # }
```
