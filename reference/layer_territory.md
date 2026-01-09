# layer_territory generates layer from the outside to the inside of a territory

layer_territory generates layer from the outside to the inside of a
territory

## Usage

``` r
layer_territory(
  vesalius_assay,
  territory = NULL,
  trial = "last",
  layer_depth = NULL,
  morphology_factor = 0,
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  vesalius_assay object

- territory:

  integer or vector of integers desrcining territories to morph.

- trial:

  character string - which territory trial that should be used to select
  territorires. Default is last one computed

- layer_depth:

  integer describing the number of final layers.

- morphology_factor:

  integer or vector of integers describing growth and/or shrink extent.

- verbose:

  logical - progress message output.

## Value

a vesalius_assay

## Details

Each territory can be subdivided into a series of layers. Each layer
will be considered a seperate territory and can be treated as such for
functions such as
[`identify_markers`](wonlab-cs.github.io/Vesalius/reference/identify_markers.md)
and
[`territory_plot`](wonlab-cs.github.io/Vesalius/reference/territory_plot.md).

However, all other territories present will be labled as "out". This
means that for the time being you can only work with a single territory
at a time.

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

# morph territory

ves <- layer_territory(ves)

# view territory morphing
territory_plot(ves)
} # }
```
