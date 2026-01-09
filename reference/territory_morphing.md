# territory_morphing applies morphological operators to a set of territoriees

territory_morphing applies morphological operators to a set of
territoriees

## Usage

``` r
territory_morphing(
  vesalius_assay,
  territory = NULL,
  trial = "last",
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

- morphology_factor:

  integer or vector of integers describing growth and/or shrink extent.

- verbose:

  logical - progress message output.

## Value

a vesalius_assay

## Details

Territory morphing can manipulate territories by growing, shrinking,
filling, and cleaning territories. Growing = Positive integers -
Territory will be dilated by x number of pixels Shrinking = Negative
integers - Territory will be contracted by x number of pixels Filling =
grow followed by shrink. Cleaning = shrink followed by grow.

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

ves <- territory_morphing(ves, 8, morphology_factor = 30)
ves <- terriotry_morphing(ves, 1, morpholgy_factor = c(-15, 15))

# view territory morphing
territory_plot(ves)
} # }
```
