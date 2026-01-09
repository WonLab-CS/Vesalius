# territory_plot - plotting Vesalius territories

territory_plot - plotting Vesalius territories

## Usage

``` r
territory_plot(
  vesalius_assay,
  trial = "last",
  split = FALSE,
  highlight = NULL,
  contour = "None",
  randomise = TRUE,
  cex = 10,
  cex_pt = 1,
  alpha = 0.65,
  use_image = FALSE
)
```

## Arguments

- vesalius_assay:

  a vesalius_Assay object

- trial:

  character string describing which segmentation trial to use. Default
  is "last" which is the last segmentation trial used.

- split:

  logical - If TRUE, territories will be plotted in separate panels

- highlight:

  numeric vector describing which territories should be highlighted.

- contour:

  if territory contours should be added. Availble: "None", "convex",
  "concave"

- randomise:

  logical - If TRUE, colour palette will be randomised.

- cex:

  numeric describing font size multiplier.

- cex_pt:

  numeric describing point size multiplier.

- alpha:

  opacity factor \]0,1\[

- use_image:

  logical - whether to use image background

## Value

a ggplot object

## Details

Territory plots show all territories in false colour after they have
been isolated from a Vesalius image.

Note that this function can be applied to image segments, territories,
and layers.

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

# Plot Territories
p <- territory_plot(ves)
} # }
```
