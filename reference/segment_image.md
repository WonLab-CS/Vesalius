# segment image

segment vesalius images to find initial territories

## Usage

``` r
segment_image(
  vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "kmeans",
  col_resolution = 10,
  compactness = 1,
  scaling = 0.5,
  threshold = 0.9,
  index_selection = "bubble",
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- dimensions:

  numeric vector of latent space dimensions to use.

- embedding:

  character string describing which embedding should be used.

- method:

  character string for which method should be used for segmentation.
  Select from "kmeans", "louvain", "leiden", "slic",
  "leiden_slic","louvain_slic","som"

- col_resolution:

  numeric colour resolution used for segmentation. (see details)

- compactness:

  numeric - factor defining super pixel compaction.

- scaling:

  numeric - scaling image ration during super pixel segmentation.

- threshold:

  numeric \[0,1\] - correlation threshold between nearest neighbors when
  generating segments from super pixels.

- index_selection:

  character - method for selecting initial indices

- verbose:

  logical - progress message output.

## Value

a vesalius_assay object

## Details

Applying image segmentation ensures a reduction in colour complexity.

Vesalius provides 7 different methods for clustering colours and
reducing color complexity: \*\*Kmeans\*\*, \*\*Louvain\*\*,
\*\*Leiden\*\*, \*\*slic\*\*, \*\*leiden_slic\*\*, \*\*louvain_slic\*\*,
and \*\*som\*\*

In the case of kmeans clustering the `col_resolution` argument shows the
number of colours that the images should be reduced to. In this case,
`col_resolution` should be an integer and we suggest first looking at
values between 3 and 20.

In the case of \*\*leiden\*\* and \*\*louvain\*\* clustering, the
`col_resolution` is the granularity of the clustering. In this case, we
suggest using values between 0.01 and 1 to start with. We recommned
uisng \*\*louvain\*\* clustering over \*\*leiden\*\* in this context.

In the case of slic, the col_resolution define the number of starting
points used to generate super pixels. Depending on the number of points
there are in the assay, we suggested using 10 number of points as
starting point. For example, if you have 1000 spatial indices, you can
set col_resolution to 100.

The optimal `col_resolution` will depend on your interest and biological
question at hand. You might be interested in more or less granular
territories. Along with smoothing, the number of segments is one way to
control this granularity.

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
} # }
```
