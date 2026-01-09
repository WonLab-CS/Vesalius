# image_plot - plotting vesalius embeddings

image_plot - plotting vesalius embeddings

## Usage

``` r
image_plot(
  vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  cex = 10
)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- dimensions:

  which dimensions to use to generate image (see details)

- embedding:

  character string descrining which embedding should be used for image
  generation.

- cex:

  numeric - font and point size resizing factor.

## Value

ggplot object

## Details

Once you have generated your embeddings, you can visualise these
embeddings using the `image_plot` function. This function will generate
a ggplot object representing the embedding image containing all pixels.
You can select any dimension and in any combination you desire, however
you can only select 1 or 3 dimensions to visualise at a time. This will
either generate grey scale image or RGB images.

This function is always applied to the active embedding. By default,
this is the last you generated. This means that you can also use this
function to visualise the effect if smoothing, equalization,
regularisation or segmentation.

Please note that if you are using "louvain" or "leiden" for image
segmentation, this will always generate grey scale images even if you
select multiple dimensions. "Kmeans" on the other hand will still
produce RGB color segments.

The usage of this function remains the same in after any processing
steps.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vesalius)
# First we build a simple object
ves <- build_vesalius_object(coordinates, counts)
# We can do a simple run 
ves <- build_vesalius_embeddings(ves)
# Plot 1st 3 PCs
p <- image_plot(ves)

} # }
```
