# equalise image histogram

equalizeHistogram image enhancement via colour histogram equalization.

## Usage

``` r
equalize_image(
  vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  method = "BalanceSimplest",
  N = 1,
  smax = 1,
  sleft = 1,
  sright = 1,
  lambda = 0.1,
  up = 100,
  down = 10,
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

  character - histogram EQ type. Select from: BalanceSimplest,
  EqualizePiecewise, SPE, EqualizeDP, EqualizeADP, ECDF (see details)

- N:

  numeric describing how each colour channel will be mapped back to the
  image (Higher N = Higher greyscale contrast). Used with
  EqualizePiecewise

- smax:

  numeric - upper limit if contrast stretching. Used with
  EqualizePiecewise

- sleft:

  numeric - Range 0 - 100. Percentage of pixel to be saturated on the
  left side of the histogram. Used with BalanceSimplest

- sright:

  numeric - Range 0 - 100. Percentage of pixel to be saturated on the
  right side of the histogram. Used with BalanceSimplest

- lambda:

  numeric - strength of background correction. Used with SPE (Screened
  Poisson Equation).

- up:

  numeric - color value threshold in the upper limit. Used with
  EqualizeDP.

- down:

  numeric color value threshold in the lower limit. Used with
  EqualizeDP.

- verbose:

  logical - progress message output.

## Value

a vesalius_assay object

## Details

Histogram equalization ensures that image details are amplified. In
turn, territories may be extract with greater precision. We recommend
balancing the histogram prior to smoothing.

For further details on each method described here, please refer to
[https://cran.r-project.org/web/packages/imagerExtra/vignettes/gettingstarted.html](wonlab-cs.github.io/Vesalius/reference/imagerExtra%20Vignette)

## Examples

``` r
if (FALSE) { # \dontrun{
data(vesalius)
# First we build a simple object
ves <- build_vesalius_object(coordinates, counts)
# We can do a simple run
ves <- build_vesalius_embeddings(ves)

# simple EQ
ves <- equalisz_image(ves, embedding = "PCA")
} # }
```
