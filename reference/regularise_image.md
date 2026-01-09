# regularise image

regularise_image denoise Vesalius images via variance regularization

## Usage

``` r
regularise_image(
  vesalius_assay,
  dimensions = seq(1, 3),
  embedding = "last",
  lambda = 1,
  niter = 100,
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

- lambda:

  numeric - positive real numbers describing regularization parameter
  (see details)

- niter:

  numeric - number of variance regularization iterations (Default = 100)

- verbose:

  logical - progress message output.

## Value

a vesalius_assay

## Details

Image regularization can be seen as a form of image denoising. Details
on each method can be found in the tvR package under the denoise2
function
[https://cran.r-project.org/web/packages/tvR/tvR.pdf](wonlab-cs.github.io/Vesalius/reference/tvR).

A higher value for lambda will results in a smoother image. It should be
noted that in the context of spatial omics the more sparse the points in
the data (the more space between coordinates), the more you will need to
increase the value of lambda to obtain better denoising.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vesalius)
# First we build a simple object
ves <- build_vesalius_object(coordinates, counts)
# We can do a simple run
ves <- build_vesalius_embeddings(ves)

# simple regularisation
ves <- regularise_image(ves, embedding = "PCA")
} # }
```
