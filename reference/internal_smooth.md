# internal smoothing function Function that does the smoothing on individual image array

internal smoothing function Function that does the smoothing on
individual image array

## Usage

``` r
internal_smooth(
  image,
  method = c("median", "iso", "box"),
  iter = 1,
  sigma = 1,
  box = 20,
  threshold = 0,
  neuman = TRUE,
  gaussian = TRUE,
  na.rm = FALSE,
  across_levels = "min"
)
```

## Arguments

- image:

  cimg array

- method:

  character describing smoothing method to use "median" , "iso" or "box"
  or a combination of them.

- iter:

  numeric - number of smoothing iteration

- sigma:

  numeric - standard deviation associated with isoblur (Gaussian)

- box:

  numeric describing box size (centered around center pixel) for
  smoothing

- threshold:

  numeric - discard pixels that are too low in value (cutoff threshold
  only applied in box/median blurs).

- neuman:

  logical describing If Neumann boundary conditions should be used,
  Dirichlet otherwise (default true, Neumann)

- gaussian:

  logical - use gaussian filter

- na.rm:

  logical describing if NA values should be removed

- across_levels:

  character - method used to account for multiple smoothing levels (see
  details). Select from: "min","mean", "max"
