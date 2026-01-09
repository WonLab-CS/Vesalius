# pixel normalisation dispatch function

pixel normalisation dispatch function

## Usage

``` r
norm_pixel(embeds, type = c("minmax", "quantile_norm", "z_norm"))
```

## Arguments

- embeds:

  a embedding vector

- type:

  string "minmax" or "quantileNorm"

## Details

how pixels should be normalised At the moment only miman is used.
Quantile needs to be tested.
