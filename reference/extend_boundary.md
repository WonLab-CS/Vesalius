# extend image boundary around territory

extend image boundary around territory

## Usage

``` r
extend_boundary(territories, morphology_factor)
```

## Arguments

- territories:

  data frame containing x/y and color value of territories to extend

- morphology_factor:

  integer or vector of integers describing growth and/or shrink extent.

## Details

we want to avoid clipping territory if they sit at the edge of the
image. To avoid this we simply extend the image boundary.
