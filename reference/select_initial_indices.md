# select inital indices

select inital indices

## Usage

``` r
select_initial_indices(
  coordinates,
  embeddings,
  type = "bubble",
  n_centers = 500,
  max_iter = 500
)
```

## Arguments

- coordinates:

  data frame containing spatial coordinates of beads

- embeddings:

  matrix containing embedding values - full pixel image

- type:

  character - method for selecting initial indices

- n_centers:

  numeric number of beads to select as super pixel centers

- max_iter:

  numeric number of iteration before returning result if no coveregnce.

## Value

barcodes of starting coordinates
