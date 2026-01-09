# Method dispatch function for neighborhood selection

Method dispatch function for neighborhood selection

## Usage

``` r
get_neighborhood_signal(coord, signal, method, k = 20, depth = 3, radius = 20)
```

## Arguments

- coord:

  data.frame - coordinates of assay (barcodes, x, y)

- signal:

  matrix - matrix or sparse matrix containing assay signal for all
  spatial indices contained in coord

- method:

  character - which method should be use to collect neighborhood -
  switch matches

- k:

  int - how many nearest neighbors from KNN algorithm

- depth:

  int - graph path depth to define neighborhood 0 = self, 1 = direct
  neigbors, 2 = neighbors of neighbors, etc

- radius:

  \- numeric - radius around center cell

## Value

matrix of average signals for each spatial index and its neighborhood.
