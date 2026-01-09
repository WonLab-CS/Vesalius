# Graph depth method based method to select niche

Graph depth method based method to select niche

## Usage

``` r
graph_neighborhood(coord, depth)
```

## Arguments

- coord:

  data.frame - coordinates of spatial indices in assay

- depth:

  int - graph path depth to define neighborhood 0 = self, 1 = direct
  neigbors, 2 = neighbors of neighbors, etc

## Value

list containing barcodes of neighbors for each spatial index.
