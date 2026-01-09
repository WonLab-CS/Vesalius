# Method dispatch function for neighborhood selection - added flavor specific to composition

Method dispatch function for neighborhood selection - added flavor
specific to composition

## Usage

``` r
niche_composition(
  coord,
  vesalius_assay,
  method,
  cell_label = NULL,
  k = 20,
  depth = 3,
  radius = 20
)
```

## Arguments

- coord:

  data.frame - coordinates of assay (barcodes, x, y)

- vesalius_assay:

  vesalius_assay object

- method:

  character - which method should be use to collect neighborhood -
  switch matches

- cell_label:

  character - cell label column name

- k:

  int - how many nearest neighbors from KNN algorithm

- depth:

  int - graph path depth to define neighborhood 0 = self, 1 = direct
  neigbors, 2 = neighbors of neighbors, etc

- radius:

  \- numeric - radius around center cell

## Value

list of all cells and cell types for each niche
