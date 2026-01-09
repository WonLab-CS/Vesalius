# optimize matching scores through batching

optimize matching scores through batching

## Usage

``` r
optimize_matching(cost_matrix, batch_size = 10000, epochs = 1, verbose = TRUE)
```

## Arguments

- cost_matrix:

  matrix containing mapping cost for each cell

- batch_size:

  int - number of cells to be assigned to each batch

- epochs:

  number of epochs to run the optimization

- verbose:

  logical - output progress messages

## Value

list with best matching cell pairs (data.frame) and total cost at each
epoch
