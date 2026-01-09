# Dispatch cells into batches

Dispatch cells into batches

## Usage

``` r
dispatch_batch(cost_matrix, batch_size = 5000)
```

## Arguments

- cost_matrix:

  matrix containing cost for each cell pair

- batch_size:

  int size of batch

## Value

Nested list. Each element of the list will contain a batched cost matrix
and the mapping pairs

## Details

Create cell batches that will dynamically adapt to the size of the data
set with respect to batch size. Smalled data sets, cells will be sampled
to match the size of the larger data set. This allows for multiple to
multiple matching. All cells will be selected at least once.
