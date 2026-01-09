# get cell signal from vesalius assays

get cell signal from vesalius assays

## Usage

``` r
get_signal(
  seed_assay,
  query_assay,
  signal,
  dimensions = seq(1:30),
  use_norm = "raw",
  scale = FALSE,
  verbose = TRUE
)
```

## Arguments

- seed_assay:

  vesalius_assay object

- query_assay:

  vesalius_assay object

- signal:

  character string where the signal should be taken from

- dimensions:

  int vector - if signal is embeddings which embeddings should be
  selected

- use_norm:

  charcater string which counts should be use when extracting signal

- scale:

  logical - should signal be scaled

- verbose:

  logical - should progress messages be outputed.

## Value

list contain seed signal, query signal and features used
