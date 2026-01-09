# Internal differantial gene expression function.

This function dispatches the groups to the various methods That are
avialbale

## Usage

``` r
vesalius_deg(
  seed,
  query,
  seed_id,
  query_id,
  method,
  log_fc,
  pval,
  min_pct,
  min_spatial_index,
  verbose = TRUE,
  args
)
```

## Arguments

- seed:

  = group1 count data

- query:

  = group 2 count data

- seed_id:

  = territory ID's for group 1

- query_id:

  = territory ID's for group 2

- method:

  = DEG stat method

- log_fc:

  = fold change threshold

- pval:

  = p value threshold

- min_pct:

  = minimum percentage of barcodes that should contain a given gene

- min_spatial_index:

  = minimum number of barcodes present in a territory

- verbose:

  = progress message output

- args:

  arguments parse to (...) in upper level function (not functional)
