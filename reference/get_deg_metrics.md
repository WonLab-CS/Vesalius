# get DEG metrics from groups

get DEG metrics from groups

## Usage

``` r
get_deg_metrics(seed, query, params)
```

## Arguments

- seed:

  count matrix for group 1

- query:

  count matrix for group 2

- params:

  parameter list (pval,log_fc, min_pct)

## Details

Computes basic metrics such as fold change and percent of cells
containing each gene. This functions is used to filter out genes so we
don't run differential expression on all genes. We also compute an effct
size estimate by running a power analysis with 2 uneven groups
importFrom pwr pwr.2p2n.test
