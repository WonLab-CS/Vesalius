# edgeR functions for DEG

edgeR functions for DEG

## Usage

``` r
vesalius_deg_edger(seed, seed_id, query, query_id, params, type)
```

## Arguments

- seed:

  count matrix for group 1

- seed_id:

  territory used in group 1

- query:

  count matrix for group 2

- query_id:

  territory used in group 2

- params:

  parameter value list (pval, log_fc, min_pct)

- type:

  either "QLF" for quasi-likelihood F-test or "LRT" fpr likelihood ratio
  test
