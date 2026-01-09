# compute_effect_size

compute_effect_size

## Usage

``` r
compute_effect_size(pval, seed, query)
```

## Arguments

- pval:

  pvalue for a given gene

- seed:

  number of cells in seed

- query:

  number of cells in query

## Value

efect size estimate for an unbalenced power analysis. importFrom pwr
pwr.2p2n.test

## Details

If pval is 0 we convert that to a very small number pwr does not take 0
as input values. DESeq might return NA pvals - just skip the effect size
and return NA
