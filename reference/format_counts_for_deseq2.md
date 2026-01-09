# format counts for DESeq

format counts for DESeq

## Usage

``` r
format_counts_for_deseq2(seed, query)
```

## Arguments

- seed:

  seed count matrix

- query:

  query count matrix

## Value

DESeq2 object

## Details

Always adding a pseudocount of 1 to avoid issues with 0 counts. We also
force coercion to int since DESeq does not handle numerics nor does it
do internal coersion.
