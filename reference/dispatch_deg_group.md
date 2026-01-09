# dispatch barcodes to seed and query groups

dispatch barcodes to seed and query groups

## Usage

``` r
dispatch_deg_group(ter, seed, query, cells, sample, verbose)
```

## Arguments

- ter:

  territories data frame from territoires slot in vesalius_assay object

- seed:

  interger vector indicating which territories should be included in
  seed group

- query:

  interger vector indicating which territories should be included in
  query group

- cells:

  cell barcodes

- sample:

  barcodes for each sample type

- verbose:

  logical if progress messages should be outputed

## Value

list with seed group and seed id as well as query group and query id

## Details

This function generates groups for DEG analysis. It creates different
groups depending on what is being parse to both seed and query. This
could probably be cleaned up and simplified. Always need to return a
list though.
