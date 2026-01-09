# Compute score based on cell type labels

Compute score based on cell type labels

## Usage

``` r
cell_type_match(seed_labels, query_labels)
```

## Arguments

- seed_labels:

  cell type labels in seed (reference) data

- query_labels:

  cell type labels in query data

## Value

Score matrix based on cell label similarity

## Details

Return 1 if same label and 0 if different.
