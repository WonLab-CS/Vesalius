# no norm

no normalisation applied simply return raw counts

## Usage

``` r
no_norm(counts, use_count = "raw")
```

## Arguments

- counts:

  seurat object containing counts

- use_count:

  string describing name that needs to be added to list element. This
  list will be appended to the count slot in the vesalius_assay.

## Value

list with seurat object used later and raw counts to be stored in the
vesalius objects

## Details

Here, either the user doesn't want to normalise the data or they provide
their custom count matrix. In this case, we parse it as "none" to avoid
writing another function and add the custom name.
