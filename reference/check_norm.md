# check if select norm method is an available option

check if select norm method is an available option

## Usage

``` r
check_norm(vesalius_assay, norm_method, method = NULL, verbose = TRUE)
```

## Arguments

- vesalius_assay:

  a vesalius_assay

- norm_method:

  string - selected normalisation parsed by user

- method:

  DEG method

- verbose:

  logical if progress message should be ouputed

## Value

count matrix

## Details

Here we check if the normalisation method has been computed in the count
slot. We also check if this is used in the context of differential gene
expression analysis. DESeq and edgeR require raw counts so if the user
parses a valid norm method but is running DEG with one those methods, we
ignore request and parse raw count instead.
