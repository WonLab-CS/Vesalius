# Add counts to vesalius assay Adding custom count matrix to a vesalius assay object.

Add counts to vesalius assay Adding custom count matrix to a vesalius
assay object.

## Usage

``` r
add_counts(
  vesalius_assay,
  counts,
  raw_counts = NULL,
  add_name = NULL,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- vesalius_assay:

  a vesalius assay object

- counts:

  matrix or sparse matrix containing normalised counts

- raw_counts:

  matrix or sparse matrix containing raw counts

- add_name:

  character string defining the name of the count matrix being added.

- force:

  logical indicating if count matrix provided should also be used as raw
  count matrix.

- verbose:

  logical indicating if progress messages should be outputed.

## Value

a vesalius_assay object

## Details

In some case, you might wish to use your own normalisation method. In
this case, you can add your own matrix and specify the name you want to
give to that count matrix.

Differential gene expression tools such as DESeq2 and edgeR require raw
counts. As such, we recommend providing raw counts as well. If you do
not have raw counts, or do not wish to provide them, you can set the
force argument to TRUE. This will force vesalius to generate a copy of
your count matrix and use this as "raw" count matrix.

Since coordinates need to be parsed to the vesalius_assay contructor,
this function will also compared the your count matrix to the coordinate
data. It will trim the count matrix based on barcodes shared between
both.
