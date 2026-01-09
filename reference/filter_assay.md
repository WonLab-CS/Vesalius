# filter vesalius assay by cells, territories, or genes

filter vesalius assay by cells, territories, or genes

## Usage

``` r
filter_assay(
  vesalius_assay,
  cells = NULL,
  territories = NULL,
  trial = "last",
  genes = NULL
)
```

## Arguments

- vesalius_assay:

  vesalius_assay object

- cells:

  character vector of cell barcodes to keep

- territories:

  character vector of territory labels to keep

- trial:

  character string defining which territory trial to use

- genes:

  character vector of genes to keep

## Value

filtered vesalius_assay object
