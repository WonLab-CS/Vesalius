# renaming counts to remain consistent with vesalius nomenclature

renaming counts to remain consistent with vesalius nomenclature

## Usage

``` r
rename_counts(integrated, seed_cells, seed_genes, query_cells, query_genes)
```

## Arguments

- integrated:

  seurat object containing integrated count data

- seed_cells:

  character vector contain seed barcodes

- seed_genes:

  character vector contain seed genes

- query_cells:

  character vector contain query barcodes

- query_genes:

  character vector contain query genes

## Value

list of count matrices with cell and gene names added and each matrix is
renamed to follow the vesalius nomencalture
