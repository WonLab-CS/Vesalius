# visualize mapping metrics for vesalius assay

visualize mapping metrics for vesalius assay

## Usage

``` r
view_mapping_metrics(
  vesalius_assay,
  trial = "cost",
  cex_pt = 1,
  cex = 15,
  randomize = TRUE
)
```

## Arguments

- vesalius_assay:

  vesalius_assay object

- trial:

  character string defining which metric trial to visualize

- cex_pt:

  numeric - point size for plotting

- cex:

  numeric - text size for plotting

- randomize:

  logical - whether to randomize point order for better visualization

## Value

ggplot object showing mapping metrics
