# Integrate counts using Seurat

Integrate counts using Seurat

## Usage

``` r
integrate_counts(
  matched,
  reference,
  method = "HarmonyIntegration",
  nfeatures = 2000,
  dimensions = 30,
  infer = FALSE,
  signal = "variable_features",
  verbose
)
```

## Arguments

- matched:

  matrix - matrix containing counts from matched/mapped assay

- reference:

  matrix - matrix containing counts from reference assay

- method:

  character - Seurat integration method to use

- nfeatures:

  integer - number of features to use during integration

- dimensions:

  interger - number of dimensions integrated latent space dimensions.

- infer:

  logical - back infer original counts by reversing reduced dimensional
  space roations.

- signal:

  character - defining which signal should be returned:
  variable_features, all_features or custom gene list.

- verbose:

  logical - print output messages
