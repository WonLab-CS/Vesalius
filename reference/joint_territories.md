# Jointly measured spatial omic assays territories

Jointly measured spatial omic assays territories

## Usage

``` r
joint_territories(
  mod1,
  mod2,
  dimensions = seq(1, 30),
  embedding = "last",
  method = "interlace",
  norm_method = "log_norm",
  dim_reduction = "PCA",
  signal = "variable_features",
  verbose = TRUE
)
```

## Arguments

- mod1:

  vesalius_assay object containing first modality

- mod2:

  vesalius_assay objecty containing second modality

- dimensions:

  numeric vector describing latent space dimensions to use during
  intergration

- embedding:

  character - embedding to use

- method:

  character - integration method. interlace - mean - concat are
  available options

- norm_method:

  character - which count values should be use for integration when
  using concat method

- dim_reduction:

  characater - which dim reduction methods should be used for concat
  integration (PCA,PCA_L,UMAP,LSI,LSI_UMAP,NMF)

- signal:

  character - signal type to use

- verbose:

  logical - should progress message be outputed to the console.

## Value

vesalius object containing new image embeddings
