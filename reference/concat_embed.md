# create new embedding from jointly measure spatial omics

create new embedding from jointly measure spatial omics

## Usage

``` r
concat_embed(seed, query, dimensions, norm_method, dim_reduction, signal)
```

## Arguments

- seed:

  vesalius_assay object of the first modality

- query:

  vesalius_assay object of the second modality

- dimensions:

  int - number of gray scale images to create

- norm_method:

  string describing which normalisation method to use. One of the
  following "log_norm", "SCT", "TFIDF", "raw"

- dim_reduction:

  string describing which dimensionality reduction method should be
  used. One of the following: "PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP"

- signal:

  character - signal type to use

## Value

embedding matrix used as grey scale image stack

## Details

Create new latent space using feature from both modalities Creates a new
feature matrix than in nromalized and converted to latent space image
stack.
