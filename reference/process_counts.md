# process counts

pre-process count matrices

## Usage

``` r
process_counts(
  counts,
  assay,
  method = "log_norm",
  use_count = "raw",
  nfeatures = 2000,
  min_cutoff = "q5",
  verbose = TRUE
)
```

## Arguments

- counts:

  count matrix in the form of a sparse matrix

- assay:

  character string describing the assay that is being pre-processed in
  the vesaliusObject or vesalius_assay

- method:

  character string describing which normalisation method to use. One of
  the following "log_norm", "SCT", "TFIDF", "none".

- use_count:

  string describing which counts should be used for the generating
  emebddings. Default = "raw".

- nfeatures:

  numeric describing the number of variable features to use.

- min_cutoff:

  only used when dimensionality reduction method is LSI or LSI_UMAP
  cutoff for feature to be included in the VariableFeatures for the
  object.

- verbose:

  logical - progress messages outputed or not

## Details

The \`use_count\` argument specifies which count matrix should be used
for normalization. This argument is only necessary if you use a custom
normalised count matrix. In this case, set this argument to the name you
gave your count matrix (see
[`add_counts`](wonlab-cs.github.io/Vesalius/reference/add_counts.md))
and \`generate_embeddings\` will skip the normalization and use your
custom count matrix to generate image embeddings.
