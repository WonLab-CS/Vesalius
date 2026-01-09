# Add embeddings Add custom embeddings to vesalius objects

Add embeddings Add custom embeddings to vesalius objects

## Usage

``` r
add_embeddings(vesalius_assay, embeddings, add_name = NULL, verbose = TRUE)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- embeddings:

  a matrix containing embedding values (see details)

- add_name:

  character string to be used as embedding name.

- verbose:

  logical indicating if progress message should be outputed.

## Value

a vesalius_assay object

## Details

Vesalius objects accepts custom embedding values that will be used to
generate images.

The embedding matrix should be in the form of a matrix with columns
being the latent space dimensiions and rows representing the spatial
indices present in the data set.

The intersection between spatial indices present in the tiles and custom
embeddings will be used. This intersection will be applied to the custom
embeddings and filter out all tiles that are not present in the tiles.
The tiles will not be filtered.

Rownames should also be present and should represent the barcode name.
These rownames are use to match the latent space embedding to its tile
in the image.
