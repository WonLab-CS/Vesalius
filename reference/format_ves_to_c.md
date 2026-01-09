# convert vesalius_assay to cimg images

convert vesalius_assay to cimg images

## Usage

``` r
format_ves_to_c(vesalius_assay, dims, embed = "last", verbose = TRUE)
```

## Arguments

- vesalius_assay:

  a vesalius_Assay object

- dims:

  integer vector indicating the number of dimensions to select from
  embeddings.

- embed:

  character indicating which embedding should be selected. Default uses
  last embedding produced

- verbose:

  logical if progress message should be outputed

## Value

list of cimg images
