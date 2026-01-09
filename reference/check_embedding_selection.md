# check embedding selection

check embedding selection

## Usage

``` r
check_embedding_selection(vesalius_assay, embed, dims)
```

## Arguments

- vesalius_assay:

  a vesalius_assay

- embed:

  string embedding selection choice

- dims:

  integer vector containing embedding dimension to extract

## Value

embedding data frame

## Details

we want to check if the embedding that the user requests is present in
the assay. If not return error. If more than one with that name return
last entry of that name and warning. Default is last that will just take
the last embedding created which should be stored in the active slot in
the vesalius_Assay object. We also want to be able to select which
dimensions we want to return. We make sure that those dimensions can be
extracted from the embedding data frame.
