# assign centroid values to active embedding

assign centroid values to active embedding

## Usage

``` r
assign_centers(vesalius_assay, clusters, kcenters, dimensions, ratio = NULL)
```

## Arguments

- vesalius_assay:

  a vesalius assy object

- clusters:

  data.frame containing cluster values

- kcenters:

  matrix containing centroid values for each dimension

- dimensions:

  vector (nummeric / int) describin which latent space

- ratio:

  if used in the context of super pixel - spatial ration dimensiuons
  shouls be used.

## Value

matrix for the active embedding usiong color segementation
