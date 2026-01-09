# build vesalius assay object

build a simple vesalius assay object from single count matrix and
spatial coordinate pair.

## Usage

``` r
build_vesalius_assay(
  coordinates,
  counts = NULL,
  image = NULL,
  assay = "spatial_omics",
  scale = "auto",
  unit = "um",
  verbose = TRUE
)
```

## Arguments

- coordinates:

  data.frame in the format barcodes, x, y. Default is NULL. See details.

- counts:

  matrix, sparse matrix containing counts. Default is NULL. See details.

- image:

  connection string or image array

- assay:

  character vector containing names of the assays (see details).

- scale:

  character \| numeric - if "auto", vesalius will compute 99 percentile
  of inter barcodes distance else provide a numeric value describing
  distance between barcodes.

- unit:

  character - units of scale

- verbose:

  logical indicating if progress message should be outputed or not.

## Value

A vesalius_assay objecy

## Details

The vesalius_assay constructor allows you to create partial or complete
vesalius_assay objects.

Partial objects contain only the coordinates.

Complete objects contain both the counts and the coordinates.

The main purpose of using partial objects (or empty objects) is for you
to be able to provide your own count matrix. This will be useful if you
want to normalise your data in a way that is not provided by vesalius.

This approach of using your own data can also apply to embeddings. If
you have generated a set of latent space embeddings that you wish to use
instead of those provided by vesalius, you can use the
[`add_embeddings`](wonlab-cs.github.io/Vesalius/reference/add_embeddings.md)
function.

Along side this input data, you can provide a name to your assay. If
none are provided, Vesalius will generate a set of names based on the
default assay name "spatial_omics".

You can decide if you want to adjust the coordinates to the origin i.e
remove unnecessary space or normalise the coordinates. Norm is not
recommened at the moment

## Examples

``` r
data(vesalius)
# Single assay object
ves <- build_vesalius_assay(coordinates, counts)
#> #--------------------------------------------------------------------------------# 
#> 2026-17-01/09/26 22:17:40  Checking Coordinates in spatial_omics 
#> 2026-17-01/09/26 22:17:40  Checking Counts in spatial_omics 
#> 2026-17-01/09/26 22:17:40  Calculating Assay scale from coordinates
#> #--------------------------------------------------------------------------------# 
```
