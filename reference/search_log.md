# search through log for parameter values or names

search through log for parameter values or names

## Usage

``` r
search_log(vesalius_assay, arg, return_assay = TRUE)
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- arg:

  string indicating which parameter value or argument name should be
  searched for

- return_assay:

  logical indicating if the log list should returned or only the value.

## Value

either a list containing log calls or values found in call

## Details

You may search through the log to see if you have used certain
parameters and if so which values did you use when running a certain
trial. If \`return_assay\` is \`TRUE\` then the entire log call will be
returned. This will include all parameter values including defaults
parse to function.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vesalius)
# First we build a simple object
ves <- build_vesalius_object(coordinates, counts)
# We can do a simple run 
ves <- build_vesalius_embeddings(ves)
# maybe we want to try a different method 
# both will be stored in the object
ves <- build_vesalius_embeddings(ves, dim_reduction = "UMAP")

# search log 
search_log(ves, "tensor_resolution")
} # }
```
