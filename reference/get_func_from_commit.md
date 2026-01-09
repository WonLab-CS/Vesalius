# get function name from commit list

get function name from commit list

## Usage

``` r
get_func_from_commit(commit)
```

## Arguments

- commit:

  commit list

## Value

function name

## Details

Using tail since if you make an explicit function call using pkg::func
you get pkg as well. Neat. We only want the function
