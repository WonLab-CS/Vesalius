# check if segment selection is a valid option

check if segment selection is a valid option

## Usage

``` r
check_segment_trial(vesalius_assay, trial = "last")
```

## Arguments

- vesalius_assay:

  a vesalius_assay object

- trial:

  string - trial selection parse by user

## Value

data frame contain selected trial

## Details

check if the trial selection exists in territory slot Default is last
that will take the last entry. This function will also reformat to only
include the necessay information.
