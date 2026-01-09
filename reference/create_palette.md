# create color palette from predefine scheme

create color palette from predefine scheme

## Usage

``` r
create_palette(territories, randomise)
```

## Arguments

- territories:

  vesalius territories taken from a vesalius_assay

- randomise:

  logical describing if colour palette should be randomised.

## Value

color vector

## Details

We use a predefined palette that use colour blind friendly base colours.
We generate a color palette based on the number of territories present.
If required the colours will be randomly assinged to each territory.
Note that as the territory plot return a ggplot object, you can easily
override the color scheme.
