# convexify

order coordinates based on their angle around a central point

## Usage

``` r
convexify(xside, yside, indx, indy, order = FALSE)
```

## Arguments

- xside:

  vector of x coordinates

- yside:

  vector of y coordinates

- indx:

  central point x coordinate

- indy:

  central point y coordinate

- order:

  logical - if TRUE return order only and not the coordinates

## Value

a data frame of ordered x and y coordinates.

## Details

For rasterisation, the shape of the polygon must be as convex as
possible. To ensure that all points are in some sense in a convex form
we order them based on their polar coordinate angle.
