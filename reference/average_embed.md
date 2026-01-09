# average image stack between seed and query

average image stack between seed and query

## Usage

``` r
average_embed(seed, query, dimensions)
```

## Arguments

- seed:

  matrix - seed embedding image stack

- query:

  matrix - query embedding image stack

- dimensions:

  int vector describing which embeddings should be selected

## Value

embedding matrix containing average pixel value for both seed query

## Details

Takes select embedding from seed and query and creates and avarage the
grey scale pixel values for each spatial location
