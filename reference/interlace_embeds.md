# interlace image stack between seed and query

interlace image stack between seed and query

## Usage

``` r
interlace_embeds(seed, query, dimensions)
```

## Arguments

- seed:

  matrix - seed embedding image stack

- query:

  matrix - query embedding image stack

- dimensions:

  int vector describing which embeddings should be selected

## Value

embedding matrix containing seed embeddings + query embeddings.

## Details

Takes selected embedding from seed and query and creates and interlaced
embedding matrix starting with the first seed embedding
