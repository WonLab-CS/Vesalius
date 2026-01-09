# compute the similarity between seed and query signals

compute the similarity between seed and query signals

## Usage

``` r
signal_similarity(seed, query, method = "pearson", digits = 4)
```

## Arguments

- seed:

  seed signal

- query:

  query signal

- method:

  character - correlation method to use

- digits:

  numeric - number of digits for rounding

## Value

matrix with query as rows and seed as colmuns

## Details

Chunking cost and signal into smaller chunks to run the correlation
score in paralell. There is room for improvement here. First we could
dispatch the longer list to future_lapply but cannot know which one it
is and we need to know that so we can subset the cost. Also the
functions calls feature_cost which is a R wrapper for a c++ function
(cost.cpp).
