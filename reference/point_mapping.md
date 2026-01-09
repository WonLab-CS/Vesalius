# mapping points between data sets

mapping points between data sets

## Usage

``` r
point_mapping(
  query_signal,
  query_assay,
  cost,
  seed_signal,
  seed_assay,
  method = "pearson",
  neighborhood = "knn",
  k = 20,
  radius = 0.05,
  depth = 1,
  batch_size = 10000,
  epochs = 1,
  use_cost = c("feature", "niche"),
  threshold = 0.5,
  filter_cells = FALSE,
  seed_territory_labels = "Territory",
  query_territory_labels = "Territory",
  seed_meta_labels = NULL,
  query_meta_labels = NULL,
  digits = 4,
  verbose = TRUE
)
```

## Arguments

- query_signal:

  processed query signal from query assay

- query_assay:

  vesalius_assay object

- cost:

  matrix - matrix of size n (query cells) by p (seed cells) containing
  custom cost matrix.

- seed_signal:

  processed seed signal from seed assay

- seed_assay:

  vesalius_assay object

- method:

  character - correlation method

- neighborhood:

  character - neighborhood method

- k:

  int size of niche (knn)

- radius:

  0.05 proportion of max distance to use as radius for neighborhood

- depth:

  graph path depth to condsider for neighborhood.

- batch_size:

  int number of points in each query batch

- epochs:

  numeric - number of epochs

- use_cost:

  character string defining how should total cost be computer Available:
  feature, niche, territory, composition (See details for combinations

- threshold:

  score threshold below which indicices should be removed. Scores will
  always be between 0 and 1

- filter_cells:

  logical - filter cells

- seed_territory_labels:

  character - seed territory labels

- query_territory_labels:

  character - query territory labels

- seed_meta_labels:

  character - seed meta labels

- query_meta_labels:

  character - query meta labels

- digits:

  numeric - number of digits

- verbose:

  logical - out progress to console

## Value

list of matched and aligned coordinates in query
