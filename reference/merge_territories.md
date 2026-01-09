# Merge territories using specified labels

Merge territories using specified labels

## Usage

``` r
merge_territories(
  matched,
  reference,
  coordinates,
  labels_mapped,
  labels_reference
)
```

## Arguments

- matched:

  matched (query) vesalius_assay object

- reference:

  reference vesalius_assay

- coordinates:

  coord data frame after merging (see merge_coordinates)

- labels_mapped:

  Territory columns to merge in mapped object

- labels_reference:

  Territory columns to merge in reference object

## Value

territory data frame with merged columns

## Details

Territory data frame generally contain a lot of information depending on
how much analysis the user has run. (Segments, Cells, territories,
morphology). Using pariwise matching of labels parse in labels_mapped
and labels_reference the data will be merged into a single column. A new
column will be added to distinguish samples. All other columns that are
not present in the labels argument will retain their name and NA will be
added to cell present in the other data set (similar to left and right
joins).
