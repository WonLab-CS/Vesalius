# load vesalius data 
data(vesalius)
set.seed(1453)

# Create Vesalius object for processing
vesalius <- build_vesalius_assay(coordinates, counts)
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)

vesalius <- generate_embeddings(vesalius,
    filter_threshold = 1,
    filter_grid = 1)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 10)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)
cells <- sample(LETTERS[1:6], size = nrow(coordinates), replace = TRUE)
names(cells) <- coordinates$barcodes
vesalius <- add_cells(vesalius, cells, add_name = "Cells")

jitter_ves <- generate_embeddings(jitter_ves,
    filter_threshold = 1,
    filter_grid = 1)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 10)
jitter_ves <- equalize_image(jitter_ves, sleft = 5, sright = 5)
jitter_ves <- segment_image(jitter_ves, col_resolution = 2)
jitter_ves <- isolate_territories(jitter_ves)
cells <- sample(LETTERS[1:6], size = nrow(jitter_coord), replace = TRUE)
names(cells) <- jitter_coord$barcodes
jitter_ves <- add_cells(jitter_ves, cells, add_name = "Cells")

matched <- map_assays(vesalius,
    jitter_ves,
    batch_size = 1000,
    epochs = 10,
    threshold = -1, # noisy data will generate poor matches
    signal = "variable_features")


test_that("simple integration", {
    # check that we can get exact match
    expect_s4_class(integrate_assays(matched,
        vesalius), "vesalius_assay")
    # checks?
    expect_error(integrate_assays(matched,
        vesalius,
        labels_mapped = "Cells",
        labels_reference = c("Cells","Something_else")))

})

test_that("Duplicated Spatial Indices", {
    matched_loc <- map_assays(jitter_ves,
        vesalius,
        batch_size = 1000,
        epochs = 10,
        threshold = -1,
        signal = "variable_features")
    expect_warning(integrate_assays(matched_loc,
        jitter_ves), "vesalius_assay")
    

})

test_that("simple integration - DEGs", {
    # check that we can get exact match
    inter <- integrate_assays(matched, vesalius)
    # checks?
    expect_s4_class(identify_markers(inter, sample = TRUE),
        "vesalius_assay")

})