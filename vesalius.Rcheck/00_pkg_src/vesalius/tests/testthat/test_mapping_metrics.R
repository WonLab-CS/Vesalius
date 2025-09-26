# load vesalius data 
data(vesalius)

# Create Vesalius object for processing
vesalius <- build_vesalius_assay(coordinates, counts)
cells <- sample(LETTERS[1:6], size = nrow(vesalius@tiles),replace =T)
names(cells) <- vesalius@tiles$barcodes
vesalius <- add_cells(vesalius, cells = cells, add_name = "Cells")
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)
cells <- sample(LETTERS[1:6], size = nrow(jitter_ves@tiles),replace =T)
names(cells) <- jitter_ves@tiles$barcodes
jitter_ves <- add_cells(jitter_ves, cells = cells, add_name = "Cells")

vesalius <- generate_embeddings(vesalius,
    filter_threshold = 1,
    filter_grid = 1)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 10)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)


jitter_ves <- generate_embeddings(jitter_ves,
    filter_threshold = 1,
    filter_grid = 1)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 10)
jitter_ves <- equalize_image(jitter_ves, sleft = 5, sright = 5)
jitter_ves <- segment_image(jitter_ves, col_resolution = 2)
jitter_ves <- isolate_territories(jitter_ves)

matched <- map_assays(vesalius,
    jitter_ves,
    threshold = 0,
    use_cost = c("feature","niche","territory","composition"),
    batch_size = 500,
    epoch = 5,
    jitter = 1)


test_that("cost contrubution", {
    expect_s4_class(get_cost_contribution(matched, method = "range"),"vesalius_assay")
})


test_that("metric clusters", {
    expect_s4_class(get_metric_clusters(matched),"vesalius_assay")
})

