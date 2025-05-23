# load vesalius data 
data(vesalius)

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

# noise_coord <- data.frame("barcodes" = paste0("bar_", 1:1200),
#     "x" = runif(1200, min = 1, max =  1000),
#     "y" = runif(1200, min = 1, max =  1000))
# noise_counts <- round(runif(1200 * nrow(counts), min = 0, max = 100))
# dim(noise_counts) <- c(nrow(counts), 1200)
# colnames(noise_counts) <- noise_coord$barcodes
# rownames(noise_counts) <- rownames(counts)
# noise_ves <- build_vesalius_assay(noise_coord, noise_counts)
# noise_ves <- generate_embeddings(noise_ves)

# gene_vec <- sample(rownames(counts), 200)


test_that("input sanity checks", {
    # signal sanity 
    expect_error(map_assays(vesalius,
        jitter_ves,
        signal = "funky",
        map = "exact"))
    # if custom genes 
    expect_error(map_assays(vesalius,
        jitter_ves,
        signal = c("Never", "Gonna", "Give", "You", "Up"),
        map = "exact"))
    # if custom matrix
    custom_matrix <- matrix(0.5, ncol = 500,
        nrow = 500)
    rownames(custom_matrix) <- sample(colnames(jitter_counts), 500)
    colnames(custom_matrix) <- sample(colnames(counts), 500)
    expect_warning(map_assays(vesalius,
        jitter_ves,
        custom_cost = custom_matrix
        ))
    rownames(custom_matrix) <- make.unique(sample(LETTERS, 500, replace = TRUE))
    colnames(custom_matrix) <- make.unique(sample(LETTERS, 500, replace = TRUE))
    expect_error(map_assays(vesalius,
        jitter_ves,
        custom_cost = custom_matrix
        ))
    custom_matrix <- matrix(0.5, ncol = ncol(counts),
        nrow = ncol(jitter_counts))
    colnames(custom_matrix) <- colnames(counts)
    rownames(custom_matrix) <- colnames(jitter_counts)
    custom_matrix <- list(custom_matrix)
    expect_warning(map_assays(vesalius,
        jitter_ves,
        custom_cost = custom_matrix
        ))
    names(custom_matrix) <- "custom_cost"
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        custom_cost = custom_matrix
        ),"vesalius_assay")
})

test_that("horizontal - neighborhood", {
    # check that we can get exact match
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        signal = "variable_features",
        neighborhood = "knn",
        k = 10), "vesalius_assay")
    # checks?
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        signal = "variable_features",
        neighborhood = "radius",
        radius = 20), "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        signal = "variable_features",
        neighborhood = "graph",
        depth = 2), "vesalius_assay")

})


test_that("batching case", {
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 600,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1001,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1052,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1053,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        batch_size = 1000,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 100,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        batch_size = 500,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1030,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        batch_size = 1030,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 5001,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        batch_size = 2000,
        signal = "variable_features"),
        "vesalius_assay")
})



test_that("batching case with multi epoch", {
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        threshold = 0,
        batch_size = 1000,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        batch_size = 600,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        batch_size = 1001,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        batch_size = 1051,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        batch_size = 1053,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        epochs = 3,
        batch_size = 1000,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 100,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        epochs = 3,
        batch_size = 500,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        batch_size = 1030,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        epochs = 3,
        batch_size = 1030,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        epochs = 3,
        batch_size = 5001,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        epochs = 3,
        batch_size = 2000,
        signal = "variable_features"),
        "vesalius_assay")
    expect_s4_class(map_assays(jitter_ves,
        vesalius,
        epochs = 50,
        batch_size = 300,
        signal = "variable_features"),
        "vesalius_assay")
})

test_that("filtering", {
   expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 2000,
        threshold = 0.6,
        signal = "variable_features"),
        "vesalius_assay")
    tmp <- map_assays(vesalius,
        jitter_ves,
        batch_size = 2000,
        threshold = 0.6,
        signal = "variable_features")
    #expect_true(nrow(get_coordinates(tmp)) == 840)
    expect_error(map_assays(jitter_ves,
        vesalius,
        batch_size = 2000,
        threshold = 1,
        signal = "variable_features"))

})

# Split this into sections for each type and composite selections


test_that("feature score", {
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = "feature"),
        "vesalius_assay")
})

test_that("niche score", {
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = "niche"),
        "vesalius_assay")
})

test_that("territory score", {
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = "territory"),
        "vesalius_assay")
})


test_that("niche composition", {

    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = "composition"),
        "vesalius_assay")
})


test_that("cell type labels", {

    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = "cell_type"),
        "vesalius_assay")
    # checking if wrong cell labels parsed throws warning

})

test_that("multi_cost_types", {
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 5000,
        threshold = 0.1,
        use_cost = c("feature", "niche", "territory","cell_type")),
        "vesalius_assay")

    custom_matrix <- matrix(0.5, ncol = ncol(counts),
        nrow = ncol(jitter_counts))
    rownames(custom_matrix) <- colnames(jitter_counts)
    colnames(custom_matrix) <- colnames(counts)
    expect_warning(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = c("feature", "niche"),
        custom_cost = custom_matrix))
    custom_matrix <- list(custom_matrix)
    names(custom_matrix) <- "funky"
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = c("feature", "niche", "funky"),
        custom_cost = custom_matrix),
        "vesalius_assay")
    expect_error(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = c("feature", "niche","Not_so_funky"),
        custom_cost = custom_matrix))
    custom_matrix <- c(custom_matrix, custom_matrix, custom_matrix)
    names(custom_matrix) <- LETTERS[1:3]
    expect_s4_class(map_assays(vesalius,
        jitter_ves,
        batch_size = 1000,
        threshold = 0.1,
        use_cost = c("feature", "niche", "A", "B", "C"),
        custom_cost = custom_matrix),
        "vesalius_assay")
})


