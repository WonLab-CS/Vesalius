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

jitter_ves <- generate_embeddings(jitter_ves,
    filter_threshold = 1,
    filter_grid = 1)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 10)
jitter_ves <- equalize_image(jitter_ves, sleft = 5, sright = 5)
jitter_ves <- segment_image(jitter_ves, col_resolution = 2)
jitter_ves <- isolate_territories(jitter_ves)

gene_vec <- sample(rownames(counts), 200)


test_that("input sanity checks", {
    # signal sanity 
    expect_error(integrate_horizontally(vesalius,
        jitter_ves,
        signal = "funky",
        map = "exact"))
    # if custom genes 
    expect_error(integrate_horizontally(vesalius,
        jitter_ves,
        signal = c("Never", "Gonna", "Give", "You", "Up"),
        map = "exact"))
    # if custom matrix
    custom_matrix <- matrix(0.5, ncol = 500,
        nrow = 500)
    rownames(custom_matrix) <- sample(colnames(jitter_counts), 500)
    colnames(custom_matrix) <- sample(colnames(counts), 500)
    expect_warning(integrate_horizontally(vesalius,
        jitter_ves,
        custom_cost = custom_matrix,
        overwrite = FALSE
        ))
    rownames(custom_matrix) <- make.unique(sample(LETTERS, 500, replace = TRUE))
    colnames(custom_matrix) <- make.unique(sample(LETTERS, 500, replace = TRUE))
    expect_error(integrate_horizontally(vesalius,
        jitter_ves,
        custom_cost = custom_matrix,
        overwrite = FALSE
        ))
})
test_that("horizontal - exact", {
    # check that we can get exact match
    expect_s4_class(integrate_horizontally(vesalius,
        jitter_ves,
        signal = "variable_features",
        map = "exact"), "vesalius_assay")
    # checks?
    expect_s4_class(integrate_horizontally(vesalius,
        jitter_ves,
        signal = gene_vec,
        map = "exact"), "vesalius_assay")

})

test_that("horizontal - div", {
    # check that we can get exact match
    expect_s4_class(integrate_horizontally(vesalius,
        jitter_ves,
        signal = "variable_features",
        map = "div"), "vesalius_assay")
    # checks?
    expect_s4_class(integrate_horizontally(vesalius,
        jitter_ves,
        signal = gene_vec,
        map = "div"), "vesalius_assay")
})


