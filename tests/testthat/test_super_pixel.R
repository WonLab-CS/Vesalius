# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts)
vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, sigma = 5, iter = 10)


test_that("Initialise super pixel centers", {
    expect_type(vesalius:::select_initial_indices(coordinates, k = 50), "list")
})

