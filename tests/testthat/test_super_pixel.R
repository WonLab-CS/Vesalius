# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts)
vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 5)
vesalius <- segment_image(vesalius,
    method = "slic",
    dimensions = 1:3,
    col_resolution = 50,
    compactness = 1,
    index_selection = "bubble",
    scaling = 0.2)


vesalius <- segment_image(vesalius, method = "kmeans", col_resolution = 10)
image_plot(vesalius, dimensions = 1:3)

test_that("Initialise super pixel centers", {
    expect_type(vesalius:::select_initial_indices(coordinates, k = 50), "list")
})

par(mfrow =c(1,2))
plot(images)
plot(images *ratio)