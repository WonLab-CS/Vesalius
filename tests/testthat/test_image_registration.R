# load vesalius data 
# We will assume that the embeddings have been produced 

data(vesalius)
image_file <- file.path(system.file("extdata", package = "vesalius"),
    "dentate.png")

test_that("vesalius assay with image", {
    img <- load.image(image_file)
    expect_s4_class(build_vesalius_assay(coordinates, counts, img),
        "vesalius_assay")
    img <- image_file
    expect_s4_class(build_vesalius_assay(coordinates, counts, img),
        "vesalius_assay")
    img <- matrix(0, nrow = 100, ncol = 100)
    expect_s4_class(build_vesalius_assay(coordinates, counts, img),
        "vesalius_assay")
    img <- array(0, dim = c(100, 100, 3))
    expect_s4_class(build_vesalius_assay(coordinates, counts, img),
        "vesalius_assay")
    img <- array(0, dim = c(100, 100, 5, 4, 3))
    expect_error(build_vesalius_assay(coordinates, counts, img))
    expect_error(build_vesalius_assay(coordinates, counts, "funky"))
})


vesalius <- build_vesalius_assay(coordinates, counts, image_file)
vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, sigma = 5, iter = 15)
vesalius <- segment_image(vesalius, col_resolution = 3)
vesalius <- isolate_territories(vesalius)



test_that("generate a image template", {
    expect_true(is(vesalius:::generate_image_template(vesalius), "cimg"))
    expect_true(is(vesalius:::generate_image_template(vesalius, "Segment"),
        "cimg"))
    expect_true(is(vesalius:::generate_image_template(vesalius, "PCA"),
        "cimg"))
    expect_true(is(vesalius:::generate_image_template(vesalius, "active"),
        "cimg"))
})

test_that("Register source to image", {
    test <- register_image(vesalius)
    img <- load.image(image_file)
    
    img_v <- vesalius:::generate_image_template(vesalius, "active")
    
})

