# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)

test_that("Checking and converting count type - Vesalius build", {
    expect_s4_class(build_vesalius_object(coordinates, counts))
    expect_s4_class(build_vesalius_object(coordinates,
        as.data.frame(as.matrix(counts))))
    expect_s4_class(build_vesalius_object(coordinates,
       as.matrix(counts)))
    expect_error(build_vesalius_object(coordinates, "counts"))
    expect_error(build_vesalius_object(coordinates, 2))
})

test_that("Checking and converting coordinate type - Vesalius build", {
    expect_s4_class(build_vesalius_object(coordinates, counts))
    expect_s4_class(build_vesalius_object(as.matrix(coordinates), counts))
    expect_s4_class(build_vesalius_object(coordinates, counts,
        adjust_coordinates = "norm"))
    expect_error(build_vesalius_object("coordinates", counts))
    expect_error(build_vesalius_object(20, counts))
    expect_error(build_vesalius_object(coordinates, counts,
        adjust_coordinates = "something"))
    colnames(coordinates) <- c("tatty", "bow", "jangles")
    expect_error(build_vesalius_object(coordinates, counts))
})