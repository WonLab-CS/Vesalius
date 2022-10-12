# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)

test_that("Vesalius build single assay", {
    # checking normal builds with format conversion
    expect_s4_class(build_vesalius_object(coordinates, counts))
    expect_s4_class(build_vesalius_object(coordinates,
        as.data.frame(as.matrix(counts))))
    expect_s4_class(build_vesalius_object(coordinates,
       as.matrix(counts)))
    # checking normal builds with unsuported formats
    expect_error(build_vesalius_object(coordinates, "counts"))
    expect_error(build_vesalius_object(coordinates, 2))
    expect_error(build_vesalius_object("coordinates", counts))
    # checking normal builds with too many assays provided
    expect_error(build_vesalius_object(coordinates, counts, c("one", "two")))
})

test_that("Vesalius build multiple assay", {
    # first we build lsit of assays
    coord <- list(coordinates, coordinates)
    co <- list(counts, counts)
    # checking build with multiple assays
    expect_s4_class(build_vesalius_object(coord, co))
    expect_s4_class(build_vesalius_object(coord,
        co,
        c("my_assay", "my_assay1")))
    # checking unsupported formats over multiple assays
    coord <- list(coordinates, "coordinates")
    expect_error(build_vesalius_object(coord, co))
    co <- list("counts", counts)
    coord <- list(coordinates, coordinates)
    expect_error(build_vesalius_object(coord, co))
    # checking normal builds with too many assays provided
    expect_error(build_vesalius_object(coordinates,
        counts,
        c("one", "two", "three")))
})
