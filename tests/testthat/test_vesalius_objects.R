# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)

test_that("Vesalius build single assay", {
    # checking normal builds with format conversion
    expect_s4_class(build_vesalius_assay(coordinates, counts),
        "vesalius_assay")
    expect_s4_class(build_vesalius_assay(coordinates,
        as.data.frame(as.matrix(counts))),
        "vesalius_assay")
    expect_s4_class(build_vesalius_assay(coordinates,
       as.matrix(counts)),
        "vesalius_assay")
    # checking normal builds with unsuported formats
    expect_error(build_vesalius_assay(coordinates, "counts"))
    expect_error(build_vesalius_assay(coordinates, 2))
    expect_error(build_vesalius_assay("coordinates", counts))
})