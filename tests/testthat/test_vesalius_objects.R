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

    # checking that you can build partial objects
    # coordinates must always be parsed 
    expect_s4_class(build_vesalius_assay(coordinates))
    expect_error(build_vesalius_assay())
    expect_error(build_vesalius_assay(counts = counts))

})

test_that("Adding counts to empty vesalius_assay",{
    # checking if count matrix can be added
    vesalius <- build_vesalius_assay(coordinates)
    ## When no counts are present
    expect_s4_class(add_counts(vesalius,
        counts = counts,
        raw = counts))
    expect_error(add_counts(vesalius,
        counts = counts))
    expect_s4_class(add_counts(vesalius,
        counts = counts,
        force = TRUE))
    ## checking that active count matrix is being tagged
    tmp <- add_counts(vesalius,
        counts = counts,
        raw = counts)
    expect_identical()
    
   
})

 ## When counts are already present 
    vesalius <- build_vesalius_assay(coordinates, counts)
    expect_s4_class(add_counts(vesalius, counts))
    ## Checking default count matrix is being parsed to comment
    vesalius <- add_counts(vesalius, counts)