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

    ## Check barcode matching 

})

test_that("Adding counts to vesalius_assay", {
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
    expect_identical(comment(tmp@counts), "custom_counts")
    tmp <- add_counts(vesalius,
        counts = counts,
        raw = counts,
        count_type = "my_counts")
    expect_identical(comment(tmp@counts), "my_counts")
    ## checking if count have been properly addjusted
    ## Note: this is check again when tiles have already 
    ## been generated.
    barcodes_in_tiles <- vesalius@tiles$barcodes
    barcodes_in_counts <- colnames(tmp@counts$my_counts)
    expect_identical(barcodes_in_tiles, barcodes_in_counts)
})

test_that("Adding counts to vesalius_assay with counts", {
    # generating base object
    vesalius <- build_vesalius_assay(coordinates, counts)
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
    expect_identical(comment(tmp@counts), "custom_counts")
    expect_true(grepl(pattern = "raw_custom_counts",
        x = names(tmp@counts)))
})

test_that("Adding custom embeddings", {
    # from coordinates only 
    vesalius <- build_vesalius_assay(coordinates)
    vesalius <- generate_tiles(vesalius)
    # Generate custom embedding matrix of identical size
    embeds <- matrix(0,
        ncol = 30,
        nrow = length(unique(vesalius@tiles$barcodes)))
    # Expect error if no rownames are provided 
    expect_error(add_embeddings(vesalius, embeds))
    # Adding rownames
    rownames(embeds) <- unique(vesalius@tiles$barcodes)
    # expect valid object now that you have rownames
    expect_s4_class(add_embeddings(vesalius, emebds))
    # We expect the same thing even if it is truncated
    expect_s4_class(add_embeddings(vesalius, emebds[-(1:20), ]))
    # Add sanity checks to make sure that barcodes overlap
})