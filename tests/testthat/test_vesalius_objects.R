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
    expect_s4_class(build_vesalius_assay(as.matrix(coordinates),
       counts),
        "vesalius_assay")
    tmp <- coordinates
    colnames(tmp) <- LETTERS[1:3]
    expect_error(build_vesalius_assay(tmp, counts))
    # checking normal builds with unsuported formats
    expect_error(build_vesalius_assay(coordinates, "counts"))
    expect_error(build_vesalius_assay(coordinates, 2))
    expect_error(build_vesalius_assay("coordinates", counts))

    # checking that you can build partial objects
    # coordinates must always be parsed
    expect_s4_class(build_vesalius_assay(coordinates),
        "vesalius_assay")
    expect_error(build_vesalius_assay())
    expect_error(build_vesalius_assay(counts = counts))

    
    #testing output
    tmp <- build_vesalius_assay(coordinates, counts)
    expect_output(show(tmp))
})



test_that("Adding counts to vesalius_assay", {
    # checking if count matrix can be added
    vesalius <- build_vesalius_assay(coordinates)
    ## When no counts are present
    expect_s4_class(add_counts(vesalius,
        counts = counts,
        raw = counts),
        "vesalius_assay")
    expect_error(add_counts(vesalius,
        counts = counts))
    expect_s4_class(add_counts(vesalius,
        counts = counts,
        force = TRUE),
        "vesalius_assay")
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
        raw = counts),
        "vesalius_assay")
    expect_error(add_counts(vesalius,
        counts = counts))
    expect_s4_class(add_counts(vesalius,
        counts = counts,
        force = TRUE),
        "vesalius_assay")
    ## checking that active count matrix is being tagged
    tmp <- add_counts(vesalius,
        counts = counts,
        raw = counts)
    expect_identical(comment(tmp@counts), "custom_counts")
    expect_true(grep(pattern = "raw_custom_counts",
        x = names(tmp@counts),
        value = TRUE) == "raw_custom_counts")
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
    expect_s4_class(add_embeddings(vesalius, embeds),
        "vesalius_assay")
    expect_s4_class(add_embeddings(vesalius, as.data.frame(embeds)),
        "vesalius_assay")
    # adding duplicated colnames
    colnames(embeds) <- rep("A", 30)
    #expect_error(add_embeddings(vesalius, embeds))
    # We expect the same thing even if it is truncated
    expect_warning(add_embeddings(vesalius, embeds[-(1:20), ]))
    # Add sanity checks to make sure that barcodes overlap
})

test_that("Scale of coordinates",{
    expect_s4_class(build_vesalius_assay(coordinates,
        counts), "vesalius_assay")
    expect_s4_class(build_vesalius_assay(coordinates,
        counts,
        scale = 15), "vesalius_assay")
    vesalius <- build_vesalius_assay(coordinates,
        counts, scale = "auto")
    expect_equal(vesalius@meta$scale$scale, 23.0623695)
    vesalius <- build_vesalius_assay(coordinates,
        counts, scale = 15)
    expect_equal(vesalius@meta$scale$scale, 15)
    expect_equal(vesalius@meta$unit$unit, "um")
})

test_that("Adding Cells",{
    # generating base object
    vesalius <- build_vesalius_assay(coordinates, counts)
    cells <- sample(LETTERS[1:10], size = nrow(coordinates), replace = TRUE)
    # expect error if no names assigned to cell type
    expect_error(add_cells(vesalius, cells))
    names(cells) <- make.unique(sample(LETTERS, nrow(coordinates), replace = TRUE))
    expect_error(add_cells(vesalius, cells))
    # adding new territory slot
    names(cells) <- coordinates$barcodes
    expect_s4_class(add_cells(vesalius, cells), "vesalius_assay")
    cells <- sample(LETTERS[1:13], size = nrow(coordinates), replace = TRUE)
    names(cells) <- coordinates$barcodes
    expect_s4_class(add_cells(vesalius, cells, add_name = "funky"),
        "vesalius_assay")

})

