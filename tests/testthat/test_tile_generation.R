# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)

test_that("Vesalius skips tile generation if already computed", {
    vesalius <- build_vesalius_assay(coordinates)
    vesalius <- generate_tiles(vesalius)
    expect_warning(generate_tiles(vesalius))
})

test_that("Vesalius reduce tensor resolution works", {
    # testing that reduced tensor works as expected 
    new_coord <- data.frame("barcodes" = paste0("barcode", seq(1, 256)),
        "xcoord" = rep(seq(1, l = 16, by = 0.5), times = 16),
        "ycoord" = rep(seq(1, l = 16, by = 0.5), each = 16))
    # building base object - no filtering
    vesalius <- build_vesalius_assay(new_coord)
    expect_s4_class(generate_tiles(vesalius,
        filter_grid = 1,
        filter_threshold = 1),
        "vesalius_assay")
    # reduce tensor resolution
    tmp <- generate_tiles(vesalius,
        tensor_resolution = 0.75,
        filter_grid = 1,
        filter_threshold = 1)
    expect_s4_class(tmp,
        "vesalius_assay")
    tiles <- get_tiles(tmp)
    # check that you only have 25 barcodes 
    expect_equal(length(unique(tiles$barcodes)), 36)
    # check that you have 25 unique origin coordinates 
    expect_equal(sum(tiles$origin), 36)
    # check that unique barcodes all have an "origin" coordinates 
    expect_equal(length(unique(tiles$barcodes)), sum(tiles$origin))
    # check that you have no negative coordinates
    expect_true(all(tiles$x > 0) && all(tiles$y > 0))
    # create a new tiles with unreasonable tensor res
    expect_warning(generate_tiles(vesalius,
        tensor_resolution = 0.5,
        filter_grid = 1,
        filter_threshold = 1))
    expect_error(generate_tiles(vesalius,
        tensor_resolution = 0.01,
        filter_grid = 1,
        filter_threshold = 1))

})

test_that("Count adjustement works as intended", {
    # testing that reduced tensor works as expected 
    new_coord <- data.frame("barcodes" = paste0("barcode", seq(1, 256)),
        "xcoord" = rep(seq(1, l = 16, by = 0.5), times = 16),
        "ycoord" = rep(seq(1, l = 16, by = 0.5), each = 16))
    # building base object - no filtering
    vesalius <- build_vesalius_assay(new_coord)
    vesalius <- generate_tiles(vesalius,
        tensor_resolution = 0.75,
        filter_grid = 1,
        filter_threshold = 1)
    tiles  <- get_tiles(vesalius)
    # test that barcodes names are compressed with et
    expect_true(length(grep("_et_", tiles$barcodes)) > 0)
    # test if all barcodes present
    uncompressed_barcodes <- unlist(strsplit(tiles$barcodes, "_et_"))
    expect_true(all(uncompressed_barcodes %in% paste0("barcode", seq(1, 256))))
})

test_that("Adjusting counts when counts added manually", {
     # generating base object
    vesalius <- build_vesalius_assay(coordinates, counts)
    # Adding tiles when tiles have been created
    vesalius <- add_counts(vesalius, counts, counts)
    vesalius <- generate_tiles(vesalius,
        tensor_resolution = 0.4)

})