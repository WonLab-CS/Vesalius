# load data and process it to get territories 
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts)
vesalius <- generate_embeddings(vesalius)
vesalius <- generate_embeddings(vesalius, normalization = "SCTransform")
vesalius <- smooth_image(vesalius, iter = 15)
vesalius <- segment_image(vesalius, col_resolution = 3)
vesalius <- isolate_territories(vesalius)
cells <- sample(get_territories(vesalius)$barcodes, 1040, replace = FALSE)
cell_error <- get_territories(vesalius)[
    get_territories(vesalius)$Territory == 1, ]$barcodes
cell_error <- sample(cell_error, 100, replace = FALSE)

test_that("Vesalius works as expected with identify markers", {
    # expect error if un suported format
    expect_error(identify_markers("vesalius"))
    # norm method does not exists
    expect_error(identify_markers(vesalius, norm_method = "funky"))
    # territory trial does not exists
    expect_error(identify_markers(vesalius, trial = "funky"))
    # if deg method does not exist 
    expect_error(identify_markers(vesalius, method = "funky"))
    # if seed territory not present
    expect_error(identify_markers(vesalius, seed = 42))
    # if query territory is not present
    expect_error(identify_markers(vesalius, query = 42))
    # If partial overlap in seed
    expect_warning(identify_markers(vesalius, seed = 4:42))
    # if partial overlap in query 
    expect_warning(identify_markers(vesalius, query = 4:42))
    # if min spatial index is too restrictive
    expect_warning(identify_markers(vesalius,
        seed = 1,
        min_spatial_index = 500))
    # expect warning if no cells present in a territory
    # Note that this will throw 2 warnings if cells are 
    # are not present in either - thestthat doesn't like that
    expect_warning(identify_markers(vesalius,
        seed = 1,
        cells = cell_error))
})

test_that("wilcox works as expected", {
    expect_s4_class(identify_markers(vesalius), "vesalius_assay")
    expect_s4_class(identify_markers(vesalius, seed = 1), "vesalius_assay")
    expect_s4_class(identify_markers(vesalius, query = 1), "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        cells = cells),
        "vesalius_assay")
    tmp <- identify_markers(vesalius,
        seed = 1,
        query = 2)
    expect_true(is(get_markers(tmp), "data.frame"))
    expect_true(is(get_markers(tmp, "DEG"), "data.frame"))
    expect_error(get_markers(tmp, "funky"))
    expect_output(show(tmp))
    # need to be expanded 
})

test_that("t.test works as expected", {
    expect_s4_class(identify_markers(vesalius,
        method = "t.test"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        method = "t.test"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        cells = cells,
        method = "t.test"),
        "vesalius_assay")
    # need to be expanded 
})

test_that("chisq works as expected", {
    # # we do not expect to have 2 by 2 table from
    # # transcriptomic count data 
    # # throw error - will need to be test with binary data
    # expect_error(identify_markers(vesalius,
    #     method = "chisq"))
    # expect_error(identify_markers(vesalius,
    #     seed = 1,
    #     query = 2,
    #     method = "chisq"))
    # expect_error(identify_markers(vesalius,
    #     seed = 1,
    #     query = 2,
    #     cells = cells,
    #     method = "chisq"))
    # # need to be expanded 
})


test_that("fisher.exact works as expected", {
    # transcriptomic data does not have binary output 
    # hense we expect an error to be thrown 
    # needs to be tested more in depeth
    # throw error - will need to be test with binary data
    # expect_error(identify_markers(vesalius,
    #     method = "fisher.exact"))
    # expect_error(identify_markers(vesalius,
    #     seed = 1,
    #     query = 2,
    #     method = "fisher.exact"))
    # expect_error(identify_markers(vesalius,
    #     seed = 1,
    #     query = 2,
    #     cells = cells,
    #     method = "fisher.exact"))
    # # need to be expanded 
})

test_that("DESeq2 works as expected", {
    expect_s4_class(identify_markers(vesalius,
        method = "DESeq2"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        method = "DESeq2"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        cells = cells,
        method = "DESeq2"),
        "vesalius_assay")
    # need to be expanded 
})


test_that("QLF works as expected", {
    # one of those territories will not work
    # not sure why but interal erro from edgR
    # message is not meaningfull
    #expect_warning(identify_markers(vesalius,
    #    method = "QLF"),
    #    "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        method = "QLF"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        cells = cells,
        method = "QLF"),
        "vesalius_assay")
    # need to be expanded 
})

test_that("LRT works as expected", {
    # one of those territories will not work
    # not sure why but interal erro from edgR
    # message is not meaningfull
    #expect_warning(identify_markers(vesalius,
    #    method = "LRT"),
    #    "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        method = "LRT"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        cells = cells,
        method = "LRT"),
        "vesalius_assay")
    # need to be expanded 
})

test_that("logit works as expected", {
    #expect_s4_class(identify_markers(vesalius,
    #    method = "logit"),
    #    "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        method = "logit"),
        "vesalius_assay")
    expect_s4_class(identify_markers(vesalius,
        seed = 1,
        query = 2,
        cells = cells,
        method = "logit"),
        "vesalius_assay")
    # need to be expanded 
})