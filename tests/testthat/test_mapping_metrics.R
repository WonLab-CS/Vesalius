# load vesalius data 
data(vesalius)

# Create Vesalius object for processing
vesalius <- build_vesalius_assay(coordinates, counts)
cells <- sample(LETTERS[1:6], size = nrow(vesalius@tiles),replace =T)
names(cells) <- vesalius@tiles$barcodes
vesalius <- add_cells(vesalius, cells = cells, add_name = "Cells")
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)
cells <- sample(LETTERS[1:6], size = nrow(jitter_ves@tiles),replace =T)
names(cells) <- jitter_ves@tiles$barcodes
jitter_ves <- add_cells(jitter_ves, cells = cells, add_name = "Cells")

vesalius <- generate_embeddings(vesalius,
    filter_threshold = 1,
    filter_grid = 1)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 10)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)


jitter_ves <- generate_embeddings(jitter_ves,
    filter_threshold = 1,
    filter_grid = 1)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 10)
jitter_ves <- equalize_image(jitter_ves, sleft = 5, sright = 5)
jitter_ves <- segment_image(jitter_ves, col_resolution = 2)
jitter_ves <- isolate_territories(jitter_ves)

matched <- map_assays(vesalius,
    jitter_ves,
    threshold = 0,
    use_cost = c("feature","niche","territory","composition"),
    batch_size = 500,
    epoch = 5,
    jitter = 1)




test_that("metric clusters", {
    expect_s4_class(get_metric_clusters(matched),"vesalius_assay")
})


test_that("comparing niches", {
    # this should return some data
    # expect_s3_class(vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "niche",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     pval = 0.5,
    #     log_fc = 0.0),"data.frame")
    # # this should return an empty data frame 
    # local <- vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "niche",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     pval = 0.0001,
    #     log_fc = 10)
    # expect_true(is.null(local))
})


test_that("comparing niches aggregate", {
    # this should return some data
    # expect_s3_class(vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "niche",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     aggregate = TRUE,
    #     pval = 0.5,
    #     log_fc = 0.0),"data.frame")
    # # this should return an empty data frame 
    # local <- vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "niche",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     aggregate = TRUE,
    #     pval = 0.0001,
    #     log_fc = 10)
    # expect_true(is.null(local))
})



test_that("comparing territories", {
    # this should return some data
    # expect_s3_class(vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "territory",
    #     group_1 = ter_ves,
    #     group_2 = ter_jitter,
    #     pval = 0.5,
    #     log_fc = 0.0),"data.frame")
    # # this should return an empty data frame 
    # local <- vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "territory",
    #     group_1 = ter_ves,
    #     group_2 = ter_jitter,
    #     pval = 0.0001,
    #     log_fc = 10)
    
    # expect_true(is.null(local))
})


test_that("comparing territories aggregate", {
    # this should return some data
    # expect_s3_class(vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "territory",
    #     group_1 = ter_ves,
    #     group_2 = ter_jitter,
    #     aggregate = TRUE,
    #     pval = 0.5,
    #     log_fc = 0.0),"data.frame")
    # # this should return an empty data frame 
    # local <- vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "territory",
    #     group_1 = ter_ves,
    #     group_2 = ter_jitter,
    #     aggregate = TRUE,
    #     pval = 0.0001,
    #     log_fc = 10)
    
    # expect_true(is.null(local))
})



test_that("comparing niche composition", {
    # this should return some data
    # expect_s3_class(vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "composition",
    #     method = "chisq",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     pval = 1,
    #     log_fc = 0.0),"data.frame")
    # # this should return an empty data frame 
    # local <- vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "composition",
    #     method = "chisq",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     pval = 0.0001,
    #     log_fc = 10)
    
    # expect_true(is.null(local))
})


test_that("comparing niche composition aggregate", {
    # this should return some data
    # expect_s3_class(vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "composition",
    #     method = "chisq",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     aggregate = TRUE,
    #     pval = 0.5,
    #     log_fc = 0.0),"data.frame")
    # # this should return an empty data frame 
    # local <- vesalius::compare_assays(vesalius,
    #     matched,
    #     compare = "composition",
    #     nethod = "chisq",
    #     group_1 = sub_ves,
    #     group_2 = sub_jitter,
    #     aggregate = TRUE,
    #     pval = 0.0001,
    #     log_fc = 10)
    
    # expect_true(is.null(local))
})