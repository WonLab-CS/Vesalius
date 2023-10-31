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
        tensor_resolution = 0.000001,
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
    expect_s4_class(add_counts(vesalius, counts, counts),
        "vesalius_assay")
    vesalius <- add_counts(vesalius, counts, counts)
    expect_s4_class(generate_tiles(vesalius,
        tensor_resolution = 0.7),
        "vesalius_assay")

})
# Loading data from the packages
data(vesalius)

test_that("Unknown embedding type request", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    expect_error(generate_embeddings(vesalius, dim_reduction = "funky"))
})

test_that("Unknown norm type request", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    expect_error(generate_embeddings(vesalius, normalisation = "funky"))
})

test_that("Vesalius PCA embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "PCA")
    expect_identical(get_embeddings(vesalius, active = FALSE)$PCA,
        get_embeddings(vesalius, active = TRUE))
    # it should work as well if run it more than once
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    expect_true(get_active_embedding_tag(vesalius) == "PCA.1")
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("PCA", "PCA.1")))
})

test_that("Vesalius NMF embeds", {
  vesalius <- build_vesalius_assay(coordinates, counts)
  vesalius <- generate_embeddings(vesalius,
                                  dim_reduction = "NMF",
                                  normalisation = "log_norm")
  # simple embedding build and tile build
  expect_s4_class(vesalius, "vesalius_assay")
  # check if active embedding is set correctly
  expect_true(get_active_embedding_tag(vesalius) == "NMF")
  expect_identical(get_embeddings(vesalius, active = FALSE)$NMF,
                   get_embeddings(vesalius, active = TRUE))
  # it should work as well if run it more than once
# NMF is slow so only run it once since this is the same as 
# other methods 
#   vesalius <- generate_embeddings(vesalius,
#                                   dim_reduction = "NMF",
#                                   normalisation = "log_norm")
#   # simple embedding build and tile build
#   expect_s4_class(vesalius, "vesalius_assay")
#   expect_true(get_active_embedding_tag(vesalius) == "NMF.1")
#   expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
#                     c("NMF", "NMF.1")))
})

test_that("Vesalius PCA_L embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA_L",
        dimensions = 3,
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "PCA_L")
    expect_identical(get_embeddings(vesalius, active = FALSE)$PCA_L,
        get_embeddings(vesalius, active = TRUE))
   
    # it should work as well if run it more than once
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA_L",
        dimensions = 3,
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    expect_true(get_active_embedding_tag(vesalius) == "PCA_L.1")
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("PCA_L", "PCA_L.1")))
})

test_that("Vesalius UMAP embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    # warning on UMAP start up from seurat
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "UMAP",
        normalisation = "log_norm")
    # simple embedding build and tile build

    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "UMAP")
    expect_identical(get_embeddings(vesalius, active = FALSE)$UMAP,
        get_embeddings(vesalius, active = TRUE))
    
    # it should work as well if run it more than once
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "UMAP",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    expect_true(get_active_embedding_tag(vesalius) == "UMAP.1")
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("UMAP", "UMAP.1")))
})

test_that("Vesalius LSI embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "LSI",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "LSI")
    expect_identical(get_embeddings(vesalius, active = FALSE)$LSI,
        get_embeddings(vesalius, active = TRUE))
   
    # it should work as well if run it more than once
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "LSI",
        normalisation = "log_norm",
        remove_lsi_1 = FALSE)
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    expect_true(get_active_embedding_tag(vesalius) == "LSI.1")
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("LSI", "LSI.1")))
})

test_that("Vesalius LSI UMAP embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "LSI_UMAP",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "LSI_UMAP")
    expect_identical(get_embeddings(vesalius, active = FALSE)$LSI_UMAP,
        get_embeddings(vesalius, active = TRUE))
    
    # it should work as well if run it more than once
    # we will also run it with removal of first embedding 
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "LSI_UMAP",
        normalisation = "log_norm",
        remove_lsi_1 = FALSE)
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    expect_true(get_active_embedding_tag(vesalius) == "LSI_UMAP.1")
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("LSI_UMAP", "LSI_UMAP.1")))

    
})

test_that("Vesalius mixed embeds", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "PCA")
    expect_identical(get_embeddings(vesalius, active = FALSE)$PCA,
        get_embeddings(vesalius, active = TRUE))
    # it should work as well if run it more than once
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "UMAP",
        normalisation = "log_norm")
    expect_true(get_active_embedding_tag(vesalius) == "UMAP")
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("PCA", "UMAP")))
    # Running the same embedding twice with multiple embeddings present 
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "PCA",
        normalisation = "log_norm")
    # simple embedding build and tile build
    expect_s4_class(vesalius, "vesalius_assay")
    # check if active embedding is set correctly
    expect_true(get_active_embedding_tag(vesalius) == "PCA.1")
    expect_identical(get_embeddings(vesalius, active = FALSE)$PCA.1,
        get_embeddings(vesalius, active = TRUE))
    expect_true(all(names(get_embeddings(vesalius, active = FALSE)) %in%
        c("PCA", "UMAP", "PCA.1")))
})

test_that("Vealius log_norm", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    # simple embeds
    vesalius <- generate_embeddings(vesalius,
        normalisation = "log_norm")
    expect_s4_class(vesalius, "vesalius_assay")
})

test_that("Vealius SCTransform", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    # simple embeds
    vesalius <- generate_embeddings(vesalius,
        normalisation = "SCTransform")
    expect_s4_class(vesalius, "vesalius_assay")
})

test_that("Vealius TFIDF", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    # signac throw warnin for zero counts
    expect_s4_class(generate_embeddings(vesalius,
        normalisation = "TFIDF"), "vesalius_assay")
})

test_that("Vealius raw", {
    vesalius <- build_vesalius_assay(coordinates, counts)
    # simple embeds - this is equivalent to no norm
    vesalius <- generate_embeddings(vesalius,
        normalisation = "raw")
    expect_s4_class(vesalius, "vesalius_assay")
})

test_that("Vealius Custom", {
    vesalius <- build_vesalius_assay(coordinates)
    vesalius <- add_counts(vesalius, counts, counts)
    vesalius <- generate_tiles(vesalius)
    # simple embeds
    vesalius <- generate_embeddings(vesalius,
        use_count = "custom_counts")
    expect_s4_class(vesalius, "vesalius_assay")
})


