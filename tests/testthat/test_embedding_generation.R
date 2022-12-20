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
        normalisation = "log_norm")
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
    vesalius <- generate_embeddings(vesalius,
        dim_reduction = "LSI_UMAP",
        normalisation = "log_norm")
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
    expect_warning(vesalius <- generate_embeddings(vesalius,
        normalisation = "TFIDF"))
    expect_s4_class(vesalius, "vesalius_assay")
})

# test_that("Vealius raw", {
#     vesalius <- build_vesalius_assay(coordinates, counts)
#     # simple embeds - this is equivalent to no norm
#     vesalius <- generate_embeddings(vesalius,
#         normalisation = "raw")
#     expect_s4_class(vesalius, "vesalius_assay")
# })

# test_that("Vealius Custom", {
#     vesalius <- build_vesalius_assay(coordinates)
#     vesalius <- add_counts(vesalius, counts, counts)
#     vesalius <- generate_tiles(vesalius)
#     # simple embeds
#     vesalius <- generate_embeddings(vesalius,
#         use_count = "custom_counts")
#     expect_s4_class(vesalius, "vesalius_assay")
# })