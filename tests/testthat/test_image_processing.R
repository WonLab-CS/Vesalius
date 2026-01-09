# load vesalius data 
# We will assume that the embeddings have been produced 

data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts)
vesalius <- generate_embeddings(vesalius)

test_that("Vesalius assay work as expected with reg", {
    expect_s4_class(regularise_image(vesalius), "vesalius_assay")
    expect_error(regularise_image("vesalius"))
    # expect s4 even with a single dim 
    expect_s4_class(regularise_image(vesalius, dimensions = 1),
        "vesalius_assay")
    # expect error if dims are out of bounds 
    expect_error(regularise_image(vesalius, dimensions = seq(1, 40)))
    # expect s4 if you make call to another embedding
    expect_s4_class(regularise_image(vesalius, embedding = "PCA"),
        "vesalius_assay")
    # expect error if unkown input to embedding
    expect_error(smooth_image(vesalius, embedding = "funky"))
    # expect warning if more than one embedding contain the same name
    tmp <- generate_embeddings(vesalius)
    expect_s4_class(regularise_image(tmp, embedding = "PCA"),
        "vesalius_assay")
    # expect s4 class if you call the correct name
    expect_s4_class(regularise_image(tmp, embedding = "PCA.1"),
        "vesalius_assay")
})

test_that("Regularisation works as expected", {
    expect_error(regularise_image(vesalius, lambda = 0))
    expect_error(regularise_image(vesalius, niter = 0))
    # checking if values are nornalised correctly
    tmp <- regularise_image(vesalius)
    expect_true(max(tmp@active) == 1)
    expect_true(min(tmp@active) == 0)
})

test_that("Vesalius assay work as expected with smooth", {
    expect_s4_class(smooth_image(vesalius), "vesalius_assay")
    expect_error(smooth_image("vesalius"))
    # expect s4 even with a single dim 
    expect_s4_class(smooth_image(vesalius, dimensions = 1),
        "vesalius_assay")
    # expect error if smooth_image are out of bounds 
    expect_error(smooth_image(vesalius, dimensions = seq(1, 40)))
    # expect s4 if you make call to another embedding
    expect_s4_class(smooth_image(vesalius, embedding = "PCA"),
        "vesalius_assay")
    # expect error if unkown input to embedding
    expect_error(smooth_image(vesalius, embedding = "funky"))
    # expect warning if more than one embedding contain the same name
    tmp <- generate_embeddings(vesalius)
    expect_s4_class(smooth_image(tmp, embedding = "PCA"),
        "vesalius_assay")
    # expect s4 class if you call the correct name
    expect_s4_class(smooth_image(tmp, embedding = "PCA.1"),
        "vesalius_assay")
    # wrong values to smooth types 
    expect_error(smooth_image(vesalius, iter = 2, method = "funky"))
})

test_that("smoothing works as expected", {
    # rotate if more than one smoothing iteration
    expect_s4_class(smooth_image(vesalius, iter = 2),
        "vesalius_assay")
    #running more that one method
    expect_s4_class(smooth_image(vesalius, method = c("iso", "box")),
        "vesalius_assay")
    # running smoothing with vector of values 
    expect_s4_class(smooth_image(vesalius, method = c("iso"),
        sigma = seq(1, 3)),
        "vesalius_assay")
    # running smoothing with vector of values 
    expect_s4_class(smooth_image(vesalius, method = c("box"),
        box = seq(3, 12)),
        "vesalius_assay")
    # running smoothing with vector of values 
    expect_s4_class(smooth_image(vesalius, method = c("median"),
        box = seq(3, 12)),
        "vesalius_assay")
    # running smoothing with vector of values and multiple methods 
    expect_s4_class(smooth_image(vesalius, method = c("iso", "box"),
        box = seq(3, 12),
        sigma = seq(1, 3)),
        "vesalius_assay")
    # running smoothing with vector of values and multiple methods 
    expect_s4_class(smooth_image(vesalius, method = c("iso", "box"),
        box = seq(3, 12),
        sigma = seq(1, 3),
        across_levels = "max"),
        "vesalius_assay")
    # running smoothing with vector of values and multiple methods 
    expect_s4_class(smooth_image(vesalius, method = c("iso", "box"),
        box = seq(3, 12),
        sigma = seq(1, 3),
        across_levels = "mean"),
        "vesalius_assay")
    # expect error if incorrect input to across levels
    expect_error(smooth_image(vesalius, method = c("iso", "box"),
        box = seq(3, 12),
        sigma = seq(1, 3),
        across_levels = "funky"))
})

test_that("Vesalius assay work as expected with eq", {
    expect_s4_class(equalize_image(vesalius), "vesalius_assay")
    expect_error(equalize_image("vesalius"))
    # expect s4 even with a single dim 
    expect_s4_class(equalize_image(vesalius, dimensions = 1),
        "vesalius_assay")
    # expect error if smooth_image are out of bounds
    expect_error(equalize_image(vesalius, dimensions = seq(1, 40)))
    # expect s4 if you make call to another embedding
    expect_s4_class(equalize_image(vesalius, embedding = "PCA"),
        "vesalius_assay")
    # expect error if unkown input to embedding
    expect_error(smooth_image(vesalius, embedding = "funky"))
    # expect warning if more than one embedding contain the same name
    tmp <- generate_embeddings(vesalius)
    expect_s4_class(equalize_image(tmp, embedding = "PCA"),
        "vesalius_assay")
    # expect s4 class if you call the correct name
    expect_s4_class(equalize_image(tmp, embedding = "PCA.1"),
        "vesalius_assay")
    # wrong input to method 
    expect_error(equalize_image(vesalius, method = "funky"))
})

test_that("eq wroks as expected", {
    expect_s4_class(equalize_image(vesalius, method = "BalanceSimplest"),
        "vesalius_assay")
    expect_s4_class(equalize_image(vesalius, method = "EqualizePiecewise"),
        "vesalius_assay")
    expect_s4_class(equalize_image(vesalius, method = "SPE"),
        "vesalius_assay")
    expect_s4_class(equalize_image(vesalius, method = "EqualizeDP"),
        "vesalius_assay")
    expect_s4_class(equalize_image(vesalius, method = "EqualizeADP"),
        "vesalius_assay")
    expect_s4_class(equalize_image(vesalius, method = "ECDF"),
        "vesalius_assay")
})

test_that("Vesalius assay work as expected in image seg", {
    expect_s4_class(segment_image(vesalius), "vesalius_assay")
    expect_error(segment_image("vesalius"))
    # expect s4 even with a single dim 
    expect_s4_class(segment_image(vesalius, dimensions = 1),
        "vesalius_assay")
    # expect error if smooth_image are out of bounds
    expect_error(segment_image(vesalius, dimensions = seq(1, 40)))
    # expect s4 if you make call to another embedding
    expect_s4_class(segment_image(vesalius, embedding = "PCA"),
        "vesalius_assay")
    # expect error if unkown input to embedding
    expect_error(smooth_image(vesalius, embedding = "funky"))
    # expect warning if more than one embedding contain the same name
    tmp <- generate_embeddings(vesalius)
    expect_s4_class(segment_image(tmp, embedding = "PCA"),
        "vesalius_assay")
    # expect s4 class if you call the correct name
    expect_s4_class(segment_image(tmp, embedding = "PCA.1"),
        "vesalius_assay")
    # wrong input to method 
    expect_error(segment_image(vesalius, method = "funky"))
})

test_that("Image segmentation works as expected", {
    expect_s4_class(segment_image(vesalius, method = "kmeans"),
        "vesalius_assay")
    tmp <- segment_image(vesalius, method = "kmeans")
    # expect 10 cluster from kmeans 
    expect_true(length(unique(get_territories(tmp)$Segment)) == 10)
    # If adding a new run expect new segment trial
    expect_s4_class(segment_image(vesalius, method = "louvain"),
        "vesalius_assay")
    tmp <- get_territories(segment_image(vesalius, method = "louvain"))
    expect_true(all(colnames(names) %in%
        c("barcodes", "x", "y", "Segment", "Segment.1")))
    # expect s4 with leiden
    expect_s4_class(segment_image(vesalius, method = "leiden"),
        "vesalius_assay")
    # expect s4 with slic
    expect_s4_class(generate_spix(vesalius, method = "slic"),
        "vesalius_assay")
})

test_that("Vesalius assay works as expect with isolate ter", {
    #throws error if no segment have been computed
    expect_error(isolate_territories(vesalius))
    tmp <- segment_image(vesalius)
    expect_s4_class(isolate_territories(tmp), "vesalius_assay")
    tmp <- segment_image(tmp)
    #expect warning if more than one trial exists and same name
    expect_s4_class(isolate_territories(tmp, trial = "Segment"),
        "vesalius_assay")
    expect_error(isolate_territories(tmp, trial = "funky"))
    # expect error if input to method is not available 
    expect_error(isolate_territories(tmp, method = "funky"))
    tmp <- isolate_territories(tmp)
    expect_output(show(tmp))
})

