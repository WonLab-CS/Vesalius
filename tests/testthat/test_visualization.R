# load data and process it to get territories 
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts)
vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, sigma = 5, iter = 15)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)
cells <- sample(get_territories(vesalius)$barcodes, 1040, replace = FALSE)
cell_error <- get_territories(vesalius)[
    get_territories(vesalius)$Territory == 1, ]$barcodes
cell_error <- sample(cell_error, 100, replace = FALSE)

test_that("image plot", {
    expect_true(is(image_plot(vesalius), "gg"))
    expect_true(is(image_plot(vesalius, dimensions = 1), "gg"))
    expect_true(is(image_plot(vesalius, embedding = "PCA"), "gg"))
    expect_error(image_plot(vesalius, dimensions = 1:4))
    expect_error(image_plot(vesalius, dimensions = 1:2))
    tmp <- segment_image(vesalius, dimensions = seq(1, 30),
        method = "louvain",
        col_resolution = 0.01)
    # related to minmax misc function - could be test elsewhere as well.
    expect_warning(image_plot(vesalius, dimensions = 1))
})

test_that("territory plot", {
    expect_true(is(territory_plot(vesalius), "gg"))
    expect_true(is(territory_plot(vesalius, split = TRUE), "gg"))
    expect_true(is(territory_plot(vesalius, randomise = FALSE), "gg"))
})

test_that("view gene expression", {
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1"),
        "gg"))
    expect_true(is(view_gene_expression(vesalius,
        genes = c("Malat1", "Pcp4")),
        "ggarrange"))
    expect_error(view_gene_expression(vesalius))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        territory_1 = 1),
        "gg"))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        territory_2 = 1),
        "gg"))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        territory_1 = 1,
        territory_2 = 8),
        "gg"))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        as_layer = TRUE),
        "gg"))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        territory_1 = 1,
        cells = cells),
        "gg"))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        territory_2 = 1,
        cells = cells),
        "gg"))
    expect_true(is(view_gene_expression(vesalius, genes = "Malat1",
        territory_1 = 1,
        territory_2 = 8,
        cells = cells),
        "gg"))
    expect_true(is(view_gene_expression(vesalius,
        genes = "Malat1",
        norm = FALSE),
        "gg"))
    expect_error(view_gene_expression(vesalius,
        genes = "funky"))
})