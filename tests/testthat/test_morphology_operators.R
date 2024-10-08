# load data and process it to get territories 
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts, scale = "auto")
vesalius <- generate_embeddings(vesalius,
    filter_threshold = 1,
    filter_grid = 1)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 10)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)
vesalius <- isolate_territories(vesalius)

test_that("Vesalius assay works as expect with morphing", {
    # expect error if no territory is specified
    expect_error(territory_morphing(vesalius))
    # expect error if wrong input
    expect_error(territory_morphing("vesalius"))
    expect_s4_class(territory_morphing(vesalius,
        territory = 1,
        trial = "Territory"),
        "vesalius_assay")
    # expect s4 if territory is specified
    expect_s4_class(territory_morphing(vesalius, territory = 1),
        "vesalius_assay")
    # expect error if territory does not exist 
    expect_error(territory_morphing(vesalius, territory = 42))
    # expect warning if not all territories exists
    expect_warning(territory_morphing(vesalius, territory = 4:42))
})

test_that("mophology operations work as expected", {
    # growing
    expect_s4_class(territory_morphing(vesalius,
        territory = 1,
        morphology_factor = 5),
        "vesalius_assay")
    # shrink
    expect_s4_class(territory_morphing(vesalius,
        territory = 1,
        morphology_factor = -5),
        "vesalius_assay")
    # clean
    expect_s4_class(territory_morphing(vesalius,
        territory = 1,
        morphology_factor = c(-5, 5)),
        "vesalius_assay")
    # fill
    expect_s4_class(territory_morphing(vesalius,
        territory = 1,
        morphology_factor = c(5, -5)),
        "vesalius_assay")
})

test_that("Vesalius assay works as expect with layering", {
    # expect error if no territory is specified
    expect_error(layer_territory(vesalius))
    # expect error if wrong input
    expect_error(layer_territory("vesalius"))
    # expect s4 if territory is specified
    expect_s4_class(layer_territory(vesalius, territory = 1),
        "vesalius_assay")
    # expect error if territory does not exist 
    expect_error(layer_territory(vesalius, territory = 42))
    # expect warning if not all territories exists
    expect_warning(layer_territory(vesalius, territory = 4:42))
})

test_that("terriotry layering works as expected", {
    # expect s4 if mophing is applied
    expect_s4_class(layer_territory(vesalius,
        morphology_factor = 5,
        territory = 1),
        "vesalius_assay")
    # checking if different depths work
    tmp <- layer_territory(vesalius,
        trial = "Territory",
        layer_depth = 2,
        territory = 1)
    expect_s4_class(tmp,
        "vesalius_assay")
    layers <- unique(get_territories(tmp)$Layer) %in% c("out", "1", "2")
    expect_true(all(layers))
    # If requesting to many layer
    expect_warning(layer_territory(vesalius,
        layer_depth = 42,
        territory = 1))
})