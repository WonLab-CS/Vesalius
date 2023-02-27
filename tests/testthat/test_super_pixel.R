# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)

test_that("Initialise super pixel centers", {
    expect_type(vesalius:::select_initial_indices(coordinates, k = 50), "list")
})