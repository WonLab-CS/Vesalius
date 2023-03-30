# load vesalius data 
# We will assume that the embeddings have been produced 
#library(anndata)
#spat_data <- read_h5ad("/common/wonklab/CosMX/PembroRT_cosmx_RUBY.h5ad")
#sc_data <- read_h5ad("/common/wonklab/CosMX/")
#data(vesalius)

#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(vesalius, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")
library(dplyr)
library(future)
set.seed(1547)


#-----------------------------------------------------------------------------#
# create output directories
#-----------------------------------------------------------------------------#
input <- "/common/martinp4/SSv2"


#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#
plan(multicore, workers = 4)
max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)

#-----------------------------------------------------------------------------#
# Generate Brain Images
#-----------------------------------------------------------------------------#
coordinates <- list.files(path = input,
    pattern = "location|coord|Locations", full.names = TRUE)
counts <- list.files(path = input,
    pattern = "expression", full.names = TRUE)
tag <- list.files(path = input,
    pattern = "expression", full.names = FALSE)
tag <- gsub(".digital_expression.txt.gz|_expression_matrix.mtx.gz|.sparse_expression.txt",
    "", tag)
f = 51
coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
count_mat <- read.table(counts[f], header = TRUE, row.names = 1)

# vesalius <- build_vesalius_assay(coord, count_mat) %>%
#     generate_embeddings(tensor_resolution = 0.3) %>%
#     regularise_image(dimensions = 1:30, lambda = 5) %>%
#     equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
#     smooth_image(dimensions = 1:30, sigma = 5, iter = 10) %>%
#     segment_image(dimensions = 1:30, col_resolution = 12) %>%
#     isolate_territories()

vesalius <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    #smooth_image(dimensions = 1:30, method =c("iso", "box"), sigma = 2, box = 10, iter = 10) %>%
    segment_image(dimensions = 1:30, col_resolution = 1000, method = "slic",threshold = "85%") %>%
    segment_image(dimensions = 1:30, col_resolution = 12, method = "kmeans") 
pdf("test.pdf")
image_plot(vesalius)
territory_plot(vesalius)
dev.off()
f = 44
coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
count_mat <- read.table(counts[f], header = TRUE, row.names = 1)
# vesalius_query <- build_vesalius_assay(coord, count_mat) %>%
#     generate_embeddings(tensor_resolution = 0.3) %>%
#     regularise_image(dimensions = 1:30, lambda = 5) %>%
#     equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
#     smooth_image(dimensions = 1:30, sigma = 5, iter = 10) %>%
#     segment_image(dimensions = 1:30, col_resolution = 12) %>%
#     isolate_territories()

vesalius_query <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, method =c("iso", "box"), sigma = 2, box = 10, iter = 10) %>%
    segment_image(dimensions = 1:30, col_resolution = 15) %>%
    isolate_territories(min_spatial_index = 50)

test <- integrate_by_territory(vesalius,
    vesalius_query,
    method = "coherence",
    use_counts = TRUE,
    k = 5,
    use_norm = "log_norm")



# coherence_x <- test$sim$x
# coherence_y <- test$sim$y
# colnames(coherence_x) <- paste0("seed_", colnames(coherence_x))
# rownames(coherence_x) <- paste0("query_", rownames(coherence_x))

# colnames(coherence_y) <- paste0("seed_", colnames(coherence_y))
# rownames(coherence_y) <- paste0("query_", rownames(coherence_y))

# seed <- lapply(test$seed, function(path) {
#     return(list("x" = fft(path$x),
#       "y" = fft(path$y)))
# })

# query <- lapply(test$query, function(path) {
#     return(list("x" = fft(path$x),
#       "y" = fft(path$y)))
# })
coh <- integrate_by_territory(vesalius,
    vesalius,
    method = "coherence",
    use_counts = TRUE,
    use_norm = "log_norm")
territories <- vesalius:::check_territory_trial(vesalius, "last") %>%
    filter(trial != "isolated")
graph <- vesalius:::generate_territory_graph(territories, k = 5)
score <- vesalius:::score_neighbor_graph(graph, graph,coh)

########
library(vesalius)
library(ggplot2)
library(dplyr)
data(vesalius)
load("Scenes/super_pixel/jitter_ves.Rda")
vesalius <- build_vesalius_assay(coordinates, counts)
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)

vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 5)
vesalius <- segment_image(vesalius,
    method = "slic",
    dimensions = 1:3,
    col_resolution = 50,
    compactness = 1,
    index_selection = "bubble",
    scaling = 0.2)


jitter_ves <- generate_embeddings(jitter_ves)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 5)
jitter_ves <- segment_image(jitter_ves,
    method = "slic",
    dimensions = 1:3,
    col_resolution = 50,
    compactness = 1,
    index_selection = "bubble",
    scaling = 0.2)


test <- integrate_assays(vesalius, jitter_ves, n_centers = 100)




# Running test for Vesalius objects
# Loading data from the packages
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts)
vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 5)
vesalius <- segment_image(vesalius,
    method = "slic",
    dimensions = 1:3,
    col_resolution = 50,
    compactness = 1,
    index_selection = "bubble",
    scaling = 0.2)


vesalius <- segment_image(vesalius, method = "kmeans", col_resolution = 10)
image_plot(vesalius, dimensions = 1:3)

test_that("Initialise super pixel centers", {
    expect_type(vesalius:::select_initial_indices(coordinates, k = 50), "list")
})

par(mfrow =c(1,2))
plot(images)
plot(images *ratio)