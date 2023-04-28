# load vesalius data 
# We will assume that the embeddings have been produced 
#library(anndata)
#spat_data <- read_h5ad("/common/wonklab/CosMX/PembroRT_cosmx_RUBY.h5ad")
#sc_data <- read_h5ad("/common/wonklab/CosMX/")
#data(vesalius)

#-----------------------------------------------------------------------------#
# Loading Libraries and setting seed
#-----------------------------------------------------------------------------#
library(dplyr)
library(future)
library(ggplot2)
library(dplyr)
library(patchwork)
library(vesalius, lib.loc = "/common/martinp4/R")
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")


set.seed(1547)


#-----------------------------------------------------------------------------#
# create output directories
#-----------------------------------------------------------------------------#
input <- "/common/wonklab/SSv2"


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



vesalius <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, method =c("iso", "box"), sigma = 2, box = 10, iter = 10) 

f = 44
coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
count_mat <- read.table(counts[f], header = TRUE, row.names = 1)


vesalius_query <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, method =c("iso", "box"), sigma = 2, box = 10, iter = 10)



test <- integrate_assays(vesalius,
    vesalius_query,
    n_centers = 50,
    compactness = 100,
    index_selection = "random",
    signal = "features",
    threshold = 0.7)

test <- generate_embeddings(test)



g <- image_plot(vesalius, embedding = "PCA")
g1 <- image_plot(test, embedding = "PCA")


pdf("test_counts.pdf", width = 16, height = 8)
print(g + g1)
dev.off()



########
set.seed(145)
library(vesalius)
library(ggplot2)
library(dplyr)
data(vesalius)
load("Scenes/super_pixel/jitter_ves.Rda")
vesalius <- build_vesalius_assay(coordinates, counts)
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)

vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 5)
# vesalius <- segment_image(vesalius,
#     method = "slic",
#     dimensions = 1:3,
#     col_resolution = 50,
#     compactness = 1,
#     index_selection = "random",
#     scaling = 0.2)


jitter_ves <- generate_embeddings(jitter_ves)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 5)
# jitter_ves <- segment_image(jitter_ves,
#     method = "slic",
#     dimensions = 1:3,
#     col_resolution = 50,
#     compactness = 1,
#     index_selection = "bubble",
#     scaling = 0.2)


test <- integrate_assays(vesalius,
    jitter_ves,
    compactness = 10,
    index_selection = "random",
    signal = "features",
    n_centers = 50)

test <- generate_embeddings(test)
test <- equalize_image(test, embedding = "PCA",sleft = 5, sright = 5)
test <- smooth_image(test, embedding = "PCA", sigma = 5, iter = 5)
test <- segment_image(test, col_resolution = 5)

g <- image_plot(test$seed)
g1 <- image_plot(test$query)
g2 <- image_plot(test$integrate)

print(g + g1 + g2)


seed_score <- igraph::graph_from_data_frame(test$seed_score, directed = FALSE)
E(seed_score)$weight <- test$seed_score$score

query_score <- igraph::graph_from_data_frame(test$query_score, directed = FALSE)
E(seed_score)$weight <- test$query_score$score

iso <- igraph::isomorphisms(graph1 = seed_score, graph2 = query_score)



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
    scaling = 0.5)


vesalius <- segment_image(vesalius, method = "kmeans", col_resolution = 10)
image_plot(vesalius, dimensions = 1:3)

test_that("Initialise super pixel centers", {
    expect_type(vesalius:::select_initial_indices(coordinates, k = 50), "list")
})

par(mfrow =c(1,2))
plot(images)
plot(images *ratio)