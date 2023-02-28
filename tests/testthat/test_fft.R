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
library(Matrix)
set.seed(1547)


#-----------------------------------------------------------------------------#
# create output directories
#-----------------------------------------------------------------------------#
input <- "/common/martinp4/SSv2"


#-----------------------------------------------------------------------------#
# Set future global for multicore processing
#-----------------------------------------------------------------------------#
plan(multicore, workers = 4)
max_size <- 1000 * 1024^2
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
f = 50
coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
count_mat <- read.table(counts[f], header = TRUE, row.names = 1)

vesalius <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, sigma = 5, iter = 10) %>%
    segment_image(dimensions = 1:30, col_resolution = 12) %>%
    isolate_territories()


f = 43
coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
count_mat <- read.table(counts[f], header = TRUE, row.names = 1)
vesalius_query <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, sigma = 5, iter = 10) %>%
    segment_image(dimensions = 1:30, col_resolution = 12) %>%
    isolate_territories()

test <- integrate_territories(vesalius,
    vesalius_query,
    global = TRUE,
    method = "coherence",
    start = "connected")

coherence_x <- test$sim$x
coherence_y <- test$sim$y
colnames(coherence_x) <- paste0("seed_", colnames(coherence_x))
rownames(coherence_x) <- paste0("query_", rownames(coherence_x))

colnames(coherence_y) <- paste0("seed_", colnames(coherence_y))
rownames(coherence_y) <- paste0("query_", rownames(coherence_y))

seed <- lapply(test$seed, function(path) {
    return(list("x" = fft(path$x),
      "y" = fft(path$y)))
})

query <- lapply(test$query, function(path) {
    return(list("x" = fft(path$x),
      "y" = fft(path$y)))
})
