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
library(Matrix)
library(deldir)
library(imager)
library(imagerExtra)
library(pwr, lib.loc = "/common/martinp4/R")
library(gsignal, lib.loc = "/common/martinp4/R")
library(kohonen, lib.loc = "/common/martinp4/R")
library(registry, lib.loc = "/common/martinp4/R")
library(rngtools, lib.loc = "/common/martinp4/R")
library(NMF, lib.loc = "/common/martinp4/R")
library(RcppHungarian, lib.loc = "/common/martinp4/R")
library(spatstat.utils, lib.loc = "/common/martinp4/R")
library(geometry, lib.loc = "/common/martinp4/R")
library(vesalius, lib.loc = "/common/martinp4/R")



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
    generate_embeddings(dim_reduction = "PCA", tensor_resolution = 0.3) %>%
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



test <- integrate_horizontally(vesalius,
    vesalius_query,
    n_centers = 15,
    depth = 4,
    signal = "features")

test <- generate_embeddings(test,tensor_resolution = 0.9999)
# test <- #regularise_image(test, dimensions = 1:30, lambda = 5) %>%
#     equalize_image(test,dimensions = 1:30, sleft = 5, sright = 5)# %>%
    #smooth_image(dimensions = 1:30, method =c("iso", "box"), sigma = 2, box = 10, iter = 10)
pdf("test_bary.pdf")
image_plot(test)
dev.off()

g <- ggplot(query, aes(x,y,col = as.factor(Segment))) + geom_point(size = 0.3)
g1 <- ggplot(query, aes(x_new,y_new, col = as.factor(Segment))) + geom_point(size=0.3)
pdf("test.pdf", width = 10, height = 5)
g+g1
dev.off()

pdf("test.pdf")
plot(0, type = "n", xlim = c(0, 650), ylim = c(0, 650))
points(seed_centers$x, seed_centers$y,col="black", cex = 3)
points(query_centers$x, query_centers$y,pch = 20, cex = 3, col="red")

for(i in seq_len(nrow(query_centers))){
    x <- c(query_centers$x[anchors$to[i]],
        seed_centers$x[anchors$from[i]])
    y <- c(query_centers$y[anchors$to[i]],
        seed_centers$y[anchors$from[i]])
    x_new <- c(query_centers$x[anchors$to[i]], query_centers$x_new[anchors$to[i]])
    y_new <- c(query_centers$y[anchors$to[i]], query_centers$y_new[anchors$to[i]])
    lines(x,y, col = "green", lwd = 2)
    lines(x_new,y_new, col = "blue", lwd = 2)
}
dev.off()

g <- image_plot(vesalius, embedding = "PCA")
g1 <- image_plot(test, embedding = "PCA")
g2 <- image_plot(vesalius_query, embedding = "PCA")

pdf("test_counts.pdf", width = 24, height = 8)
print(g + g2 + g1)
dev.off()


plot(0, type = "n", xlim = c(0, 650), ylim = c(0, 650))
points(seed_spix$x, seed_spix$y,col="black", cex = 3)
points(query_spix$x, query_spix$y,pch = 20, cex = 3, col="red")
for(i in seq_len(nrow(anchor_map))){
    x <- c(anchor_map$x[i], query_spix$x[anchor_map$from[i]])
    y <- c(anchor_map$y[i], query_spix$y[anchor_map$from[i]])
    lines(x,y, lwd = 2, col = "blue")
}
dev.off()



## Vertical 
# data prep
rna <- "/common/wonklab/Spatial_CITE/GSM6578061_mousekidney_RNA.tsv"
prot <- "/common/wonklab/Spatial_CITE/GSM6578070_mousekidney_protein.tsv"

rna <- read.table(rna, header = TRUE, sep = "\t")
prot <- read.table(prot, header = TRUE, sep = "\t")

CITE_coordinates <- function(assay) {
    locations <- assay[, 1]
    coordinates <- strsplit(locations,"x")
    x <- sapply(coordinates,"[[", 1)
    y <- sapply(coordinates,"[[", 2)
    return(data.frame("barcodes" = locations,
        "x" = x,
        "y" = y))
}

CITE_counts <- function(assay) {
    locations <- assay[, 1]
    genes <- colnames(assay)[-1]
    counts <- t(assay[, -1])
    colnames(counts) <- locations
    rownames(counts) <- genes
    return(counts)
}

rna_coord <- CITE_coordinates(rna)
rna_counts <- CITE_counts(rna)

prot_coord <- CITE_coordinates(prot)
prot_counts <- CITE_counts(prot)

ves_rna <- build_vesalius_assay(rna_coord, rna_counts)
ves_rna <- generate_embeddings(ves_rna, filter_grid = 1) %>%
    equalize_image(dimensions = 1:5, sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:5,sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:5,col_resolution = 10) %>%
    isolate_territories()

embed_plot <- vector("list", 30)
for(i in seq_along(embed_plot)){
    embed_plot[[i]] <- image_plot(ves_rna, dimensions = i,
        embedding = "PCA")
}
pdf("rna_embed.pdf", width = 30, height = 30)
print(wrap_plots(embed_plot, nrow = 6, ncol = 5))
dev.off()
pdf("rna_territories.pdf", width = 10, height = 8)
print(territory_plot(ves_rna, cex_pt = 5))
dev.off()

ves_prot <- build_vesalius_assay(prot_coord,prot_counts)
ves_prot <- generate_embeddings(ves_prot, filter_grid = 1) %>%
    equalize_image(dimensions = 1:5,sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:5,sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:5,col_resolution = 10) %>%
    isolate_territories()

embed_plot <- vector("list", 30)
for(i in seq_along(embed_plot)){
    embed_plot[[i]] <- image_plot(ves_prot, dimensions = i,
        embedding = "PCA")
}
pdf("prot_embed.pdf", width = 30, height = 30)
print(wrap_plots(embed_plot, nrow = 6, ncol = 5))
dev.off()
pdf("prot_territories.pdf", width = 10, height = 8)
print(territory_plot(ves_prot, cex_pt = 5))
dev.off()

integrated <- integrate_vertically(ves_rna,
    ves_prot,
    dimensions = 1:10,
    method = "interlace")
integrated <- smooth_image(integrated, sigma = 2, iter = 5)
integrated <- segment_image(integrated, dimensions = 1:10)
integrated <- isolate_territories(integrated)
pdf("test_vertical_int_interlace.pdf", width = 12, height = 10)
territory_plot(integrated, cex_pt= 5)
dev.off()

integrated <- integrate_vertically(ves_rna,
    ves_prot,
    dimensions = 1:10,
    method = "mean")
integrated <- smooth_image(integrated, sigma = 2, iter = 5)
integrated <- segment_image(integrated, dimensions = 1:10)
integrated <- isolate_territories(integrated)
pdf("test_vertical_int_mean.pdf", width = 12, height = 10)
territory_plot(integrated, cex_pt= 5)
dev.off()

integrated <- integrate_vertically(ves_rna,
    ves_prot,
    dimensions = 1:10,
    method = "concat")
integrated <- equalize_image(integrated, dimensions = 1:10,sleft = 2, sright = 2)
integrated <- smooth_image(integrated, dimensions = 1:10,sigma = 1, iter = 5)
integrated <- segment_image(integrated, dimensions = 1:10, col_resolution = 10)
integrated <- isolate_territories(integrated)
pdf("test_vertical_int_concat.pdf", width = 12, height = 10)
territory_plot(integrated, cex_pt= 5)
dev.off()




########
set.seed(145)
library(vesalius)
library(ggplot2)
library(dplyr)
library(kohonen)
data(vesalius)
load("Scenes/super_pixel/jitter_ves.Rda")
vesalius <- build_vesalius_assay(coordinates, counts)
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)

vesalius <- generate_embeddings(vesalius)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 5)
# vesalius <- segment_image(vesalius,
#     method = "som",
#     dimensions = 1:3,
#     col_resolution = 10,
#     compactness = 1,
#     scaling = 0.2)


jitter_ves <- generate_embeddings(jitter_ves)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 5)
# jitter_ves <- segment_image(jitter_ves,
#     method = "slic",
#     dimensions = 1:3,
#     col_resolution = 10,
#     compactness = 1,
#     scaling = 0.2)


test <- integrate_horizontally(vesalius,
    jitter_ves,
    signal = "variable_features",
    mapping = "mesh",
    index_selection = "bubble",
    n_landmarks = 50,
    threshold = 0.2,
    grid = 80)


test <- generate_embeddings(test)
test <- equalize_image(test, sleft = 5, sright = 5)
test <- smooth_image(test, sigma = 5, iter = 5)
#test <- segment_image(test, col_resolution = 5)


image_plot(test) + image_plot(jitter_ves)

plot(0, type = "n", xlim = c(0, 650), ylim = c(0, 650))
points(seed_centers$x, seed_centers$y,col="black", cex = 3)
points(query_centers$x, query_centers$y,pch = 20, cex = 3, col="red")

for(i in seq_len(nrow(query_centers))){
    # x <- c(query_centers$x[anchors$from[i]],
    #     seed_centers$x[anchors$to[i]])
    # y <- c(query_centers$y[anchors$from[i]],
    #     seed_centers$y[anchors$to[i]])
    x_new <- c(query_centers$x[i], query_centers$x_new[i])
    y_new <- c(query_centers$y[i], query_centers$y_new[i])
    # lines(x,y, col = "green", lwd = 2)
    lines(x_new,y_new, col = "blue", lwd = 2)
}



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