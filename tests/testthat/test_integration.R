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
library(Morpho)
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

sub_sample <- sample(seq_len(nrow(coord)), 3000, replace = FALSE)
sub_barcodes <- coord$barcodes[sub_sample]
coord <- coord[coord$barcodes %in% sub_barcodes, ]
count_mat <- count_mat[, colnames(count_mat) %in% sub_barcodes]

vesalius <- build_vesalius_assay(coord, count_mat) %>%
    generate_embeddings(dim_reduction = "PCA",filter_threshold =1, filter_grid=1, tensor_resolution = 0.15) %>%
    equalize_image(dimensions = 1:15, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:15, method =c("iso", "box"), sigma = 2, box = 15, iter = 10) %>%
    segment_image(dimensions = 1:15, col_resolution = 15)%>%
    isolate_territories()

f = 44
coord <- read.csv(coordinates[f], header = FALSE, skip = 1)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
count_mat <- read.table(counts[f], header = TRUE, row.names = 1)
sub_sample <- sample(seq_len(nrow(coord)), 3000, replace = FALSE)
sub_barcodes <- coord$barcodes[sub_sample]
coord <- coord[coord$barcodes %in% sub_barcodes, ]
count_mat <- count_mat[, colnames(count_mat) %in% sub_barcodes]

vesalius_query <- build_vesalius_assay(coord, count_mat) %>%
     generate_embeddings(dim_reduction = "PCA",filter_threshold =1, filter_grid=1, tensor_resolution = 0.15) %>%
    equalize_image(dimensions = 1:15, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:15, method =c("iso", "box"), sigma = 2, box = 15, iter = 10) %>%
    segment_image(dimensions = 1:15, col_resolution = 15)%>%
    isolate_territories()


test <- integrate_horizontally(vesalius,
    vesalius_query,
    k = 20,
    mapping = "exact",
    signal = "variable_features")

test <- generate_embeddings(test, tensor_resolution = 1, filter_grid = 1, filter_threshold =1)%>%
    equalize_image(dimensions = 1:15, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:15, method =c("iso", "box"), sigma = 2, box = 15, iter = 10)%>%
    segment_image(dimensions = 1:15, col_resolution = 15)%>%
    isolate_territories()



g1 <- image_plot(vesalius, embedding = "PCA") + 
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Reference")
g2 <- image_plot(vesalius_query, embedding = "PCA")+
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Query")
g3 <- image_plot(test, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Aligned")

pdf("horizontal_embeddings.pdf", width = 24, height = 8)
print(g1 + g2 + g3)
dev.off()


library(RColorBrewer)
cols <- c("#a77447ff","#c98c57ff","#e3ab80ff","#f1dbceff","#cf826eff",
    "#9e5652ff","#cd8581ff","#d7a7a5ff","#dde2f4ff",
    "#9babe1ff","#5372c0ff","#304476ff","#111b34ff")
palette <- colorRampPalette(cols)


colors <- palette(length(unique(vesalius@territories$Territory)))
#colors <- colors[sample(seq_along(colors), size = length(colors))]
ter_1 <- territory_plot(vesalius, cex_pt = 3, alpha =1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Reference")
colors <- palette(length(unique(vesalius_query@territories$Territory)))
#colors <- colors[sample(seq_along(colors), size = length(colors))]
ter_2 <- territory_plot(vesalius_query, cex_pt = 3, alpha=1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Query")
colors <- palette(length(unique(test@territories$Territory)))
#colors <- colors[sample(seq_along(colors), size = length(colors))]
ter_3 <- territory_plot(test, cex_pt = 3, alpha=1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Aligned")

pdf("horizontal_territories.pdf", width = 24, height = 8)
print(ter_1 + ter_2 + ter_3)
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
gene_rna <- toupper(rownames(rna_counts))

prot_coord <- CITE_coordinates(prot)
prot_counts <- CITE_counts(prot)
prot_counts <- prot_counts[, -ncol(prot_counts)]
gene_prot <- toupper(rownames(prot_counts))



ves_rna <- build_vesalius_assay(rna_coord, rna_counts)
ves_rna <- generate_embeddings(ves_rna, nfeatures = 2000, filter_grid = 1) %>%
    equalize_image(dimensions = 1:5, sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:5,sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:5,col_resolution = 10) %>%
    isolate_territories()



ves_prot <- build_vesalius_assay(prot_coord,prot_counts)
ves_prot <- generate_embeddings(ves_prot, nfeatures = 199, filter_grid = 1) %>%
    equalize_image(dimensions = 1:5,sleft = 2, sright = 2) %>%
    smooth_image(dimensions = 1:5,sigma = 1, iter = 5) %>%
    segment_image(dimensions = 1:5,col_resolution = 10) %>%
    isolate_territories()



integrated <- integrate_vertically(ves_rna,
    ves_prot,
    norm_method = "log_norm",
    dimensions = 1:5,
    method = "concat")
integrated <- equalize_image(integrated, dimensions = 1:5,sleft = 2, sright = 2)
integrated <- smooth_image(integrated, dimensions = 1:5,sigma = 1, iter = 5)
integrated <- segment_image(integrated, dimensions = 1:5, col_resolution = 10)
integrated <- isolate_territories(integrated)




library(RColorBrewer)
# cols <- c("#a77447ff","#c98c57ff","#e3ab80ff","#f1dbceff","#cf826eff",
#     "#9e5652ff","#cd8581ff","#d7a7a5ff","#dde2f4ff",
#     "#9babe1ff","#5372c0ff","#304476ff","#111b34ff")
# palette <- colorRampPalette(cols)

cols <- c("#a77447ff","#c98c57ff","#e3ab80ff","#f1dbceff","#cf826eff",
    "#9e5652ff","#cd8581ff","#d7a7a5ff","#dde2f4ff",
    "#9babe1ff","#5372c0ff","#304476ff","#111b34ff")
palette <- colorRampPalette(cols)

colors <- palette(length(unique(ves_rna@territories$Territory)))
#colors <- colors[sample(seq_along(colors), size = length(colors))]
v_1 <- territory_plot(ves_rna, cex_pt = 4, alpha =1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "RNA Territories")

colors <- palette(length(unique(ves_prot@territories$Territory)))
#colors <- colors[sample(seq_along(colors), size = length(colors))]
v_2 <- territory_plot(ves_prot, cex_pt = 4, alpha =1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Protein Territories")

colors <- palette(length(unique(integrated@territories$Territory)))
#colors <- colors[sample(seq_along(colors), size = length(colors))]
v_3 <- territory_plot(integrated, cex_pt = 4, alpha =1) +
    #scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Joint Territories")

pdf("vertical_territories_clean.pdf", width = 24, height = 8)
print(v_1 + v_2 + v_3)
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

########
set.seed(145)
library(vesalius)
library(ggplot2)
library(dplyr)
data(vesalius)
load("Scenes/super_pixel/jitter_ves.Rda")
vesalius <- build_vesalius_assay(coordinates, counts)
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)

vesalius <- generate_embeddings(vesalius, filter_threshold = 1, filter_grid = 1)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 10)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)

jitter_ves <- generate_embeddings(jitter_ves, filter_threshold = 1, filter_grid = 1)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 10)
jitter_ves <- equalize_image(jitter_ves, sleft = 5, sright = 5)
jitter_ves <- segment_image(jitter_ves, col_resolution = 2)
jitter_ves <- isolate_territories(jitter_ves)

test <- integrate_horizontally(vesalius,
    jitter_ves,
    signal = "variable_features",
    map = "exact")



test <- generate_embeddings(test, filter_threshold = 1, filter_grid = 1)
test <- smooth_image(test, embedding = "PCA", sigma = 5, iter = 10)
test <- equalize_image(test, sleft = 5, sright = 5)
test <- segment_image(test, col_resolution = 2)
test <- isolate_territories(test)


e1 <- image_plot(vesalius, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Reference")

e2 <- image_plot(jitter_ves, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Query")
e3 <- image_plot(test, embedding = "PCA") +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Aligned")
pdf("Scenes/super_pixel/integrated_embeddings.pdf", width= 21, height = 7)
print(e1 + e2 + e3)
dev.off()

library(RColorBrewer)
cols <- c("#9e5652ff","#efb4a6","#abbeed","#111b34ff")
palette <- colorRampPalette(cols)


vesalius@territories$Territory <- gsub("isolated", 1, vesalius@territories$Territory)
colors <- palette(length(unique(vesalius@territories$Territory)))
ter_1 <- territory_plot(vesalius, cex_pt = 5, alpha =1) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Reference Territories")


colors <- palette(length(unique(jitter_ves@territories$Territory)))

ter_2 <- territory_plot(jitter_ves, cex_pt = 5,alpha =1) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Query Territories")

colors <- palette(length(unique(test@territories$Territory)))

ter_3 <- territory_plot(test, cex_pt = 5,alpha =1) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.position = "none") +
    labs(color = "", title = "Aligned Territories")

pdf("Scenes/super_pixel/integrated_territories.pdf", width= 21, height = 7)
print(ter_1 + ter_2 + ter_3)
dev.off()