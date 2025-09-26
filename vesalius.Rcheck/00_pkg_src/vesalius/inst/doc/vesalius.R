## ----loading, eval = TRUE, echo = TRUE----------------------------------------
library(vesalius)
library(ggplot2)
library(patchwork)
library(ggpubr)
data(vesalius, package = "vesalius")


## ----str_counts, eval = TRUE, echo = TRUE-------------------------------------
str(counts)

## ----str_coordinates, eval = TRUE, echo = TRUE--------------------------------
str(coordinates)

## ----build_vesalius_arrary, eval = TRUE, echo = TRUE--------------------------
vesalius <- build_vesalius_assay(
  coordinates = coordinates, # spatial coordinates
  counts  = counts, # count matrix
  assay = "spatial_omics", # name you wish to give your assay
  verbose = FALSE # Do you want progress messages?
)


## ----show_ves, eval = TRUE, echo = TRUE---------------------------------------
vesalius

## ----build_embeddings, eval = TRUE, echo = TRUE-------------------------------
vesalius <- generate_embeddings(vesalius,
  dim_reduction = "PCA",
  normalization = "log_norm",
  nfeatures = 100, # Setting number of features low for low run time
  verbose = FALSE)


## ----embedding_view, eval = TRUE, echo = FALSE--------------------------------
vesalius

## ----adding_embeddings, eval = TRUE, echo = TRUE------------------------------
vesalius <- generate_embeddings(vesalius,
  dim_reduction = "UMAP",
  nfeatures = 100, # Setting number of features low for low run time
  verbose = FALSE)

vesalius <- generate_embeddings(vesalius,
  dim_reduction = "PCA",
  nfeatures = 200, # Setting number of features low for low run time
  verbose = FALSE)

vesalius <- generate_embeddings(vesalius,
  dim_reduction = "PCA",
  verbose = FALSE)

## ----checkick_embed_names, eval = TRUE, echo = TRUE---------------------------
vesalius

## ----grey_scale_view, eval = TRUE, echo = TRUE--------------------------------
p1 <- image_plot(vesalius, dimensions = 1) + labs(title = "Grey PCA dim 1")
p2 <- image_plot(vesalius, dimensions = seq(1, 3)) + labs(title = "RGB PCA")

## ----grey_scale_view_PCA, eval = TRUE, echo = TRUE----------------------------
p3 <- image_plot(vesalius, dimensions = 1, embedding = "UMAP") +
  labs(title = "Grey UMAP dim 1")
p4 <- image_plot(vesalius, dimensions = c(1, 2, 3), embedding = "UMAP") +
  labs(title = "RGB UMAP dim 1, 2, and 3")

## ----out_plot_grey, eval = TRUE, echo = TRUE, fig.width = 10, fig.height = 10----
(p1 + p2) / (p3 + p4)

## ----acive_embeds, eval = TRUE, echo = TRUE-----------------------------------
get_active_embedding_tag(vesalius)

## ----image_processing, eval = TRUE, echo = TRUE-------------------------------
vesalius <- regularise_image(vesalius, lambda = 1)
vesalius <- smooth_image(vesalius, sigma = 5, iter = 10)
vesalius <- equalize_image(vesalius, sleft = 5, sright = 5)

## ----selecting_PCS, eval = FALSE, echo = TRUE---------------------------------
# # Selecting a subset of PCs
# dims <- c(1, 3, 4, 5, 7:11)
# 
# # running smoothing o
# vesalius <- regularise_image(vesalius,
#   dimensions = dims,
#   embedding = "PCA",
#   verbose = FALSE)
# vesalius <- equalize_image(vesalius,
#   dimensions = dims,
#   verbose = FALSE)
# vesalius <- smooth_image(vesalius,
#   dimensions = dims,
#   iter = 10,
#   sigma = 1,
#   verbose = FALSE)
# 

## ----segment_image, eval = TRUE, echo = TRUE----------------------------------
vesalius <- segment_image(vesalius,
  method = "kmeans",
  col_resolution = 2,
  verbose = FALSE)

## ----seg_check, eval = TRUE, echo = TRUE--------------------------------------
vesalius

## ----plot_segments_only, eval = TRUE, echo = TRUE, fig.width = 10, fig.height = 8----
p5 <- image_plot(vesalius) + labs(title = "Segments only")
print(p5)

## ----iso_territories, eval = TRUE, echo = TRUE--------------------------------
vesalius <- isolate_territories(vesalius, capture_radius = 0.05)

## ----vesalius_with_ter, eval = TRUE, echo = TRUE------------------------------
vesalius

## ----plot_territories, eval = TRUE, echo = TRUE, fig.width = 10, fig.height = 8----
p6 <- territory_plot(vesalius, cex_pt = 3.5)
p6

## ----compare_territories, eval = TRUE, echo = TRUE----------------------------
vesalius <- identify_markers(vesalius, seed = 1, query = 2)
deg <- get_markers(vesalius)

## ----view_DEG, eval = TRUE, echo = TRUE, fig.width = 10,fig.height = 4--------

p7 <- view_gene_expression(vesalius, genes = "Malat1")
p8 <- view_gene_expression(vesalius, genes = "Malat1", as_layer = TRUE)
p7 + p8

## ----mapping_build, eval = TRUE, echo = TRUE----------------------------------
# load vesalius data 
data(vesalius)

# Create Vesalius object for processing
vesalius <- build_vesalius_assay(coordinates, counts)
cells <- sample(LETTERS[1:6], size = nrow(vesalius@tiles),replace =T)
names(cells) <- vesalius@tiles$barcodes
vesalius <- add_cells(vesalius, cells = cells, add_name = "Cells")
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)
cells <- sample(LETTERS[1:6], size = nrow(jitter_ves@tiles),replace =T)
names(cells) <- jitter_ves@tiles$barcodes
jitter_ves <- add_cells(jitter_ves, cells = cells, add_name = "Cells")

vesalius <- generate_embeddings(vesalius,
    filter_threshold = 1,
    filter_grid = 1)
vesalius <- smooth_image(vesalius, embedding = "PCA", sigma = 5, iter = 10)
vesalius <- equalize_image(vesalius, sleft = 5, sright = 5)
vesalius <- segment_image(vesalius, col_resolution = 2)
vesalius <- isolate_territories(vesalius)


jitter_ves <- generate_embeddings(jitter_ves,
    filter_threshold = 1,
    filter_grid = 1)
jitter_ves <- smooth_image(jitter_ves, embedding = "PCA", sigma = 5, iter = 10)
jitter_ves <- equalize_image(jitter_ves, sleft = 5, sright = 5)
jitter_ves <- segment_image(jitter_ves, col_resolution = 2)
jitter_ves <- isolate_territories(jitter_ves)




## ----mapping, eval = TRUE, echo = TRUE----------------------------------------
matched <- map_assays(vesalius,
    jitter_ves,
    threshold = 0,
    use_cost = c("feature","niche","territory","composition"),
    batch_size = 500,
    epoch = 5,
    jitter = 1)

## ----map_embeddings, eval = TRUE, echo = TRUE---------------------------------
matched <- generate_embeddings(matched,
    dim_reduction = "PCA",
  nfeatures = 200, 
  verbose = FALSE)
matched <- smooth_image(matched, sigma = 5, iter = 10)
matched <- equalize_image(matched, sleft = 5, sright = 5)
matched <- segment_image(matched,col_resolution = 2)

## ----map_image_plot, eval = TRUE, echo = TRUE,fig.width = 12------------------
m <- image_plot(vesalius)
m1 <- image_plot(jitter_ves)
m2 <- image_plot(matched)

m + m1 + m2


## ----map_territory, eval = TRUE, echo = TRUE,fig.width = 12-------------------
t <- territory_plot(vesalius)
t1 <- territory_plot(jitter_ves)
t2 <- territory_plot(matched)
t + t1 + t2

## ----map_inter, eval = TRUE, echo = TRUE--------------------------------------
inter <- integrate_assays(matched,
    vesalius)



## ----inter_territory, eval = TRUE, echo = TRUE--------------------------------
inter <- smooth_image(inter, sigma = 5, iter = 10)
inter <- equalize_image(inter, sleft = 5, sright = 5)
inter <- segment_image(inter, col_resolution = 2)

## ----inter_plot, eval = TRUE, echo = TRUE, fig.width = 9----------------------
i <- image_plot(inter)
i1 <- territory_plot(inter)

i + i1

## ----inter_DEG, eval = TRUE, echo  = TRUE-------------------------------------
inter <- identify_markers(inter, sample = TRUE)
degs <- get_markers(inter)
head(degs)

## ----session, eval = TRUE, echo = TRUE----------------------------------------
sessionInfo()

