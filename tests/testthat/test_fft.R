# load vesalius data 
# We will assume that the embeddings have been produced 
library(anndata)
spat_data <- read_h5ad("/common/wonklab/CosMX/PembroRT_cosmx_RUBY.h5ad")
sc_data <- read_h5ad("/common/wonklab/CosMX/")
data(vesalius)
vesalius <- build_vesalius_assay(coordinates, counts) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, sigma = 5, iter = 10) %>%
    segment_image(dimensions = 1:30, col_resolution = 12) %>%
    isolate_territories()
vesalius_query <- build_vesalius_assay(coordinates, counts) %>%
    generate_embeddings(tensor_resolution = 0.3) %>%
    regularise_image(dimensions = 1:30, lambda = 5) %>%
    equalize_image(dimensions = 1:30, sleft = 5, sright = 5) %>%
    smooth_image(dimensions = 1:30, sigma = 5, iter = 10) %>%
    segment_image(dimensions = 1:30, col_resolution = 12) %>%
    isolate_territories()

test <- integrate_territories(vesalius, vesalius_query, method = "coherence")

testx <- test$x
colnames(testx) <- paste("s_", colnames(testx))
rownames(testx) <- paste("q_", rownames(testx))
testy <- test$y
colnames(testy) <- paste("s_", colnames(testy))
rownames(testy) <- paste("q_", rownames(testy))
test_b <- testx + testy
pdf("sim.pdf", width = 10, height = 10)
test_b %>% 
  as.data.frame() %>%
  rownames_to_column("q_id") %>%
  pivot_longer(-c(q_id), names_to = "samples", values_to = "counts") %>%
  ggplot(aes(x=samples, y=q_id, fill=counts)) + 
  geom_tile() +
  scale_fill_viridis_c()
dev.off()
# pdf("path_seed.pdf")
# partials <- 50000
# tmp <- test[[1]]
# x <- lapply(tmp, function(x, partials) {
#     partials <- min(partials, length(x$x))
#     return(Re(fft(x$x[1:partials], inverse = TRUE))/ (length(x$x)))
# }, partials)
# y <- lapply(tmp, function(y, partials) {
#     partials <- min(partials, length(y$y))
#     return(Re(fft(y$y[1:partials], inverse = TRUE))/ (length(y$y)))
# }, partials)
# range_x <- c(min(sapply(x, min)), max(sapply(x, max)))
# range_y <- c(min(sapply(y, min)), max(sapply(y, max)))
# plot(0, type = "n", xlim = range_x, ylim = range_y)
# for (i in seq_along(tmp)) {
#     lines(x[[i]], y[[i]], pch = 19, col = rainbow(length(tmp))[i])
# }
# dev.off()

# cor <- test$cor
# cov <- test$cov

# pdf("cor.pdf")
# image(cor)
# image(cov)
# dev.off()