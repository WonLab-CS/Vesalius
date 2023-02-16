# load vesalius data 
# We will assume that the embeddings have been produced 
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

test <- integrate_territories(vesalius, vesalius_query)


pdf("path_seed.pdf")
partials <- 50000
tmp <- test[[1]]
x <- lapply(tmp, function(x, partials) {
    partials <- min(partials, length(x$x))
    return(Re(fft(x$x[1:partials], inverse = TRUE))/ (length(x$x)))
}, partials)
y <- lapply(tmp, function(y, partials) {
    partials <- min(partials, length(y$y))
    return(Re(fft(y$y[1:partials], inverse = TRUE))/ (length(y$y)))
}, partials)
range_x <- c(min(sapply(x, min)), max(sapply(x, max)))
range_y <- c(min(sapply(y, min)), max(sapply(y, max)))
plot(0, type = "n", xlim = range_x, ylim = range_y)
for (i in seq_along(tmp)) {
    lines(x[[i]], y[[i]], pch = 19, col = rainbow(length(tmp))[i])
}
dev.off()

cor <- test$cor
cov <- test$cov

pdf("cor.pdf")
image(cor)
image(cov)
dev.off()