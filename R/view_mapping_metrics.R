
#' visualize mapping metrics for vesalius assay
#' @param vesalius_assay vesalius_assay object
#' @param trial character string defining which metric trial to visualize
#' @param cex_pt numeric - point size for plotting
#' @param cex numeric - text size for plotting
#' @param randomize logical - whether to randomize point order for better visualization
#' @return ggplot object showing mapping metrics
#' @export
#' @importFrom tools toTitleCase
view_mapping_metrics <- function(vesalius_assay,
    trial = "cost",
    cex_pt = 1,
    cex = 15,
    randomize = TRUE) {
    coord <- get_coordinates(vesalius_assay)
    trial <- check_metric_trial(vesalius_assay, trial)
    trial_name <- tail(colnames(trial))
    colnames(trial)[which(colnames(trial) == "from")] <- "barcodes"
    trial <- right_join(coord, trial, by = "barcodes")
    if (grepl("Map_cluster", tail(colnames(trial),1))) {
        colnames(trial)[ncol(trial)] <- "trial"
        trial$trial <- as.factor(trial$trial)
        colors <- create_palette(trial, randomize)
        g <- ggplot(trial, aes(x,y, col = trial)) +
            geom_point(size = cex_pt) +
            scale_color_manual(values = colors) +
            theme_classic() +
            theme(legend.title = element_text(size = cex),
                legend.text = element_text(size = cex),
                plot.margin = margin(1, 1, 1, 1, "cm")) +
            labs(col = "Mapping Clusters", title = tools::toTitleCase(trial_name)) +
            guides(color = guide_legend(override.aes = list(size = cex * 0.5)))
        return(g)
    } else {
        colnames(trial)[ncol(trial)] <- "trial"
        trial$trial <- as.numeric(trial$trial)
        g <- ggplot(trial, aes(x,y, col = trial)) +
            geom_point(size = cex_pt) +
            scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
            theme_classic() +
            theme(legend.title = element_text(size = cex),
                legend.text = element_text(size = cex),
                plot.margin = margin(1, 1, 1, 1, "cm")) +
            labs(col = "Mapping score", title = tools::toTitleCase(trial_name)) +
            guides(shape = guide_legend(override.aes = list(size = cex * 0.5)))
        return(g)
    }
}
