#-----------------------------/Simulated Data/---------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to data simumlation in high res data sets
# We used slide-seq data as basis for data simulation.
#------------------------------------------------------------------------------#
library(dplyr)
library(splatter)



#------------------------------/ Functions /-----------------------------------#



simulate_spatial <- function(territories = 3,
    n_cell_types = 3,
    n_cells = 6000,
    n_genes = 5000,
    type = "pure",
    xlim = c(1, 3000),
    ylim = c(1, 3000),
    file = "Simulated") {
    #-------------------------------------------------------------------------#
    # First let's set n cells random position
    #-------------------------------------------------------------------------#
    simulated <- data.frame(
        "x" = runif(n = n_cells, min = min(xlim) ,max = max(xlim)),
        "y" = runif(n = n_cells, min = min(ylim), max = max(ylim)),
        "territory" = rep(0, n_cells),
        "cells" = rep(0, n_cells))
    #-------------------------------------------------------------------------#
    # now we deine how many cell types we actually need 
    #-------------------------------------------------------------------------#
    if (type != "pure") {
        n_cell_types <- n_cell_types * territories
    }
    #-------------------------------------------------------------------------#
    # we can assign territories 
    #-------------------------------------------------------------------------#
    simulated <- switch(type,
        "pure" = pure_dist(simulated, territories),
        "exp" = exp_dist(simulated, territories, n_cell_types),
        "uni" = uni_dist(simulated, territories, n_cell_types),
        "dot" = dot_dist(simulated, territories, n_cell_types))
    #-------------------------------------------------------------------------#
    # get relative cell type proportion 
    #-------------------------------------------------------------------------#
    relative_proportion <- as.numeric(table(simulated$cells) / nrow(simulated))
    params <- newSplatParams(batchCells = n_cells,
        nGenes = n_genes)
    counts <- splatSimulate(params,
        group.prob = relative_proportion,
        method = "groups",
        verbose = FALSE)
    
    
    simulated <- assign_cells_to_group(simulated, counts)
    #write.csv(simulated, file = file_coord)
    #write.table(counts, file = file_count)
}

assign_cells_to_group <- function(simulated, counts) {
    group_prob <- table(colData(counts)$Group) / nrow(colData(counts))
    cell_prob <- table(simulated$cells) / nrow(simulated)
    browser()
    locs <- names(group_prob)[sapply(group_prob, function(g, c) {
        return(which(abs(g - c) == min(abs(g - c))))
    }, c = cell_prob)]


}

pure_dist <- function(simulated, territories) {
    cut <- seq(1, max(simulated$x) + 1, length.out = territories + 1)
    for (i in seq(1, territories)) {
        simulated$territory[simulated$x >= cut[i] &
            simulated$x < cut[i + 1]] <- i
        simulated$cells[simulated$territory == i] <- i
    }
    return(simulated)
}

uni_dist <- function(simulated, territories, n_cell_types) {
    cut <- seq(1, max(simulated$x), length.out = territories + 1)
    cell_cut <- seq(1, n_cell_types, by = territories)
    for (i in seq(1, territories)) {
        simulated$territory[simulated$x >= cut[i] &
            simulated$x < cut[i + 1]] <- i
        cells <- sample(seq(cell_cut[i], l = n_cell_types),
            sum(simulated$territory == i),
            replace = TRUE,
            prob = rep(1 / n_cell_types, times = n_cell_types))
        simulated$cells[simulated$territory == i] <- cells
    }
    return(simulated)
}

exp_dist <- function(simulated, territories, n_cell_types) {
    cut <- seq(1, max(simulated$x), length.out = territories + 1)
    cell_cut <- seq(1, n_cell_types, by = territories)
    divs <- exp(1)^seq(1, n_cell_types)
    divs <- round(divs / sum(divs), digits = 2)
    mat <- matrix(0, ncol = territories, nrow = n_cell_types)
    for(i in seq_len(n_cell_types)){
        mat[i,] <- seq(divs[i],
            divs[1 + n_cell_types - i],
            length.out = territories)

    }
    for (i in seq(1, territories)) {
        simulated$territory[simulated$x >= cut[i] &
            simulated$x < cut[i + 1]] <- i
        cells <- sample(seq(cell_cut[i], l = n_cell_types),
            sum(simulated$territory == i),
            replace = TRUE,
            prob = mat[i,])
        simulated$cells[simulated$territory == i] <- cells
    }
    return(simulated)
}

dot_dist <- function(simulated,
    territories,
    n_cell_types,
    radius = c(100, 500),
    rare = FALSE) {
    territories <- territories - 1
    seed <- sample(seq(1, nrow(simulated)), territories, replace = FALSE)
    radius <- sample(seq(radius[1L],radius[2L]), territories, replace = FALSE)
    cell_cut <- seq(1, n_cell_types, by = territories)
    for (i in seq_along(seed)) {
         simulated$territory[
            sqrt(((simulated$x - simulated$x[seed[i]])^2 +
            (simulated$y - simulated$y[seed[i]])^2)) <= radius[i]] <- i
        cells <- sample(seq(cell_cut[i], l = n_cell_types),
            sum(simulated$territory == i),
            replace = TRUE,
            prob = rep(1 / n_cell_types, times = n_cell_types))
        simulated$cells[simulated$territory == i] <- cells
    }
    # replacing back groud
    cells <- sample(seq(tail(cell_cut, 1) + 1, l = n_cell_types),
            sum(simulated$territory == 0),
            replace = TRUE,
            prob = rep(1 / n_cell_types, times = n_cell_types))
    simulated$cells[simulated$territory == i] <- cells
    if (rare) {
        cells_by_territory <- split(simulated, simulated$territory)
        cells_by_territory <- names(cells_by_territory)[
            which(sapply(cells_by_territory, nrow) ==
            min(sapply(cells_by_territory, nrow)))
        ]
        simulated$cells[simulated$territory == cells_by_territory] <- cells
    }
    return(simulated)
}




# PARAMS
n_runs <- 10
n_territories <- c(3,3,3,3,3,3,6,3)
n_cells <- c(3,3,4,4,5,5,3,3)
method <- c("uni","exp","uni","exp","uni","exp","dot","pure")

