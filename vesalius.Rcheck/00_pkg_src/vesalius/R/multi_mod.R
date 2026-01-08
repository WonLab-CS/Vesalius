

#-----------------------------------------------------------------------------#
############################### Joint Measure #################################
#-----------------------------------------------------------------------------#

#' Jointly measured spatial omic assays territories
#' @param mod1 vesalius_assay object containing first modality
#' @param mod2 vesalius_assay objecty containing second modality
#' @param dimensions numeric vector describing latent space dimensions 
#' to use during intergration
#' @param method character - integration method. interlace - mean - concat 
#' are available options
#' @param norm_method character - which count values should be use 
#' for integration when using concat method
#' @param dim_reduction characater - which dim reduction methods should be 
#' used for concat integration (PCA,PCA_L,UMAP,LSI,LSI_UMAP,NMF)
#' @param verbose logical - should progress message be outputed to the 
#' console.
#' @return vesalius object containing new image embeddings
#' @export
joint_territories <- function(mod1,
    mod2,
    dimensions = seq(1, 30),
    embedding = "last",
    method = "interlace",
    norm_method = "log_norm",
    dim_reduction = "PCA",
    signal = "variable_features",
    verbose = TRUE) {
    simple_bar(verbose)
    #-------------------------------------------------------------------------#
    # check place holder 
    #-------------------------------------------------------------------------#

    #-------------------------------------------------------------------------#
    # Get embeddings - need make some changes here
    # for now we assume we get the active 
    #-------------------------------------------------------------------------#
    mod1_embed <- check_embedding_selection(mod1, embedding, dimensions)
    mod2_embed <- check_embedding_selection(mod2, embedding, dimensions)
    #-------------------------------------------------------------------------#
    # Get counts from each vesalius object
    #-------------------------------------------------------------------------#
    mod1_counts <- get_counts(mod1, type = "all")
    names(mod1_counts) <- paste0("mod1_", names(mod1_counts))
    mod2_counts <- get_counts(mod2, type = "all")
    names(mod2_counts) <- paste0("mod2_", names(mod2_counts))
    #-------------------------------------------------------------------------#
    # method switch - which method is best 
    #-------------------------------------------------------------------------#
    integrated_embeds <- switch(EXPR = method,
        "interlace" = interlace_embeds(mod1_embed, mod2_embed, dimensions),
        "mean" = average_embed(mod1_embed, mod2_embed, dimensions),
        "concat" = concat_embed(mod1,
            mod2,
            dimensions,
            norm_method,
            dim_reduction,
            signal = signal))
    integrated_embeds <- list(integrated_embeds)
    names(integrated_embeds) <- method
    #-------------------------------------------------------------------------#
    # get tile overlap - it is possible that there is not a perfect overlap
    # so we will filter
    #-------------------------------------------------------------------------#
    tiles <- mod1@tiles
    tiles <- tiles[tiles$barcodes %in% rownames(integrated_embeds[[1]]), ]
    integrated <- new("vesalius_assay",
        assay = "integrated",
        embeddings = integrated_embeds,
        active = integrated_embeds[[1]],
        tiles = tiles)
    integrated_counts <- c(mod1_counts, mod2_counts)
    integrated <- update_vesalius_assay(vesalius_assay = integrated,
      data = integrated_counts,
      slot = "counts",
      append = FALSE)
    #--------------------------------------------------------------------------#
    # we can update the comment on the count slot list 
    # this comment will indicate which count matrix is set as default
    #--------------------------------------------------------------------------#
    integrated <- add_active_count_tag(integrated, norm = "joint")
    #--------------------------------------------------------------------------#
    # Finally we update the vesalius commit log
    #--------------------------------------------------------------------------#
    commit <- create_commit_log(arg_match = as.list(match.call()),
      default = formals(joint_territories))
    integrated <- commit_log(integrated,
      commit,
      "integrated")
    simple_bar(verbose)
    return(integrated)

}

#' interlace image stack between seed and query
#' @param seed matrix - seed embedding image stack
#' @param query matrix - query embedding image stack
#' @param dimensions int vector describing which embeddings
#' should be selected
#' @details Takes selected embedding from seed and query and 
#' creates and interlaced embedding matrix starting with the 
#' first seed embedding
#' @return embedding matrix containing seed embeddings + query 
#' embeddings.
interlace_embeds <- function(seed, query, dimensions) {
    locs <- intersect(rownames(seed), rownames(query))
    seed <- seed[rownames(seed) %in% locs,
        dimensions]
    query <- query[rownames(query) %in% locs,
        dimensions]
    seed <- seed[order(rownames(seed)), ]
    query <- query[order(rownames(query)), ]

    interlaced_embed <- matrix(0,
        ncol = ncol(seed) + ncol(query),
        nrow = nrow(seed))
    rownames(interlaced_embed) <- rownames(seed)
    dimensions <- rep(dimensions, each = 2)
    for (i in seq(1, ncol(interlaced_embed), by = 2)) {
        interlaced_embed[, i] <- seed[, dimensions[i]]
        interlaced_embed[, i + 1] <- query[, dimensions[i + 1]]
    }
    return(interlaced_embed)
}

#' average image stack between seed and query
#' @param seed matrix - seed embedding image stack
#' @param query matrix - query embedding image stack
#' @param dimensions int vector describing which embeddings
#' should be selected
#' @details Takes select embedding from seed and query and 
#' creates and avarage the grey scale pixel values for each spatial
#' location
#' @return embedding matrix containing average pixel value for both seed
#' query
average_embed <- function(seed, query, dimensions) {
    locs <- intersect(rownames(seed), rownames(query))
    seed <- seed[rownames(seed) %in% locs,
        dimensions]
    query <- query[rownames(query) %in% locs,
        dimensions]
    seed <- seed[order(rownames(seed)), ]
    query <- query[order(rownames(query)), ]
    averaged_embed <- matrix(0,
        ncol = length(dimensions),
        nrow = nrow(seed))
    rownames(averaged_embed) <- rownames(seed)
    for (i in seq(1, ncol(averaged_embed))) {
        averaged_embed[, i] <- apply(cbind(seed[, dimensions[i]],
            query[, dimensions[i]]),
            MARGIN = 1,
            mean)
    }
    return(averaged_embed)
}


#' create new embedding from jointly measure spatial omics 
#' @param seed vesalius_assay object of the first modality 
#' @param query vesalius_assay object of the second modality
#' @param dimensions int - number of gray scale images to create
#' @param norm_method string describing which normalisation 
#' method to use. One of the following "log_norm", "SCT", "TFIDF", "raw"
#' @param dim_reduction string describing which dimensionality
#' reduction method should be used. One of the following:
#' "PCA", "PCA_L", "UMAP", "LSI", "LSI_UMAP"
#' @details Create new latent space using feature from both modalities
#' Creates a new feature matrix than in nromalized and converted to latent
#' space image stack.
#' @return embedding matrix used as grey scale image stack
concat_embed <- function(seed,
    query,
    dimensions,
    norm_method,
    dim_reduction,
    signal) {
    #-------------------------------------------------------------------------#
    # first we get signal - i.e. these features that should be used
    #-------------------------------------------------------------------------#
    seed_features <- check_signal(seed, signal, type = norm_method)
    query_features <- check_signal(query, signal, type = norm_method)
    #-------------------------------------------------------------------------#
    # next we extract counts and make sure that they are int he right order
    # before rbidning the whole thing 
    #-------------------------------------------------------------------------#
    seed_counts <- get_counts(seed, type = "raw")
    query_counts <- get_counts(query, type = "raw")
    locs <- intersect(colnames(seed_counts), colnames(query_counts))
    seed_counts <- seed_counts[rownames(seed_counts) %in% seed_features,
        colnames(seed_counts) %in% locs]
    query_counts <- query_counts[rownames(query_counts) %in% query_features,
        colnames(query_counts) %in% locs]
    seed_counts <- seed_counts[, order(colnames(seed_counts))]
    query_counts <- query_counts[, order(colnames(query_counts))]
    integrated_counts <- rbind(seed_counts, query_counts)
    #-------------------------------------------------------------------------#
    # Now we can process the counts and create a new embedding
    #-------------------------------------------------------------------------#
    integrated_counts <- process_counts(integrated_counts,
        assay = "integrated",
        method = norm_method,
        use_count = "raw",
        nfeatures = sum(c(length(seed_features), length(query_features))))
    integrated_embeds <- embed_latent_space(integrated_counts$SO,
        assay = "integrated",
        dim_reduction = dim_reduction,
        dimensions = max(dimensions),
        remove_lsi_1 = FALSE,
        verbose = FALSE)
    
    return(integrated_embeds[[1]])
}

