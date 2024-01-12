# load vesalius data 
data(vesalius)
library(Seurat)
library(dplyr)


# Building object to test integration
vesalius <- vesalius:::process_counts(counts = counts,
    assay = "batch1")
features <- vesalius:::check_features(vesalius$SO)
vesalius <- vesalius:::embed_latent_space(vesalius$SO,
    assay = "batch1",
    "UMAP",
    dimensions = 30,
    features = features,
    verbose = FALSE)
vesalius <- as.data.frame(vesalius$UMAP)
vesalius$batch <- "batch_1"
jitter <- vesalius:::process_counts(counts = jitter_counts,
    assay = "batch2")
features <- vesalius:::check_features(jitter$SO)
jitter <- vesalius:::embed_latent_space(jitter$SO,
    assay = "batch2",
    "UMAP",
    dimensions = 30,
    features = features,
    verbose = FALSE)
jitter <- as.data.frame(jitter$UMAP)
jitter$batch <- "batch_2"

# Directly combining
combined_counts <- cbind(counts,jitter_counts)
combined <- vesalius:::process_counts(counts = combined_counts,
    assay = "batch2")
features <- vesalius:::check_features(combined$SO)
combined <- vesalius:::embed_latent_space(combined$SO,
    assay = "batch2",
    "UMAP",
    dimensions = 30,
    features = features,
    verbose = FALSE)
combined <- as.data.frame(combined$UMAP)
combined$batch <- "batch"
combined$batch[!grepl("q_",rownames(combined))] <- "batch_1"
combined$batch[grepl("q_",rownames(combined))] <- "batch_2"

# Vesalius Integration
vesalius <- build_vesalius_assay(coordinates, counts)
jitter_ves <- build_vesalius_assay(jitter_coord, jitter_counts)

vesalius <- generate_embeddings(vesalius,
    filter_threshold = 1,
    filter_grid = 1)

jitter_ves <- generate_embeddings(jitter_ves,
    filter_threshold = 1,
    filter_grid = 1)

adjusted_jitter <- vesalius:::intergrate_counts(vesalius,
    jitter_ves,
    "embeddings",
    use_norm = "log_norm")
combined_counts <- cbind(counts[rownames(counts) %in% rownames(adjusted_jitter),],adjusted_jitter)
combined <- vesalius:::process_counts(counts = combined_counts,
    assay = "batch2")
features <- vesalius:::check_features(combined$SO)
combined <- vesalius:::embed_latent_space(combined$SO,
    assay = "batch2",
    "UMAP",
    dimensions = 30,
    features = features,
    verbose = FALSE)
combined <- as.data.frame(combined$UMAP)
combined$batch <- "batch"
combined$batch[!grepl("q_",rownames(combined))] <- "batch_1"
combined$batch[grepl("q_",rownames(combined))] <- "batch_2"

# Seurat Integration
seurat_base <- Seurat::CreateSeuratObject(counts)
seurat_base <- Seurat::AddMetaData(seurat_base,"batch_1", col.name = "batch")
seurat_jitter <- Seurat::CreateSeuratObject(jitter_counts)
seurat_jitter <- Seurat::AddMetaData(seurat_jitter,"batch_2", col.name = "batch")
seurat_combined <- merge(seurat_base,seurat_jitter)
#seurat_combined[["RNA"]] <- split(seurat_combined,f = seurat_combined$batch)
seurat_combined <- Seurat::NormalizeData(seurat_combined) %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA()
seurat_combined <- IntegrateLayers(object = seurat_combined,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

seurat_combined[["RNA"]] <- JoinLayers(seurat_combined[["RNA"]])
seurat_combined <- RunUMAP(seurat_combined, dims = 1:30, reduction = "integrated.cca")
