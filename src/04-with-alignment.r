library(Seurat)
library(parallel)
library(clustree)
library(viridis)

seurat_objects <- readRDS("data/processed/output-of-seurat-path.rds")

# Normalization is required but in this case is performed in script "02-seurat-path.r"

features <- SelectIntegrationFeatures(
    object.list = seurat_objects,
    nfeatures = 2000,
    verbose = TRUE
)

seurat_objects <- PrepSCTIntegration(
    object.list = seurat_objects,
    anchor.features = features,
    verbose = TRUE
)


seurat_objects <- mclapply(
    X = seurat_objects,
    FUN = RunPCA,
    features = features,
    mc.cores = 32
)

anchors <- FindIntegrationAnchors(
    object.list = seurat_objects,
    normalization.method = "SCT",
    anchor.features = features,
    dims = 1:30,
    reduction = "rpca",
    k.anchor = 5,
    verbose = TRUE
)

fre_combined_sct <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    dims = 1:30,
    verbose = TRUE
)

fre_combined_sct <- RunPCA(fre_combined_sct, verbose = FALSE)
fre_combined_sct <- RunUMAP(fre_combined_sct, reduction = "pca", dims = 1:30)

fre_combined_sct <- FindNeighbors(
    fre_combined_sct,
    reduction = "pca",
    dims = 1:30
)

for (res in c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.40, 0.60, 0.80, 1.0, 1.2, 1.4)) {
    fre_combined_sct <- FindClusters(
        fre_combined_sct,
        resolution = res,
        verbose = TRUE
    )
}

### clustree image ###
png(
    "output/plots/with-alignment-clustree.png",
    width = 750,
    height = 750
)
clustree(fre_combined_sct, prefix = "integrated_snn_res.")
dev.off()

### set the "best" cluster resolution ###
fre_combined_sct@meta.data$seurat_clusters <- fre_combined_sct@meta.data$integrated_snn_res.0.1

### clusters by cluster id ###
png(
    "output/plots/with-alignment-clusters-by-cluster-id.png",
    width = 750,
    height = 750
)
DimPlot(
    fre_combined_sct,
    pt.size = 0.1,
    label = TRUE,
    shuffle = TRUE,
    cols = viridis(length(unique(fre_combined_sct@meta.data$seurat_clusters))),
    group.by = "seurat_clusters"
) + NoLegend()
dev.off()


### clusters by cell type ###
png(
    "output/plots/with-alignment-clusters-by-cell-type.png",
    width = 750,
    height = 750
)
DimPlot(
    fre_combined_sct,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    # cols = viridis(length(unique(fre@meta.data$seurat_clusters))),
    group.by = "SingleR_labels"
)
dev.off()

### clusters by sample ###
png(
    "output/plots/with-alignment-clusters-by-sample.png",
    width = 750,
    height = 750
)
DimPlot(
    fre_combined_sct,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    cols = viridis(length(unique(fre_combined_sct@meta.data$orig.ident))),
    group.by = "orig.ident"
)
dev.off()

### clusters by sample type ###
fre_combined_sct@meta.data$sample_type <- ifelse(
    fre_combined_sct@meta.data$orig.ident %in% c("RQU001", "RQU002"),
    "control",
    "treatment"
)
png(
    "output/plots/with-alignment-clusters-by-sample-type.png",
    width = 750,
    height = 750
)
DimPlot(
    fre_combined_sct,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    cols = viridis(length(unique(fre_combined_sct@meta.data$sample_type))),
    group.by = "sample_type"
)
dev.off()


saveRDS(fre_combined_sct, "data/processed/fre-aligned.rds")
