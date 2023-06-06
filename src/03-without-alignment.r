# What does this script do?
#

# Output:
#

library(Seurat)
library(scCustomize)
library(patchwork)
# library(viridis)
library(clustree)

set.seed(42)

setwd("/data/high-level-analyses")

# pal <- viridis(n = 10, option = "D")

seurat_objects <- readRDS("data/processed/output-of-seurat-path.rds")

fre <- Reduce(
    f = function(x, y) {
        merge(x, y, merge.data = FALSE)
    },
    x = seurat_objects
)

# remove the meta data columns associated with "sample-wise normalization/clustering"
# fre@meta.data <- fre@meta.data[, 1:18]
unwanted_cols <- c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.7", "seurat_clusters")
fre@meta.data <- fre@meta.data[, !names(fre@meta.data) %in% unwanted_cols]

fre@meta.data$sample_type <- ifelse(
    fre@meta.data$orig.ident %in% c("RQU001", "RQU002"),
    "control",
    "treatment"
)

fre <- SetIdent(fre, value = "orig.ident")

# sample_cols <- viridis(4) # this is fed into scCustomize's plotting functions

png(
    "output/plots/QC-post-filtering.png",
    width = 750,
    height = 750
)

# wrap_plots(p1, p2, p3, p4, ncol = 4)
QC_Plots_Combined_Vln(
    seurat_object = fre,
    feature_cutoffs = c(250, 4500),
    # UMI_cutoffs = c(1200,45000),
    mito_cutoffs = 5,
    pt.size = 0
    #    colors_use = sample_cols
)
dev.off()


fre <- SCTransform(
    fre,
    # vars.to.regress = "percent_ribo",
    vst.flavor = "v2",
    verbose = TRUE,
    conserve.memory = FALSE,
    return.only.var.genes = FALSE,
    min_cells = 0
)

fre <- RunPCA(fre, verbose = TRUE)
fre <- RunUMAP(fre, dims = 1:30, verbose = TRUE)
fre <- FindNeighbors(fre, dims = 1:30, verbose = TRUE)

for (res in c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.40, 0.60, 0.80, 1.0, 1.2, 1.4)) {
    fre <- FindClusters(
        fre,
        resolution = res,
        verbose = TRUE
    )
}

### clustree image ###
png(
    "output/plots/wo-alignment-clustree.png",
    width = 750,
    height = 750
)
clustree(fre, prefix = "SCT_snn_res.")
dev.off()

### set the "best" cluster resolution ###
fre@meta.data$seurat_clusters <- fre@meta.data$SCT_snn_res.0.4

### clusters by cluster id ###
png(
    "output/plots/wo-alignment-clusters-by-cluster-id.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(
    fre,
    pt.size = 0.1,
    label = TRUE,
    shuffle = TRUE,
    group.by = "seurat_clusters"
) + NoLegend()
dev.off()


### clusters by cell type: loo et al 2019 ###
png(
    "output/plots/wo-alignment-clusters-by-loo-et-al-2019.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(
    fre,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    group.by = "loo_et_al_2019_annots"
)
dev.off()


### clusters by cell type: dilbella et al 2021 ###
png(
    "output/plots/wo-alignment-clusters-by-dilbella-et-al-2021.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(
    fre,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    group.by = "dibella_et_al_2021_annots"
)
dev.off()


### clusters by cell phase
png(
    "output/plots/wo-alignment-clusters-by-cell-phase.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(
    fre,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    group.by = "cyclone_phase"
)
dev.off()







### clusters by sample ###
png(
    "output/plots/wo-alignment-clusters-by-sample.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(
    fre,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    group.by = "orig.ident"
)
dev.off()

### clusters by sample type ###
png(
    "output/plots/wo-alignment-clusters-by-sample-type.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(
    fre,
    pt.size = 0.1,
    label = FALSE,
    shuffle = TRUE,
    group.by = "sample_type"
)
dev.off()

### meta data variables ###
numeric_features <- colnames(fre@meta.data)[sapply(fre@meta.data, is.numeric)]

png(
    "output/plots/wo-alignment-clusters-by-meta-data-vars.png",
    width = 750,
    height = 750
)
FeaturePlot_scCustom(
    fre,
    features = numeric_features,
    colors_use = viridis_light_high
)
dev.off()

# saveRDS(fre, "data/processed/fre-no-alignment-percent-ribo-regressed-out.rds")
saveRDS(fre, "data/processed/fre-no-alignment.rds")

library(Seurat)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(tidyr)

pal <- polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")

summary_df <- fre@meta.data |>
    group_by(seurat_clusters) |> 
    summarise(
        num_cells = n(),
        mean_G1 = mean(G1),
        mean_S = mean(S),
        mean_G2M = mean(G2M),
        median_G1 = median(G1),
        median_S = median(S),
        median_G2M = median(G2M),
        .groups = 'drop') |>
  as.data.frame()

png(
    "output/plots/wo-alignment-cell-counts.png",
    width = 750,
    height = 750
)

ggplot(
    summary_df,
    aes(
        x = seurat_clusters,
        y = num_cells,
        fill = seurat_clusters)
    ) +
    geom_bar(stat="identity") +
    theme_bw() +
    scale_fill_manual(values = pal)
dev.off()  

summary <- pivot_longer(
    summary_df,
    cols = c("mean_G1", "mean_S", "mean_G2M", "median_G1", "median_S", "median_G2M")
)
summary$group <- ifelse(
    startsWith(summary$name, "mean_" ),
    "mean",
    "median"
)
png(
    "output/plots/wo-alignment-mean-cell-scores.png",
    width = 750,
    height = 750
)
ggplot(
    summary,
    aes(
        x = seurat_clusters,
        y = value,
        fill = name
    )
) +
    geom_bar(stat = "identity") +
    facet_wrap(~ group)

dev.off()

write.csv(
    summary_df,
    "cell-cycle-scores-per-cluster",
    row.names = FALSE
)
