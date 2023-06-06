library(Seurat)
library(scCustomize)
library(viridis)

set.seed(42)

setwd("/data/high-level-analyses")

pal <- viridis(n = 10, option = "D")

fre <- readRDS("data/processed/fre-no-alignment.rds")

Idents(fre) <- fre@meta.data$seurat_clusters
all_markers <- FindAllMarkers(object = fre)

write.table(
    all_markers,
    file = "data/processed/cluster-markers.csv",
    append = FALSE)

top_markers <- Extract_Top_Markers(
    marker_dataframe = all_markers,
    num_genes = 10,
    named_vector = FALSE,
    make_unique = TRUE
)

png(
    "output/plots/cluster-markers.png",
    width = 1500,
    height = 1500,
    res = 100
)
Clustered_DotPlot(
    seurat_object = fre,
    features = top_markers,
    colors_use_exp = viridis_light_high,
    k = 17
    )
dev.off()

