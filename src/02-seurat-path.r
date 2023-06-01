# What does this script do?
# 1) Creates Seurat objects
# 2) Other functionality explained within process_seurat_object()
#

# Output:
#

library(Seurat)
library(parallel)
library(scCustomize)
library(ggplot2)
library(viridis)

set.seed(42)

setwd("/data/high-level-analyses")

pal <- viridis(n = 10, option = "D")

# this will be used to transfer the doublet/multiplet predictions of SCE objects
# to the corresponding Seurat objects
sce_objects_coldata <- readRDS("data/processed/sce-objects-colData.rds") |>
    as.data.frame()

### READ DATA AND CREATE SEURAT OBJECT ###

sample_names <- paste0(
    "RQU00",
    1:4
)

paths <- paste0(
    "../nextflow-run/nfcore-scrnaseq-results/cellranger/count/sample-",
    sample_names,
    "_cDNA/outs/raw_feature_bc_matrix"
)
names(paths) <- sample_names

create_seurat_from_10x <- function(path, sample) {
    tenx_data <- Read10X(data.dir = path)
    seurat_object <- CreateSeuratObject(
        counts = tenx_data,
        project = sample,
        min.cells = 0,
        min.features = 250
    )
    seurat_object <- RenameCells(
        seurat_object,
        add.cell.id = sample
    )
}

seurat_objects <- mcmapply(
    create_seurat_from_10x,
    paths,
    names(paths),
    mc.cores = 32
)

process_seurat_object <- function(seurat_object, seurat_object_name) {
    ### PROCESS EACH SEURAT OBJECT INDIVIDUALLY ###

    # add QC metrics
    # add doublet/multiplet information from the corresponding SCE object
    # perform QC: max and min features per cell, mito, doublet, hemo
    # normalize
    # cluster

    # scCustomize for %mitochondrial and %ribosomal genes
    seurat_object <- Add_Mito_Ribo_Seurat(seurat_object, species = "Mouse")
    seurat_object <- Add_Cell_Complexity_Seurat(seurat_object)

    # hemoglobin genes
    seurat_object[["hemo_pct"]] <- PercentageFeatureSet(seurat_object, pattern = "^Hb[ab]-")
    # adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    max_feature_cutoff <- 4000
    min_feature_cutoff <- 500
    mito_cutoff <- 5
    hemo_pct_cutoff <- 5

    more_meta_data <- merge(
        seurat_object@meta.data,
        sce_objects_coldata,
        by = "row.names",
        all.x = TRUE
    )
    rownames(more_meta_data) <- more_meta_data[, "Row.names"]
    more_meta_data[, "Row.names"] <- NULL

    # extra filtering to align the SCE and Seurat objects:
    # NAs are introduced after the "merging" as Seurat and SCE does not use the
    # same set of cells
    common_cells <- more_meta_data[!is.na(more_meta_data$Sample), ] |>
        rownames()
    seurat_object <- subset(
        seurat_object,
        cells = common_cells
    )
    seurat_object <- AddMetaData(seurat_object, metadata = more_meta_data)

    seurat_object <- subset(
        seurat_object,
        subset =
            nFeature_RNA < max_feature_cutoff &
                nFeature_RNA > min_feature_cutoff &
                percent_mito < mito_cutoff &
                scDblFinder.class == "singlet" &
                hemo_pct < hemo_pct_cutoff
    )

    print(paste("Starting SCTransform for", seurat_object_name))

    seurat_object <- SCTransform(
        seurat_object,
        vst.flavor = "v2",
        verbose = FALSE,
        min_cells = 0
    ) |>
        RunPCA(npcs = 30, verbose = FALSE) |>
        RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) |>
        FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) |>
        FindClusters(resolution = 0.7, verbose = FALSE)

    is_numeric <- unlist(lapply(seurat_object@meta.data, is.numeric))
    numeric_cols <- colnames(seurat_object@meta.data)[is_numeric]
    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-feature-plots",
            ".png"
        ),
        width = 1000,
        height = 1000
    )
    FeaturePlot_scCustom(
        seurat_object,
        features = numeric_cols,
        colors_use = pal,
        # ncol = 3,
        pt.size = 0.1
    ) |>
        print()
    dev.off()

    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-doublet-class",
            ".png"
        ),
        width = 750,
        height = 750
    )
    DimPlot_scCustom(
        seurat_object,
        pt.size = 0.1,
        label = FALSE,
        shuffle = TRUE,
        group.by = "scDblFinder.class"
    ) |> print()
    dev.off()

    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-clusters",
            ".png"
        ),
        width = 750,
        height = 750
    )
    DimPlot_scCustom(
        seurat_object,
        pt.size = 0.1,
        label = FALSE,
        shuffle = TRUE,
        group.by = "seurat_clusters"
    ) |> print()
    dev.off()

    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-loo-et-al-2019-annots",
            ".png"
        ),
        width = 750,
        height = 750
    )
    DimPlot_scCustom(
        seurat_object,
        pt.size = 0.1,
        label = FALSE,
        shuffle = TRUE,
        group.by = "loo_et_al_2019_annots"
    ) |> print()
    dev.off()

    cell_counts <- as.data.frame(
        table(seurat_object@meta.data$dibella_et_al_2021_annots)
    )
    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-cell-counts",
            ".png"
        ),
        width = 750,
        height = 750
    )
    print(
        ggplot(
            cell_counts,
            aes(x = Var1, y = Freq, fill = Var1)
        ) +
            geom_bar(stat = "identity") +
            theme(
                axis.text.x = element_text(angle = 45)
            )
    )
    dev.off()

    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-dibella-et-al-2021-annots",
            ".png"
        ),
        width = 750,
        height = 750
    )
    DimPlot_scCustom(
        seurat_object,
        pt.size = 0.1,
        label = FALSE,
        shuffle = TRUE,
        group.by = "dibella_et_al_2021_annots"
    ) |> print()
    dev.off()

    png(
        paste0(
            "output/plots/",
            seurat_object_name,
            "-hist-hemo-pct",
            ".png"
        ),
        width = 750,
        height = 750
    )

    hist(seurat_object@meta.data$hemo_pct)

    dev.off()

    return(seurat_object)
}

processed_seurat_objects <- mcmapply(
    process_seurat_object,
    seurat_objects,
    names(seurat_objects),
    mc.cores = 32
)

print("Saving object list...")
saveRDS(processed_seurat_objects, "data/processed/output-of-seurat-path.rds")
