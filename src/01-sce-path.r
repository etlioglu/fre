# What does this script do?
# 1) read 10x data as SCE objects
# 2) remove empty droplets
# 3) remove cells with low coverage
# 4) find doublets/multiplets
# 5) annotate cells with SingleR with Loo et al., 2019 as reference

# Output:
# 1) individual SingleCellExperiment objects
# 2) Aggregated "colData" to be transferred to "the Seurat path"

library(DropletUtils) # read 10x data
library(parallel)
library(scuttle) # for QC metric calculations
library(scDblFinder) # doublet detection
library(scran) # cell cycle score
library(SingleR)

set.seed(42)

setwd("/data/high-level-analyses")

### sample paths ###

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

### read data ###

sce_objects <- mcmapply(
    read10xCounts,
    samples = paths,
    sample.names = sample_names,
    col.names = TRUE,
    mc.cores = 32
)

### remove empty droplets ###

remove_empty_droplets <- function(sce_object) {
    e_out <- emptyDrops(assay(sce_object, "counts"))
    sce_object <- sce_object[, which(e_out$FDR <= 0.001)]
    return(sce_object)
}

sce_objects <- mclapply(
    X = sce_objects,
    FUN = remove_empty_droplets,
    mc.cores = 32
)

### add per cell qc metrics, to be used to filter out cells with low coverage ###

sce_objects <- mclapply(
    X = sce_objects,
    FUN = addPerCellQCMetrics,
    mc.cores = 32
)

remove_low_cov_cells <- function(sce) {
    sce <- sce[, sce$sum >= 200]
    return(sce)
}

sce_objects <- mclapply(
    X = sce_objects,
    FUN = remove_low_cov_cells,
    mc.cores = 32
)

### doublet/multiplet detection ###

sce_objects <- mclapply(
    X = sce_objects,
    FUN = scDblFinder,
    mc.cores = 32
)

### cell cycle ###

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
    package="scran"))

add_cell_cycle_score <- function(sce, mm.pairs) {
    assignments <- cyclone(
        sce,
        mm.pairs,
        gene.names = rownames(rowData(sce))
    )

    colData(sce)$cyclone_phase <- assignments$phases
    colData(sce) <- cbind(colData(sce), assignments$scores)

    return(sce)
    
}

sce_objects <- mclapply(
    X = sce_objects,
    FUN = add_cell_cycle_score,
    mm.pairs = mm.pairs,
    mc.cores = 32
)


### annotation ###

e10_e11 <- readRDS("data/processed/dilbella_et_al_2021.rds")
e14 <- readRDS("data/processed/loo_et_al_2019_e14.rds")

annotate_sce_wit_ref <- function(sce_object, ref, label_col_ref, label_col_target) {
    rownames(sce_object) <- rowData(sce_object)$Symbol

    annot_df <- SingleR(
        test = sce_object,
        ref = ref,
        labels = ref[[label_col_ref]],
        de.method = "wilcox",
        assay.type.test = 1
    )

    if (all(colnames(sce_object) == rownames(annot_df))) {
        colData(sce_object)[[label_col_target]] <- annot_df$pruned.labels
    }

    return(sce_object)
}

sce_objects <- mclapply(
    X = sce_objects,
    FUN = annotate_sce_wit_ref,
    ref = e14,
    label_col_ref = "label",
    label_col_target = "loo_et_al_2019_annots",
    mc.cores = 32
)

sce_objects <- mclapply(
    X = sce_objects,
    FUN = annotate_sce_wit_ref,
    label_col_ref = "dibella_et_al_2021_annots",
    label_col_target = "dibella_et_al_2021_annots",
    ref = e10_e11,
    mc.cores = 32
)

saveRDS(sce_objects, "data/processed/output-of-sce-path.rds")

fetch_colData <- function(sce_object) {
    sce_colData <- colData(sce_object)
    rownames(sce_colData) <- paste0(
        sce_colData$Sample,
        "_",
        sce_colData$Barcode
    )
    return(sce_colData)
}

sce_objects_colData <- mclapply(
    sce_objects,
    fetch_colData,
    mc.cores = 32
)

sce_objects_colData <- Reduce(rbind, sce_objects_colData)

saveRDS(sce_objects_colData, "data/processed/sce-objects-colData.rds")
