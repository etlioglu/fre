# TODO
# use the function for GSE123335 as well

# What does this script do?
# 1) Reads count and annotation data as a SingleCellExperiment object
# 2) Performs QC
# 3) Log normalizes

# Output:
# SCE objects with cell annotations

library(Seurat)
library(SingleCellExperiment)
library(scater)

set.seed(42)

setwd("/data/high-level-analyses")

### accessory functions ###

sce_create_qc_lognorm <- function(count_matrix) {

    sce <- SingleCellExperiment(
        list(counts = count_matrix)
    ) 
    
    sce <- addPerCellQC(
        sce,
        subsets = list(
            mito = grep("^mt-", rownames(sce))
        )
    )

    print("Mitochondrial genes found: ")
    print(
        rownames(sce)[grep("^mt-", rownames(sce))]
    )

    qc <- quickPerCellQC(
        colData(
            sce,
            percent_subsets = c("subsets_mito_percent")
        )
    )
    
    print(paste("Object pre-qc:", dim(sce)))
    sce <- sce[, which(!qc$discard)]
    print(paste("Object post-qc:", dim(sce)))

    sce <- logNormCounts(sce)

    return(sce)

}

### GSE153164 ###

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153164

# annotation file from https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
e10_e11_annot <- read.table(
    "data/external/GSE153164-metaData_scDevSC.txt",
    header = TRUE,
    sep = "\t",
    quote = ""
)

# removing the extra row stating column types
e10_e11_annot <- e10_e11_annot[-1,]

### E10 of GSE153164 ###

e10_count <- Read10X_h5(
    "data/external/GSM5277843_E10_v1_filtered_feature_bc_matrix.h5"
)
# count data do not have "sample prefixes" but plain cell IDs
colnames(e10_count) <- paste0("E10_v1_", colnames(e10_count))

e10_sce <- sce_create_qc_lognorm(e10_count)

# merge the count and annotation data
common_cells <- intersect(
    colnames(e10_sce),
    e10_e11_annot$NAME
)
e10_sce <- e10_sce[, common_cells]
e10_annot <- e10_e11_annot[match(common_cells, e10_e11_annot$NAME), ] 
if (all(rownames(colData(e10_sce)) == e10_annot$NAME)) {
    print("Cell IDs OK!")
    colData(e10_sce)$dibella_et_al_2021_annots <- e10_annot$New_cellType
} else {
    print("Cell IDs do not match!")
}

### E11 of GSE153164 ###

e11_count <- Read10X_h5(
    "data/external/GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5"
)
# count data do not have "sample prefixes" but plain cell IDs
colnames(e11_count) <- paste0("E11_5_", colnames(e11_count))
# count data has a suffix ("-1") that is not present in the meta data
colnames(e11_count) <- sub("-1", "", colnames(e11_count))

e11_sce <- sce_create_qc_lognorm(e11_count)

# merge the count and annotation data
common_cells <- intersect(
    colnames(e11_sce),
    e10_e11_annot$NAME
)
e11_sce <- e11_sce[, common_cells]
e11_annot <- e10_e11_annot[match(common_cells, e10_e11_annot$NAME), ] 
if (all(rownames(colData(e11_sce)) == e11_annot$NAME)) {
    print("Cell IDs OK!")
    colData(e11_sce)$dibella_et_al_2021_annots <- e11_annot$New_cellType
} else {
    print("Cell IDs do not match!")
}

dilbella_et_al_2021 <- cbind(e10_sce, e11_sce)
saveRDS(dilbella_et_al_2021, "data/processed/dilbella_et_al_2021.rds")

### GSE123335 ###
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123335

e14_count <- read.table(
    "data/external/GSE123335_E14_combined_matrix.txt.gz",
    header = TRUE,
    row.names = 1
)

e14_annot <- read.table(
    "data/external/GSE123335_E14_combined_matrix_ClusterAnnotations.txt.gz",
    header = TRUE,
    sep = "\t"
)

e14_count_mtx <- e14_count[, e14_annot$CellID] |>
    data.matrix()

if (all(colnames(e14_count_mtx) == e14_annot$CellID)) {
    e14 <- SingleCellExperiment(
        list(counts = e14_count_mtx),
        colData = DataFrame(label = e14_annot$Cluster)
    )
} else {
    print("Cell names do not match!")
}

e14 <- addPerCellQC(
    e14,
    subsets = list(
        mito = grep("^mt-", rownames(e14))
    )
)

qc <- quickPerCellQC(
    colData(e14),
    percent_subsets = c("subsets_mito_percent")
)

e14 <- e14[, which(!qc$discard)]

e14 <- logNormCounts(e14)

saveRDS(e14, "data/processed/loo_et_al_2019_e14.rds")
