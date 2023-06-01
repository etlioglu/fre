# TODO
# use the function sce_create_qc_lognorm for GSE123335 as well

# What does this script do?
# 1) Reads count and annotation data as a SingleCellExperiment object
# 2) Performs QC
# 3) Log normalizes
# Studies included: GSE153164, GSE123335

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

    # print("Mitochondrial genes found: ")
    # print(
    #     rownames(sce)[grep("^mt-", rownames(sce))]
    # )

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
annots <- read.table(
    "data/external/GSE153164-metaData_scDevSC.txt",
    header = TRUE,
    sep = "\t",
    quote = ""
)
# removing the extra row stating column types
annots <- annots[-1, ]

# not all cell names in the annotations end with "-1", however, counts have this suffix
annots$NAME[!(endsWith(annots$NAME, "-1"))] <- paste0(annots$NAME[!(endsWith(annots$NAME, "-1"))], "-1")


file_names <- c(
    "GSM5277843_E10_v1_filtered_feature_bc_matrix.h5",
    "GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5",
    "GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5",
    "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5",
    "GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5",
    "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5",
    "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
    "GSM5277844_E17_5_filtered_feature_bc_matrix.h5",
    "GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5",
    "GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5",
    "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
    "GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5",
    "GSM5277845_P4_filtered_feature_bc_matrix.h5"
)

file_paths <- paste0("data/external/", file_names)

sample_prefixes <- c(
    "E10_v1",
    "E11_5",
    "E12_5",
    "E13_5",
    "E14_5",
    "E15_5",
    "E16_5",
    "E17_5",
    "E18_5_S1",
    "E18_S3",
    "P1_S1",
    "P1_S2",
    "P4"
)


# works for GSE153164, needs to be tweaked for other data sets
read_and_sanitize_count_file <- function(file_path, prefix) {
    count_data <- Read10X_h5(file_path)

    # adding "sample prefixes" to count data
    colnames(count_data) <- paste0(prefix, "_", colnames(count_data))

    sce <- sce_create_qc_lognorm(count_data)

    # merge the count and annotation data
    common_cells <- intersect(
        colnames(sce),
        annots$NAME
    )
    sce <- sce[, common_cells]

    print("Dimension of the SCE object")
    print(prefix)
    print(dim(sce))

    annot <- annots[match(common_cells, annots$NAME), ]
    if (all(rownames(colData(sce)) == annot$NAME)) {
        print("Cell IDs OK!")
        colData(sce)$dibella_et_al_2021_annots <- annot$New_cellType
    } else {
        print("Cell IDs do not match!")
    }
    return(sce)
}

dilbella_et_al_2021_sce_objects <- mapply(read_and_sanitize_count_file, file_paths, sample_prefixes)
dilbella_et_al_2021 <- Reduce(cbind, dilbella_et_al_2021_sce_objects)
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
