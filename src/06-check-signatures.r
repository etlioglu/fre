library(Seurat)
library(AUCell)
library(parallel)
library(geneSynonym)
library(ggplot2)

set.seed(42)

setwd("/data/high-level-analyses")

# fre <- readRDS("data/processed/fre-aligned.rds")
# fre <- readRDS("data/processed/fre-no-alignment.rds")
fre <- readRDS("data/processed/fre-no-alignment-percent-ribo-regressed-out.rds")

exp_mtx <- fre@assays$SCT@counts
seurat_genes <- rownames(exp_mtx)

# signatures
ir_uniq <- read.csv("other/IR_unique_RNAseq_Mm.csv", header = FALSE)
overlap <- read.csv("other/Overlap_RNAseq_Mm.csv", header = FALSE)
apoptosis <- read.csv("other/Apoptosis_Mm.csv", header = FALSE)
differentiation <- read.csv("other/Differentiation_Mm.csv", header = FALSE)

check_genes <- function(gene) {
    if (gene %in% seurat_genes) {
        return(gene)
    } else {
        syns <- mouseSyno(gene) |> unlist()
        if (length(syns) == 1) {
            print(
                paste(gene, "could  not be found...")
            )
            return(gene)
        } else {
            for (syn in syns) {
                if (syn %in% seurat_genes) {
                    return(syn)
                } else {
                    print(
                        paste(gene, "could  not be found...")
                    )
                    return(gene)
                }
            }
        }
    }
}

ir_uniq$seurat_genes <- apply(ir_uniq, 1, check_genes)
ir_uniq$seurat_genes[ir_uniq$seurat_genes == "AW555464"] <- "Cep170b"
ir_uniq$seurat_genes[ir_uniq$seurat_genes == "BC037703"] <- "Mab21l3"
ir_uniq$seurat_genes[ir_uniq$seurat_genes == "C77370"] <- "Nexmif"

overlap$seurat_genes <- apply(overlap, 1, check_genes)

apoptosis$seurat_genes <- apply(apoptosis, 1, check_genes)

differentiation$seurat_genes <- apply(differentiation, 1, check_genes)


# gs_ir_uniq <- list(ir_uniq = ir_uniq$seurat_genes)
# gs_overlap <- list(overlap = overlap$seurat_genes)

# gene_sets <- c(
#     GeneSet(ir_uniq$seurat_genes, setName="ir_uniq"),
#     GeneSet(overlap$seurat_genes, setName = "overlap")
# )

gene_sets <- list(
    ir_uniq = ir_uniq$seurat_genes,
    overlap = overlap$seurat_genes
)

cells_AUC <- AUCell_run(
    exp_mtx,
    gene_sets,
    BPPARAM=BiocParallel::MulticoreParam(32)
)

png("aucell-test.png")
cells_assignment <- AUCell_exploreThresholds(
    cells_AUC,
    plotHist=TRUE,
    assign=TRUE
) |> print()
dev.off()

ir_uniq_cells <- cells_assignment$ir_uniq$assignment
overlap_cells <- cells_assignment$overlap$assignment

fre@meta.data$ir_uniq <- ifelse(
    rownames(fre@meta.data) %in% ir_uniq_cells,
    "SIGNAL IR_UNIQ",
    "no_ir_uniq"
)

png(
    "output/plots/aucell-fre-not-aligned-dimplot-ir-uniq.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(fre,
    group.by = "ir_uniq",
    pt.size = 0.1,
    shuffle = TRUE
)
dev.off()

fre@meta.data$overlap <- ifelse(
    rownames(fre@meta.data) %in% overlap_cells,
    "SIGNAL OVERLAP",
    "no_overlap"
)

png(
    "output/plots/aucell-fre-not-aligned-dimplot-overlap.png",
    width = 750,
    height = 750
)
DimPlot_scCustom(fre,
    group.by = "overlap",
    pt.size = 0.1,
    shuffle = TRUE
)
dev.off()

# gene_sets <- subsetGeneSets(
#     gene_sets,
#     rownames(exp_mtx)
# )

# cbind(nGenes(gene_sets))

# check_signatures <- function(seurat_object, seurat_object_name) {
#     counts <- GetAssayData(
#         seurat_object,
#         slot = "counts"
#     )
#     cell_rankings <- AUCell_buildRankings(counts)

#     ### IR UNIQUE ###

#     cells_AUC_ir_uniq <- AUCell_calcAUC(
#         geneSets = gs_ir_uniq,
#         rankings = cell_rankings
#     )

#     png(
#         paste0(
#             "output/plots/aucell-",
#             seurat_object_name,
#             "-hist-ir-uniq",
#             ".png"
#         ),
#         width = 750,
#         height = 750
#     )
#     cells_assignment_ir_uniq <- AUCell_exploreThresholds(
#         cells_AUC_ir_uniq,
#         plotHist = TRUE,
#         assign = TRUE
#     ) |> print()
#     dev.off()

#     seurat_object$ir_uniq <- ifelse(
#         colnames(seurat_object) %in% cells_assignment_ir_uniq$ir_uniq$assignment,
#         "IR-UNIQUE",
#         "not-ir-unique"
#     )

#     png(
#         paste0(
#             "output/plots/aucell-",
#             seurat_object_name,
#             "-dimplot-ir-uniq",
#             ".png"
#         ),
#         width = 750,
#         height = 750
#     )
#     DimPlot(seurat_object, group.by = "ir_uniq", pt.size = 0.1, shuffle = TRUE) |> print()
#     dev.off()

#     ### OVERLAP ###

#     cells_AUC_overlap <- AUCell_calcAUC(
#         geneSets = gs_overlap,
#         rankings = cell_rankings
#     )

#     png(
#         paste0(
#             "output/plots/aucell-",
#             seurat_object_name,
#             "-hist-overlap",
#             ".png"
#         ),
#         width = 750,
#         height = 750
#     )
#     cells_assignment_overlap <- AUCell_exploreThresholds(
#         cells_AUC_overlap,
#         plotHist = TRUE,
#         assign = TRUE
#     ) |> print()
#     dev.off()

#     seurat_object$overlap <- ifelse(
#         colnames(seurat_object) %in% cells_assignment_overlap$overlap$assignment,
#         "OVERLAP",
#         "not-overlap"
#     )

#     png(
#         paste0(
#             "output/plots/aucell-",
#             seurat_object_name,
#             "-dimplot-overlap",
#             ".png"
#         ),
#         width = 750,
#         height = 750
#     )
#     DimPlot(seurat_object, group.by = "overlap", pt.size = 0.1, shuffle = TRUE) |> print()
#     dev.off()

#     return(seurat_object)
# }

# seurat_objects <- mcmapply(
#     FUN = check_signatures,
#     seurat_objects,
#     names(seurat_objects),
#     mc.cores = 32
# )

saveRDS(fre, "data/processed/fre-not-aligned-aucell.rds")


# scaled_exp <- fre@assays$SCT@scale.data
# scaled_exp_ir_uniq <- scaled_exp[rownames(scaled_exp) %in% ir_uniq$seurat_genes] 

exp <- fre@assays$SCT@data

### IR UNIQUE ###

exp_ir_uniq <- exp[rownames(exp) %in% ir_uniq$seurat_genes,]
ir_uniq_means <- apply(exp_ir_uniq, 2, mean)

if ( all( names(ir_uniq_means) == rownames(fre@meta.data) )) {
    fre@meta.data$ir_uniq_mean <- ir_uniq_means
} else {
    print("Something is wrong...")
}

png(
    "output/plots/mean-ir-uniq-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
FeaturePlot(fre,
    features = "ir_uniq_mean",
    pt.size = 0.1)
dev.off()

png(
    "output/plots/hist-mean-ir-uniq-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
ggplot(fre@meta.data, aes(x = ir_uniq_mean, fill = orig.ident)) + geom_histogram(alpha = 0.5)
dev.off()

### OVERLAP ###

exp_overlap <- exp[rownames(exp) %in% overlap$seurat_genes,]
overlap_means <- apply(exp_overlap, 2, mean)

if ( all( names(overlap_means) == rownames(fre@meta.data) )) {
    fre@meta.data$overlap_mean <- overlap_means
} else {
    print("Something is wrong...")
}

png(
    "output/plots/mean-overlap-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
FeaturePlot(fre,
    features = "overlap_mean",
    pt.size = 0.1)
dev.off()

png(
    "output/plots/hist-mean-overlap-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
ggplot(fre@meta.data, aes(x = overlap_mean, fill = orig.ident)) + geom_histogram(alpha = 0.5)
dev.off()

### APOPTOSIS ###

exp_apoptosis <- exp[rownames(exp) %in% apoptosis$seurat_genes,]
apoptosis_means <- apply(exp_apoptosis, 2, mean)

if ( all( names(apoptosis_means) == rownames(fre@meta.data) )) {
    fre@meta.data$apoptosis_mean <- apoptosis_means
} else {
    print("Something is wrong...")
}

png(
    "output/plots/mean-apoptosis-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
FeaturePlot(fre,
    features = "apoptosis_mean",
    pt.size = 0.1)
dev.off()

png(
    "output/plots/hist-mean-apoptosis-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
ggplot(fre@meta.data, aes(x = apoptosis_mean, fill = orig.ident)) + geom_histogram(alpha = 0.5)
dev.off()

### DIFFERENTIATION ###

exp_differentiation <- exp[rownames(exp) %in% differentiation$seurat_genes,]
differentiation_means <- apply(exp_differentiation, 2, mean)

if ( all( names(differentiation_means) == rownames(fre@meta.data) )) {
    fre@meta.data$differentiation_mean <- differentiation_means
} else {
    print("Something is wrong...")
}

png(
    "output/plots/mean-differentiation-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
FeaturePlot(fre,
    features = "differentiation_mean",
    pt.size = 0.1)
dev.off()

png(
    "output/plots/hist-mean-differentiation-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
ggplot(fre@meta.data, aes(x = differentiation_mean, fill = orig.ident)) + geom_histogram(alpha = 0.5)
dev.off()

## cor ir vs unique ###

png(
    "output/plots/mean-overlap-vs-mean-ir-uniq-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
ggplot(fre@meta.data, aes(x = ir_uniq_mean, y = overlap_mean)) + geom_point(size = 0.1)
dev.off()

## CORRELATION APOPTOSIS VS DIFFERENTIATION ###

png(
    "output/plots/mean-apoptosis-vs-mean-differentiation-with-no-alignment-with-regression.png",
    width = 750,
    height = 750
)
ggplot(fre@meta.data, aes(y = apoptosis_mean, x = differentiation_mean)) + geom_point(size = 0.1)
dev.off()
