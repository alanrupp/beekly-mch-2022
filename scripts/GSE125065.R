# - Analyzing Mickelsen et al. dataset ----------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/GSE125065")) dir.create("results/GSE125065")

# load annotated LHA dataset from Jax single cell portal
load("data/lha-2019_neuronal.rds")
log2cpm <- as.matrix(log2cpm)

# identify cluster that expresses the highest Pmch
mch_ensembl <- rownames(featuredata)[featuredata$Associated.Gene.Name == "Pmch"]
mch <- split(log2cpm[mch_ensembl, ], tsne.data$dbCluster) %>% sapply(mean)
mch_cluster <- names(mch)[mch == max(mch)]
mch_cells <- rownames(tsne.data)[tsne.data$dbCluster == mch_cluster]

# read in raw data from GEO
read_data <- function(sample) {
  barcodes <- read.table(paste0("data/GSE125605/", sample, "_barcodes.tsv.gz"))
  genes <- read.table(paste0("data/GSE125605/", sample, "_genes.tsv.gz"))
  mtx <- Matrix::readMM(paste0("data/GSE125605/", sample, "_matrix.mtx.gz"))
  colnames(mtx) <- barcodes$V1
  rownames(mtx) <- genes$V1
  # remove all 0 rows
  mtx <- mtx[apply(mtx, 1, function(x) any(x > 0)), ]
  genes <- tibble::deframe(genes) %>% .[rownames(mtx)]
  # remove duplicate gene names
  if (any(duplicated(genes))) {
    # keep the highest-expressing version of a gene name
    dups <- unique(genes[duplicated(genes)])
    drop <- vector("character")
    for (dup in dups) {
      expr <- Matrix::rowSums(mtx[names(genes)[genes == dup], ])
      drop <- c(drop, names(expr)[expr != max(expr)])
    }
    mtx <- mtx[!rownames(mtx) %in% drop, ]
  }
  rownames(mtx) <- genes[rownames(mtx)]
  return(mtx)
}
samples <- list.files("data/GSE125605")
samples <- unique(str_extract(samples, "GSM[0-9]+_AJ[0-9]+"))
mtx <- lapply(samples, read_data)

# rename second sample barcodes to avoid clashes
colnames(mtx[[2]]) <- paste0(str_extract(colnames(mtx[[2]]), "[A-Z]+"), "-2")

# combine matrices
common_genes <- intersect(rownames(mtx[[1]]), rownames(mtx[[2]]))
mtx <- lapply(mtx, function(x) x[common_genes, ])
mtx <- do.call(cbind, mtx)

# - Make seurat object and process --------------------------------------------
srt <- CreateSeuratObject(mtx) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:15)

# add sample names and Jax annotation
srt$Sample <- str_extract(Cells(srt), "[1,2]$")
srt$Sample <- recode(srt$Sample, !!!setNames(samples, c(1, 2)))
srt$Jax <- NA
srt$Jax[rownames(tsne.data)] <- tsne.data$dbCluster

# plot all cells by sample and jax cluster ID
plot_settings <- list(
  theme_void(),
  theme(legend.text = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.key.width = unit(0.02, "in"),
        legend.title = element_text(size = 7))
)
p <- DimPlot(srt, group.by = "Sample", cells = sample(Cells(srt))) +
  plot_settings
ggsave("results/GSE125065/umap_sample.png", p, width = 5, height = 4,
       units = "in", dpi = 400, bg = "white")
p <- DimPlot(srt, group.by = "Jax", label = TRUE, label.size = 2) +
  ggtitle("Jax neuron cluster label") +
  plot_settings +
  theme(legend.position = "none")
ggsave("results/GSE125065/umap_jax.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")

# pull neurons and cluster to identify MCH cells
neurons <- subset(srt, cells = rownames(tsne.data)) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 1)

# plot neuorn clusters and MCH expression
p <- DimPlot(neurons, label = TRUE, label.size = 2) + plot_settings +
  ggtitle("Neuron clusters") +
  theme(legend.position = "none")
ggsave("results/GSE125065/umap_neuron_clusters.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")
p <- FeaturePlot(neurons, "Pmch") + plot_settings
ggsave("results/GSE125065/umap_neurons_Pmch.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")

# keep cluster with most matching MCH neurons
mch_table <- table(Idents(neurons), neurons$Jax == mch_cluster)
mch_cluster <- sort(mch_table[, "TRUE"], decreasing = TRUE)[1] %>% names()

# save MCH neurons
mch <- subset(neurons, idents = mch_cluster)
saveRDS(mch, "results/GSE125605.rds")
