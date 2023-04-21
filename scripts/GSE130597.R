# - Analyzing Rossi et al. dataset --------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/GSE130597")) dir.create("results/GSE130597")

# read in data
samples <- list.files("data/GSE130597")
mtx <- lapply(samples, function(x) {
  filename <- paste0("data/GSE130597/", x)
  read.table(filename, header = TRUE, sep = "\t", row.names = 1)
})

# keep genes that are in all datasets
gene_frequency <- table(unlist(sapply(mtx, rownames)))
common_genes <- names(gene_frequency)[gene_frequency == length(mtx)]
mtx <- lapply(mtx, function(x) x[common_genes, ])

# make sure there are no duplicate cell names
for (i in 1:length(mtx)) colnames(mtx[[i]]) <- paste0(colnames(mtx[[i]]), "-", i)

# get sample lengths
samples <- str_extract(samples, "GSM[0-9]+")
lengths <- sapply(mtx, ncol)

# combine matrices and create Seurat object
mtx <- do.call(cbind, mtx)
srt <- CreateSeuratObject(mtx)
srt$Sample <- rep(samples, lengths)

# process all cells
srt <- NormalizeData(srt) %>%
  FindVariableFeatures() %>%
  ScaleData(split.by = "Sample", vars.to.regress = "nFeature_RNA") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:15)

# plot by sample
plot_settings <- list(
  theme_void(),
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5))
)
p <- DimPlot(srt, group.by = "Sample", cells = sample(Cells(srt))) + plot_settings
ggsave("results/GSE130597/umap_sample.png", p, width = 5, height = 4,
       units = "in", dpi = 400, bg = "white")

# neuron markers
p <- FeaturePlot(srt, "Rbfox1") + plot_settings +
  theme(legend.key.width = unit(0.02, "in"))
ggsave("results/GSE130597/umap_Rbfox1.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")

# identify neurons by expression of known neuron markers
srt <- FindNeighbors(srt, dims = 1:15) %>% FindClusters()
markers <- c("Syn1", "Syp", "Rbfox1", "Rbfox3")
mtx <- sapply(levels(Idents(srt)), function(x) {
  Matrix::rowMeans(srt[["RNA"]]@scale.data[markers, Idents(srt) == x])
})
neuron_clusters <- apply(mtx, 2, function(x) any(x > 0)) %>% .[.] %>% names()

# process neurons and remove low depth cells
neurons <- subset(srt, idents = neuron_clusters) %>%
  subset(nFeature_RNA > 700) %>%
  FindVariableFeatures() %>%
  ScaleData(split.by = "Sample", vars.to.regress = "nFeature_RNA") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:25) %>%
  FindNeighbors(dims = 1:25) %>%
  FindClusters()

# plot neuron clusters and MCH expression
p <- DimPlot(neurons, label = TRUE, label.size = 2) +
  ggtitle("Neuron clusters") +
  plot_settings +
  theme(legend.position = "none")
ggsave("results/GSE130597/umap_neuron_clusters.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")

# plot MCH expression
p <- FeaturePlot(neurons, "Pmch") +
  plot_settings +
  theme(legend.key.width = unit(0.02, "in"))
ggsave("results/GSE130597/umap_neurons_Pmch.png", p, width = 4, height = 4,
              units = "in", dpi = 400, bg = "white")

# identify MCH-expressing population
mch_expr <- sapply(split(neurons[["RNA"]]@data["Pmch", ], Idents(neurons)), mean)
mch_cluster <- names(mch_expr)[mch_expr == max(mch_expr)]
mch <- subset(neurons, idents = mch_cluster)

# save object
saveRDS(mch, "results/GSE130597.rds")
