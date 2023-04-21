# - Analyzing Jiang et al. dataset --------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/GSE146020")) dir.create("results/GSE146020")

# read in data
df <- read.csv("data/GSE146020/GSM4354862_count_matrix.csv.gz", row.names = 1)

# process and embed
srt <- CreateSeuratObject(df) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters()

# plot cells
plot_settings <- list(
  theme_void(),
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5))
)
p <- DimPlot(srt, label = TRUE, label.size = 2) +
  ggtitle("MCH clusters") +
  plot_settings +
  theme(legend.position = "none")
ggsave("results/GSE146020/umap_clusters.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")

# plot MCH
p <- FeaturePlot(srt, "Pmch") +
  plot_settings +
  theme(legend.key.width = unit(0.02, "in"))
ggsave("results/GSE146020/umap_Pmch.png", p, width = 4, height = 4,
       units = "in", dpi = 400, bg = "white")

# save object
saveRDS(srt, "results/GSE146020.rds")
