# - Merging datasets ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
source("scripts/functions.R")
if (!dir.exists("results/plots")) dir.create("results/plots")

# read in data
samples <- list.files("results", pattern = "GSE[0-9]+.rds")
srt <- lapply(samples, function(x) readRDS(paste0("results/", x)))
samples <- str_extract(samples, "GSE[0-9]+")
names(srt) <- samples
for (i in 1:length(srt)) srt[[i]]$Dataset <- samples[i]

# run Seurat CCA to integrate
anchors <- FindIntegrationAnchors(srt)
k <- min(sapply(srt, ncol))-1
srt <- IntegrateData(anchors, k.weight = k)

# process integrated data
srt <- ScaleData(srt) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:8) %>%
  FindNeighbors(dims = 1:8)

# run clustering
srt <- cluster(srt)
mtx <- sapply(unique(Idents(srt)), function(x) {
  Matrix::rowMeans(srt[["integrated"]]@scale.data[, Idents(srt) == x])
})
colnames(mtx) <- unique(Idents(srt))
tree <- hclust(dist(t(mtx)))
cluster_order <- tree$labels[tree$order]
Idents(srt) <- factor(Idents(srt), levels = cluster_order, labels = 1:length(cluster_order))

# plot datasets and clusters in UMAP space
dataset_colors <- RColorBrewer::brewer.pal(n = length(unique(srt$Dataset)), "Accent")
cluster_colors <- RColorBrewer::brewer.pal(n = length(cluster_order), "Set2")
names(cluster_colors) <- levels(Idents(srt))
p <- cowplot::plot_grid(
  DimPlot(srt, group.by = "Dataset", cells = sample(Cells(srt))) +
    ggtitle("Dataset") +
    scale_color_manual(values = dataset_colors) +
    theme_void() +
    theme(legend.text = element_text(size = 6),
          plot.title = element_text(size = 8, hjust = 0.5)),
  DimPlot(srt, label = TRUE, label.size = 2) +
    ggtitle("Clusters") +
    theme_void() +
    scale_color_manual(values = cluster_colors) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, hjust = 0.5)),
  rel_widths = c(1.1, 1)
)
ggsave("results/plots/umap.png", p, width = 6.5, height = 3, units = "in", dpi = 400)

# proportions
p <- DimPlot(srt, split.by = "Dataset", label = TRUE, label.size = 2) +
  scale_color_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 8))
ggsave("results/plots/umap_by_dataset.png", p, width = 6, height = 2, units = "in",
       dpi = 400, bg = "white")
df <- table(Idents(srt), srt$Dataset) %>% as.data.frame() %>%
  group_by(Var2) %>%
  mutate("Percent" = Freq / sum(Freq) * 100)
p <- ggplot(df, aes(x = Var2, y = Percent, fill = Var1)) +
  geom_col(color = "black") +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   color = "black", size = 8),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.width = unit(0.05, "in"),
        legend.key.height = unit(0.05, "in"))
ggsave("results/plots/proportions.png", p, width = 2, height = 3, units = "in", dpi = 400)

# Pmch levels
df <- bind_cols(
  as.data.frame(srt[["umap"]]@cell.embeddings),
  data.frame("Pmch" = srt[["RNA"]]@data["Pmch", ], "Dataset" = srt$Dataset)
)
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Pmch)) +
  geom_point(stroke = 0, size = 1) +
  theme_void() +
  scale_color_gradient(low = "gray90", high = "firebrick4") +
  facet_wrap(~Dataset, ncol = 1) +
  theme(legend.key.width = unit(0.02, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        strip.text = element_text(size = 8))
ggsave("results/plots/umap_Pmch_by_dataset.png", p, width = 2.5, height = 6.4,
       units = "in", dpi = 400, bg = "white")

# find markers
markers <- FindAllMarkers(srt, only.pos = TRUE)
write.csv(markers, "results/markers.csv", row.names = FALSE)
top_markers <- filter(markers, p_val_adj < 0.05) %>%
  group_by(gene) %>%
  arrange(desc(pct.1 - pct.2)) %>%
  slice(1) %>%
  group_by(cluster) %>%
  arrange(desc(pct.1 - pct.2)) %>%
  slice(1:20) %>%
  pull(gene)

if (!all(top_markers %in% rownames(srt[["integrated"]]@scale.data))) {
  to_scale <- union(rownames(srt[["integrated"]]@scale.data), top_markers)
  srt <- ScaleData(srt, features = to_scale)
}
p <- plot_heatmap(srt, top_markers, assay = "integrated",
                  colors = cluster_colors, flipped = TRUE)
ggsave("results/plots/heatmap.png", p, width = 4, height = 4, units = "in", dpi = 400, bg = "white")

# plot a summary heatmap
DefaultAssay(srt) <- "RNA"
srt <- ScaleData(srt, top_markers, split.by = "Dataset")
mtx <- sapply(levels(Idents(srt)), function(x) {
  cells <- Idents(srt) == x
  Matrix::rowMeans(srt[["RNA"]]@scale.data[, cells])
})
mtx <- apply(mtx, c(1, 2), function(x) ifelse(x > 2, 2, ifelse(x < -2, -2, x)))

# order genes & clusters
tree1 <- hclust(dist(mtx))
tree2 <- hclust(dist(t(mtx)))
df <- mtx %>% as.data.frame() %>% tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Cluster", values_to = "Z") %>%
  mutate(Gene = factor(Gene, levels = tree1$labels[tree1$order]),
         Cluster = factor(Cluster, levels = tree2$labels[tree2$order]))
p <- ggplot(df, aes(x = Gene, y = Cluster, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080",
                       name = "z-score") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 8, face = "italic"),
        axis.text.y = element_text(size = 8),
        legend.key.width = unit(0.02, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
ggsave("results/plots/marker_average_heatmap.png", p, width = 10, height = 2.5,
       units = "in", dpi = 400, bg = "white")

# cluster 3 vs. 4 markers
markers <- FindConservedMarkers(
  srt,
  ident.1 = 3,
  ident.2 = 4,
  grouping.var = "Dataset",
  meta.method = metap::logitp
)
# reformat for saving
stats <- apply(markers, 1, function(x) {
  c("avg_log2FC" = mean(x[str_detect(names(x), "avg_log2FC")]),
    "pct.3" = mean(x[str_detect(names(x), "pct.1")]),
    "pct.4" = mean(x[str_detect(names(x), "pct.2")]),
    "p_val_adj" = unname(x["minimump_p_val"])
  )
})
markers <- t(stats) %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
  mutate_at(vars(avg_log2FC, pct.3, pct.4), ~ round(.x, 2)) %>%
  mutate(p_val_adj = format(p_val_adj, digits = 2, scientific = TRUE))
write.csv(markers, "results/markers_3_v_4.csv", row.names = FALSE, na = "")

# save integrated object
saveRDS(srt, "results/integrated.rds")
