# - Figures -------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("figures/panels")) dir.create("figures/panels")

# read in data
srt <- readRDS("results/integrated.rds")
dataset_colors <- c(
  "GSE125065" = "#117733",
  "GSE130597" = "#6699CC",
  "GSE146020" = "#888888"
)
cluster_colors <- c("#88CCEE", "#CC6677", "#DDCC77", "#AA4499", "#44AA99", "#999933")
names(cluster_colors) <- levels(Idents(srt))

# - A - plot dataset in UMAP --------------------------------------------------
p <- DimPlot(srt, group.by = "Dataset", cells = sample(Cells(srt))) +
  scale_color_manual(values = dataset_colors) +
  theme_void() +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")
ggsave("figures/panels/umap_dataset.png", p, width = 2.3, height = 2.5,
       units = "in", dpi = 400, bg = "white")

# - B - plot clusters with markers and proportions ----------------------------
# UMAP plot of clusters
p1 <- DimPlot(srt, label = TRUE, label.box = TRUE, label.size = 2,
              cols = rep(alpha("white", 0.5), 6)) +
  scale_color_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = "none")
ggsave("figures/panels/umap_clusters.png", p1, width = 2, height = 2,
       units = "in", dpi = 400, bg = "white")

# heatmap plot of top markers
markers <- read.csv("results/markers.csv")
top_markers <- filter(markers, p_val_adj < 0.05) %>%
  group_by(gene) %>%
  arrange(desc(pct.1 - pct.2)) %>%
  slice(1) %>%
  group_by(cluster) %>%
  arrange(desc(pct.1 - pct.2)) %>%
  slice(1:20) %>%
  pull(gene)
if (!all(top_markers %in% rownames(srt[["integrated"]]@scale.data))) {
  srt <- ScaleData(srt, features = top_markers)
}
p2 <- plot_heatmap(srt, top_markers, assay = "integrated",
                   colors = cluster_colors, flipped = TRUE)
ggsave("figures/panels/heatmap.png", p2, width = 1.7, height = 1.9, units = "in",
       dpi = 400, bg = "white")

# plot proportions
df <- table(Idents(srt), srt$Dataset) %>% as.data.frame() %>%
  group_by(Var2) %>%
  mutate("Percent" = Freq / sum(Freq) * 100)
p3 <- ggplot(df, aes(x = Var2, y = Percent, fill = Var1)) +
  geom_col(color = "black") +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   color = "black", size = 6),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 7),
        legend.position = "none")
ggsave("figures/panels/proportions.png", p3, width = 1.1, height = 1.7, units = "in",
       dpi = 400, bg = "white")

# - C - neurotransmitter markers ----------------------------------------------
genes <- c("Slc17a6", "Slc32a1", "Gad1", "Gad2")
plots <- lapply(genes, function(x) plot_gene(srt, x))
p <- cowplot::plot_grid(plotlist = plots, ncol = length(genes))
ggsave("figures/panels/neurotransmitters.png", p, width = 5, height = 1.5,
       units = "in", dpi = 400, bg = "white")
