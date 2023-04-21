# - Find neurotransmitter markers across MCH subtypes -------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

# read in data
srt <- readRDS("results/integrated.rds")

genes <- c("Gria3", "Grin2b", "Grm2", "Grm3", "Grm4", "Grm5", "Cdk5r1", "Avp",
           "P2rx7", "Aldh5a1", "Adora2a", "Gabbr2", "Gabra4", "Slc6a12",
           "Slc6a13", "Abat", "Gls", "Slc17a6", "Slc17a7", "Slc17a8", "Slc32a1",
           "Gad1", "Gad2", "Slc18a2")

DefaultAssay(srt) <- "RNA"
srt <- ScaleData(srt, genes, split.by = "Dataset")
mtx <- sapply(levels(Idents(srt)), function(x) {
  cells <- Idents(srt) == x
  Matrix::rowMeans(srt[["RNA"]]@scale.data[, cells])
})

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
ggsave("results/plots/neurotransmitter_heatmap.png", p, width = 4, height = 2.5,
       units = "in", dpi = 400, bg = "white")

#' plot genes
plot_gene <- function(gene, dataset = FALSE, only = NULL) {
  df <- bind_cols(
    as.data.frame(srt[["umap"]]@cell.embeddings),
    data.frame("Gene" = srt[["RNA"]]@data[gene, ], "Dataset" = srt$Dataset)
  )
  if (!is.null(only)) {
    df <- filter(df, Dataset == only)
  }
  if (all(df$Gene == 0)) {
    colors <- c("gray90", "gray90")
  } else {
    colors <- c("gray90", "firebrick4")
  }
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Gene)) +
    geom_point(stroke = 0, size = 1) +
    ggtitle(gene) +
    scale_color_gradient(low = colors[1], high = colors[2]) +
    theme_void() +
    theme(legend.key.width = unit(0.02, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9, face = "italic"),
          strip.text = element_text(size = 8))
  if (dataset) {
    p <- p + facet_wrap(~Dataset)
  }
  return(p)
}

plots <- lapply(tree1$labels[tree1$order], plot_gene)
p <- cowplot::plot_grid(plotlist = plots, ncol = 6)
ggsave("results/plots/umap_neurotransmitters.png", p, width = 10, height = 6.7,
       units = "in", dpi = 400, bg = "white")

# plot by dataset
by_dataset <- function(dataset) {
  plots <- lapply(tree1$labels[tree1$order], function(x) plot_gene(x, only = dataset))
  p <- cowplot::plot_grid(plotlist = plots, ncol = 6)
  filename <- paste0("results/plots/umap_neurotransmitters_", dataset, ".png")
  ggsave(filename, p, width = 10, height = 6.7, units = "in", dpi = 400, bg = "white")
}
for (dataset in unique(srt$Dataset)) {
  by_dataset(dataset)
}
