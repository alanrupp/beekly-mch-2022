# - Helper functions ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

#' choose optimal clustering based on maximizing silhouette width
cluster <- function(object, resolutions = seq(0.2, 2, by = 0.2)) {
  # cluster at all resolutions until reaching maximum silhouette width
  repeat {
    clusters <- lapply(resolutions, function(x) {
      FindClusters(object, resolution = x, verbose = FALSE)@active.ident
    })
    clusters <- clusters[!sapply(clusters, is.null)]
    clusters <- do.call(cbind, clusters)
    if (length(clusters) == 0) break()
    if (ncol(clusters) == 1) break()
    distances <- cluster::daisy(object[["pca"]]@cell.embeddings[, 1:8])
    # get mean silhouette width for a given clustering
    widths <- apply(clusters, 2, function(x) {
      if (length(unique(x)) > 1) {
        sil <- cluster::silhouette(as.numeric(x), distances)
        mean(sil[, 3])
      } else {
        NA
      }
    })
    if (all(is.na(widths))) {
      best_width <- 1
    } else {
      best_width <- which(widths == max(widths, na.rm = TRUE))[1]
    }
    if (best_width == length(widths)) {
      resolutions <- seq(max(resolutions), max(resolutions) + 2, by = 0.2)
    } else {
      break()
    }
  }
  # assign maximum silhouette width clustering
  if (length(clusters) == 0) {
    clusters <- rep("x", ncol(object))
  } else if (ncol(clusters) == 1) {
    clusters <- clusters[, 1]
  } else {
    clusters <- clusters[, best_width]
  }
  # return object with clusters
  Idents(object) <- clusters
  return(object)
}

#' Plot gene expression across cells
plot_heatmap <- function(object, genes, limit = 2,
                         assay = "RNA", slot = "scale.data",
                         n_cells = 1000, flipped = FALSE, bars = TRUE,
                         colors = NULL, balanced = TRUE) {
  mtx <- slot(object[[assay]], slot)
  if (!all(genes %in% rownames(mtx))) {
    stop("Not all genes are in the ", slot, " slot")
  }
  if (ncol(mtx) > n_cells) {
    if (balanced) {
    each <- round(n_cells / length(unique(Idents(object))), 0)
    cells <- lapply(levels(Idents(object)), function(x) {
      options <- Cells(object)[Idents(object) == x]
      sample(options, min(each, length(options)))
    }) %>% unlist()
    mtx <- mtx[, cells]
    } else {
      mtx <- mtx[, sample(colnames(mtx), n_cells)]
    }
  }
  cell_order <- Idents(object)[colnames(mtx)] %>% sort() %>% names()
  if (any(abs(mtx) > limit)) {
    mtx <- apply(mtx, c(1,2), function(x) ifelse(x < -limit, -limit, ifelse(x > limit, limit, x)))
  }
  # make tidy
  df <- as.data.frame(mtx) %>% tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Cell", values_to = "Expr")
  if (flipped) {
    df <- df %>%
      mutate(Gene = factor(Gene, levels = genes[length(genes):1])) %>%
      mutate(Cell = factor(Cell, levels = cell_order))
  } else {
    df <- df %>%
      mutate(Gene = factor(Gene, levels = genes)) %>%
      mutate(Cell = factor(Cell, levels = cell_order[length(cell_order):1]))
  }
  # plot
  p <- ggplot(df, aes(x = Gene, y = Cell, fill = Expr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080",
                         breaks = seq(-2, 2), labels = seq(-2, 2)) +
    theme_void() +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.key.height = unit(0.02, "in"),
          legend.key.width = unit(0.1, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_blank())
  if (bars) {
    if (is.null(colors)) {
      colors <- sample(grDevices::colors(), length(unique(Idents(object))))
      names(colors) <- levels(Idents(object))
    }
    height <- length(unique(df$Gene)) * 0.03
    for (i in levels(Idents(object))) {
      p <- p + annotate(
        "rect",
        xmin = 0.5 - height,
        xmax = 0.5,
        ymin = min(which(levels(df$Cell) %in% Cells(object)[Idents(object) == i])),
        ymax = max(which(levels(df$Cell) %in% Cells(object)[Idents(object) == i])),
        fill = colors[i],
        color = "black",
        size = 0.5
      )
    }
  }
  if (flipped) p <- p + coord_flip()
  p
}

#' plot gene in UMAP space
plot_gene <- function(object, gene) {
  df <- bind_cols(
    as.data.frame(object[["umap"]]@cell.embeddings),
    data.frame("Gene" = object[["RNA"]]@data[gene, ])
  )
  df <- slice_sample(df, n = nrow(df))
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Gene)) +
    geom_point(stroke = 0, size = 1) +
    theme_void() +
    ggtitle(gene) +
    scale_color_gradient(low = "gray90", high = "#332288") +
    theme(legend.key.height = unit(0.02, "in"),
          legend.key.width = unit(0.1, "in"),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          plot.title = element_text(size = 8, face = "italic", hjust = 0.5))
}
