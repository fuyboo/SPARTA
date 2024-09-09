#' Quality control
#'
#' This function for quality control of cells and to determine the thresholds for sgRNA identities in a Seurat object.
#'
#' @param SUM159 A Seurat object containing GDO and motif data.
#' @param count_threshold
#' @param feature_threshold
#' @param percent_mt_threshold
#' @param save_path
#' @param assay used to determine the identity of gRNAs in cells.
cell_quality <- function(SUM159,
                         count_threshold = 100,
                         feature_threshold = 200,
                         percent_mt_threshold = 10,
                         assay = "sgRNA",
                         save_path = NULL) {
  # Filter the data
  DefaultAssay(object = SUM159) <- "RNA"
  SUM159[["percent.mt"]] <- PercentageFeatureSet(SUM159, pattern = "^MT-")
  SUM159 <- subset(SUM159, subset = nCount_aptamer > count_threshold)
  SUM159 <- subset(SUM159,
    subset = nFeature_RNA > feature_threshold &
      percent.mt < percent_mt_threshold
  )

  # sgRNA first classification
  cell_sgrna <- as.data.frame(SUM159$nCount_sgRNA)
  colnames(cell_sgrna) <- "count"

  p1 <- ggplot(cell_sgrna) +
    geom_histogram(aes(x = log10(count))) +
    theme(
      legend.position = "none",
      text = element_text(size = 6),
      legend.key.width = unit(0.2, "line"),
      legend.key.height = unit(0.2, "line"),
      legend.text.align = 0
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.7),
      plot.title = element_text(size = 20, family = "sans"),
      axis.text = element_text(size = 10, family = "sans"),
      axis.title.x = element_text(size = 22, family = "sans"),
      axis.title.y = element_text(size = 22, family = "sans")
    ) +
    xlab("log10(sgRNA UMIs)") +
    ylab("Frequency")
  if (!is.null(save_path)) {
    ggsave(p1, filename = paste0(save_path, "/sgRNA_qc1.pdf"), width = 6, height = 6)
  }
  p1

  # Calculate top gRNA ratio
  top_gRNA_ratio <- apply(SUM159[[assay]]@counts, 2, function(t) {
    sort(t, decreasing = TRUE)[1] / sum(t)
  })
  top_gRNA_ratio <- as.data.frame(top_gRNA_ratio)
  colnames(top_gRNA_ratio) <- "top_gRNA_ratio"
  p2 <- ggplot(top_gRNA_ratio) +
    geom_histogram(aes(x = top_gRNA_ratio)) +
    theme(
      legend.position = "none",
      text = element_text(size = 6),
      legend.key.width = unit(0.2, "line"),
      legend.key.height = unit(0.2, "line"),
      legend.text.align = 0
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.7),
      plot.title = element_text(size = 20, family = "sans"),
      axis.text = element_text(size = 10, family = "sans"),
      axis.title.x = element_text(size = 22, family = "sans"),
      axis.title.y = element_text(size = 22, family = "sans")
    ) +
    xlab("Enrichment Ratio") +
    ylab("Frequency")
  p2
  if (!is.null(save_path)) {
    ggsave(p2, filename = paste0(save_path, "/sgrna_qc2.pdf"), width = 6, height = 6)
  }

  grid.arrange(p1, p2, ncol = 2)

  return(SUM159)
}
