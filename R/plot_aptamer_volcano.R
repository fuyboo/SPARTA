# ------------------------------------------------------------------------------
# Aptamer-level specificity validation: volcano plots
#
# This script provides a reusable function `plot_aptamer_volcano()` that:
#   - takes a Seurat object and a target protein name,
#   - performs differential aptamer binding analysis (KO vs Control) using limma,
#   - and generates a volcano plot showing log2 fold change vs -log10(P value).
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# ------------------------------------------------------------------------------
# Helper: select cells closest to the median value in each group
# ------------------------------------------------------------------------------

select_near_median <- function(sub_df, n = 150) {
  median_value <- median(sub_df$value)
  sub_df <- sub_df %>%
    mutate(dist_to_median = abs(value - median_value)) %>%
    arrange(dist_to_median) %>%
    head(n)
  sub_df$dist_to_median <- NULL
  sub_df
}

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

#' Plot volcano for aptamer-level specificity validation
#'
#' @param seurat_obj    A Seurat object containing aptamer and motif assays.
#' @param protein       Target protein name (used to define KO group).
#' @param aptamer_assay Name of aptamer assay in the Seurat object (default: "aptamer").
#' @param motif_assay   Name of motif assay in the Seurat object (default: "motif1w").
#' @param need_motif    Character vector of motif cluster IDs (e.g. "Clust-1", ...).
#' @param exclude_targets Character vector of sgRNA targets to exclude from the analysis.
#' @param protein_colors Named vector of base colors for each protein.
#' @param protein_clusters Named list: each protein -> vector of motif clusters to use.
#' @param protein_highlights Named list: each protein -> vector of aptamer IDs to highlight.
#' @param n_cells_per_group Number of cells per group (near median) per cluster.
#' @param output_dir    Directory to save volcano plots (PDF).
#'
#' @return A ggplot object of the volcano plot.
#'
plot_aptamer_volcano <- function(
  seurat_obj,
  protein,
  aptamer_assay      = "aptamer",
  motif_assay        = "motif1w",
  need_motif         = paste0("Clust-", c(1, 2, 3, 4, 5, 6, 7, 11, 13, 14, 15, 17)),
  exclude_targets    = c("ITGB1-2", "CD151-3", "PTPRF-3", "ITGA3-2", "CDCP1-1", "NRP2-2"),
  protein_colors     = NULL,
  protein_clusters   = NULL,
  protein_highlights = NULL,
  n_cells_per_group  = 150,
  output_dir         = "volcano"
) {
  # --------------------------------------------------------------------------
  # Set default color palette and cluster mapping if not provided
  # --------------------------------------------------------------------------
  if (is.null(protein_colors)) {
    protein_colors <- c(
      "CDCP1" = "#53ABD8",
      "ITGA3" = "#FF8C00",
      "ITGB1" = "#AB82FF",
      "NRP1"  = "#7CCD7C",
      "NRP2"  = "#DEB887",
      "PTK7"  = "#F08080",
      "PTPRF" = "#B291B5",
      "PTPRD" = "#BB6B84"
    )
  }

  if (is.null(protein_clusters)) {
    # NOTE: this mapping encodes which motif clusters are used per protein.
    # You can modify it to match your analysis exactly.
    protein_clusters <- list(
      "PTPRD" = c("Clust-1", "Clust-15"),
      "PTPRF" = c("Clust-4"),
      "ITGA3" = c("Clust-3"),
      "CDCP1" = c("Clust-5", "Clust-13", "Clust-17"),
      "NRP1"  = c("Clust-2", "Clust-6", "Clust-7"),
      "NRP2"  = c("Clust-11", "Clust-14")
      # Add more proteins here if needed
    )
  }

  if (is.null(protein_highlights)) {
    # Optional: aptamers to highlight specifically for each protein
    protein_highlights <- list(
      "PTK7"  = c("Apt-1-1",  "Apt-15-27"),
      "CDCP1" = c("Apt-13-18", "Apt-5-4",   "Apt-17-39"),
      "NRP1"  = c("Apt-2-2",   "Apt-6-6",   "Apt-7-10"),
      "NRP2"  = c("Apt-11-15", "Apt-14-26"),
      "ITGA3" = c("Apt-3-3",   "Apt-3-42"),
      "PTPRF" = c("Apt-4-5")
      # Add PTPRD-specific aptamers here if desired
    )
  }

  if (!protein %in% names(protein_clusters)) {
    stop("No cluster mapping for protein: ", protein)
  }

  if (!protein %in% names(protein_colors)) {
    stop("No color defined for protein: ", protein)
  }

  # --------------------------------------------------------------------------
  # Local copy of the Seurat object to avoid modifying the original
  # --------------------------------------------------------------------------
  obj <- seurat_obj

  # --------------------------------------------------------------------------
  # Define "group" (KO vs Control) based on metadata
  # Assumes obj$group is something like "PTPRD-1", "NRP1-2", "NC-1", etc.
  # --------------------------------------------------------------------------
  if (!"group" %in% colnames(obj@meta.data)) {
    stop("Metadata column 'group' not found in Seurat object.")
  }

  obj$group <- ifelse(
    grepl("NC", obj$group),
    "Control",
    str_split(obj$group, "\\-", simplify = TRUE)[, 1]
  )

  # --------------------------------------------------------------------------
  # Optionally subset sgRNA targets
  # --------------------------------------------------------------------------
  if ("target" %in% colnames(obj@meta.data) && length(exclude_targets) > 0) {
    obj <- subset(obj, !(target %in% exclude_targets))
  }

  # --------------------------------------------------------------------------
  # Aptamer expression matrix
  # --------------------------------------------------------------------------
  if (!aptamer_assay %in% names(obj@assays)) {
    stop("Aptamer assay '", aptamer_assay, "' not found in Seurat object.")
  }

  apt_exp <- as.data.frame(
    t(as.matrix(obj@assays[[aptamer_assay]]@counts))
  )

  # --------------------------------------------------------------------------
  # Motif expression matrix
  # --------------------------------------------------------------------------
  if (!motif_assay %in% names(obj@assays)) {
    stop("Motif assay '", motif_assay, "' not found in Seurat object.")
  }

  DefaultAssay(obj) <- motif_assay

  rownames(obj@assays[[motif_assay]]@counts) <- gsub(
    "\\.", "-",
    rownames(obj@assays[[motif_assay]]@counts)
  )

  data_motif <- FetchData(
    obj,
    vars = c(need_motif, "group"),
    slot = "counts"
  )

  # --------------------------------------------------------------------------
  # Filter to the current KO group + Control, and keep only relevant clusters
  # --------------------------------------------------------------------------
  clust_vec <- protein_clusters[[protein]]

  data_motif_f <- data_motif[
    data_motif$group %in% c(protein, "Control"),
    c(clust_vec, "group")
  ]
  data_motif_f$cell <- rownames(data_motif_f)

  # --------------------------------------------------------------------------
  # Select cells near the median in each group / cluster
  # --------------------------------------------------------------------------
  filtered_df <- data.frame()

  for (t in clust_vec) {
    data_ff <- data_motif_f[, c(t, "group", "cell")]
    colnames(data_ff)[1] <- "value"

    c_df <- data_ff %>%
      group_by(group) %>%
      group_modify(~ select_near_median(.x, n = n_cells_per_group)) %>%
      ungroup()

    filtered_df <- rbind(filtered_df, c_df)
  }

  filtered_df <- filtered_df[!duplicated(filtered_df$cell), ]

  # --------------------------------------------------------------------------
  # Build aptamer expression matrix for the selected cells
  # --------------------------------------------------------------------------
  exp_unique <- apt_exp[filtered_df$cell, ]
  exp_unique <- as.data.frame(t(exp_unique))
  exp_unique <- log2(exp_unique + 1)

  group_df <- data.frame(df = filtered_df$group, cell = filtered_df$cell)
  rownames(group_df) <- group_df$cell

  # --------------------------------------------------------------------------
  # limma differential analysis: KO vs Control
  # --------------------------------------------------------------------------
  design <- model.matrix(~ 0 + factor(group_df$df))
  colnames(design) <- ifelse(
    grepl(protein, colnames(design)),
    "ss",
    "Control"
  )
  rownames(design) <- rownames(group_df)

  contrast.matrix <- makeContrasts(ss - Control, levels = design)

  fit  <- lmFit(exp_unique, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  DEG <- topTable(fit2, coef = 1, n = Inf)
  DEG$gene <- rownames(DEG)
  DEG$gg   <- ifelse(DEG$logFC >= 0, "no_change", "change")

  # --------------------------------------------------------------------------
  # Define labels for volcano plot
  #   - by default: top 10 aptamers by |logFC|
  #   - optionally: overwrite labels for predefined highlighted aptamers
  # --------------------------------------------------------------------------
  DEG$label <- ""
  top_idx <- head(order(abs(DEG$logFC), decreasing = TRUE), 10)
  DEG$label[top_idx] <- rownames(DEG)[top_idx]

  if (protein %in% names(protein_highlights)) {
    hl <- protein_highlights[[protein]]
    idx_hl <- which(rownames(DEG) %in% hl)
    DEG$label[idx_hl] <- rownames(DEG)[idx_hl]
  }

  # --------------------------------------------------------------------------
  # Volcano plot
  # --------------------------------------------------------------------------
  base_col <- as.character(protein_colors[protein])

  p <- ggplot(DEG, aes(x = logFC, y = -log10(P.Value), color = gg)) +
    geom_point() +
    scale_color_manual(values = c(base_col, "#696969")) +
    theme_bw() +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      text             = element_text(size = 12, family = "sans"),
      axis.title       = element_text(size = 12, family = "sans"),
      axis.line        = element_line(color = "black"),
      axis.text.x      = element_text(size = 12, angle = 0,
                                      vjust = 0.85, hjust = 0.75),
      plot.title       = element_text(hjust = 0.5, size = 12, family = "sans")
    ) +
    labs(
      title = paste0(protein, " sgRNA vs Control"),
      x     = "logFC",
      y     = "log10(P.Value)"
    ) +
    guides(color = "none")

  # label color: default = base_col; highlighted aptamers can use a darker color
  data_lab <- DEG
  lab_col  <- rep(base_col, nrow(data_lab))

  if (protein %in% names(protein_highlights)) {
    hl <- protein_highlights[[protein]]
    idx_hl <- which(rownames(data_lab) %in% hl)
    lab_col[idx_hl] <- "#EC3333"  # highlight color; modify as needed
  }

  p <- p +
    geom_text_repel(
      data          = data_lab[data_lab$label != "", ],
      aes(label     = label),
      size          = 3.5,
      box.padding   = unit(0.5, "lines"),
      point.padding = unit(0.8, "lines"),
      colour        = lab_col[data_lab$label != ""],
      max.overlaps  = 4000,
      show.legend   = FALSE
    )

  # --------------------------------------------------------------------------
  # Save plot
  # --------------------------------------------------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  outfile <- file.path(output_dir, paste0(protein, "_volcano.pdf"))
  ggsave(filename = outfile, plot = p, width = 5, height = 5)

  message("Volcano plot saved to: ", outfile)

  return(p)
}
