#' Visualize Prediction Results
#'
#' This function generates a plot for prediction results and saves it if a path is provided.
#'
#' @param df A data frame containing the prediction results.
#' @param clust A character string specifying the cluster to use.
#' @param save_path A character string specifying the path to save the plot. If NULL, the plot will not be saved.
#' @return A ggplot object.
#' @export
visualize_aptamer_difference <-
  function(predict_result, clust, save_path = NULL) {
    # Set your color palettes
    default_color <- "#696969"
    colors <-
      c(
        "#53ABD8",
        "#BB7018",
        "#7A4BDE",
        "#549154",
        "#CB8F40",
        "#DA544D",
        "#B291B5",
        "#BB6B84",
        "#83A388"
      )
    protein_colors <-
      c(
        "#53ABD8",
        "#FF8C00",
        "#AB82FF",
        "#7CCD7C",
        "#DEB887",
        "#F08080",
        "#B291B5",
        "#BB6B84",
        "#83A388"
      )
    title_colors <-
      c(
        "#53ABD8",
        "#FF8C00",
        "#AB82FF",
        "#549154",
        "#CB8F40",
        "#DA544D",
        "#B291B5",
        "#B291B5",
        "#83A388"
      )

    # Example data
    df <- predict_result$gRNA_cell_dif_median_fi
    aptamer_pro <- predict_result$predict_pro
    rownames(aptamer_pro) <- aptamer_pro$col_name

    proteins <- unlist(aptamer_pro[clust, 2])
    protein_number <- length(proteins)

    dat <- df[, c(clust, "group")] %>%
      filter(!is.na(.[[1]]))
    colnames(dat)[1] <- "ss"

    # Filter significant data
    dat_sig <- dat
    threshold <- as.numeric(unlist(aptamer_pro[clust, 3]))
    dat_sig$log <-
      ifelse(
        0 < dat_sig$ss &
          dat_sig$ss < 1 | -1 < dat_sig$ss & dat_sig$ss < 0,
        dat_sig$ss,
        sign(dat_sig$ss) * log2(abs(dat_sig$ss) + 1)
      )
    dat_sig$group <-
      ifelse(dat_sig$log < threshold &
        (dat_sig$group %in% proteins),
      dat_sig$group,
      "Control"
      )
    counts <- table(dat_sig$group)
    dat_sig$group[dat_sig$group %in% names(counts[counts == 1])] <-
      "Control"

    com_dif <- compare_means(ss ~ group, dat_sig, method = "t.test")
    com_dif <- com_dif[com_dif$group1 == "Control", ]

    filter_target <- rownames(dat)[dat$group %in% proteins]
    filter_target <-
      filter_target[!(filter_target %in% rownames(dat_sig)[dat_sig$group %in% proteins])]
    if (length(filter_target) > 0) {
      dat <- dat[-which(rownames(dat) == filter_target), ]
    }

    dat %>%
      group_by(group) %>%
      summarize(
        mean_value = mean(ss),
        std_error = sd(ss) / sqrt(n())
      ) %>%
      ungroup() -> new.dat
    rownames(new.dat) <- new.dat$group

    need_pro <- new.dat[proteins, ]

    word_col <- rep(default_color, nrow(new.dat))
    errbar_col <- rep(default_color, nrow(new.dat))
    word_col[new.dat$group %in% proteins] <- colors[1:protein_number]
    title_col <- title_colors[1:protein_number]
    errbar_col[new.dat$group %in% proteins] <-
      protein_colors[1:protein_number]
    protein_col <- protein_colors[1:protein_number]

    line_position <- sign(threshold) * (2**abs(threshold) - 1)

    p <- ggplot(new.dat, aes(x = group, y = mean_value)) +
      geom_point(size = 4.5, color = "#696969") +
      geom_errorbar(
        aes(
          ymin = mean_value - std_error,
          ymax = mean_value + std_error
        ),
        width = 0.3,
        size = 1,
        color = errbar_col
      ) +
      geom_point(
        data = need_pro,
        aes(x = group, y = mean_value, fill = proteins),
        size = 5,
        color = protein_col
      ) +
      geom_hline(
        yintercept = line_position,
        linetype = "dashed",
        size = 1.5,
        color = default_color
      ) +
      theme_minimal() +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      ggtitle(clust) +
      ylab("Median difference") +
      xlab("") +
      theme(
        panel.background = element_blank(),
        text = element_text(
          size = 12,
          family = "sans",
          face = "bold"
        ),
        axis.title = element_text(size = 12, family = "sans"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(
          size = 12,
          angle = 30,
          vjust = 0.85,
          hjust = 0.75
        )
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(
        plot.title = element_text(
          size = 20,
          family = "sans",
          colour = title_col,
          vjust = -1
        ),
        axis.title.x = element_text(size = 12, family = "sans"),
        axis.title.y = element_text(size = 12, family = "sans")
      ) +
      theme(axis.text.x = element_text(colour = word_col)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      NoLegend()

    annotations <- data.frame(
      gene = new.dat$group,
      label = ifelse(new.dat$group %in% com_dif$group2, com_dif$p.signif[match(new.dat$group, com_dif$group2)], "")
    )

    y_min <- min(new.dat$mean_value - new.dat$std_error)
    y_max <- max(new.dat$mean_value + new.dat$std_error)
    y_ann <- y_min - (y_max - y_min) * 1 / 4
    ann_color <- rep(default_color, nrow(new.dat))
    ann_color[match(com_dif$group2, new.dat$group)] <-
      colors[1:protein_number]
    p <- p + coord_cartesian(clip = "off", ylim = c(y_min, y_max)) +
      geom_text(
        data = annotations,
        aes(x = gene, y = y_ann, label = label),
        color = ann_color,
        size = 6,
        vjust = 1
      )

    # Save the plot if a save_path is provided
    if (!is.null(save_path)) {
      ggsave(p,
        filename = save_path,
        width = 5,
        height = 5
      )
    }

    return(p)
  }
