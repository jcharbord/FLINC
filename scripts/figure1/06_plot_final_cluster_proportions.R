suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
output_dir <- args[2]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stage_levels <- c(
  "F0",
  "F1",
  "F2",
  "F3",
  "F4"
)

stage_labels <- c(
  "F0" = "Healthy",
  "F1" = "F1",
  "F2" = "F2",
  "F3" = "F3",
  "F4" = "F4"
)

stage_palette <- c(
  "F0" = "#BDBDBD",
  "F1" = "#F8C8D8",
  "F2" = "#F06292",
  "F3" = "#FF7043",
  "F4" = "#C62828"
)

p_to_stars <- function(p) {
  ifelse(
    is.na(p),
    "ns",
    ifelse(
      p < 0.001,
      "***",
      ifelse(
        p < 0.01,
        "**",
        ifelse(
          p < 0.05,
          "*",
          "ns"
        )
      )
    )
  )
}

run_wilcoxon <- function(data, group1, group2) {

  test_data <- data %>%
    filter(
      as.character(F_score) %in% c(group1, group2)
    ) %>%
    mutate(
      F_score = factor(
        as.character(F_score),
        levels = c(group1, group2)
      )
    ) %>%
    droplevels()

  if (length(unique(test_data$F_score)) < 2) {
    return(NA_real_)
  }

  test_result <- suppressWarnings(
    wilcox.test(
      prop ~ F_score,
      data = test_data,
      exact = FALSE,
      correct = FALSE
    )
  )

  test_result$p.value
}

make_pairwise_table <- function(data) {

  pair_tests <- lapply(
    combn(stage_levels, 2, simplify = FALSE),
    function(groups) {

      group1 <- groups[1]
      group2 <- groups[2]

      p_raw <- run_wilcoxon(
        data,
        group1,
        group2
      )

      data.frame(
        group1 = group1,
        group2 = group2,
        p_raw = p_raw,
        n1 = sum(
          as.character(data$F_score) == group1
        ),
        n2 = sum(
          as.character(data$F_score) == group2
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(pair_tests) %>%
    mutate(
      p_adj = p.adjust(
        p_raw,
        method = "holm"
      ),
      label = paste0(
        "p=",
        signif(p_adj, 2),
        " (",
        p_to_stars(p_adj),
        ")"
      )
    )
}

add_bracket <- function(
    plot,
    group1,
    group2,
    y,
    label,
    stage_levels,
    y_top
) {

  x1 <- match(
    group1,
    stage_levels
  )

  x2 <- match(
    group2,
    stage_levels
  )

  bracket_height <- max(
    y_top * 0.03,
    y * 0.012
  )

  plot +
    annotate(
      "segment",
      x = x1,
      xend = x1,
      y = y - bracket_height,
      yend = y
    ) +
    annotate(
      "segment",
      x = x2,
      xend = x2,
      y = y - bracket_height,
      yend = y
    ) +
    annotate(
      "segment",
      x = x1,
      xend = x2,
      y = y,
      yend = y
    ) +
    annotate(
      "text",
      x = mean(c(x1, x2)),
      y = y + bracket_height * 1.7,
      label = label,
      size = 4
    )
}

save_plot_pair <- function(
    plot,
    output_dir,
    stem,
    width,
    height
) {

  pdf_file <- file.path(
    output_dir,
    paste0(stem, ".pdf")
  )

  png_file <- file.path(
    output_dir,
    paste0(stem, ".png")
  )

  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

plot_specs <- list(
  list(
    label = "C0_Portal_tracts",
    comparison = c("F1", "F4")
  ),
  list(
    label = "C1_Periportal_hepatocytes",
    comparison = c("F0", "F4")
  ),
  list(
    label = "C2_Centrilobular_hepatocytes",
    comparison = c("F0", "F4")
  ),
  list(
    label = "C3_Immune_cells",
    comparison = c("F0", "F3")
  ),
  list(
    label = "C4_Myofibroblasts",
    comparison = c("F0", "F4")
  ),
  list(
    label = "C5_Cholangiocytes",
    comparison = c("F0", "F4")
  ),
  list(
    label = "C6_Matrix_rich_fibroblast",
    comparison = c("F0", "F4")
  ),
  list(
    label = "C4C6_Myofibroblast_MatrixRich",
    comparison = c("F0", "F4")
  )
)

all_selected_tests <- list()
all_pairwise_tests <- list()

for (spec in plot_specs) {

  label <- spec$label

  cluster_number <- sub("_.*", "", label)
  if (cluster_number == "C4C6") cluster_number <- "C4+C6"

  input_file <- file.path(
    input_dir,
    paste0(
      label,
      "_donor_level_proportions.csv"
    )
  )

  if (!file.exists(input_file)) {
    stop(
      "Input file not found: ",
      input_file
    )
  }

  donor_data <- read.csv(
    input_file,
    stringsAsFactors = FALSE
  )

  required_columns <- c(
    "F_score",
    "donor",
    "prop"
  )

  if (!all(
    required_columns %in% colnames(donor_data)
  )) {
    stop(
      "Missing required columns in ",
      input_file
    )
  }

  donor_data <- donor_data %>%
    mutate(
      F_score = factor(
        as.character(F_score),
        levels = stage_levels
      )
    ) %>%
    filter(
      !is.na(F_score),
      is.finite(prop)
    )

  stage_summary <- donor_data %>%
    group_by(F_score) %>%
    summarise(
      mean_prop = mean(
        prop,
        na.rm = TRUE
      ),
      sd_prop = sd(
        prop,
        na.rm = TRUE
      ),
      n = n(),
      se_prop = sd_prop / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      F_score = factor(
        as.character(F_score),
        levels = stage_levels
      )
    )

  pairwise_table <- make_pairwise_table(
    donor_data
  )

  pairwise_table$cluster_label <- label

  all_pairwise_tests[[label]] <- pairwise_table

  tops <- stage_summary %>%
    mutate(
      top = mean_prop +
        ifelse(
          is.na(se_prop),
          0,
          se_prop
        )
    )

  y_top <- max(
    tops$top,
    na.rm = TRUE
  )

  if (!is.finite(y_top) || y_top <= 0) {
    y_top <- 0.01
  }

  if (identical(
    label,
    "C4_Myofibroblasts"
  )) {

    selected_tests <- data.frame(
      group1 = c("F0", "F0"),
      group2 = c("F1", "F2", "F3", "F4"),
      stringsAsFactors = FALSE
    )

    selected_tests$p_raw <- mapply(
      function(group1, group2) {
        run_wilcoxon(
          donor_data,
          group1,
          group2
        )
      },
      selected_tests$group1,
      selected_tests$group2
    )

    selected_tests$p_adj <- p.adjust(
      selected_tests$p_raw,
      method = "holm"
    )

    selected_tests$significant <- (
      selected_tests$p_adj < 0.05
    )

    selected_tests$p_label <- paste0(
      "p=",
      signif(
        selected_tests$p_adj,
        2
      ),
      " (",
      p_to_stars(
        selected_tests$p_adj
      ),
      ")"
    )

    selected_tests$correction <- paste0(
      "Holm correction across ",
      nrow(selected_tests),
      " displayed C4 comparisons"
    )

    selected_tests$cluster_label <- label

    comparisons_to_plot <- selected_tests %>%
      filter(significant)

    output_stem <- paste0(
      label,
      "_Healthy_vs_F1_F2_F3_and_F4_wilcox"
    )

  } else {

    selected_tests <- data.frame(
      group1 = spec$comparison[1],
      group2 = spec$comparison[2],
      stringsAsFactors = FALSE
    )

    selected_tests$p_raw <- run_wilcoxon(
      donor_data,
      selected_tests$group1,
      selected_tests$group2
    )

    selected_tests$p_adj <- NA_real_

    selected_tests$significant <- (
      selected_tests$p_raw < 0.05
    )

    selected_tests$p_label <- paste0(
      "p=",
      signif(
        selected_tests$p_raw,
        2
      ),
      " (",
      p_to_stars(
        selected_tests$p_raw
      ),
      ")"
    )

    selected_tests$correction <- (
      "Original selected comparison: unadjusted"
    )

    selected_tests$cluster_label <- label

    comparisons_to_plot <- selected_tests

    comparison_filename <- ifelse(
      spec$comparison == "F0",
      "Healthy",
      spec$comparison
    )

    output_stem <- paste0(
      label,
      "_",
      paste(
        comparison_filename,
        collapse = "_vs_"
      ),
      "_wilcox"
    )
  }

  all_selected_tests[[label]] <- selected_tests

  if (nrow(comparisons_to_plot) > 0) {

    comparisons_to_plot$y_position <- (
      y_top *
        (
          1.12 +
            0.22 *
            (
              seq_len(
                nrow(comparisons_to_plot)
              ) -
                1
            )
        )
    )

    label_top <- max(
      comparisons_to_plot$y_position
    ) + y_top * 0.10

    y_max <- label_top + y_top * 0.10

  } else {

    y_max <- y_top * 1.18
  }

  y_max <- min(
    1,
    y_max
  )

  p_main <- ggplot(
    stage_summary,
    aes(
      x = F_score,
      y = mean_prop,
      fill = F_score
    )
  ) +
    geom_col(
      width = 0.65
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(
          mean_prop -
            ifelse(
              is.na(se_prop),
              0,
              se_prop
            ),
          0
        ),
        ymax = pmin(
          mean_prop +
            ifelse(
              is.na(se_prop),
              0,
              se_prop
            ),
          1
        )
      ),
      width = 0.18
    ) +
    scale_fill_manual(
      values = stage_palette,
      guide = "none",
      drop = FALSE
    ) +
    scale_x_discrete(
      labels = stage_labels,
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      limits = c(
        0,
        y_max
      ),
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = paste0(
        "n(",
        cluster_number,
        " spots) / n(total spots)"
      )
    ) +
    theme_minimal(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.title = element_text(
        face = "bold"
      ),
      panel.grid.major.y = element_line(
        color = "gray90"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(
        size = 12
      ),
      axis.text.y = element_text(
        size = 11
      )
    )

  if (nrow(comparisons_to_plot) > 0) {

    for (
      comparison_index in seq_len(
        nrow(comparisons_to_plot)
      )
    ) {

      p_main <- add_bracket(
        plot = p_main,
        group1 = comparisons_to_plot$group1[
          comparison_index
        ],
        group2 = comparisons_to_plot$group2[
          comparison_index
        ],
        y = comparisons_to_plot$y_position[
          comparison_index
        ],
        label = comparisons_to_plot$p_label[
          comparison_index
        ],
        stage_levels = stage_levels,
        y_top = y_top
      )
    }
  }

  save_plot_pair(
    plot = p_main,
    output_dir = output_dir,
    stem = output_stem,
    width = 4.6,
    height = 4.2
  )

  adjacent_table <- pairwise_table %>%
    mutate(
      idx1 = match(
        group1,
        stage_levels
      ),
      idx2 = match(
        group2,
        stage_levels
      )
    ) %>%
    filter(
      idx2 - idx1 == 1
    ) %>%
    arrange(
      idx1
    )

  if (nrow(adjacent_table) > 0) {

    tops_adjacent <- stage_summary %>%
      mutate(
        top = mean_prop +
          ifelse(
            is.na(se_prop),
            0,
            se_prop
          )
      ) %>%
      select(
        F_score,
        top
      )

    adjacent_y_base <- vapply(
      2:length(stage_levels),
      function(stage_index) {
        max(
          tops_adjacent$top[
            tops_adjacent$F_score %in%
              c(
                stage_levels[
                  stage_index - 1
                ],
                stage_levels[
                  stage_index
                ]
              )
          ],
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    adjacent_y_positions <- (
      adjacent_y_base * 1.08 +
        (
          seq_along(
            adjacent_y_base
          ) -
            1
        ) *
        0.035
    )

    adjacent_table$y_position <- pmin(
      adjacent_y_positions,
      1
    )

    y_max_adjacent <- max(
      adjacent_table$y_position,
      na.rm = TRUE
    ) * 1.08

    if (
      !is.finite(y_max_adjacent) ||
        y_max_adjacent <= 0
    ) {
      y_max_adjacent <- y_max
    }

    p_adjacent <- ggplot(
      stage_summary,
      aes(
        x = F_score,
        y = mean_prop,
        fill = F_score
      )
    ) +
      geom_col(
        width = 0.65
      ) +
      geom_errorbar(
        aes(
          ymin = pmax(
            mean_prop -
              ifelse(
                is.na(se_prop),
                0,
                se_prop
              ),
            0
          ),
          ymax = pmin(
            mean_prop +
              ifelse(
                is.na(se_prop),
                0,
                se_prop
              ),
            1
          )
        ),
        width = 0.18
      ) +
      scale_fill_manual(
        values = stage_palette,
        guide = "none",
        drop = FALSE
      ) +
      scale_x_discrete(
        labels = stage_labels,
        drop = FALSE
      ) +
      scale_y_continuous(
        labels = percent_format(
          accuracy = 1
        ),
        limits = c(
          0,
          min(
            1,
            y_max_adjacent
          )
        ),
        expand = expansion(
          mult = c(
            0,
            0.02
          )
        )
      ) +
      labs(
        title = NULL,
        x = "F score",
        y = paste0(
          cluster_number,
          " spots (proportion per donor)"
        )
      ) +
      theme_minimal(
        base_size = 14
      ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        ),
        axis.title = element_text(
          face = "bold"
        ),
        panel.grid.major.y = element_line(
          color = "gray90"
        ),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )

    for (
      bracket_index in seq_len(
        nrow(adjacent_table)
      )
    ) {

      p_adjacent <- add_bracket(
        plot = p_adjacent,
        group1 = adjacent_table$group1[
          bracket_index
        ],
        group2 = adjacent_table$group2[
          bracket_index
        ],
        y = adjacent_table$y_position[
          bracket_index
        ],
        label = adjacent_table$label[
          bracket_index
        ],
        stage_levels = stage_levels,
        y_top = y_top
      )
    }

    save_plot_pair(
      plot = p_adjacent,
      output_dir = output_dir,
      stem = paste0(
        label,
        "_adjacent_stage_wilcox_Healthy"
      ),
      width = 5.5,
      height = 4.2
    )
  }
}

selected_results <- bind_rows(
  all_selected_tests
)

pairwise_results <- bind_rows(
  all_pairwise_tests
)

write.csv(
  selected_results,
  file = file.path(
    output_dir,
    "ALL_selected_comparisons_Healthy_C4_Holm.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_results,
  file = file.path(
    output_dir,
    "ALL_pairwise_Wilcoxon_tests_Holm_corrected.csv"
  ),
  row.names = FALSE
)

c4_results <- selected_results %>%
  filter(
    cluster_label == "C4_Myofibroblasts"
  )

write.csv(
  c4_results,
  file = file.path(
    output_dir,
    "C4_Myofibroblasts_Healthy_vs_F1_F2_F3_and_F4_tests.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("C4 HEALTHY COMPARISONS\n")
cat("============================================\n")
print(c4_results)

pdf_files <- list.files(
  output_dir,
  pattern = "\\.pdf$",
  full.names = TRUE
)

png_files <- list.files(
  output_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

cat("\nPDF files generated:", length(pdf_files), "\n")
cat("PNG files generated:", length(png_files), "\n")

if (
  length(pdf_files) != 16 ||
    length(png_files) != 16
) {
  stop(
    "Expected 16 PDFs and 16 PNGs."
  )
}

cat("\nAll plots completed successfully.\n")
cat("Output folder:\n")
cat(output_dir, "\n")
