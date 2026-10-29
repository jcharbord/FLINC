# QC summaries and supplementary violin plots from the final FLINC object.
# In the R session containing the final annotated object:
#   source("scripts/figure1/07_summarize_spatial_QC.R")
#   qc_result <- run_flinc_qc(obj, output_parent = "results")
# No object slots or metadata are modified.

run_flinc_qc <- function(
    object,
    output_parent = "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code/results",
    sample_col = "donor",
    expected_spots = 26056L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    stop("This plot requires ggplot2 and scales in the current R environment.")
  }

  metrics <- c("nFeature_Spatial", "nCount_Spatial")
  md <- object[[]]
  needed <- c(sample_col, metrics)
  missing <- setdiff(needed, colnames(md))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(md), collapse = ", "))
  }
  if (nrow(md) != ncol(object) ||
      !identical(rownames(md), colnames(object))) {
    stop("Metadata rows do not match the object spot order.")
  }
  if (!is.null(expected_spots) && nrow(md) != expected_spots) {
    stop("Found ", nrow(md), " spots; the manuscript object had ",
         expected_spots, ". Check that this is the intended final object.")
  }
  for (metric in metrics) {
    if (!is.numeric(md[[metric]]) || any(!is.finite(md[[metric]])) ||
        any(md[[metric]] < 0)) {
      stop(metric, " must contain finite, non-negative numeric values.")
    }
  }

  sample_ids <- as.character(md[[sample_col]])
  sample_ids <- sub("^X(?=[1-6][A-D]$)", "", sample_ids, perl = TRUE)
  sample_order <- unlist(lapply(1:6, function(i) paste0(i, LETTERS[1:4])))
  if (anyNA(sample_ids) || !setequal(unique(sample_ids), sample_order)) {
    stop("Expected sample IDs 1A-6D in ", sample_col, ". Found: ",
         paste(sort(unique(sample_ids)), collapse = ", "),
         "\nCheck table(object[[\"", sample_col, "\"]]).")
  }
  qc <- data.frame(
    barcode = rownames(md),
    sample_id = factor(sample_ids, levels = sample_order),
    nFeature_Spatial = md$nFeature_Spatial,
    nCount_Spatial = md$nCount_Spatial
  )

  # Pooled summaries describe the distributions of all retained spatial spots.
  overall <- do.call(rbind, lapply(metrics, function(metric) {
    x <- qc[[metric]]
    data.frame(metric = metric, n_spots = length(x), mean = mean(x),
               median = median(x), q25 = unname(quantile(x, 0.25)),
               q75 = unname(quantile(x, 0.75)), min = min(x), max = max(x))
  }))
  by_sample <- do.call(rbind, lapply(sample_order, function(id) {
    d <- qc[qc$sample_id == id, , drop = FALSE]
    data.frame(sample_id = id, n_spots = nrow(d),
               mean_genes = mean(d$nFeature_Spatial),
               median_genes = median(d$nFeature_Spatial),
               mean_counts = mean(d$nCount_Spatial),
               median_counts = median(d$nCount_Spatial))
  }))

  # These alternative summaries weight samples equally; they are not pooled
  # per-spot means or medians. Keep the definitions explicit in the manuscript.
  definitions <- data.frame(
    metric = metrics,
    pooled_spot_mean = overall$mean,
    pooled_spot_median = overall$median,
    mean_of_sample_means = c(mean(by_sample$mean_genes), mean(by_sample$mean_counts)),
    mean_of_sample_medians = c(mean(by_sample$median_genes), mean(by_sample$median_counts)),
    median_of_sample_medians = c(median(by_sample$median_genes), median(by_sample$median_counts))
  )
  comparison <- data.frame(
    metric = metrics,
    current_manuscript_value = c(4845, 20064),
    pooled_spot_mean = overall$mean,
    pooled_spot_median = overall$median
  )
  cat("\nSPOTS:", nrow(qc), " SAMPLES:", nlevels(qc$sample_id), "\n")
  cat("\nCHECK AGAINST THE CURRENT MANUSCRIPT VALUES:\n")
  print(comparison, row.names = FALSE, digits = 8)
  cat("\nPOOLED PER-SPOT QC SUMMARIES:\n")
  print(overall, row.names = FALSE, digits = 8)
  cat("\nALTERNATIVE SAMPLE-WEIGHTED DEFINITIONS:\n")
  print(definitions, row.names = FALSE, digits = 8)

  output_dir <- file.path(output_parent, paste0(
    "261029_QC_final_object_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()))
  if (file.exists(output_dir)) stop("Output folder already exists: ", output_dir)
  if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create output folder.")
  write.csv(overall, file.path(output_dir, "QC_overall_spots.csv"), row.names = FALSE)
  write.csv(by_sample, file.path(output_dir, "QC_by_sample.csv"), row.names = FALSE)
  write.csv(definitions, file.path(output_dir, "QC_summary_definitions.csv"), row.names = FALSE)
  write.csv(comparison, file.path(output_dir, "QC_manuscript_value_check.csv"), row.names = FALSE)
  write.csv(qc, file.path(output_dir, "QC_spot_metadata.csv"), row.names = FALSE)

  palette <- setNames(scales::hue_pal()(length(sample_order)), sample_order)
  make_plot <- function(metric, y_label) {
    ggplot2::ggplot(qc, ggplot2::aes(x = sample_id, y = .data[[metric]], fill = sample_id)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.35) +
      ggplot2::stat_summary(fun = stats::median, geom = "point", size = 1.1, color = "black") +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
      ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
      ggplot2::labs(x = "Sample", y = y_label, title = metric) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  }
  plots <- list(
    D = make_plot("nFeature_Spatial", "Detected genes per spot"),
    E = make_plot("nCount_Spatial", "Spatial UMI counts per spot")
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  for (panel in names(plots)) {
    stem <- paste0("SupFig1", panel, "_", metrics[match(panel, names(plots))])
    ggplot2::ggsave(file.path(output_dir, paste0(stem, ".pdf")), plots[[panel]],
                    width = 10, height = 4, units = "in", device = pdf_device, bg = "white")
    ggplot2::ggsave(file.path(output_dir, paste0(stem, "_600dpi.png")), plots[[panel]],
                    width = 10, height = 4, units = "in", dpi = 600, bg = "white")
  }
  grDevices::pdf(file.path(output_dir, "SupFig1D_E_QC_combined.pdf"),
                width = 10, height = 8, useDingbats = FALSE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
  print(plots$D, vp = grid::viewport(layout.pos.row = 1))
  print(plots$E, vp = grid::viewport(layout.pos.row = 2))
  grDevices::dev.off()

  writeLines(c(
    "Supplementary Figure 1D-E. Spatial sequencing quality metrics across samples.",
    "Violins show the per-spot distributions of detected genes (D) and Spatial UMI counts (E).",
    "Black dots indicate the median within each sample.",
    "No additional spot filtering was applied by this plotting script."
  ), file.path(output_dir, "QC_plot_legend.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  cat("\nSaved tables and plots to:\n", normalizePath(output_dir), "\n", sep = "")
  invisible(list(overall = overall, by_sample = by_sample,
                 definitions = definitions, plots = plots, output_dir = output_dir))
}
