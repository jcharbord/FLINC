#!/usr/bin/env Rscript

# Figure 3A-B and Supplementary Figure 8A-E: warning-free deterministic plotting from the
# retained Step 02 analysis tables. This script performs no new biological
# analysis and does not use historical reference outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(liana)
  library(circlize)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(scales)
  library(grid)
})

options(width = 180)

started_at <- Sys.time()

repo_root <- "/Users/jianwu/Documents/FLINC"
repro_root <- file.path(repo_root, "Jeremie", "Reproducible_code")
output_dir <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "02_main_and_supplementary_figures"
)

selected_rds <- file.path(output_dir, "02_selected_interactions_top60.rds")
row_z_csv <- file.path(output_dir, "02_selected_interactions_row_z_scores.csv")
chord_rds <- file.path(output_dir, "02_chord_significant_lr_counts.rds")
cholangiocyte_counts_csv <- file.path(
  output_dir,
  "02_cholangiocyte_partner_counts.csv"
)
dotplot_paths <- setNames(
  file.path(
    output_dir,
    paste0("02_dotplot_input_", c("F1", "F2", "F3", "F4"), "_vs_F0.csv")
  ),
  c("F1", "F2", "F3", "F4")
)
step02_checksums_tsv <- file.path(output_dir, "02_OUTPUT_SHA256.tsv")
step02_completion_tsv <- file.path(output_dir, "02_COMPLETION.tsv")

stage_levels <- c("F0", "F1", "F2", "F3", "F4")
cell_types <- c(
  "Bcells", "CentralHep", "Cholangiocyte", "Erythroid",
  "InfMac", "InterHep", "NKTcell", "NonInfMac",
  "PortalEndo", "PortalHep", "Stellate", "cvLSECs"
)

expected_selected_interactions <- 60L
expected_significant_interactions <- 29205L
expected_cholangiocyte_partner_rows <- 23L
expected_plot_count <- 7L

required_packages <- c(
  "dplyr", "ggplot2", "liana", "circlize", "ComplexHeatmap",
  "RColorBrewer", "scales", "grid", "digest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

required_inputs <- c(
  selected_rds,
  row_z_csv,
  chord_rds,
  cholangiocyte_counts_csv,
  unname(dotplot_paths),
  step02_checksums_tsv,
  step02_completion_tsv
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required Step 02 inputs: ", paste(missing_inputs, collapse = ", "))
}

if (!dir.exists(output_dir)) {
  stop("Step 02 output directory does not exist: ", output_dir)
}

plot_paths <- c(
  Figure3A = file.path(
    output_dir,
    "03_Figure3A_chord_LIANA_aggregate_rank_le_0.01.pdf"
  ),
  Figure3B = file.path(
    output_dir,
    "03_Figure3B_stage_activity_heatmap_with_table.pdf"
  ),
  SupplementaryFigure8A = file.path(
    output_dir,
    "03_SuppFigure8A_cholangiocyte_partner_counts.pdf"
  ),
  SupplementaryFigure8B = file.path(
    output_dir,
    "03_SuppFigure8B_dotplot_F1_vs_F0.pdf"
  ),
  SupplementaryFigure8C = file.path(
    output_dir,
    "03_SuppFigure8C_dotplot_F2_vs_F0.pdf"
  ),
  SupplementaryFigure8D = file.path(
    output_dir,
    "03_SuppFigure8D_dotplot_F3_vs_F0.pdf"
  ),
  SupplementaryFigure8E = file.path(
    output_dir,
    "03_SuppFigure8E_dotplot_F4_vs_F0.pdf"
  )
)

metadata_paths <- c(
  plot_manifest_csv = file.path(output_dir, "03_plot_manifest.csv"),
  parameters_tsv = file.path(output_dir, "03_plot_parameters.tsv"),
  input_manifest_tsv = file.path(output_dir, "03_input_manifest.tsv"),
  session_info_txt = file.path(output_dir, "03_sessionInfo.txt")
)

planned_outputs <- c(plot_paths, metadata_paths)
blocked_outputs <- unname(planned_outputs)[
  file.exists(unname(planned_outputs)) |
    file.exists(paste0(unname(planned_outputs), ".partial"))
]

if (length(blocked_outputs)) {
  stop(
    "Refusing to replace existing Step 03 outputs: ",
    paste(blocked_outputs, collapse = ", ")
  )
}

atomic_write_csv <- function(object, path, row.names = FALSE) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  utils::write.csv(object, partial, row.names = row.names)
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
}

atomic_write_tsv <- function(object, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  utils::write.table(
    object,
    partial,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
}

atomic_write_lines <- function(lines, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
}

atomic_base_pdf <- function(path, width, height, draw_function) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing plot: ", path)
  }

  device_open <- FALSE
  tryCatch(
    {
      grDevices::pdf(partial, width = width, height = height)
      device_open <- TRUE
      draw_function()
      grDevices::dev.off()
      device_open <- FALSE
    },
    finally = {
      if (device_open) {
        try(grDevices::dev.off(), silent = TRUE)
      }
    }
  )

  if (!file.exists(partial) || file.info(partial)$size <= 0) {
    stop("PDF device did not create a non-empty plot: ", path)
  }
  if (!file.rename(partial, path)) {
    stop("Could not finalize plot: ", path)
  }
}

atomic_ggsave_pdf <- function(path, plot, width, height) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing plot: ", path)
  }

  ggplot2::ggsave(
    filename = partial,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )

  if (!file.exists(partial) || file.info(partial)$size <= 0) {
    stop("ggsave did not create a non-empty plot: ", path)
  }
  if (!file.rename(partial, path)) {
    stop("Could not finalize plot: ", path)
  }
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
} else {
  NA_character_
}

message("Figure 3 and Supplementary Figure 8 plotting started at ", format(started_at))
message("Verifying Step 02 completion and checksums...")

step02_completion <- utils::read.delim(
  step02_completion_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

step02_status <- step02_completion$value[step02_completion$field == "status"]
if (length(step02_status) != 1L || step02_status != "complete") {
  stop("Step 02 completion status is not 'complete'.")
}

step02_checksums <- utils::read.delim(
  step02_checksums_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("file", "bytes", "sha256") %in% names(step02_checksums))) {
  stop("Step 02 checksum manifest has unexpected columns.")
}

checksum_verified_inputs <- c(
  selected_rds,
  row_z_csv,
  chord_rds,
  cholangiocyte_counts_csv,
  unname(dotplot_paths)
)

for (path in checksum_verified_inputs) {
  expected <- step02_checksums$sha256[
    step02_checksums$file == basename(path)
  ]
  if (length(expected) != 1L) {
    stop("No unique Step 02 checksum found for: ", basename(path))
  }
  if (!identical(sha256_file(path), expected)) {
    stop("Step 02 checksum mismatch for: ", path)
  }
}

message("Loading retained plotting inputs...")

selected_interactions <- readRDS(selected_rds)
saved_row_z <- utils::read.csv(
  row_z_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
chord_counts <- readRDS(chord_rds)
cholangiocyte_counts <- utils::read.csv(
  cholangiocyte_counts_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dotplot_inputs <- lapply(
  dotplot_paths,
  utils::read.csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_selected_columns <- c(
  "source", "target", "lig_gene", "rec_gene",
  "ligand.complex", "receptor.complex", stage_levels, "rowlab"
)
missing_selected_columns <- setdiff(
  required_selected_columns,
  names(selected_interactions)
)
if (length(missing_selected_columns)) {
  stop(
    "Selected-interaction table lacks required columns: ",
    paste(missing_selected_columns, collapse = ", ")
  )
}

if (nrow(selected_interactions) != expected_selected_interactions) {
  stop("Expected 60 selected interactions.")
}

if (!is.matrix(chord_counts) ||
    !identical(dim(chord_counts), c(length(cell_types), length(cell_types))) ||
    !setequal(rownames(chord_counts), cell_types) ||
    !setequal(colnames(chord_counts), cell_types)) {
  stop("Chord-count matrix does not have the expected 12 x 12 cell-type structure.")
}

if (sum(chord_counts) != expected_significant_interactions) {
  stop("Chord counts do not sum to 29,205 significant interactions.")
}

if (!all(c("direction", "partner", "n_lr_pairs") %in%
         names(cholangiocyte_counts))) {
  stop("Cholangiocyte partner table has unexpected columns.")
}

if (nrow(cholangiocyte_counts) != expected_cholangiocyte_partner_rows) {
  stop("Unexpected cholangiocyte partner-table row count.")
}

for (stage in names(dotplot_inputs)) {
  required_dotplot_columns <- c(
    "source", "target", "ligand.complex", "receptor.complex",
    stage, paste0("delta_", stage, "_vs_F0"), "comparison", "delta_pos"
  )
  missing_dotplot_columns <- setdiff(
    required_dotplot_columns,
    names(dotplot_inputs[[stage]])
  )
  if (length(missing_dotplot_columns)) {
    stop(
      "Dot-plot input for ", stage, " lacks required columns: ",
      paste(missing_dotplot_columns, collapse = ", ")
    )
  }
  if (nrow(dotplot_inputs[[stage]]) != expected_selected_interactions) {
    stop("Dot-plot input for ", stage, " does not contain 60 interactions.")
  }
  if (any(dotplot_inputs[[stage]]$delta_pos < 0, na.rm = TRUE)) {
    stop("Dot-plot positive deltas contain a negative value for ", stage, ".")
  }
}

message("Preparing Figure 3B heatmap input...")

stage_matrix <- as.matrix(
  selected_interactions[, stage_levels, drop = FALSE]
)
row_z_matrix <- t(scale(t(stage_matrix)))
row_z_matrix[is.na(row_z_matrix)] <- 0

row_z_keys <- c("source", "target", "ligand", "receptor", "rowlab")
if (!all(c(row_z_keys, stage_levels) %in% names(saved_row_z))) {
  stop("Saved row-z table lacks required columns.")
}

selected_keys_for_z <- data.frame(
  source = as.character(selected_interactions$source),
  target = as.character(selected_interactions$target),
  ligand = as.character(selected_interactions$lig_gene),
  receptor = as.character(selected_interactions$rec_gene),
  rowlab = as.character(selected_interactions$rowlab),
  stringsAsFactors = FALSE
)

if (!identical(selected_keys_for_z, saved_row_z[, row_z_keys, drop = FALSE])) {
  stop("Saved row-z annotation keys are inconsistent with selected interactions.")
}

saved_row_z_matrix <- as.matrix(saved_row_z[, stage_levels, drop = FALSE])
if (max(abs(row_z_matrix - saved_row_z_matrix)) > 1e-12) {
  stop("Saved and recomputed row-z values are inconsistent.")
}

truncate9 <- function(x) {
  substr(as.character(x), 1, 9)
}

heatmap_labels <- selected_interactions %>%
  dplyr::transmute(
    Source = truncate9(source),
    Target = truncate9(target),
    Ligand = as.character(lig_gene),
    Receptor = as.character(rec_gene)
  )

right_annotation <- ComplexHeatmap::rowAnnotation(
  Source = ComplexHeatmap::anno_text(
    heatmap_labels$Source,
    just = "left",
    gp = grid::gpar(fontsize = 8)
  ),
  Target = ComplexHeatmap::anno_text(
    heatmap_labels$Target,
    just = "left",
    gp = grid::gpar(fontsize = 8)
  ),
  Ligand = ComplexHeatmap::anno_text(
    heatmap_labels$Ligand,
    just = "left",
    gp = grid::gpar(fontsize = 8)
  ),
  Receptor = ComplexHeatmap::anno_text(
    heatmap_labels$Receptor,
    just = "left",
    gp = grid::gpar(fontsize = 8)
  ),
  annotation_name_side = "top",
  annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  annotation_width = grid::unit(c(20, 20, 20, 24), "mm")
)

heatmap_colors <- grDevices::colorRampPalette(
  rev(RColorBrewer::brewer.pal(7, "RdYlBu"))
)(100)
heatmap_values <- seq(
  min(row_z_matrix, na.rm = TRUE),
  max(row_z_matrix, na.rm = TRUE),
  length.out = length(heatmap_colors)
)
heatmap_color_function <- circlize::colorRamp2(
  heatmap_values,
  heatmap_colors
)

heatmap_object <- ComplexHeatmap::Heatmap(
  row_z_matrix,
  name = "row z",
  col = heatmap_color_function,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_title = "Stage activity (row z-score)",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  column_names_gp = grid::gpar(fontsize = 10),
  right_annotation = right_annotation,
  gap = grid::unit(1.5, "mm"),
  heatmap_legend_param = list(
    legend_height = grid::unit(100, "mm")
  )
)

message("Plotting Figure 3A chord diagram...")

chord_order <- sort(cell_types)
chord_plot_matrix <- chord_counts[
  chord_order,
  chord_order,
  drop = FALSE
]
diag(chord_plot_matrix) <- pmax(diag(chord_plot_matrix), 1e-3)
chord_colors <- setNames(
  scales::hue_pal()(length(chord_order)),
  chord_order
)

atomic_base_pdf(
  plot_paths[["Figure3A"]],
  width = 9,
  height = 9,
  draw_function = function() {
    circlize::circos.clear()
    circlize::chordDiagramFromMatrix(
      chord_plot_matrix,
      order = chord_order,
      grid.col = chord_colors,
      transparency = 0.4
    )
  }
)

message("Plotting Figure 3B stage-activity heatmap...")

atomic_base_pdf(
  plot_paths[["Figure3B"]],
  width = 14,
  height = 8,
  draw_function = function() {
    ComplexHeatmap::draw(
      heatmap_object,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
  }
)

message("Plotting Supplementary Figure 8A cholangiocyte partner counts...")

cholangiocyte_partner_plot <- cholangiocyte_counts %>%
  dplyr::group_by(direction) %>%
  dplyr::arrange(desc(n_lr_pairs), .by_group = TRUE) %>%
  dplyr::mutate(
    partner = factor(partner, levels = unique(partner))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = partner, y = n_lr_pairs)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~direction, scales = "free_y") +
  ggplot2::labs(
    x = NULL,
    y = "# significant LR pairs (LIANA aggregate_rank <= 0.01)",
    title = "Cholangiocyte communication partners (cell-type resolved)"
  ) +
  ggplot2::theme_classic(base_size = 12)

atomic_base_pdf(
  plot_paths[["SupplementaryFigure8A"]],
  width = 14,
  height = 4,
  draw_function = function() {
    print(cholangiocyte_partner_plot)
  }
)

make_dotplot <- function(stage, output_path) {
  data <- dotplot_inputs[[stage]]

  plot <- data %>%
    liana::liana_dotplot(
      source_groups = sort(unique(data$source)),
      target_groups = sort(unique(data$target)),
      ntop = nrow(data),
      magnitude = stage,
      specificity = "delta_pos"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 10)
    )

  atomic_ggsave_pdf(
    output_path,
    plot = plot,
    width = 24,
    height = 7
  )
}

message("Plotting Supplementary Figure 8B-E stage-versus-F0 dot plots...")

make_dotplot("F1", plot_paths[["SupplementaryFigure8B"]])
make_dotplot("F2", plot_paths[["SupplementaryFigure8C"]])
make_dotplot("F3", plot_paths[["SupplementaryFigure8D"]])
make_dotplot("F4", plot_paths[["SupplementaryFigure8E"]])

if (length(plot_paths) != expected_plot_count ||
    any(!file.exists(unname(plot_paths))) ||
    any(file.info(unname(plot_paths))$size <= 0)) {
  stop("The expected seven non-empty PDF panels were not all generated.")
}

plot_manifest <- data.frame(
  figure = c(
    "Figure 3", "Figure 3",
    rep("Supplementary Figure 8", 5)
  ),
  panel = c("A", "B", "A", "B", "C", "D", "E"),
  content = c(
    "LIANA significant-interaction chord diagram",
    "Stage-activity row-z heatmap with LR annotation table",
    "Cholangiocyte incoming and outgoing significant LR counts",
    "Selected LR activity: F1 versus F0",
    "Selected LR activity: F2 versus F0",
    "Selected LR activity: F3 versus F0",
    "Selected LR activity: F4 versus F0"
  ),
  filename = basename(unname(plot_paths)),
  width_inches = c(9, 14, 14, 24, 24, 24, 24),
  height_inches = c(9, 8, 4, 7, 7, 7, 7),
  pages = 1L,
  stringsAsFactors = FALSE
)

parameter_table <- data.frame(
  parameter = c(
    "figure3a_chord_celltype_order",
    "figure3a_chord_palette",
    "figure3a_chord_transparency",
    "figure3a_minimum_diagonal_plot_value",
    "figure3b_row_scaling",
    "figure3b_row_clustering",
    "figure3b_column_clustering",
    "figure3b_palette",
    "figure3b_source_target_max_characters",
    "suppfigure8a_bar_color",
    "suppfigure8b_e_color",
    "suppfigure8b_e_size",
    "suppfigure8b_e_negative_deltas",
    "pdf_panels"
  ),
  value = c(
    paste(chord_order, collapse = ";"),
    "scales::hue_pal()",
    "0.4",
    "0.001",
    "z-score within each selected LR interaction across F0-F4",
    "hierarchical; ComplexHeatmap defaults",
    "FALSE; F0-F4 retained in fixed order",
    "100-color reverse RdYlBu interpolation",
    "9",
    "ggplot2 default discrete grey",
    "activity at the indicated stage",
    "positive activity change from F0",
    "set to zero by Step 02",
    length(plot_paths)
  ),
  stringsAsFactors = FALSE
)

input_manifest <- data.frame(
  role = c(
    "selected_top60_interactions",
    "selected_interaction_row_z_scores",
    "significant_lr_chord_count_matrix",
    "cholangiocyte_partner_counts",
    paste0("dotplot_input_", names(dotplot_paths), "_vs_F0"),
    "step02_output_checksum_manifest",
    "step02_completion_record"
  ),
  path = normalizePath(required_inputs, mustWork = TRUE),
  bytes = unname(file.info(required_inputs)$size),
  modified = format(file.info(required_inputs)$mtime, tz = "UTC", usetz = TRUE),
  sha256 = vapply(required_inputs, sha256_file, character(1)),
  stringsAsFactors = FALSE
)

session_lines <- capture.output(sessionInfo())

atomic_write_csv(plot_manifest, metadata_paths[["plot_manifest_csv"]])
atomic_write_tsv(parameter_table, metadata_paths[["parameters_tsv"]])
atomic_write_tsv(input_manifest, metadata_paths[["input_manifest_tsv"]])
atomic_write_lines(session_lines, metadata_paths[["session_info_txt"]])

output_checksums <- data.frame(
  file = basename(unname(c(plot_paths, metadata_paths))),
  bytes = unname(file.info(unname(c(plot_paths, metadata_paths)))$size),
  sha256 = vapply(
    unname(c(plot_paths, metadata_paths)),
    sha256_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

checksum_path <- file.path(output_dir, "03_OUTPUT_SHA256.tsv")
atomic_write_tsv(output_checksums, checksum_path)

finished_at <- Sys.time()

completion_table <- data.frame(
  field = c(
    "status",
    "started_at",
    "finished_at",
    "elapsed_minutes",
    "script_path",
    "script_sha256",
    "pdf_panels",
    "figure3_panels",
    "supplementary_figure8_panels",
    "selected_interactions",
    "chord_significant_interactions",
    "cholangiocyte_partner_rows",
    "dotplot_rows_per_stage"
  ),
  value = c(
    "complete",
    format(started_at, tz = "UTC", usetz = TRUE),
    format(finished_at, tz = "UTC", usetz = TRUE),
    sprintf("%.3f", as.numeric(difftime(finished_at, started_at, units = "mins"))),
    script_path,
    if (is.na(script_path)) NA_character_ else sha256_file(script_path),
    length(plot_paths),
    2L,
    5L,
    nrow(selected_interactions),
    sum(chord_counts),
    nrow(cholangiocyte_counts),
    paste(
      names(dotplot_inputs),
      vapply(dotplot_inputs, nrow, integer(1)),
      sep = ":",
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

completion_path <- file.path(output_dir, "03_COMPLETION.tsv")
atomic_write_tsv(completion_table, completion_path)

message("Figure 3 and Supplementary Figure 8 plotting completed at ", format(finished_at))
message("PDF panels generated: ", length(plot_paths))
message("Outputs: ", output_dir)
