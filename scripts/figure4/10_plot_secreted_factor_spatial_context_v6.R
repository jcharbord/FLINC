#!/usr/bin/env Rscript

# FLINC Figure 4: donor-balanced spatial-cluster and stage context of selected secreted factors.
#
# This script extends, but does not alter, the validated strict secreted-factor
# analysis in 07_run_secreted_factor_analysis.R. It uses the same Spatial-count
# normalization (log1p[count / spot library size * 10,000]) and reuses the
# scale-factor-based Visium geometry used for the validated PROGENy maps.

options(stringsAsFactors = FALSE, digits = 10, scipen = 999)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
set.seed(42L)

fail <- function(...) stop(paste0(...), call. = FALSE)

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table",
  "ggplot2", "cowplot", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  fail("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) fail("Arguments must be supplied as --key value pairs")
arg_list <- as.list(args[seq(2L, length(args), by = 2L)])
names(arg_list) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])

required_args <- c(
  "input", "topgenes", "output-dir", "log-dir", "expected-input-sha256"
)
missing_args <- setdiff(required_args, names(arg_list))
if (length(missing_args) > 0L) {
  fail("Missing arguments: ", paste(missing_args, collapse = ", "))
}

input_path <- normalizePath(arg_list[["input"]], mustWork = TRUE)
topgenes_path <- normalizePath(arg_list[["topgenes"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
expected_input_sha <- tolower(arg_list[["expected-input-sha256"]])
map_images <- if ("map-images" %in% names(arg_list)) {
  trimws(strsplit(arg_list[["map-images"]], ",", fixed = TRUE)[[1L]])
} else if ("map-image" %in% names(arg_list)) {
  arg_list[["map-image"]]
} else {
  c("X1A", "X6D")
}
map_images <- unique(map_images[nzchar(map_images)])
if (length(map_images) != 2L) {
  fail("Exactly two comma-separated map images are required (Healthy F0 and F4)")
}

if (dir.exists(output_dir) || file.exists(output_dir)) {
  fail("Final output path already exists: ", output_dir)
}
if (dir.exists(log_dir) || file.exists(log_dir)) {
  fail("Final log path already exists: ", log_dir)
}

dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_dir), recursive = TRUE, showWarnings = FALSE)

run_stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
result_tmp <- paste0(output_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
log_tmp <- paste0(log_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
dir.create(result_tmp, recursive = FALSE, showWarnings = FALSE)
dir.create(log_tmp, recursive = FALSE, showWarnings = FALSE)

completed <- FALSE
on.exit(
  if (!completed) {
    message("Incomplete run retained: ", result_tmp, "\n", log_tmp)
  },
  add = TRUE
)

atomic_write_table <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(
    x, temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

slot_or_null <- function(x, slot_name) {
  if (slot_name %in% slotNames(x)) slot(x, slot_name) else NULL
}

as_single_numeric <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value)) {
    fail(label, " must be one finite numeric value")
  }
  value
}

started_at <- Sys.time()
input_sha <- sha256(input_path)
topgenes_sha <- sha256(topgenes_path)
if (!identical(tolower(input_sha), expected_input_sha)) {
  fail("Input SHA256 mismatch. Expected ", expected_input_sha,
       "; observed ", input_sha)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

# Prespecified biological selection. Every gene must be present in the validated
# strict top-30 table and in the stated regenerated WGCNA module.
selected_gene_definition <- data.frame(
  gene = c(
    "IGFBP7", "DCN", "CXCL12", "BGN", "C1QC", "TGFBI", "TIMP1", "ADAMTSL2",
    "APOA1", "APOC1", "C3", "TTR", "RBP4", "VTN", "HP", "APOB"
  ),
  module = c(rep("SM_1", 8L), rep("SM_6", 8L)),
  display_group = c(
    rep("SM_1 - fibro-immune ECM remodeling", 8L),
    rep("SM_6 - hepatic plasma/secretory output", 8L)
  ),
  biological_context = c(
    "ECM/stromal remodeling", "ECM organization", "chemokine/stromal signaling",
    "ECM organization", "complement/immune", "matrix-cell interaction",
    "matrix turnover", "ECM-associated candidate marker",
    "lipoprotein output", "lipoprotein/immune-associated output", "complement",
    "transport protein", "retinoid transport", "matrix/adhesion", "acute-phase output",
    "lipoprotein output"
  ),
  stringsAsFactors = FALSE
)
map_genes <- c("TIMP1", "IGFBP7", "ADAMTSL2", "DCN")
main_heatmap_genes <- c(
  "TIMP1", "IGFBP7", "ADAMTSL2", "DCN", "CXCL12", "C1QC",
  "APOA1", "TTR", "RBP4", "HP"
)
main_heatmap_groups <- c(
  "TIMP1" = "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
  "IGFBP7" = "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
  "ADAMTSL2" = "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
  "DCN" = "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
  "CXCL12" = "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
  "C1QC" = "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
  "APOA1" = "Hepatocyte-associated\nsecreted factors (SM_6)",
  "TTR" = "Hepatocyte-associated\nsecreted factors (SM_6)",
  "RBP4" = "Hepatocyte-associated\nsecreted factors (SM_6)",
  "HP" = "Hepatocyte-associated\nsecreted factors (SM_6)"
)

message("Reading validated strict top-30 secreted-factor table...")
topgenes <- data.table::fread(topgenes_path, data.table = FALSE)
if (!all(c("gene", "module", "rank") %in% names(topgenes))) {
  fail("Top-gene table must contain gene, module, and rank columns")
}
topgenes$gene <- toupper(trimws(as.character(topgenes$gene)))
topgenes$module <- gsub("-", "_", trimws(as.character(topgenes$module)))

selection_audit <- merge(
  selected_gene_definition,
  topgenes[, c("gene", "module", "rank", "mean_normalized_expression"), drop = FALSE],
  by = c("gene", "module"), all.x = TRUE, sort = FALSE
)
missing_from_validated_top30 <- selection_audit$gene[is.na(selection_audit$rank)]
if (length(missing_from_validated_top30) > 0L) {
  fail(
    "Selected genes absent from the validated strict top-30 table/module: ",
    paste(missing_from_validated_top30, collapse = ", ")
  )
}

message("Loading final reproducibility object read-only...")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) fail("Input is not a Seurat object")
if (!("Spatial" %in% names(object@assays))) fail("Input lacks the Spatial assay")

metadata <- object[[]]
required_metadata <- c("cluster_named", "donor", "F_score")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0L) {
  fail("Input metadata lacks: ", paste(missing_metadata, collapse = ", "))
}
if (!identical(rownames(metadata), colnames(object))) {
  fail("Metadata row order is not identical to Seurat spot order")
}

cluster_levels <- c(
  "C0 Portal tracts",
  "C1 Periportal hepatocytes",
  "C2 Centrilobular hepatocytes",
  "C3 Immune cells",
  "C4 Myofibroblasts",
  "C5 Cholangiocytes",
  "C6 Matrix-rich fibroblast"
)
observed_clusters <- sort(unique(as.character(metadata$cluster_named)))
if (!setequal(observed_clusters, cluster_levels)) {
  fail(
    "cluster_named does not contain exactly the expected seven labels. Observed: ",
    paste(observed_clusters, collapse = "; ")
  )
}

stage_levels <- c("F0", "F1", "F2", "F3", "F4")
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
observed_stages <- sort(unique(as.character(metadata$F_score)))
if (!setequal(observed_stages, stage_levels)) {
  fail("F_score does not contain exactly F0-F4. Observed: ",
       paste(observed_stages, collapse = ", "))
}

counts <- SeuratObject::LayerData(object[["Spatial"]], layer = "counts")
if (!(inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix"))) {
  fail("Spatial counts are not a sparse matrix")
}
if (!identical(colnames(counts), colnames(object))) {
  fail("Spatial count columns are not identical to Seurat spot order")
}

feature_symbols <- toupper(rownames(counts))
if (anyDuplicated(feature_symbols)) fail("Upper-cased Spatial feature names are duplicated")
rownames(counts) <- feature_symbols

selected_genes <- selected_gene_definition$gene
missing_expression_genes <- setdiff(selected_genes, rownames(counts))
if (length(missing_expression_genes) > 0L) {
  fail("Selected genes absent from Spatial counts: ",
       paste(missing_expression_genes, collapse = ", "))
}

library_sizes <- Matrix::colSums(counts)
if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
  fail("Invalid spot library sizes")
}

selected_counts <- counts[selected_genes, , drop = FALSE]
normalized_selected <- selected_counts %*%
  Matrix::Diagonal(n = ncol(selected_counts), x = 10000 / library_sizes)
normalized_selected <- log1p(normalized_selected)
dimnames(normalized_selected) <- dimnames(selected_counts)

if (!identical(rownames(normalized_selected), rownames(selected_counts)) ||
    !identical(colnames(normalized_selected), colnames(selected_counts))) {
  fail("Normalized selected-gene matrix did not retain gene/spot identities")
}

message("Calculating expression and detection summaries across spatial clusters...")
cluster_rows <- lapply(cluster_levels, function(cluster_name) {
  spot_index <- which(as.character(metadata$cluster_named) == cluster_name)
  if (length(spot_index) == 0L) fail("No spots found for ", cluster_name)
  data.frame(
    gene = selected_genes,
    cluster_named = cluster_name,
    n_spots_cluster = length(spot_index),
    mean_normalized_expression = as.numeric(
      Matrix::rowMeans(normalized_selected[, spot_index, drop = FALSE])
    ),
    pct_expressed = 100 * as.numeric(
      Matrix::rowMeans(selected_counts[, spot_index, drop = FALSE] > 0)
    ),
    stringsAsFactors = FALSE
  )
})
cluster_statistics <- do.call(rbind, cluster_rows)
cluster_statistics <- merge(
  cluster_statistics,
  selected_gene_definition,
  by = "gene", all.x = TRUE, sort = FALSE
)

cluster_statistics$expression_scaled_within_gene <- ave(
  cluster_statistics$mean_normalized_expression,
  cluster_statistics$gene,
  FUN = function(x) {
    current_sd <- stats::sd(x)
    if (!is.finite(current_sd) || current_sd == 0) rep(0, length(x)) else as.numeric(scale(x))
  }
)

display_gene_order <- selected_gene_definition$gene
cluster_statistics$gene <- factor(
  cluster_statistics$gene,
  levels = rev(display_gene_order)
)
cluster_statistics$cluster_named <- factor(
  cluster_statistics$cluster_named,
  levels = cluster_levels
)
cluster_statistics$display_group <- factor(
  cluster_statistics$display_group,
  levels = unique(selected_gene_definition$display_group)
)

message("Calculating donor-balanced expression summaries across spatial clusters...")
donor_cluster_key <- interaction(
  as.character(metadata$donor),
  as.character(metadata$cluster_named),
  drop = TRUE,
  lex.order = TRUE
)
donor_cluster_rows <- lapply(levels(donor_cluster_key), function(group_name) {
  spot_index <- which(donor_cluster_key == group_name)
  donors <- unique(as.character(metadata$donor[spot_index]))
  clusters <- unique(as.character(metadata$cluster_named[spot_index]))
  if (length(donors) != 1L || length(clusters) != 1L) {
    fail("Donor/cluster group is not unique: ", group_name)
  }
  data.frame(
    donor = donors,
    cluster_named = clusters,
    gene = selected_genes,
    n_spots = length(spot_index),
    mean_normalized_expression = as.numeric(
      Matrix::rowMeans(normalized_selected[, spot_index, drop = FALSE])
    ),
    pct_expressed = 100 * as.numeric(
      Matrix::rowMeans(selected_counts[, spot_index, drop = FALSE] > 0)
    ),
    stringsAsFactors = FALSE
  )
})
donor_cluster_statistics <- do.call(rbind, donor_cluster_rows)

donor_balanced_cluster_statistics <- stats::aggregate(
  cbind(mean_normalized_expression, pct_expressed) ~ gene + cluster_named,
  data = donor_cluster_statistics,
  FUN = mean
)
donor_counts_by_cluster <- stats::aggregate(
  donor ~ gene + cluster_named,
  data = donor_cluster_statistics,
  FUN = function(x) length(unique(x))
)
names(donor_counts_by_cluster)[names(donor_counts_by_cluster) == "donor"] <- "n_donors"
donor_balanced_cluster_statistics <- merge(
  donor_balanced_cluster_statistics,
  donor_counts_by_cluster,
  by = c("gene", "cluster_named"),
  all.x = TRUE,
  sort = FALSE
)
donor_balanced_cluster_statistics <- merge(
  donor_balanced_cluster_statistics,
  selected_gene_definition,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
donor_balanced_cluster_statistics$expression_scaled_within_gene <- ave(
  donor_balanced_cluster_statistics$mean_normalized_expression,
  donor_balanced_cluster_statistics$gene,
  FUN = function(x) {
    current_sd <- stats::sd(x)
    if (!is.finite(current_sd) || current_sd == 0) {
      rep(0, length(x))
    } else {
      as.numeric(scale(x))
    }
  }
)

dot_plot <- ggplot(
  cluster_statistics,
  aes(
    x = cluster_named,
    y = gene,
    size = pct_expressed,
    colour = expression_scaled_within_gene
  )
) +
  geom_point(alpha = 0.92) +
  facet_grid(display_group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    name = "Mean expression\nscaled within gene"
  ) +
  scale_size_continuous(
    limits = c(0, 100), range = c(1.2, 8.5),
    breaks = c(25, 50, 75, 100), name = "Spots expressed (%)"
  ) +
  labs(
    title = "Extended selected secreted factors across spatial clusters",
    subtitle = "Colour: relative mean log-normalized expression within each gene; size: percentage of expressing spots",
    x = NULL,
    y = NULL
  ) +
  cowplot::theme_cowplot(font_size = 11) +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "#EEEEEE", colour = NA),
    strip.text.y = element_text(angle = 0, face = "plain"),
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.3),
    legend.position = "right",
    plot.title.position = "plot"
  )

ggsave(
  file.path(result_tmp, "SupFig_extended_secreted_factors_by_spatial_cluster.pdf"),
  dot_plot, width = 12.0, height = 8.5, units = "in", useDingbats = FALSE
)
ggsave(
  file.path(result_tmp, "SupFig_extended_secreted_factors_by_spatial_cluster.png"),
  dot_plot, width = 12.0, height = 8.5, units = "in", dpi = 300
)

message("Calculating donor-level stage summaries for key mapped factors...")
donor_stage_key <- interaction(
  as.character(metadata$donor),
  as.character(metadata$F_score),
  drop = TRUE,
  lex.order = TRUE
)
donor_stage_rows <- lapply(levels(donor_stage_key), function(group_name) {
  spot_index <- which(donor_stage_key == group_name)
  donors <- unique(as.character(metadata$donor[spot_index]))
  stages <- unique(as.character(metadata$F_score[spot_index]))
  if (length(donors) != 1L || length(stages) != 1L) {
    fail("Donor/stage group is not unique: ", group_name)
  }
  data.frame(
    donor = donors,
    F_score = stages,
    gene = selected_genes,
    n_spots = length(spot_index),
    mean_normalized_expression = as.numeric(
      Matrix::rowMeans(normalized_selected[, spot_index, drop = FALSE])
    ),
    pct_expressed = 100 * as.numeric(
      Matrix::rowMeans(selected_counts[, spot_index, drop = FALSE] > 0)
    ),
    stringsAsFactors = FALSE
  )
})
donor_stage_statistics <- do.call(rbind, donor_stage_rows)
donor_stage_statistics$F_score <- factor(
  donor_stage_statistics$F_score,
  levels = stage_levels
)

message("Testing donor-level associations with ordinal fibrosis stage...")
stage_association_rows <- lapply(map_genes, function(current_gene) {
  current <- donor_stage_statistics[
    donor_stage_statistics$gene == current_gene,
    , drop = FALSE
  ]
  ordinal_stage <- match(as.character(current$F_score), stage_levels) - 1L
  test <- suppressWarnings(stats::cor.test(
    ordinal_stage,
    current$mean_normalized_expression,
    method = "spearman",
    exact = FALSE
  ))
  data.frame(
    gene = current_gene,
    n_donors = nrow(current),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})
stage_association_statistics <- do.call(rbind, stage_association_rows)
stage_association_statistics$p_adjust_bh <- stats::p.adjust(
  stage_association_statistics$p_value,
  method = "BH"
)

format_adjusted_q <- function(x) {
  ifelse(x < 0.001, "q < 0.001", sprintf("q = %.3f", x))
}
stage_facet_labels <- setNames(
  paste0(
    stage_association_statistics$gene,
    "\nSpearman r = ", sprintf("%.2f", stage_association_statistics$spearman_rho),
    "; ", format_adjusted_q(stage_association_statistics$p_adjust_bh)
  ),
  stage_association_statistics$gene
)

stage_plot_table <- donor_stage_statistics[
  donor_stage_statistics$gene %in% map_genes,
  , drop = FALSE
]
stage_plot_table$gene <- factor(stage_plot_table$gene, levels = map_genes)
stage_plot_table$gene_panel <- factor(
  unname(stage_facet_labels[as.character(stage_plot_table$gene)]),
  levels = unname(stage_facet_labels[map_genes])
)
stage_plot_table$stage_display <- factor(
  as.character(stage_plot_table$F_score),
  levels = stage_levels,
  labels = unname(stage_labels[stage_levels])
)

stage_plot <- ggplot(
  stage_plot_table,
  aes(
    x = stage_display,
    y = mean_normalized_expression,
    fill = F_score
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    colour = "#444444",
    alpha = 0.88
  ) +
  geom_point(
    aes(group = donor),
    position = position_jitter(width = 0.10, height = 0),
    shape = 21,
    size = 2.0,
    stroke = 0.35,
    alpha = 0.88,
    colour = "#444444"
  ) +
  scale_fill_manual(values = stage_palette, guide = "none", drop = FALSE) +
  facet_wrap(~ gene_panel, scales = "free_y", nrow = 1L) +
  labs(
    title = "H. Selected secreted factors across fibrosis stages",
    subtitle = "Each point represents one donor; values are donor means across spatial spots",
    x = "Fibrosis stage",
    y = "Mean log-normalized expression"
  ) +
  cowplot::theme_cowplot(font_size = 11) +
  theme(strip.background = element_rect(fill = "#EEEEEE", colour = NA))

ggsave(
  file.path(result_tmp, "Fig4H_selected_secreted_factor_donor_stage_profiles.pdf"),
  stage_plot, width = 11.5, height = 3.8, units = "in", useDingbats = FALSE
)
ggsave(
  file.path(result_tmp, "Fig4H_selected_secreted_factor_donor_stage_profiles.png"),
  stage_plot, width = 11.5, height = 3.8, units = "in", dpi = 300
)

message("Rendering donor-balanced ten-marker spatial-cluster heatmap...")
heatmap_table <- donor_balanced_cluster_statistics[
  donor_balanced_cluster_statistics$gene %in% main_heatmap_genes,
  , drop = FALSE
]
heatmap_table$display_group_main <- factor(
  unname(main_heatmap_groups[heatmap_table$gene]),
  levels = c(
    "Fibro-immune remodeling-associated\nsecreted factors (SM_1)",
    "Hepatocyte-associated\nsecreted factors (SM_6)"
  )
)
heatmap_table$gene <- factor(
  heatmap_table$gene,
  levels = rev(main_heatmap_genes)
)
heatmap_table$cluster_named <- factor(
  as.character(heatmap_table$cluster_named), levels = cluster_levels
)

heatmap_plot <- ggplot(
  heatmap_table,
  aes(x = cluster_named, y = gene, fill = expression_scaled_within_gene)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    name = "Mean expression\nscaled within gene"
  ) +
  facet_grid(
    display_group_main ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  labs(
    title = "I. Spatial-niche distribution of selected secreted factors",
    subtitle = paste0(
      "Donor-balanced mean log-normalized expression; ",
      "values scaled separately within each gene"
    ),
    x = NULL,
    y = NULL
  ) +
  cowplot::theme_cowplot(font_size = 11) +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title.position = "plot",
    strip.placement = "outside",
    strip.background = element_rect(fill = "#EEEEEE", colour = NA),
    strip.text.y.left = element_text(
      angle = 90,
      face = "plain",
      size = 9.5,
      lineheight = 1.05,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(t = 5, r = 6, b = 5, l = 6)
    )
  )

ggsave(
  file.path(result_tmp, "Fig4I_selected_secreted_factor_spatial_cluster_heatmap_10markers.pdf"),
  heatmap_plot, width = 10.8, height = 7.0, units = "in", useDingbats = FALSE
)
ggsave(
  file.path(result_tmp, "Fig4I_selected_secreted_factor_spatial_cluster_heatmap_10markers.png"),
  heatmap_plot, width = 10.8, height = 7.0, units = "in", dpi = 300
)

message("Preparing matched Healthy F0 and F4 H&E spatial maps: ",
        paste(map_images, collapse = ", "), "...")
missing_map_images <- setdiff(map_images, names(object@images))
if (length(missing_map_images) > 0L) {
  fail("Requested map image(s) absent: ", paste(missing_map_images, collapse = ", "))
}

section_data <- lapply(map_images, function(map_image) {
  image_object <- object[[map_image]]
  image_cells <- SeuratObject::Cells(image_object)
  metadata_rows <- match(image_cells, rownames(metadata))
  if (anyNA(metadata_rows)) fail("Image cells are absent from metadata for ", map_image)
  image_donors <- unique(as.character(metadata$donor[metadata_rows]))
  image_stages <- unique(as.character(metadata$F_score[metadata_rows]))
  if (length(image_donors) != 1L || length(image_stages) != 1L) {
    fail("Image does not map uniquely to one donor and stage: ", map_image)
  }

  coordinates <- SeuratObject::GetTissueCoordinates(image_object)
  if (!identical(rownames(coordinates), image_cells)) {
    fail("Coordinate order does not match image-cell order for ", map_image)
  }
  if (!all(c("imagerow", "imagecol") %in% colnames(coordinates))) {
    fail("Expected imagerow/imagecol coordinates for ", map_image)
  }

  scale_factors <- slot_or_null(image_object, "scale.factors")
  if (is.null(scale_factors) || !all(c("spot", "lowres") %in% names(scale_factors))) {
    fail("Image ", map_image, " lacks spot and lowres scale factors")
  }
  lowres_spot_diameter <-
    as_single_numeric(scale_factors$spot, "spot scale factor") *
    as_single_numeric(scale_factors$lowres, "lowres scale factor")
  lowres_spot_radius <- lowres_spot_diameter / 2

  image_array <- slot_or_null(image_object, "image")
  if (is.null(image_array) || length(dim(image_array)) != 3L) {
    fail("Image raster unavailable for ", map_image)
  }

  list(
    image_name = map_image,
    image_object = image_object,
    image_cells = image_cells,
    donor = image_donors,
    stage = image_stages,
    coordinates = coordinates,
    scale_factors = scale_factors,
    spot_diameter = lowres_spot_diameter,
    spot_radius = lowres_spot_radius,
    image_array = image_array,
    image_height = dim(image_array)[1L],
    image_width = dim(image_array)[2L],
    values = as.matrix(normalized_selected[map_genes, image_cells, drop = FALSE])
  )
})
names(section_data) <- map_images

section_stages <- vapply(section_data, `[[`, character(1), "stage")
if (!setequal(section_stages, c("F0", "F4"))) {
  fail("The two requested images must comprise exactly one F0 and one F4 section; observed: ",
       paste(section_stages, collapse = ", "))
}
section_data <- section_data[order(match(section_stages, c("F0", "F4")))]

map_limits <- vapply(map_genes, function(gene) {
  pooled <- unlist(lapply(section_data, function(section) section$values[gene, ]),
                   use.names = FALSE)
  positive <- pooled[pooled > 0]
  if (length(positive) == 0L) fail("No positive mapped values for ", gene)
  current <- as.numeric(stats::quantile(positive, 0.99, names = FALSE, type = 8))
  if (!is.finite(current) || current <= 0) fail("Invalid map colour limit for ", gene)
  current
}, numeric(1))

map_palette <- grDevices::colorRampPalette(
  c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
)(257L)

map_colour <- function(values, limit) {
  clipped <- pmax(0, pmin(limit, values))
  index <- floor(clipped / limit * 256) + 1L
  colours <- map_palette[pmax(1L, pmin(257L, index))]
  colours[values <= 0] <- grDevices::adjustcolor("#BDBDBD", alpha.f = 0.22)
  colours[values > 0] <- grDevices::adjustcolor(colours[values > 0], alpha.f = 0.78)
  colours
}

draw_gene_map <- function(gene, section) {
  values <- as.numeric(section$values[gene, ])
  limit <- map_limits[[gene]]
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, section$image_width), ylim = c(section$image_height, 0), asp = 1,
    xaxs = "i", yaxs = "i"
  )
  graphics::rasterImage(
    as.raster(section$image_array), 0, section$image_height, section$image_width, 0,
    interpolate = TRUE
  )
  graphics::symbols(
    x = section$coordinates$imagecol,
    y = section$coordinates$imagerow,
    circles = rep(section$spot_radius, nrow(section$coordinates)),
    inches = FALSE,
    bg = map_colour(values, limit),
    fg = NA,
    add = TRUE
  )
  graphics::box(col = "#666666", lwd = 1.0)
  stage_label <- if (identical(section$stage, "F0")) "Healthy" else section$stage
  graphics::title(
    main = paste0(gene, " - ", stage_label),
    sub = paste0(section$image_name, " | donor ", section$donor),
    line = 0.4, cex.main = 1.15, cex.sub = 0.65, font.main = 1
  )
}

draw_map_legend <- function(gene) {
  limit <- map_limits[[gene]]
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  breaks <- seq(0.22, 0.55, length.out = length(map_palette) + 1L)
  positions <- seq(0.14, 0.86, length.out = length(map_palette) + 1L)
  graphics::rect(
    xleft = 0.22, ybottom = positions[-length(positions)],
    xright = 0.55, ytop = positions[-1L],
    col = map_palette, border = NA
  )
  graphics::rect(0.22, 0.14, 0.55, 0.86, border = "#555555", lwd = 0.6)
  graphics::text(0.63, c(0.14, 0.50, 0.86),
                 c("0", format(signif(limit / 2, 3), trim = TRUE),
                   format(signif(limit, 3), trim = TRUE)),
                 adj = c(0, 0.5), cex = 0.62)
  graphics::text(0.46, 0.94, gene, cex = 0.72, font = 2)
  graphics::text(0.46, 0.06, "Shared Healthy-F4\nlog-normalized scale", cex = 0.55)
}

write_spatial_panel <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(path, width = 12.0, height = 11.5, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = 3600, height = 3450, res = 300)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(
    matrix(seq_len(length(map_genes) * 3L), nrow = length(map_genes), byrow = TRUE),
    widths = c(1, 1, 0.16), heights = rep(1, length(map_genes))
  )
  graphics::par(oma = c(0.2, 0.2, 1.2, 0.2))
  for (gene in map_genes) {
    for (section in section_data) {
      graphics::par(mar = c(0.5, 0.5, 2.6, 0.5))
      draw_gene_map(gene, section)
    }
    graphics::par(mar = c(0.5, 0.1, 2.6, 0.1))
    draw_map_legend(gene)
  }
  graphics::mtext(
    "Representative spatial maps: Healthy versus F4",
    outer = TRUE, side = 3, line = 0.1, cex = 1.0, font = 2
  )
}

write_spatial_panel(
  file.path(result_tmp, "SupFig_selected_secreted_factor_spatial_maps_F0_vs_F4.pdf"), "pdf"
)
write_spatial_panel(
  file.path(result_tmp, "SupFig_selected_secreted_factor_spatial_maps_F0_vs_F4.png"), "png"
)

map_spot_values <- do.call(rbind, lapply(section_data, function(section) {
  data.frame(
    image = section$image_name,
    donor = section$donor,
    stage = section$stage,
    barcode = rep(section$image_cells, each = length(map_genes)),
    gene = rep(map_genes, times = length(section$image_cells)),
    imagerow = rep(section$coordinates$imagerow, each = length(map_genes)),
    imagecol = rep(section$coordinates$imagecol, each = length(map_genes)),
    normalized_expression = as.numeric(section$values),
    stringsAsFactors = FALSE
  )
}))

map_geometry <- do.call(rbind, lapply(section_data, function(section) {
  data.frame(
    image = section$image_name,
    donor = section$donor,
    stage = section$stage,
    n_spots = length(section$image_cells),
    image_height_lowres_px = section$image_height,
    image_width_lowres_px = section$image_width,
    spot_diameter_fullres_px = as_single_numeric(
      section$scale_factors$spot, "spot scale factor"
    ),
    lowres_scale_factor = as_single_numeric(
      section$scale_factors$lowres, "lowres scale factor"
    ),
    spot_diameter_lowres_px = section$spot_diameter,
    radius_used_lowres_px = section$spot_radius,
    stringsAsFactors = FALSE
  )
}))

selection_audit <- selection_audit[
  match(selected_gene_definition$gene, selection_audit$gene),
  , drop = FALSE
]
atomic_write_table(selected_gene_definition,
                   file.path(result_tmp, "selected_secreted_factor_definition.tsv"))
atomic_write_table(selection_audit,
                   file.path(result_tmp, "selected_secreted_factor_top30_audit.tsv"))
atomic_write_table(cluster_statistics,
                   file.path(result_tmp, "selected_secreted_factor_cluster_statistics.tsv"))
atomic_write_table(donor_cluster_statistics,
                   file.path(result_tmp, "selected_secreted_factor_donor_cluster_statistics.tsv"))
atomic_write_table(donor_balanced_cluster_statistics,
                   file.path(result_tmp, "selected_secreted_factor_donor_balanced_cluster_statistics.tsv"))
atomic_write_table(donor_stage_statistics,
                   file.path(result_tmp, "selected_secreted_factor_donor_stage_statistics.tsv"))
atomic_write_table(stage_association_statistics,
                   file.path(result_tmp, "selected_secreted_factor_stage_spearman.tsv"))
atomic_write_table(map_spot_values,
                   file.path(result_tmp, "selected_secreted_factor_map_spot_values.tsv"))
atomic_write_table(map_geometry,
                   file.path(result_tmp, "selected_secreted_factor_map_geometry.tsv"))

script_path <- {
  full_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", full_arguments, value = TRUE)
  if (length(file_argument) == 1L) sub("^--file=", "", file_argument) else NA_character_
}
script_sha <- if (!is.na(script_path) && file.exists(script_path)) sha256(script_path) else NA_character_

finished_at <- Sys.time()
run_manifest <- c(
  paste0("started_utc\t", format(started_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("finished_utc\t", format(finished_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("runtime_seconds\t", signif(as.numeric(difftime(finished_at, started_at, units = "secs")), 8)),
  paste0("input_path\t", input_path),
  paste0("input_sha256\t", input_sha),
  paste0("topgenes_path\t", topgenes_path),
  paste0("topgenes_sha256\t", topgenes_sha),
  paste0("script_sha256\t", script_sha),
  paste0("n_spots\t", ncol(object)),
  paste0("n_selected_genes\t", length(selected_genes)),
  paste0("cluster_metadata\tcluster_named"),
  paste0("stage_metadata\tF_score"),
  paste0("donor_metadata\tdonor"),
  paste0("map_images\t", paste(vapply(section_data, `[[`, character(1), "image_name"), collapse = ",")),
  paste0("map_donors\t", paste(vapply(section_data, `[[`, character(1), "donor"), collapse = ",")),
  paste0("map_stages\t", paste(vapply(section_data, `[[`, character(1), "stage"), collapse = ",")),
  paste0("normalization\tlog1p(count / spot_library_size * 10000)"),
  paste0("dot_colour\tmean expression z-scaled within gene across seven clusters"),
  paste0("dot_size\tpercentage of spots with count > 0"),
  paste0("main_heatmap_summary\tmean within donor x cluster, then equal-weight mean across donors"),
  paste0("stage_test\tdonor-level Spearman correlation with ordinal stage 0-4; BH adjustment across four mapped genes"),
  paste0("map_scale\tper-gene shared F0-F4 scale from 0 to pooled 99th percentile of positive expression")
)
atomic_write_lines(run_manifest, file.path(log_tmp, "run_manifest.tsv"))
atomic_write_lines(capture.output(sessionInfo()), file.path(log_tmp, "sessionInfo.txt"))

output_files <- sort(list.files(result_tmp, full.names = TRUE))
output_checksums <- data.frame(
  file = basename(output_files),
  sha256 = vapply(output_files, sha256, character(1)),
  stringsAsFactors = FALSE
)
atomic_write_table(output_checksums, file.path(log_tmp, "output_sha256.tsv"))

if (!file.rename(result_tmp, output_dir)) fail("Could not finalize result directory")
if (!file.rename(log_tmp, log_dir)) fail("Could not finalize log directory")
completed <- TRUE

message(
  "SECRETED_FACTOR_SPATIAL_CONTEXT_OK\n",
  "Results: ", output_dir, "\n",
  "Logs: ", log_dir
)
