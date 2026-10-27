#!/usr/bin/env Rscript

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
  "Seurat", "SeuratObject", "Matrix", "dorothea", "digest", "dplyr",
  "igraph", "ggraph", "ggplot2", "ggrepel", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  fail("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) fail("Arguments must be supplied as --key value pairs")
arg_list <- as.list(args[seq(2L, length(args), by = 2L)])
names(arg_list) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])

required_args <- c("input", "activity", "output-dir", "log-dir", "expected-input-sha256")
missing_args <- setdiff(required_args, names(arg_list))
if (length(missing_args) > 0L) {
  fail("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

input_path <- normalizePath(arg_list[["input"]], mustWork = TRUE)
activity_path <- normalizePath(arg_list[["activity"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
expected_input_sha <- tolower(arg_list[["expected-input-sha256"]])

if (!grepl("^[0-9a-f]{64}$", expected_input_sha)) {
  fail("--expected-input-sha256 must be a lowercase 64-character SHA256")
}
if (dir.exists(output_dir) || file.exists(output_dir)) fail("Output path already exists: ", output_dir)
if (dir.exists(log_dir) || file.exists(log_dir)) fail("Log path already exists: ", log_dir)

dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_dir), recursive = TRUE, showWarnings = FALSE)
run_stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
result_tmp <- paste0(output_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
log_tmp <- paste0(log_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
dir.create(result_tmp, recursive = FALSE, showWarnings = FALSE)
dir.create(log_tmp, recursive = FALSE, showWarnings = FALSE)
completed <- FALSE
on.exit({
  if (!completed) {
    message("Run did not complete; incomplete directories were retained for inspection:\n", result_tmp, "\n", log_tmp)
  }
}, add = TRUE)

atomic_write_table <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(x, temporary, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "")
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}
atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}
atomic_save_rds <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

script_path <- {
  full_args <- commandArgs(trailingOnly = FALSE)
  found <- grep("^--file=", full_args, value = TRUE)
  if (length(found) == 1L) normalizePath(sub("^--file=", "", found), mustWork = TRUE) else NA_character_
}
script_sha <- if (is.na(script_path)) NA_character_ else digest::digest(script_path, algo = "sha256", file = TRUE)
started_at <- Sys.time()

message("Computing input checksums...")
input_sha <- digest::digest(input_path, algo = "sha256", file = TRUE)
activity_sha <- digest::digest(activity_path, algo = "sha256", file = TRUE)
if (!identical(tolower(input_sha), expected_input_sha)) {
  fail("Input SHA256 mismatch. Expected ", expected_input_sha, "; observed ", input_sha)
}

message("Loading definitive DoRothEA activity object read-only...")
activity_object <- readRDS(activity_path)
required_activity_items <- c("activity", "spot_metadata", "selected_statistic", "times", "seed")
if (!all(required_activity_items %in% names(activity_object))) {
  fail("Activity object lacks: ", paste(setdiff(required_activity_items, names(activity_object)), collapse = ", "))
}
if (!identical(activity_object$selected_statistic, "norm_wmean") ||
    !identical(as.integer(activity_object$times), 1000L) ||
    !identical(as.integer(activity_object$seed), 42L)) {
  fail("Activity object is not the definitive norm_wmean/times=1000/seed=42 result")
}
if (!("RFX5" %in% rownames(activity_object$activity))) fail("RFX5 is absent from definitive activity matrix")

message("Loading regenerated Seurat object read-only...")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) fail("Input is not a Seurat object")
if (!("Spatial" %in% names(object@assays))) fail("Input lacks a Spatial assay")
counts <- SeuratObject::LayerData(object[["Spatial"]], layer = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
if (ncol(counts) != 26056L) fail("Unexpected number of Spatial spots: ", ncol(counts))

metadata <- object@meta.data
if (!identical(rownames(metadata), colnames(counts))) {
  metadata <- metadata[colnames(counts), , drop = FALSE]
  if (!identical(rownames(metadata), colnames(counts))) fail("Metadata cannot be aligned to Spatial count columns")
}
if (!all(c("donor", "F_score") %in% colnames(metadata))) fail("Input metadata lacks donor and/or F_score")
spot_metadata <- data.frame(
  barcode = colnames(counts),
  donor = as.character(metadata$donor),
  stage = as.character(metadata$F_score),
  stringsAsFactors = FALSE
)
if (anyNA(spot_metadata$donor) || anyNA(spot_metadata$stage)) fail("Missing donor or stage metadata")
spot_metadata$stage <- ifelse(grepl("^F", spot_metadata$stage), spot_metadata$stage, paste0("F", spot_metadata$stage))
spot_metadata$stage <- factor(spot_metadata$stage, levels = paste0("F", 0:4), ordered = TRUE)
if (anyNA(spot_metadata$stage)) fail("Unexpected stage labels")

donor_stage <- unique(spot_metadata[c("donor", "stage")])
donor_stage <- donor_stage[order(donor_stage$stage, donor_stage$donor), , drop = FALSE]
if (anyDuplicated(donor_stage$donor)) fail("A donor maps to more than one stage")
expected_donor_counts <- c(F0 = 6L, F1 = 6L, F2 = 3L, F3 = 3L, F4 = 6L)
observed_donor_counts <- table(donor_stage$stage)
if (!identical(as.integer(observed_donor_counts), as.integer(expected_donor_counts))) {
  fail("Unexpected donor counts by stage")
}

activity_matrix <- activity_object$activity
if (!identical(colnames(activity_matrix), spot_metadata$barcode)) {
  activity_matrix <- activity_matrix[, spot_metadata$barcode, drop = FALSE]
  if (!identical(colnames(activity_matrix), spot_metadata$barcode)) fail("Activity matrix cannot be aligned to input spots")
}
if (any(!is.finite(activity_matrix["RFX5", ]))) fail("RFX5 activity contains non-finite values")

message("Loading and validating canonical DoRothEA A/B resource...")
resource_environment <- new.env(parent = emptyenv())
utils::data("dorothea_hs", package = "dorothea", envir = resource_environment)
dorothea_hs <- get("dorothea_hs", envir = resource_environment, inherits = FALSE)
required_resource_columns <- c("tf", "confidence", "target", "mor")
if (!all(required_resource_columns %in% colnames(dorothea_hs))) fail("DoRothEA resource has unexpected columns")
resource_ab <- dorothea_hs[dorothea_hs$confidence %in% c("A", "B"), , drop = FALSE]
resource_ab <- dplyr::arrange(resource_ab, tf, target, confidence, mor)
resource_sha <- digest::digest(resource_ab, algo = "sha256", serialize = TRUE)
expected_resource_sha <- "7bf9a0895d3737354ab1d777c25165d57ba583b43e24daf0d0b0589c6bf73117"
if (!identical(resource_sha, expected_resource_sha)) fail("Canonical DoRothEA A/B resource checksum mismatch")

tf_name <- "RFX5"
requested_top_k <- 16L
regulon <- resource_ab[
  resource_ab$tf == tf_name,
  c("target", "confidence", "mor"),
  drop = FALSE
]
colnames(regulon)[colnames(regulon) == "target"] <- "gene"
if (nrow(regulon) == 0L || anyDuplicated(regulon$gene)) fail("Unexpected RFX5 A/B regulon structure")
regulon$mor <- as.numeric(regulon$mor)
regulon$mor_sign <- ifelse(regulon$mor >= 0, "activates", "represses")
regulon$present_in_input <- regulon$gene %in% rownames(counts)
regulon$total_count <- 0
present_index <- match(regulon$gene[regulon$present_in_input], rownames(counts))
regulon$total_count[regulon$present_in_input] <- Matrix::rowSums(counts[present_index, , drop = FALSE])
regulon$observed_nonzero <- regulon$present_in_input & regulon$total_count > 0
regulon_observed <- regulon[regulon$observed_nonzero, , drop = FALSE]
if (nrow(regulon_observed) < 2L) fail("Fewer than two observed RFX5 targets")
displayed_top_k <- min(requested_top_k, nrow(regulon_observed))

message("Computing donor-median target expression changes from regenerated counts...")
target_counts <- counts[regulon_observed$gene, , drop = FALSE]
library_sizes <- Matrix::colSums(counts)
if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) fail("Invalid spot library sizes")
scale_factor <- 10000
target_log <- target_counts %*% Matrix::Diagonal(x = scale_factor / library_sizes)
target_log@x <- log1p(target_log@x)
rownames(target_log) <- regulon_observed$gene
colnames(target_log) <- spot_metadata$barcode

donor_expression <- do.call(rbind, lapply(donor_stage$donor, function(donor_name) {
  donor_spots <- spot_metadata$barcode[spot_metadata$donor == donor_name]
  data.frame(
    gene = rownames(target_log),
    donor = donor_name,
    stage = as.character(donor_stage$stage[donor_stage$donor == donor_name]),
    expression_median = apply(target_log[, donor_spots, drop = FALSE], 1L, stats::median),
    stringsAsFactors = FALSE
  )
}))

stage_expression <- do.call(rbind, lapply(c("F0", "F4"), function(stage_name) {
  stage_values <- donor_expression[donor_expression$stage == stage_name, , drop = FALSE]
  data.frame(
    gene = unique(stage_values$gene),
    stage = stage_name,
    donor_median_expression = vapply(unique(stage_values$gene), function(gene_name) {
      stats::median(stage_values$expression_median[stage_values$gene == gene_name])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}))
f0_expression <- stage_expression[stage_expression$stage == "F0", c("gene", "donor_median_expression")]
f4_expression <- stage_expression[stage_expression$stage == "F4", c("gene", "donor_median_expression")]
colnames(f0_expression)[2L] <- "F0_median_donor_expression"
colnames(f4_expression)[2L] <- "F4_median_donor_expression"

edge_table <- merge(regulon_observed, f0_expression, by = "gene", all.x = TRUE, sort = FALSE)
edge_table <- merge(edge_table, f4_expression, by = "gene", all.x = TRUE, sort = FALSE)
edge_table$delta_expression_F4_minus_F0 <- edge_table$F4_median_donor_expression - edge_table$F0_median_donor_expression
edge_table$regulatory_effect <- edge_table$mor * edge_table$delta_expression_F4_minus_F0
edge_table$absolute_regulatory_effect <- abs(edge_table$regulatory_effect)
if (any(!is.finite(edge_table$delta_expression_F4_minus_F0)) || any(!is.finite(edge_table$regulatory_effect))) {
  fail("Non-finite target expression change or regulatory effect")
}
edge_table <- edge_table[order(-edge_table$absolute_regulatory_effect, edge_table$gene), , drop = FALSE]
edge_table$effect_rank <- seq_len(nrow(edge_table))
selected_edges <- edge_table[seq_len(displayed_top_k), , drop = FALSE]

edge_range <- range(selected_edges$absolute_regulatory_effect, na.rm = TRUE)
if (diff(edge_range) <= .Machine$double.eps) {
  selected_edges$edge_width <- 1.5
} else {
  selected_edges$edge_width <- scales::rescale(selected_edges$absolute_regulatory_effect, to = c(0.7, 2.8), from = edge_range)
}
selected_edges$from <- tf_name
selected_edges$to <- selected_edges$gene

graph <- igraph::graph_from_data_frame(
  selected_edges[c("from", "to", "mor_sign", "edge_width", "regulatory_effect", "absolute_regulatory_effect")],
  directed = TRUE
)
center_index <- match(tf_name, igraph::V(graph)$name)
coordinates <- igraph::layout_as_star(graph, center = center_index)
layout <- ggraph::create_layout(graph, layout = "manual", x = coordinates[, 1L], y = coordinates[, 2L])
layout$node_type <- ifelse(layout$name == tf_name, "index_TF", "target")
layout$delta_expression_F4_minus_F0 <- selected_edges$delta_expression_F4_minus_F0[
  match(layout$name, selected_edges$gene)
]
layout$label_x <- ifelse(layout$name == tf_name, 0, layout$x * 1.18)
layout$label_y <- ifelse(layout$name == tf_name, 0.20, layout$y * 1.18)

target_layout <- layout[layout$node_type == "target", , drop = FALSE]
tf_layout <- layout[layout$node_type == "index_TF", , drop = FALSE]
fill_limit <- max(abs(target_layout$delta_expression_F4_minus_F0), na.rm = TRUE)
if (!is.finite(fill_limit) || fill_limit <= 0) fill_limit <- 1

plot_object <- ggraph::ggraph(layout) +
  ggraph::geom_edge_link(
    ggplot2::aes(edge_colour = mor_sign, edge_width = edge_width),
    alpha = 0.90, show.legend = TRUE
  ) +
  ggraph::scale_edge_colour_manual(
    values = c(activates = "#D73027", represses = "#4575B4"),
    name = "DoRothEA regulation"
  ) +
  ggraph::scale_edge_width_identity(guide = "none") +
  ggraph::geom_node_point(
    data = target_layout,
    ggplot2::aes(fill = delta_expression_F4_minus_F0),
    shape = 21, size = 7.5, colour = "grey25", stroke = 0.35
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
    limits = c(-fill_limit, fill_limit), oob = scales::squish,
    name = expression(Delta~"target expression"~"(F4 - Healthy)")
  ) +
  ggraph::geom_node_point(
    data = tf_layout,
    fill = "#E66101", shape = 21, size = 10, colour = "grey20", stroke = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = target_layout,
    ggplot2::aes(x = x, y = y, label = name),
    nudge_x = target_layout$label_x - target_layout$x,
    nudge_y = target_layout$label_y - target_layout$y,
    size = 3.7, box.padding = 0.35, point.padding = 0.15,
    segment.colour = "grey70", min.segment.length = 0, max.overlaps = Inf,
    seed = 42
  ) +
  ggrepel::geom_label_repel(
    data = tf_layout,
    ggplot2::aes(x = x, y = y, label = name),
    nudge_y = 0.18, size = 4.0, fontface = "bold", fill = "white",
    label.size = 0.25, box.padding = 0.25, point.padding = 0.15,
    segment.colour = "grey60", min.segment.length = 0, max.overlaps = Inf,
    seed = 42
  ) +
  ggplot2::coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.25, 1.35), clip = "off") +
  ggplot2::labs(title = paste0("RFX5 A/B regulon: ", displayed_top_k, " observed targets")) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    plot.margin = ggplot2::margin(8, 30, 8, 8)
  ) +
  ggplot2::guides(
    edge_colour = ggplot2::guide_legend(order = 1),
    fill = ggplot2::guide_colorbar(order = 2)
  )

figure_stem <- "Fig4D_RFX5_regulon_observed_AB_targets_F4_vs_F0"
ggplot2::ggsave(file.path(result_tmp, paste0(figure_stem, ".pdf")), plot_object,
                width = 11.5, height = 7.5, device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(result_tmp, paste0(figure_stem, ".png")), plot_object,
                width = 11.5, height = 7.5, dpi = 300, bg = "white")

atomic_write_table(regulon, file.path(result_tmp, "07_rfx5_regulon_resource_targets.tsv"))
atomic_write_table(edge_table, file.path(result_tmp, "07_rfx5_all_observed_target_effects.tsv"))
atomic_write_table(selected_edges, file.path(result_tmp, "07_rfx5_displayed_regulatory_effects.tsv"))
atomic_write_table(donor_expression, file.path(result_tmp, "07_rfx5_target_donor_medians.tsv"))
atomic_write_table(donor_stage, file.path(result_tmp, "07_donor_stage_map.tsv"))

parameters <- data.frame(
  parameter = c(
    "input_object", "input_sha256", "dorothea_activity_object", "dorothea_activity_sha256",
    "assay", "input_layer", "normalization", "scale_factor", "regulon_database",
    "regulon_confidence", "regulon_resource_sha256", "index_TF", "contrast",
    "target_expression_summary", "target_ranking", "requested_top_k", "displayed_top_k", "edge_colour",
    "edge_width", "target_node_fill", "index_node_fill", "layout", "plot_seed", "blas_threads"
  ),
  value = c(
    input_path, input_sha, activity_path, activity_sha,
    "Spatial", "counts", "log1p(count / spot_library_size * scale_factor)", "10000",
    "dorothea::dorothea_hs", "A;B", resource_sha, tf_name, "F4_vs_F0",
    "median across spots per donor; median across donors within stage; F4 minus F0",
    "descending absolute(mor * target expression delta); gene-name tie break", as.character(requested_top_k), as.character(displayed_top_k),
    "activates=#D73027;represses=#4575B4", "absolute regulatory effect",
    "diverging target expression delta, F4 minus F0", "orange; not encoded as a target-expression value",
    "deterministic igraph star", "42", "1"
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(parameters, file.path(result_tmp, "07_run_parameters.tsv"))

package_versions <- data.frame(
  package = c("R", required_packages),
  version = c(R.version.string, vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))),
  stringsAsFactors = FALSE
)
atomic_write_table(package_versions, file.path(log_tmp, "07_package_versions.tsv"))

finished_at <- Sys.time()
summary_lines <- c(
  "FLINC Figure 4 RFX5 regulon execution summary",
  paste("Started UTC:", format(started_at, tz = "UTC", usetz = TRUE)),
  paste("Completed UTC:", format(finished_at, tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", as.numeric(difftime(finished_at, started_at, units = "secs"))),
  paste("Input SHA256:", input_sha),
  paste("DoRothEA activity SHA256:", activity_sha),
  paste("Script SHA256:", script_sha),
  paste("DoRothEA A/B resource SHA256:", resource_sha),
  paste("Observed RFX5 targets:", nrow(edge_table)),
  paste("Canonical RFX5 A/B targets:", nrow(regulon)),
  paste("Observed nonzero RFX5 targets:", nrow(edge_table)),
  paste("Displayed RFX5 targets:", nrow(selected_edges)),
  paste("RFX5 activity selected statistic:", activity_object$selected_statistic),
  paste("RFX5 activity times:", activity_object$times),
  paste("RFX5 activity seed:", activity_object$seed),
  paste("Command:", paste(commandArgs(), collapse = " "))
)
atomic_write_lines(summary_lines, file.path(log_tmp, "07_execution_summary.txt"))
atomic_write_lines(capture.output(sessionInfo()), file.path(log_tmp, "07_sessionInfo.txt"))

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
manifest <- data.frame(
  file = basename(result_files),
  bytes = file.info(result_files)$size,
  sha256 = vapply(result_files, digest::digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
atomic_write_table(manifest, file.path(result_tmp, "07_output_SHA256.tsv"))
atomic_write_lines("SUCCESS", file.path(result_tmp, "07_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "07_COMPLETED.txt"))

if (!file.rename(result_tmp, output_dir)) fail("Could not atomically publish results to ", output_dir)
if (!file.rename(log_tmp, log_dir)) fail("Could not atomically publish logs to ", log_dir)
completed <- TRUE
message("RFX5 regulon analysis completed successfully.")
message("Results: ", output_dir)
message("Compact logs: ", log_dir)
