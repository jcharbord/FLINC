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

required_packages <- c("Seurat", "SeuratObject", "Matrix", "digest", "dplyr", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) fail("Missing required packages: ", paste(missing_packages, collapse = ", "))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) fail("Arguments must be supplied as --key value pairs")
arg_list <- as.list(args[seq(2L, length(args), by = 2L)])
names(arg_list) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])
required_args <- c("input", "activity", "output-dir", "log-dir", "expected-input-sha256")
missing_args <- setdiff(required_args, names(arg_list))
if (length(missing_args) > 0L) fail("Missing required arguments: ", paste(missing_args, collapse = ", "))

input_path <- normalizePath(arg_list[["input"]], mustWork = TRUE)
activity_path <- normalizePath(arg_list[["activity"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
expected_input_sha <- tolower(arg_list[["expected-input-sha256"]])
if (!grepl("^[0-9a-f]{64}$", expected_input_sha)) fail("--expected-input-sha256 must be a lowercase 64-character SHA256")
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
  if (!completed) message("Run did not complete; incomplete directories were retained for inspection:\n", result_tmp, "\n", log_tmp)
}, add = TRUE)

atomic_write_table <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(x, temporary, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}
atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
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
if (!identical(tolower(input_sha), expected_input_sha)) fail("Input SHA256 mismatch. Expected ", expected_input_sha, "; observed ", input_sha)

message("Loading definitive DoRothEA activity object read-only...")
activity_object <- readRDS(activity_path)
required_activity_items <- c("activity", "spot_metadata", "selected_statistic", "times", "seed")
if (!all(required_activity_items %in% names(activity_object))) fail("Activity object lacks required items")
if (!identical(activity_object$selected_statistic, "norm_wmean") ||
    !identical(as.integer(activity_object$times), 1000L) ||
    !identical(as.integer(activity_object$seed), 42L)) {
  fail("Activity object is not the definitive norm_wmean/times=1000/seed=42 result")
}
if (!("RFX5" %in% rownames(activity_object$activity))) fail("RFX5 is absent from definitive activity matrix")

message("Loading regenerated Seurat object read-only...")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) fail("Input is not a Seurat object")
if (!all(c("Spatial", "singlecell") %in% names(object@assays))) {
  fail("Input must contain Spatial and singlecell assays")
}
counts <- SeuratObject::LayerData(object[["Spatial"]], layer = "counts")
if (ncol(counts) != 26056L) fail("Unexpected number of Spatial spots: ", ncol(counts))
metadata <- object@meta.data
if (!identical(rownames(metadata), colnames(counts))) {
  metadata <- metadata[colnames(counts), , drop = FALSE]
  if (!identical(rownames(metadata), colnames(counts))) fail("Metadata cannot be aligned to Spatial count columns")
}

celltype_order <- c(
  "CentralHep", "InterHep", "PortalHep", "InfMac", "NonInfMac", "Bcells",
  "NKTcell", "Stellate", "Cholangiocyte", "PortalEndo", "Erythroid", "cvLSECs"
)
required_metadata <- c("donor", "F_score")
if (!all(required_metadata %in% colnames(metadata))) {
  fail("Input metadata lacks: ", paste(setdiff(required_metadata, colnames(metadata)), collapse = ", "))
}

deconvolution_values <- SeuratObject::LayerData(
  object[["singlecell"]], layer = "counts"
)
if (!all(celltype_order %in% rownames(deconvolution_values))) {
  fail(
    "singlecell counts layer lacks deconvolution cell types: ",
    paste(setdiff(celltype_order, rownames(deconvolution_values)), collapse = ", ")
  )
}
if (!identical(colnames(deconvolution_values), colnames(counts))) {
  fail("singlecell deconvolution and Spatial count spot barcodes are not identically ordered")
}
proportion_matrix <- t(as.matrix(
  deconvolution_values[celltype_order, colnames(counts), drop = FALSE]
))
if (any(!is.finite(proportion_matrix)) ||
    min(proportion_matrix) < -1e-12 || max(proportion_matrix) > 1 + 1e-12) {
  fail("singlecell deconvolution fractions are non-finite or outside [0, 1]")
}
maximum_proportion_sum_error <- max(abs(rowSums(proportion_matrix) - 1))
if (maximum_proportion_sum_error > 1e-8) {
  fail(
    "singlecell deconvolution fractions do not sum to one per spot. Maximum error: ",
    format(maximum_proportion_sum_error, scientific = TRUE)
  )
}

spot_metadata <- data.frame(
  barcode = colnames(counts),
  donor = as.character(metadata$donor),
  stage = as.character(metadata$F_score),
  proportion_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
spot_metadata$stage <- ifelse(grepl("^F", spot_metadata$stage), spot_metadata$stage, paste0("F", spot_metadata$stage))
spot_metadata$stage <- factor(spot_metadata$stage, levels = paste0("F", 0:4), ordered = TRUE)
if (anyNA(spot_metadata$donor) || anyNA(spot_metadata$stage)) fail("Missing donor or stage metadata")

donor_stage <- unique(spot_metadata[c("donor", "stage")])
donor_stage <- donor_stage[order(donor_stage$stage, donor_stage$donor), , drop = FALSE]
if (anyDuplicated(donor_stage$donor)) fail("A donor maps to more than one stage")
expected_donor_counts <- c(F0 = 6L, F1 = 6L, F2 = 3L, F3 = 3L, F4 = 6L)
observed_donor_counts <- table(donor_stage$stage)
if (!identical(as.integer(observed_donor_counts), as.integer(expected_donor_counts))) fail("Unexpected donor counts by stage")

activity_matrix <- activity_object$activity
if (!identical(colnames(activity_matrix), spot_metadata$barcode)) {
  activity_matrix <- activity_matrix[, spot_metadata$barcode, drop = FALSE]
  if (!identical(colnames(activity_matrix), spot_metadata$barcode)) fail("Activity matrix cannot be aligned to input spots")
}
spot_metadata$RFX5_activity_norm_wmean <- as.numeric(activity_matrix["RFX5", spot_metadata$barcode])
if (any(!is.finite(spot_metadata$RFX5_activity_norm_wmean))) fail("RFX5 activity contains non-finite values")

for (cell_type in celltype_order) {
  spot_metadata[[cell_type]] <- as.numeric(spot_metadata[[cell_type]])
  if (any(!is.finite(spot_metadata[[cell_type]]))) fail("Non-finite values in deconvolution column: ", cell_type)
  if (any(spot_metadata[[cell_type]] < -1e-8 | spot_metadata[[cell_type]] > 1 + 1e-8)) {
    fail("Deconvolution fractions outside [0, 1] in column: ", cell_type)
  }
}

message("Computing within-stage spot-level Spearman correlations...")
stages <- paste0("F", 0:4)
correlations <- do.call(rbind, lapply(stages, function(stage_name) {
  stage_data <- spot_metadata[spot_metadata$stage == stage_name, , drop = FALSE]
  n_donors <- length(unique(stage_data$donor))
  result <- do.call(rbind, lapply(celltype_order, function(cell_type) {
    fractions <- stage_data[[cell_type]]
    activity <- stage_data$RFX5_activity_norm_wmean
    valid <- is.finite(fractions) & is.finite(activity)
    if (sum(valid) < 10L || stats::sd(fractions[valid]) <= 1e-12) {
      return(data.frame(
        stage = stage_name, cell_type = cell_type, rho = NA_real_, p_value = NA_real_,
        n_spots = sum(valid), n_donors = n_donors, stringsAsFactors = FALSE
      ))
    }
    test <- suppressWarnings(stats::cor.test(
      activity[valid], fractions[valid], method = "spearman", exact = FALSE
    ))
    data.frame(
      stage = stage_name, cell_type = cell_type, rho = unname(test$estimate), p_value = test$p.value,
      n_spots = sum(valid), n_donors = n_donors, stringsAsFactors = FALSE
    )
  }))
  result$FDR_BH_within_stage <- stats::p.adjust(result$p_value, method = "BH")
  result
}))
correlations$stage <- factor(correlations$stage, levels = stages, ordered = TRUE)
correlations$cell_type <- factor(correlations$cell_type, levels = celltype_order, ordered = TRUE)
correlations <- correlations[order(correlations$cell_type, correlations$stage), , drop = FALSE]

if (any(!is.finite(correlations$rho))) fail("At least one stage/cell-type correlation is undefined")
rho_limit <- max(abs(correlations$rho))
if (!is.finite(rho_limit) || rho_limit <= 0) rho_limit <- 1
rho_limit <- ceiling(rho_limit * 20) / 20

heatmap <- ggplot2::ggplot(
  correlations,
  ggplot2::aes(x = stage, y = cell_type, fill = rho)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
    limits = c(-rho_limit, rho_limit), oob = scales::squish,
    name = "Spearman rho"
  ) +
  ggplot2::scale_x_discrete(position = "top", labels = function(x) ifelse(x == "F0", "Healthy", x)) +
  ggplot2::scale_y_discrete(limits = rev(celltype_order)) +
  ggplot2::labs(
    x = "Stage", y = "Cell type",
    title = "RFX5 activity and deconvolution fractions"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 9),
    axis.title = ggplot2::element_text(face = "bold")
  )

figure_stem <- "Fig4E_RFX5_activity_celltype_spearman_by_stage"
ggplot2::ggsave(file.path(result_tmp, paste0(figure_stem, ".pdf")), heatmap,
                width = 5.4, height = 6.3, device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(result_tmp, paste0(figure_stem, ".png")), heatmap,
                width = 5.4, height = 6.3, dpi = 300, bg = "white")

deconvolution_summary <- do.call(rbind, lapply(celltype_order, function(cell_type) {
  data.frame(
    cell_type = cell_type,
    min_fraction = min(spot_metadata[[cell_type]]),
    median_fraction = stats::median(spot_metadata[[cell_type]]),
    max_fraction = max(spot_metadata[[cell_type]]),
    stringsAsFactors = FALSE
  )
}))
atomic_write_table(correlations, file.path(result_tmp, "08_rfx5_celltype_spearman_by_stage.tsv"))
atomic_write_table(deconvolution_summary, file.path(result_tmp, "08_deconvolution_fraction_summary.tsv"))
atomic_write_table(donor_stage, file.path(result_tmp, "08_donor_stage_map.tsv"))

parameters <- data.frame(
  parameter = c(
    "input_object", "input_sha256", "dorothea_activity_object", "dorothea_activity_sha256",
    "selected_statistic", "times", "seed", "index_TF", "activity_level", "association_level",
    "deconvolution_assay", "deconvolution_layer", "maximum_deconvolution_sum_error",
    "correlation_method", "correlation_test", "multiple_testing", "stage_order", "cell_type_order",
    "colour_scale", "symmetric_rho_limit", "main_figure_interpretation", "blas_threads"
  ),
  value = c(
    input_path, input_sha, activity_path, activity_sha,
    activity_object$selected_statistic, as.character(activity_object$times), as.character(activity_object$seed),
    "RFX5", "spot", "within-stage spot-level descriptive association",
    "singlecell", "counts", format(maximum_proportion_sum_error, scientific = TRUE),
    "Spearman rho", "two-sided cor.test(exact=FALSE)", "Benjamini-Hochberg within each stage across 12 cell types",
    paste(stages, collapse = ";"), paste(celltype_order, collapse = ";"),
    "#4575B4=negative;white=zero;#D73027=positive", as.character(rho_limit),
    "rho heatmap; p-values/FDR retained in table but not used for a claim of independent spot-level inference", "1"
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(parameters, file.path(result_tmp, "08_run_parameters.tsv"))

package_versions <- data.frame(
  package = c("R", required_packages),
  version = c(R.version.string, vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))),
  stringsAsFactors = FALSE
)
atomic_write_table(package_versions, file.path(log_tmp, "08_package_versions.tsv"))

finished_at <- Sys.time()
summary_lines <- c(
  "FLINC Figure 4 RFX5 cell-type association execution summary",
  paste("Started UTC:", format(started_at, tz = "UTC", usetz = TRUE)),
  paste("Completed UTC:", format(finished_at, tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", as.numeric(difftime(finished_at, started_at, units = "secs"))),
  paste("Input SHA256:", input_sha),
  paste("DoRothEA activity SHA256:", activity_sha),
  paste("Script SHA256:", script_sha),
  paste("RFX5 statistic:", activity_object$selected_statistic),
  paste("RFX5 scoring times:", activity_object$times),
  paste("RFX5 scoring seed:", activity_object$seed),
  paste("Deconvolution source:", "singlecell assay, counts layer"),
  paste("Maximum deconvolution per-spot sum error:", format(maximum_proportion_sum_error, scientific = TRUE)),
  paste("Correlation level: within-stage spot-level"),
  paste("Spearman rho colour limit:", rho_limit),
  paste("Command:", paste(commandArgs(), collapse = " "))
)
atomic_write_lines(summary_lines, file.path(log_tmp, "08_execution_summary.txt"))
atomic_write_lines(capture.output(sessionInfo()), file.path(log_tmp, "08_sessionInfo.txt"))

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
manifest <- data.frame(
  file = basename(result_files), bytes = file.info(result_files)$size,
  sha256 = vapply(result_files, digest::digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
atomic_write_table(manifest, file.path(result_tmp, "08_output_SHA256.tsv"))
atomic_write_lines("SUCCESS", file.path(result_tmp, "08_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "08_COMPLETED.txt"))

if (!file.rename(result_tmp, output_dir)) fail("Could not atomically publish results to ", output_dir)
if (!file.rename(log_tmp, log_dir)) fail("Could not atomically publish logs to ", log_dir)
completed <- TRUE
message("RFX5 cell-type association analysis completed successfully.")
message("Results: ", output_dir)
message("Compact logs: ", log_dir)
