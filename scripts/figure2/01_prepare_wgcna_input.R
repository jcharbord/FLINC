#!/usr/bin/env Rscript

# FLINC Figure 2 — prepare the spatial-expression object for hdWGCNA.
# This script stops before SetupForWGCNA and metaspot construction.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

require_columns <- function(metadata, columns) {
  missing_columns <- setdiff(columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stopf(
      "Required metadata columns are missing: %s",
      paste(missing_columns, collapse = ", ")
    )
  }
}

write_lines_atomic <- function(lines, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) {
    stopf("Partial file already exists; inspect it before retrying: %s", partial)
  }
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) {
    stopf("Could not move partial file into place: %s", path)
  }
}

write_table_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) {
    stopf("Partial file already exists; inspect it before retrying: %s", partial)
  }
  write.table(
    x,
    file = partial,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  if (!file.rename(partial, path)) {
    stopf("Could not move partial file into place: %s", path)
  }
}

save_rds_atomic <- function(object, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) {
    stopf("Partial file already exists; inspect it before retrying: %s", partial)
  }
  say("Saving prepared object to", partial)
  saveRDS(object, partial, compress = "gzip", version = 3)
  if (!file.rename(partial, path)) {
    stopf("Could not move partial RDS into place: %s", path)
  }
}

repro_root <- Sys.getenv(
  "FLINC_REPRO_ROOT",
  unset = "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)

input_file <- Sys.getenv(
  "FIG1_INPUT",
  unset = file.path(
    repro_root,
    "results",
    "260701_approach1_cluster_go_marker_deconv",
    paste0(
      "object_full_with_final_annotation_markers_deconv_proportions_",
      "ready_C3green_C5yellow.rds"
    )
  )
)

output_dir <- file.path(repro_root, "objects", "Fig2_WGCNA")
repro_log_dir <- file.path(repro_root, "reproducibility_logs", "Fig2_WGCNA")
output_file <- Sys.getenv(
  "FIG2_PREPARED_OUTPUT",
  unset = file.path(output_dir, "01_wgcna_input_prepared.rds")
)
summary_file <- file.path(repro_log_dir, "01_wgcna_input_summary.tsv")
session_file <- file.path(repro_log_dir, "01_wgcna_input_sessionInfo.txt")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_log_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(output_file, summary_file, session_file)) {
  if (file.exists(path)) {
    stopf("Refusing to overwrite existing output: %s", path)
  }
}

if (!file.exists(input_file)) {
  stopf("Input object does not exist: %s", input_file)
}

say("R version:", R.version.string)
say("Seurat version:", as.character(packageVersion("Seurat")))
say("SeuratObject version:", as.character(packageVersion("SeuratObject")))
say("Input:", input_file)
say("Input size (GiB):", sprintf("%.2f", file.info(input_file)$size / 1024^3))
say("Loading input object")

load_start <- Sys.time()
object <- readRDS(input_file)
say(
  "Input loaded in minutes:",
  sprintf("%.2f", as.numeric(difftime(Sys.time(), load_start, units = "mins")))
)

if (!inherits(object, "Seurat")) {
  stopf("Input is not a Seurat object; class is: %s", paste(class(object), collapse = ", "))
}

stored_seurat_version <- as.character(object@version)
initial_assays <- Assays(object)
initial_reductions <- Reductions(object)
initial_images <- Images(object)

if (!("Spatial" %in% initial_assays)) {
  stopf("Input object does not contain the Spatial assay")
}

required_metadata <- c(
  "donor",
  "F_score",
  "cluster",
  "annotation",
  "cluster_named"
)
require_columns(object[[]], required_metadata)

for (column_name in required_metadata) {
  if (anyNA(object[[column_name, drop = TRUE]])) {
    stopf("Required metadata column contains missing values: %s", column_name)
  }
}

observed_f_scores <- sort(unique(as.character(object$F_score)))
expected_f_scores <- paste0("F", 0:4)
if (!identical(observed_f_scores, expected_f_scores)) {
  stopf(
    "Unexpected F_score values. Expected [%s], observed [%s]",
    paste(expected_f_scores, collapse = ", "),
    paste(observed_f_scores, collapse = ", ")
  )
}

donor_to_stage <- tapply(
  as.character(object$F_score),
  as.character(object$donor),
  function(x) length(unique(x))
)
if (any(donor_to_stage != 1L)) {
  stopf("At least one donor maps to more than one F_score")
}

sample_id <- paste0(as.character(object$F_score), "_", as.character(object$donor))
sample_levels <- sort(unique(sample_id))
if (length(sample_levels) != 24L) {
  stopf("Expected 24 F_score_donor samples, observed %d", length(sample_levels))
}
object$F_score_donor <- factor(sample_id, levels = sample_levels)

say("Validated", ncol(object), "spots across", length(sample_levels), "samples")
say("Stored Seurat object version:", stored_seurat_version)
say("Initial assays:", paste(initial_assays, collapse = ", "))
say("Initial reductions:", paste(initial_reductions, collapse = ", "))
say("Initial images:", paste(initial_images, collapse = ", "))

# Import spatial coordinates from every Visium image.
if (length(initial_images) == 0L) {
  stopf("Input object contains no spatial images")
}

coordinate_list <- lapply(initial_images, function(slice_name) {
  coordinates <- as.data.frame(object@images[[slice_name]]@coordinates)
  if (nrow(coordinates) == 0L) {
    stopf("Spatial image has no coordinates: %s", slice_name)
  }
  coordinates$barcode <- rownames(coordinates)
  coordinates$slice <- slice_name
  coordinates
})

coordinates_all <- do.call(rbind, coordinate_list)
if (anyDuplicated(coordinates_all$barcode)) {
  duplicated_barcodes <- unique(
    coordinates_all$barcode[duplicated(coordinates_all$barcode)]
  )
  stopf(
    "Duplicated spatial barcodes found across images (first examples: %s)",
    paste(head(duplicated_barcodes, 10L), collapse = ", ")
  )
}

rownames(coordinates_all) <- coordinates_all$barcode
missing_coordinates <- setdiff(colnames(object), rownames(coordinates_all))
if (length(missing_coordinates) > 0L) {
  stopf(
    "%d object spots have no image coordinates (first examples: %s)",
    length(missing_coordinates),
    paste(head(missing_coordinates, 10L), collapse = ", ")
  )
}

coordinate_columns <- intersect(
  c("row", "col", "imagerow", "imagecol", "barcode"),
  colnames(coordinates_all)
)
if (length(coordinate_columns) < 3L) {
  stopf(
    "Too few expected coordinate columns were found: %s",
    paste(coordinate_columns, collapse = ", ")
  )
}

object <- AddMetaData(
  object,
  metadata = coordinates_all[colnames(object), coordinate_columns, drop = FALSE]
)

# Preserve images, metadata, and reductions while retaining only the raw
# Spatial assay used for network construction.
DefaultAssay(object) <- "Spatial"
object@assays <- object@assays["Spatial"]

initial_spatial_layers <- Layers(object[["Spatial"]])
if (!identical(initial_spatial_layers, "counts")) {
  stopf(
    "The validated input should contain only the Spatial counts layer; observed: %s",
    paste(initial_spatial_layers, collapse = ", ")
  )
}

say("Splitting Spatial counts by F_score_donor before normalization")
object[["Spatial"]] <- split(object[["Spatial"]], f = object$F_score_donor)

split_layers <- Layers(object[["Spatial"]])
count_layers <- grep("^counts\\.", split_layers, value = TRUE)
if (length(count_layers) != length(sample_levels)) {
  stopf(
    "Expected %d split count layers, observed %d",
    length(sample_levels),
    length(count_layers)
  )
}

say("Running NormalizeData per split sample with Seurat defaults")
object <- NormalizeData(object, assay = "Spatial", verbose = TRUE)

say("Running FindVariableFeatures with Seurat defaults")
object <- FindVariableFeatures(object, assay = "Spatial", verbose = TRUE)
if (length(VariableFeatures(object[["Spatial"]])) == 0L) {
  stopf("FindVariableFeatures returned no variable features")
}

say("Running ScaleData with Seurat defaults")
object <- ScaleData(object, assay = "Spatial", verbose = TRUE)

say("Joining the per-sample Spatial layers")
object <- JoinLayers(object, assay = "Spatial")
object$F_score_donor <- factor(as.character(object$F_score_donor), levels = sample_levels)
Idents(object) <- "cluster_named"

final_layers <- Layers(object[["Spatial"]])
required_final_layers <- c("counts", "data", "scale.data")
missing_final_layers <- setdiff(required_final_layers, final_layers)
if (length(missing_final_layers) > 0L) {
  stopf(
    "Prepared Spatial assay is missing expected layers: %s",
    paste(missing_final_layers, collapse = ", ")
  )
}

if (ncol(object) != 26056L || nrow(object[["Spatial"]]) != 17922L) {
  stopf(
    "Unexpected prepared dimensions: Spatial=%d features x %d spots",
    nrow(object[["Spatial"]]),
    ncol(object)
  )
}

sample_counts <- as.data.frame(table(object$F_score_donor), stringsAsFactors = FALSE)
colnames(sample_counts) <- c("F_score_donor", "n_spots")

summary_table <- rbind(
  data.frame(metric = "input_file", value = input_file),
  data.frame(metric = "input_size_GiB", value = sprintf("%.4f", file.info(input_file)$size / 1024^3)),
  data.frame(metric = "stored_Seurat_version", value = stored_seurat_version),
  data.frame(metric = "current_Seurat_version", value = as.character(packageVersion("Seurat"))),
  data.frame(metric = "current_SeuratObject_version", value = as.character(packageVersion("SeuratObject"))),
  data.frame(metric = "n_features_Spatial", value = nrow(object[["Spatial"]])),
  data.frame(metric = "n_spots", value = ncol(object)),
  data.frame(metric = "n_samples", value = length(sample_levels)),
  data.frame(metric = "n_images", value = length(Images(object))),
  data.frame(metric = "n_variable_features", value = length(VariableFeatures(object[["Spatial"]]))),
  data.frame(metric = "final_layers", value = paste(final_layers, collapse = ",")),
  data.frame(
    metric = paste0("spots_", sample_counts$F_score_donor),
    value = sample_counts$n_spots
  )
)

write_table_atomic(summary_table, summary_file)
write_lines_atomic(capture.output(sessionInfo()), session_file)
save_rds_atomic(object, output_file)

say("Prepared object size (GiB):", sprintf("%.2f", file.info(output_file)$size / 1024^3))
say("Summary:", summary_file)
say("Session information:", session_file)
say("Fig. 2 WGCNA input preparation completed successfully")
