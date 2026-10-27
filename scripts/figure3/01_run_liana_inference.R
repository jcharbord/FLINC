#!/usr/bin/env Rscript

# Figure 3: LIANA inference from the validated annotated FLINC object.
# This script never writes to or overwrites the input Seurat RDS.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(liana)
  library(SingleCellExperiment)
})

options(width = 180)

started_at <- Sys.time()

repo_root <- "/Users/jianwu/Documents/FLINC"
repro_root <- file.path(repo_root, "Jeremie", "Reproducible_code")

input_rds <- file.path(
  repro_root,
  "results",
  "260701_approach1_cluster_go_marker_deconv",
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)

output_dir <- file.path(repro_root, "results", "Fig3_CCC", "01_liana_inference")

species_use <- "human"
assay_use <- "Spatial"
input_layer <- "counts"
deconv_assay <- "singlecell"
deconv_layer <- "counts"

cell_types <- c(
  "Bcells", "CentralHep", "Cholangiocyte", "Erythroid",
  "InfMac", "InterHep", "NKTcell", "NonInfMac",
  "PortalEndo", "PortalHep", "Stellate", "cvLSECs"
)

sample_column <- "donor"
condition_column <- "F_score"
stage_levels <- c("F0", "F1", "F2", "F3", "F4")

resource_use <- "OmniPath"
methods_use <- c("cellchat", "natmi", "sca", "logfc", "connectome")
pseudobulk_scale_factor <- 1e6
significance_rank_threshold <- 0.01

expected_spots <- 26056L
expected_input_spatial_genes <- 17922L
expected_retained_spatial_genes <- 17775L
expected_zero_genes_dropped <- 147L
expected_donors <- 24L
expected_pseudobulks <- length(cell_types) * expected_donors

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr",
  "tibble", "liana", "SingleCellExperiment", "S4Vectors", "digest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

required_liana_functions <- c("liana_wrap", "liana_aggregate")
missing_liana_functions <- required_liana_functions[
  !vapply(
    required_liana_functions,
    exists,
    logical(1),
    envir = asNamespace("liana"),
    inherits = FALSE
  )
]

if (length(missing_liana_functions)) {
  stop(
    "Required LIANA functions are unavailable: ",
    paste(missing_liana_functions, collapse = ", ")
  )
}

if (!file.exists(input_rds)) {
  stop("Input RDS not found: ", input_rds)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

existing_outputs <- list.files(
  output_dir,
  all.files = TRUE,
  no.. = TRUE,
  full.names = TRUE
)

if (length(existing_outputs)) {
  stop(
    "Output directory is not empty. Preserve or back it up before rerunning: ",
    output_dir
  )
}

atomic_save_rds <- function(object, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  saveRDS(object, partial)
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
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

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
} else {
  NA_character_
}

message("Figure 3 LIANA inference started at ", format(started_at))
message("Input: ", input_rds)
message("Output: ", output_dir)
message("Calculating input SHA256 before loading the object...")

input_info <- file.info(input_rds)
input_sha256 <- sha256_file(input_rds)

message("Input SHA256: ", input_sha256)
message("Loading validated annotated Seurat object...")

sub_obj <- readRDS(input_rds)

if (!inherits(sub_obj, "Seurat")) {
  stop("Input is not a Seurat object.")
}

if (!all(c(assay_use, deconv_assay) %in% names(sub_obj@assays))) {
  stop("Required assays Spatial and singlecell were not both found.")
}

if (!all(c(sample_column, condition_column) %in% colnames(sub_obj@meta.data))) {
  stop("Required donor and F_score metadata columns were not found.")
}

spatial_counts <- SeuratObject::LayerData(
  sub_obj[[assay_use]],
  layer = input_layer
)

deconv_values <- SeuratObject::GetAssayData(
  sub_obj,
  assay = deconv_assay,
  layer = deconv_layer
)

metadata <- sub_obj@meta.data[, c(sample_column, condition_column), drop = FALSE]

if (!identical(colnames(spatial_counts), rownames(metadata))) {
  stop("Spatial-count columns and metadata rows are not identically ordered.")
}

if (!all(cell_types %in% rownames(deconv_values))) {
  stop(
    "Missing deconvolution cell types: ",
    paste(setdiff(cell_types, rownames(deconv_values)), collapse = ", ")
  )
}

if (!identical(colnames(spatial_counts), colnames(deconv_values))) {
  stop("Spatial-count and deconvolution spot names are not identically ordered.")
}

if (ncol(spatial_counts) != expected_spots) {
  stop("Unexpected spot count: ", ncol(spatial_counts), " (expected ", expected_spots, ").")
}

if (nrow(spatial_counts) != expected_input_spatial_genes) {
  stop(
    "Unexpected input Spatial-gene count: ",
    nrow(spatial_counts),
    " (expected ", expected_input_spatial_genes, ")."
  )
}

proportions <- t(as.matrix(
  deconv_values[cell_types, colnames(spatial_counts), drop = FALSE]
))

if (any(!is.finite(proportions))) {
  stop("Deconvolution proportions contain non-finite values.")
}

if (min(proportions) < -1e-12 || max(proportions) > 1 + 1e-12) {
  stop("Deconvolution proportions fall outside the expected [0, 1] interval.")
}

proportion_sum_error <- max(abs(rowSums(proportions) - 1))
if (proportion_sum_error > 1e-8) {
  stop(
    "Deconvolution proportions do not sum to one per spot. Maximum error: ",
    format(proportion_sum_error, scientific = TRUE)
  )
}

donors <- as.character(metadata[[sample_column]])
conditions <- as.character(metadata[[condition_column]])

if (!setequal(unique(conditions), stage_levels)) {
  stop(
    "Observed stages do not match F0-F4. Observed: ",
    paste(sort(unique(conditions)), collapse = ", ")
  )
}

donor_stage_check <- unique(data.frame(
  donor = donors,
  condition = conditions,
  stringsAsFactors = FALSE
))

if (anyDuplicated(donor_stage_check$donor)) {
  stop("At least one donor maps to more than one fibrosis stage.")
}

if (nrow(donor_stage_check) != expected_donors) {
  stop(
    "Unexpected donor count: ",
    nrow(donor_stage_check),
    " (expected ", expected_donors, ")."
  )
}

rm(sub_obj, deconv_values)
invisible(gc())

message("Dropping zero-count Spatial genes in memory to match the historical Figure 3 input...")

keep_gene <- Matrix::rowSums(spatial_counts) > 0
zero_genes_dropped <- sum(!keep_gene)
spatial_counts <- spatial_counts[keep_gene, , drop = FALSE]

if (zero_genes_dropped != expected_zero_genes_dropped) {
  stop(
    "Unexpected zero-gene count: ",
    zero_genes_dropped,
    " (expected ", expected_zero_genes_dropped, ")."
  )
}

if (nrow(spatial_counts) != expected_retained_spatial_genes) {
  stop(
    "Unexpected retained Spatial-gene count: ",
    nrow(spatial_counts),
    " (expected ", expected_retained_spatial_genes, ")."
  )
}

donor_levels <- unique(donors)

spot_to_donor <- Matrix::Matrix(
  0,
  nrow = length(donors),
  ncol = length(donor_levels),
  dimnames = list(rownames(metadata), donor_levels),
  sparse = TRUE
)

spot_to_donor[cbind(
  seq_along(donors),
  match(donors, donor_levels)
)] <- 1

message("Constructing 12 cell-type x 24 donor weighted pseudobulks...")

weight_list <- lapply(
  seq_along(cell_types),
  function(cell_type_index) {
    Matrix::Diagonal(x = proportions[, cell_type_index]) %*% spot_to_donor
  }
)

weights <- do.call(cbind, weight_list)
colnames(weights) <- as.vector(vapply(
  cell_types,
  function(cell_type) paste(cell_type, donor_levels, sep = "||"),
  character(length(donor_levels))
))

pseudobulk_counts <- Matrix::drop0(spatial_counts %*% weights)

if (ncol(pseudobulk_counts) != expected_pseudobulks) {
  stop(
    "Unexpected pseudobulk count: ",
    ncol(pseudobulk_counts),
    " (expected ", expected_pseudobulks, ")."
  )
}

pseudobulk_metadata <- tidyr::separate(
  data.frame(
    pb_id = colnames(pseudobulk_counts),
    stringsAsFactors = FALSE
  ),
  col = "pb_id",
  into = c("celltype", "donor"),
  sep = "\\|\\|",
  remove = FALSE
)

donor_stage_map <- metadata %>%
  tibble::rownames_to_column("spot") %>%
  dplyr::transmute(
    donor = as.character(.data[[sample_column]]),
    condition = as.character(.data[[condition_column]])
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    condition = factor(condition, levels = stage_levels)
  ) %>%
  dplyr::arrange(condition, donor) %>%
  dplyr::mutate(condition = as.character(condition))

pseudobulk_metadata <- pseudobulk_metadata %>%
  dplyr::left_join(donor_stage_map, by = "donor") %>%
  dplyr::mutate(condition_stage = condition)

if (anyNA(pseudobulk_metadata$condition_stage)) {
  stop("At least one pseudobulk could not be assigned to a fibrosis stage.")
}

rownames(pseudobulk_metadata) <- pseudobulk_metadata$pb_id

pseudobulk_object <- CreateSeuratObject(
  counts = pseudobulk_counts,
  assay = assay_use,
  project = "FLINC_Fig3_CCC_pseudobulk"
)

pseudobulk_object <- AddMetaData(
  pseudobulk_object,
  pseudobulk_metadata[colnames(pseudobulk_object), , drop = FALSE]
)

Idents(pseudobulk_object) <- "celltype"
DefaultAssay(pseudobulk_object) <- assay_use

pseudobulk_object <- NormalizeData(
  pseudobulk_object,
  normalization.method = "LogNormalize",
  scale.factor = pseudobulk_scale_factor,
  verbose = FALSE
)

rm(
  spatial_counts,
  proportions,
  spot_to_donor,
  weight_list,
  weights,
  pseudobulk_counts
)
invisible(gc())

message("Running LIANA with OmniPath and methods: ", paste(methods_use, collapse = ", "))

liana_coldata <- S4Vectors::DataFrame(
  celltype = factor(
    pseudobulk_metadata[colnames(pseudobulk_object), "celltype"],
    levels = cell_types
  ),
  row.names = colnames(pseudobulk_object)
)

pseudobulk_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts = SeuratObject::LayerData(
      pseudobulk_object[[assay_use]],
      layer = "counts"
    ),
    logcounts = SeuratObject::LayerData(
      pseudobulk_object[[assay_use]],
      layer = "data"
    )
  ),
  colData = liana_coldata
)

liana_raw <- liana::liana_wrap(
  pseudobulk_sce,
  idents_col = "celltype",
  resource = resource_use,
  species = species_use,
  method = methods_use,
  base = exp(1)
)

message("Aggregating LIANA method results...")

liana_aggregate <- liana::liana_aggregate(liana_raw)

rm(pseudobulk_sce, liana_coldata)
invisible(gc())

required_aggregate_columns <- c(
  "source", "target", "ligand.complex", "receptor.complex", "aggregate_rank"
)

missing_aggregate_columns <- setdiff(
  required_aggregate_columns,
  colnames(liana_aggregate)
)

if (length(missing_aggregate_columns)) {
  stop(
    "LIANA aggregate output lacks required columns: ",
    paste(missing_aggregate_columns, collapse = ", ")
  )
}

liana_significant <- liana_aggregate %>%
  dplyr::filter(
    is.finite(aggregate_rank),
    aggregate_rank <= significance_rank_threshold
  )

if (!nrow(liana_aggregate)) {
  stop("LIANA aggregate output is empty.")
}

if (!nrow(liana_significant)) {
  stop("No LIANA interactions passed aggregate_rank <= 0.01.")
}

finished_at <- Sys.time()

parameter_table <- data.frame(
  parameter = c(
    "species",
    "expression_assay",
    "expression_input_layer",
    "deconvolution_assay",
    "deconvolution_layer",
    "cell_types",
    "sample_column",
    "condition_column",
    "stage_levels",
    "resource",
    "methods",
    "liana_input_class",
    "liana_logarithm_base",
    "pseudobulk_normalization_method",
    "pseudobulk_scale_factor",
    "significance_rank_threshold",
    "input_spots",
    "input_spatial_genes",
    "zero_count_genes_dropped",
    "retained_spatial_genes",
    "donors",
    "pseudobulks",
    "maximum_deconvolution_sum_error"
  ),
  value = c(
    species_use,
    assay_use,
    input_layer,
    deconv_assay,
    deconv_layer,
    paste(cell_types, collapse = ";"),
    sample_column,
    condition_column,
    paste(stage_levels, collapse = ";"),
    resource_use,
    paste(methods_use, collapse = ";"),
    "SingleCellExperiment",
    format(exp(1), digits = 17),
    "LogNormalize",
    format(pseudobulk_scale_factor, scientific = FALSE),
    format(significance_rank_threshold, scientific = FALSE),
    nrow(metadata),
    expected_input_spatial_genes,
    zero_genes_dropped,
    expected_retained_spatial_genes,
    length(donor_levels),
    nrow(pseudobulk_metadata),
    format(proportion_sum_error, scientific = TRUE)
  ),
  stringsAsFactors = FALSE
)

input_manifest <- data.frame(
  role = "validated_annotated_figure1_object",
  path = normalizePath(input_rds, mustWork = TRUE),
  bytes = input_info$size,
  modified = format(input_info$mtime, tz = "UTC", usetz = TRUE),
  sha256 = input_sha256,
  stringsAsFactors = FALSE
)

session_lines <- capture.output(sessionInfo())

message("Saving retained Figure 3 inference outputs...")

output_paths <- c(
  pseudobulk_rds = file.path(output_dir, "01_pseudobulk_normalized.rds"),
  liana_raw_rds = file.path(output_dir, "01_liana_raw.rds"),
  aggregate_rds = file.path(output_dir, "01_liana_aggregate.rds"),
  aggregate_csv = file.path(output_dir, "01_liana_aggregate.csv"),
  significant_rds = file.path(output_dir, "01_liana_significant_rank_le_0.01.rds"),
  significant_csv = file.path(output_dir, "01_liana_significant_rank_le_0.01.csv"),
  pseudobulk_metadata_csv = file.path(output_dir, "01_pseudobulk_metadata.csv"),
  donor_stage_csv = file.path(output_dir, "01_donor_stage_map.csv"),
  parameters_tsv = file.path(output_dir, "01_analysis_parameters.tsv"),
  input_manifest_tsv = file.path(output_dir, "01_input_manifest.tsv"),
  session_info_txt = file.path(output_dir, "01_sessionInfo.txt")
)

atomic_save_rds(pseudobulk_object, output_paths[["pseudobulk_rds"]])
atomic_save_rds(liana_raw, output_paths[["liana_raw_rds"]])
atomic_save_rds(liana_aggregate, output_paths[["aggregate_rds"]])
atomic_write_csv(liana_aggregate, output_paths[["aggregate_csv"]])
atomic_save_rds(liana_significant, output_paths[["significant_rds"]])
atomic_write_csv(liana_significant, output_paths[["significant_csv"]])
atomic_write_csv(pseudobulk_metadata, output_paths[["pseudobulk_metadata_csv"]])
atomic_write_csv(donor_stage_map, output_paths[["donor_stage_csv"]])
atomic_write_tsv(parameter_table, output_paths[["parameters_tsv"]])
atomic_write_tsv(input_manifest, output_paths[["input_manifest_tsv"]])
atomic_write_lines(session_lines, output_paths[["session_info_txt"]])

output_checksums <- data.frame(
  file = basename(unname(output_paths)),
  bytes = unname(file.info(unname(output_paths))$size),
  sha256 = vapply(unname(output_paths), sha256_file, character(1)),
  stringsAsFactors = FALSE
)

checksum_path <- file.path(output_dir, "01_OUTPUT_SHA256.tsv")
atomic_write_tsv(output_checksums, checksum_path)

completion_table <- data.frame(
  field = c(
    "status",
    "started_at",
    "finished_at",
    "elapsed_minutes",
    "script_path",
    "script_sha256",
    "input_sha256",
    "aggregate_interactions",
    "significant_interactions",
    "pseudobulks",
    "retained_spatial_genes"
  ),
  value = c(
    "complete",
    format(started_at, tz = "UTC", usetz = TRUE),
    format(finished_at, tz = "UTC", usetz = TRUE),
    sprintf("%.3f", as.numeric(difftime(finished_at, started_at, units = "mins"))),
    script_path,
    if (is.na(script_path)) NA_character_ else sha256_file(script_path),
    input_sha256,
    nrow(liana_aggregate),
    nrow(liana_significant),
    nrow(pseudobulk_metadata),
    nrow(pseudobulk_object)
  ),
  stringsAsFactors = FALSE
)

completion_path <- file.path(output_dir, "01_COMPLETION.tsv")
atomic_write_tsv(completion_table, completion_path)

message("LIANA inference completed successfully at ", format(finished_at))
message("Aggregate interactions: ", nrow(liana_aggregate))
message("Significant interactions (aggregate_rank <= 0.01): ", nrow(liana_significant))
message("Outputs: ", output_dir)
