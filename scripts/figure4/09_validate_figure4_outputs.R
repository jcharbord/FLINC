args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L || any(!startsWith(args[seq(1L, length(args), 2L)], "--"))) {
  stop("Arguments must be supplied as --key value pairs", call. = FALSE)
}
arg_list <- setNames(
  args[seq(2L, length(args), 2L)],
  sub("^--", "", args[seq(1L, length(args), 2L)])
)
required <- c("script-dir", "resource-dir", "manifest", "record-dir")
missing <- setdiff(required, names(arg_list))
if (length(missing) > 0L) {
  stop("Missing argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package digest is required", call. = FALSE)
}

script_dir <- normalizePath(arg_list[["script-dir"]], mustWork = TRUE)
resource_dir <- normalizePath(arg_list[["resource-dir"]], mustWork = TRUE)
manifest_path <- normalizePath(arg_list[["manifest"]], mustWork = FALSE)
record_dir <- normalizePath(arg_list[["record-dir"]], mustWork = FALSE)
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)

expected_scripts <- c(
  "01_run_progeny_activity.R",
  "02_plot_progeny_spatial_maps.R",
  "03_run_dorothea_activity.R",
  "04_build_rfx5_regulon.R",
  "05_rfx5_celltype_associations.R",
  "06_build_ego_networks.R",
  "07_run_secreted_factor_analysis.R",
  "08_plot_secreted_enrichment.R",
  "09_validate_figure4_outputs.R",
  "10_plot_secreted_factor_spatial_context_v6.R",
  "11_plot_SM1_SM6_selected4_BH30.R",
  "README.md",
  "ego_network_functional_context_annotations.tsv"
)
expected_resources <- c(
  "uniprot_human_secreted_only_KW-0964_reviewed.tsv",
  "reference_secreted_counts_by_module.csv",
  "reference_top30_secreted_factors_by_module.csv",
  "reference_secreted_enrichment_top2.csv"
)

script_paths <- file.path(script_dir, expected_scripts)
resource_paths <- file.path(resource_dir, expected_resources)
missing_paths <- c(script_paths[!file.exists(script_paths)], resource_paths[!file.exists(resource_paths)])
if (length(missing_paths) > 0L) {
  stop("Public Figure 4 package is missing: ", paste(missing_paths, collapse = "; "), call. = FALSE)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256")
file_record <- function(path, kind) {
  data.frame(
    kind = kind,
    path = path,
    bytes = as.numeric(file.info(path)$size),
    sha256 = sha256(path),
    stringsAsFactors = FALSE
  )
}
manifest <- rbind(
  do.call(rbind, lapply(script_paths, file_record, kind = "script_or_annotation")),
  do.call(rbind, lapply(resource_paths, file_record, kind = "versioned_resource"))
)
repro_root <- normalizePath(dirname(dirname(script_dir)), mustWork = TRUE)
manifest$path <- sub(paste0("^", repro_root, "/"), "", manifest$path)

atomic_write <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.table(x, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
}

atomic_write_lines <- function(lines, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
}

atomic_write_lines(
  paste(manifest$sha256, manifest$path, sep = "  "),
  manifest_path
)
atomic_write(manifest, file.path(record_dir, "09_public_figure4_package_manifest.tsv"))
writeLines(
  c(
    "FLINC Figure 4 public-package validation",
    "Status: PASS",
    paste("Validated files:", nrow(manifest)),
    paste("Manifest:", normalizePath(manifest_path, mustWork = TRUE)),
    paste("R:", R.version.string)
  ),
  file.path(record_dir, "09_public_figure4_package_validation.txt")
)
cat("FIGURE4_PUBLIC_PACKAGE_VALIDATION_OK\n")

