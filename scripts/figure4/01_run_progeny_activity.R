#!/usr/bin/env Rscript

# Figure 4: deterministic spot-level PROGENy activity and donor-level tests.
#
# Primary input: the regenerated FLINC Seurat object (Spatial/counts).
# Historical scored objects are deliberately not accepted by this script.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(progeny)
  library(decoupleR)
  library(digest)
})

fail <- function(...) stop(paste0(...), call. = FALSE)

parse_args <- function(x) {
  allowed <- c("input", "output-dir", "log-dir", "expected-input-sha256")
  out <- setNames(as.list(rep(NA_character_, length(allowed))), allowed)
  i <- 1L
  while (i <= length(x)) {
    token <- x[[i]]
    if (!startsWith(token, "--")) fail("Unexpected argument: ", token)
    key <- sub("^--", "", token)
    if (!key %in% allowed) fail("Unknown option: --", key)
    if (i == length(x)) fail("Missing value for --", key)
    value <- x[[i + 1L]]
    if (startsWith(value, "--")) fail("Missing value for --", key)
    out[[key]] <- value
    i <- i + 2L
  }
  required <- c("input", "output-dir", "log-dir")
  missing <- required[is.na(unlist(out[required])) | !nzchar(unlist(out[required]))]
  if (length(missing)) {
    fail(
      "Missing required option(s): ", paste0("--", missing, collapse = ", "),
      "\nUsage: Rscript 01_run_progeny_activity.R ",
      "--input OBJECT.rds --output-dir RESULT_DIR --log-dir LOG_DIR ",
      "[--expected-input-sha256 HEX]"
    )
  }
  out
}

assert_scalar_path <- function(x, label) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) fail(label, " is empty")
  if (grepl("[\r\n]", x)) fail(label, " contains a newline")
  invisible(x)
}

atomic_write_table <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(tmp, "wt") else file(tmp, "wt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  write.table(
    x, file = con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  close(con)
  if (!file.rename(tmp, path)) fail("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_write_lines <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) fail("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_save_rds <- function(x, path, compress = "xz") {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = compress, version = 3)
  if (!file.rename(tmp, path)) fail("Atomic rename failed for: ", path)
  invisible(path)
}

sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE)

package_table <- function(packages) {
  do.call(rbind, lapply(packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) fail("Required package is missing: ", pkg)
    data.frame(
      package = pkg,
      version = as.character(packageVersion(pkg)),
      library_path = find.package(pkg),
      stringsAsFactors = FALSE
    )
  }))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_path <- normalizePath(args[["input"]], mustWork = TRUE)
result_final <- normalizePath(args[["output-dir"]], mustWork = FALSE)
log_final <- normalizePath(args[["log-dir"]], mustWork = FALSE)
expected_input_sha <- args[["expected-input-sha256"]]

assert_scalar_path(input_path, "input")
assert_scalar_path(result_final, "output-dir")
assert_scalar_path(log_final, "log-dir")
if (identical(result_final, log_final)) fail("output-dir and log-dir must differ")
if (dir.exists(result_final) || file.exists(result_final)) {
  fail("Refusing to overwrite existing output: ", result_final)
}
if (dir.exists(log_final) || file.exists(log_final)) {
  fail("Refusing to overwrite existing log directory: ", log_final)
}

dir.create(dirname(result_final), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_final), recursive = TRUE, showWarnings = FALSE)

run_tag <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-pid", Sys.getpid())
result_tmp <- paste0(result_final, ".incomplete-", run_tag)
log_tmp <- paste0(log_final, ".incomplete-", run_tag)
if (file.exists(result_tmp) || file.exists(log_tmp)) fail("Temporary run path already exists")
dir.create(result_tmp, recursive = FALSE, showWarnings = FALSE)
dir.create(log_tmp, recursive = FALSE, showWarnings = FALSE)

completed <- FALSE
on.exit({
  if (!completed) {
    message("Run did not publish final directories.")
    message("Inspect incomplete result directory: ", result_tmp)
    message("Inspect incomplete log directory: ", log_tmp)
  }
}, add = TRUE)

start_time <- Sys.time()
start_utc <- format(start_time, tz = "UTC", usetz = TRUE)
message("Computing SHA256 for the regenerated input object...")
input_sha <- sha256_file(input_path)
if (!is.na(expected_input_sha) && nzchar(expected_input_sha)) {
  expected_input_sha <- tolower(expected_input_sha)
  if (!grepl("^[0-9a-f]{64}$", expected_input_sha)) {
    fail("expected-input-sha256 is not a 64-character hexadecimal SHA256")
  }
  if (!identical(tolower(input_sha), expected_input_sha)) {
    fail("Input SHA256 mismatch. Expected ", expected_input_sha, "; observed ", input_sha)
  }
}

parameters <- data.frame(
  parameter = c(
    "input_object", "input_sha256", "assay", "input_layer",
    "zero_count_gene_rule", "normalization", "scale_factor",
    "gene_standardization", "progeny_organism", "progeny_top_targets",
    "decoupler_method", "selected_statistic", "minsize", "times", "seed",
    "randomize_type", "sparse_argument", "donor_summary",
    "primary_test", "multiple_testing", "effect_size",
    "main_contrast", "supplementary_contrast",
    "validation_only_contrasts", "blas_threads"
  ),
  value = c(
    input_path, input_sha, "Spatial", "counts",
    "remove genes with rowSums(counts) == 0",
    "log1p(count / spot_library_size * scale_factor)", "10000",
    "per-gene z-score across all spots; exclude zero-variance targets",
    "Human", "100", "decoupleR::run_wmean", "norm_wmean",
    "5", "1000", "42", "rows", "FALSE", "median across spots per donor",
    "two-sided exact Wilcoxon rank-sum when no ties; asymptotic fallback recorded",
    "Benjamini-Hochberg within each 14-pathway contrast",
    "mean(stage donor medians) - mean(F0 donor medians)",
    "F4_vs_F0", "F1_vs_F0", "F2_vs_F0;F3_vs_F0", "1"
  ),
  stringsAsFactors = FALSE
)

packages <- package_table(c(
  "Seurat", "SeuratObject", "Matrix", "progeny", "decoupleR", "digest"
))
expected_package_versions <- c(
  Seurat = "5.4.0",
  SeuratObject = "5.3.0",
  Matrix = "1.7.3",
  progeny = "1.28.0",
  decoupleR = "2.12.0",
  digest = "0.6.39"
)
observed_package_versions <- setNames(packages$version, packages$package)
package_mismatch <- names(expected_package_versions)[
  observed_package_versions[names(expected_package_versions)] != expected_package_versions
]
if (length(package_mismatch)) {
  fail(
    "Package-version mismatch: ",
    paste(
      paste0(
        package_mismatch, " expected=", expected_package_versions[package_mismatch],
        " observed=", observed_package_versions[package_mismatch]
      ),
      collapse = "; "
    )
  )
}

message("Loading regenerated Seurat object read-only...")
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) fail("Input is not a Seurat object")
if (!"Spatial" %in% Assays(obj)) fail("Spatial assay is missing")
if (!"counts" %in% Layers(obj[["Spatial"]])) fail("Spatial/counts layer is missing")

meta <- obj[[]]
required_meta <- c("donor", "F_score")
if (!all(required_meta %in% colnames(meta))) {
  fail("Missing metadata: ", paste(setdiff(required_meta, colnames(meta)), collapse = ", "))
}
counts <- LayerData(obj[["Spatial"]], layer = "counts")
if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
if (!identical(colnames(counts), rownames(meta))) fail("Counts and metadata spot order differ")
if (nrow(counts) != 17922L || ncol(counts) != 26056L) {
  fail(
    "Unexpected Spatial/counts dimensions: ",
    paste(dim(counts), collapse = " x "), " (expected 17922 x 26056)"
  )
}
if (length(obj@images) != 24L) fail("Expected 24 Visium images")

object_summary <- data.frame(
  metric = c(
    "object_class", "spots", "spatial_features_before_filter",
    "spatial_nonzero_entries", "images", "donors", "stages"
  ),
  value = c(
    paste(class(obj), collapse = ";"), ncol(counts), nrow(counts),
    Matrix::nnzero(counts), length(obj@images),
    length(unique(meta$donor)), paste(sort(unique(as.character(meta$F_score))), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
rm(obj)
invisible(gc())

donor_stage <- unique(meta[, c("donor", "F_score"), drop = FALSE])
if (anyDuplicated(donor_stage$donor)) fail("At least one donor maps to multiple stages")
stage_levels <- c("F0", "F1", "F2", "F3", "F4")
if (!setequal(as.character(donor_stage$F_score), stage_levels)) {
  fail("Unexpected stage set: ", paste(sort(unique(donor_stage$F_score)), collapse = ", "))
}
donor_stage$F_score <- factor(donor_stage$F_score, levels = stage_levels, ordered = TRUE)
donor_stage <- donor_stage[order(donor_stage$F_score, donor_stage$donor), , drop = FALSE]

gene_totals <- Matrix::rowSums(counts)
keep_gene <- is.finite(gene_totals) & gene_totals > 0
zero_count_genes <- rownames(counts)[!keep_gene]
counts <- counts[keep_gene, , drop = FALSE]
if (length(zero_count_genes) != 147L || nrow(counts) != 17775L) {
  fail(
    "Unexpected zero-count filtering result: removed ", length(zero_count_genes),
    " and retained ", nrow(counts), " genes (expected 147 and 17775)"
  )
}
library_sizes <- Matrix::colSums(counts)
if (any(!is.finite(library_sizes) | library_sizes <= 0)) fail("Non-positive spot library size")

model <- as.data.frame(progeny::getModel(organism = "Human", top = 100L))
model <- model[order(rownames(model)), sort(colnames(model)), drop = FALSE]
model_sha_serialized <- digest(model, algo = "sha256", serialize = TRUE)
expected_model_sha_serialized <- "75a780f64a38e169d26db873a6a23e84d96b5ec9c5aed815783d827bd06fe494"
if (!identical(model_sha_serialized, expected_model_sha_serialized)) {
  fail(
    "PROGENy model checksum mismatch. Expected ", expected_model_sha_serialized,
    "; observed ", model_sha_serialized
  )
}
model_matrix <- as.matrix(model)
model_index <- which(is.finite(model_matrix) & model_matrix != 0, arr.ind = TRUE)
network <- data.frame(
  source = colnames(model_matrix)[model_index[, "col"]],
  target = rownames(model_matrix)[model_index[, "row"]],
  weight = as.numeric(model_matrix[model_index]),
  stringsAsFactors = FALSE
)
network <- network[order(network$source, network$target, network$weight), , drop = FALSE]
rownames(network) <- NULL
if (nrow(network) != 1400L || length(unique(network$source)) != 14L) {
  fail("Unexpected PROGENy top-100 model dimensions")
}

available_targets <- sort(intersect(unique(network$target), rownames(counts)))
if (length(available_targets) < 5L) fail("Fewer than five PROGENy targets are present")
message("Normalizing ", length(available_targets), " observed PROGENy target genes...")
target_counts <- as.matrix(counts[available_targets, , drop = FALSE])
target_normalized <- log1p(sweep(target_counts, 2L, library_sizes, "/") * 10000)
rm(target_counts, counts)
invisible(gc())

target_means <- rowMeans(target_normalized)
centered <- sweep(target_normalized, 1L, target_means, "-")
target_sd <- sqrt(rowSums(centered * centered) / (ncol(centered) - 1L))
informative <- is.finite(target_sd) & target_sd > 0
zero_variance_targets <- rownames(centered)[!informative]
target_z <- sweep(centered[informative, , drop = FALSE], 1L, target_sd[informative], "/")
rm(target_normalized, centered)
invisible(gc())
if (any(!is.finite(target_z))) fail("Non-finite values remain after z-scoring")

network_scored <- network[network$target %in% rownames(target_z), , drop = FALSE]
coverage <- do.call(rbind, lapply(sort(unique(network$source)), function(pathway) {
  expected <- unique(network$target[network$source == pathway])
  available <- intersect(expected, available_targets)
  scored <- intersect(expected, rownames(target_z))
  data.frame(
    pathway = pathway,
    expected_top_targets = length(expected),
    observed_targets = length(available),
    zero_variance_targets = length(setdiff(available, scored)),
    scored_targets = length(scored),
    stringsAsFactors = FALSE
  )
}))
if (any(coverage$scored_targets < 5L)) {
  fail("At least one pathway has fewer than five scored targets")
}

message("Running decoupleR::run_wmean with explicit norm_wmean selection...")
activity_long_all <- decoupleR::run_wmean(
  mat = target_z,
  network = network_scored,
  .source = "source",
  .target = "target",
  .mor = "weight",
  times = 1000L,
  seed = 42L,
  sparse = FALSE,
  randomize_type = "rows",
  minsize = 5L
)
expected_statistics <- c("corr_wmean", "norm_wmean", "wmean")
if (!setequal(unique(activity_long_all$statistic), expected_statistics)) {
  fail("Unexpected run_wmean statistic set: ", paste(unique(activity_long_all$statistic), collapse = ", "))
}
activity_long <- activity_long_all[activity_long_all$statistic == "norm_wmean", , drop = FALSE]
pathways <- sort(unique(network$source))
spots <- colnames(target_z)
if (nrow(activity_long) != length(pathways) * length(spots)) {
  fail("norm_wmean result does not contain exactly one score per pathway and spot")
}
if (anyDuplicated(paste(activity_long$source, activity_long$condition, sep = "\r"))) {
  fail("Duplicate pathway/spot scores returned")
}

score_matrix <- matrix(
  NA_real_, nrow = length(pathways), ncol = length(spots),
  dimnames = list(pathways, spots)
)
row_index <- match(activity_long$source, pathways)
col_index <- match(activity_long$condition, spots)
if (anyNA(row_index) || anyNA(col_index)) fail("Could not align returned activity scores")
score_matrix[cbind(row_index, col_index)] <- activity_long$score
if (any(!is.finite(score_matrix))) fail("Activity matrix contains missing or non-finite scores")

meta <- meta[spots, , drop = FALSE]
spot_key <- data.frame(
  barcode = spots,
  donor = as.character(meta$donor),
  F_score = as.character(meta$F_score),
  stringsAsFactors = FALSE
)

donor_rows <- vector("list", length(pathways) * nrow(donor_stage))
k <- 0L
for (pathway in pathways) {
  values <- score_matrix[pathway, ]
  for (donor in donor_stage$donor) {
    k <- k + 1L
    selected <- spot_key$donor == donor
    donor_rows[[k]] <- data.frame(
      pathway = pathway,
      donor = donor,
      F_stage = spot_key$F_score[which(selected)[1L]],
      n_spots = sum(selected),
      median_activity = median(values[selected]),
      stringsAsFactors = FALSE
    )
  }
}
donor_medians <- do.call(rbind, donor_rows)
if (nrow(donor_medians) != 14L * 24L) fail("Unexpected donor-median table size")

contrast_stage <- c("F1", "F2", "F3", "F4")
contrast_rows <- list()
k <- 0L
for (stage_k in contrast_stage) {
  for (pathway in pathways) {
    x <- donor_medians$median_activity[
      donor_medians$pathway == pathway & donor_medians$F_stage == stage_k
    ]
    y <- donor_medians$median_activity[
      donor_medians$pathway == pathway & donor_medians$F_stage == "F0"
    ]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    has_ties <- anyDuplicated(c(x, y)) > 0L
    exact_available <- !has_ties
    p_primary <- if (exact_available) {
      wilcox.test(x, y, alternative = "two.sided", exact = TRUE)$p.value
    } else {
      suppressWarnings(wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value)
    }
    p_asymptotic <- suppressWarnings(
      wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value
    )
    k <- k + 1L
    contrast_rows[[k]] <- data.frame(
      contrast = paste0(stage_k, "_vs_F0"),
      figure_role = if (stage_k == "F4") "main" else if (stage_k == "F1") "supplementary" else "validation_only",
      pathway = pathway,
      n_F0_donors = length(y),
      n_FK_donors = length(x),
      mean_F0_donor_medians = mean(y),
      mean_FK_donor_medians = mean(x),
      delta_mean_donor_medians = mean(x) - mean(y),
      ties_present = has_ties,
      primary_test = if (exact_available) "exact_wilcoxon_rank_sum" else "asymptotic_wilcoxon_rank_sum_ties_fallback",
      p_primary = p_primary,
      p_asymptotic_sensitivity = p_asymptotic,
      stringsAsFactors = FALSE
    )
  }
}
contrast_stats <- do.call(rbind, contrast_rows)
contrast_stats$fdr_BH_primary <- ave(
  contrast_stats$p_primary,
  contrast_stats$contrast,
  FUN = function(p) p.adjust(p, method = "BH")
)
contrast_stats$fdr_BH_asymptotic_sensitivity <- ave(
  contrast_stats$p_asymptotic_sensitivity,
  contrast_stats$contrast,
  FUN = function(p) p.adjust(p, method = "BH")
)
contrast_stats <- contrast_stats[
  order(match(contrast_stats$contrast, paste0(contrast_stage, "_vs_F0")), contrast_stats$fdr_BH_primary, contrast_stats$pathway),
  , drop = FALSE
]

spot_activity_long <- data.frame(
  pathway = rep(rownames(score_matrix), each = ncol(score_matrix)),
  barcode = rep(colnames(score_matrix), times = nrow(score_matrix)),
  activity_norm_wmean = as.vector(t(score_matrix)),
  stringsAsFactors = FALSE
)
spot_match <- match(spot_activity_long$barcode, spot_key$barcode)
spot_activity_long$donor <- spot_key$donor[spot_match]
spot_activity_long$F_stage <- spot_key$F_score[spot_match]

model_table <- data.frame(target = rownames(model), model, check.names = FALSE)
run_object <- list(
  activity = score_matrix,
  spot_metadata = spot_key,
  parameters = setNames(parameters$value, parameters$parameter),
  model_serialized_sha256 = model_sha_serialized,
  model_coverage = coverage
)

atomic_save_rds(run_object, file.path(result_tmp, "01_progeny_spot_activity.rds"))
atomic_write_table(spot_activity_long, file.path(result_tmp, "01_progeny_spot_activity.tsv.gz"))
atomic_write_table(donor_medians, file.path(result_tmp, "01_progeny_donor_medians.tsv"))
atomic_write_table(contrast_stats, file.path(result_tmp, "01_progeny_contrast_statistics.tsv"))
for (contrast_i in unique(contrast_stats$contrast)) {
  atomic_write_table(
    contrast_stats[contrast_stats$contrast == contrast_i, , drop = FALSE],
    file.path(result_tmp, paste0("01_PROGENy_", contrast_i, ".tsv"))
  )
}
atomic_write_table(model_table, file.path(result_tmp, "01_progeny_Human_top100_model.tsv"))
atomic_write_table(network, file.path(result_tmp, "01_progeny_Human_top100_network.tsv"))
atomic_write_table(coverage, file.path(result_tmp, "01_progeny_model_coverage.tsv"))
atomic_write_table(
  data.frame(
    gene = zero_count_genes,
    reason = rep("zero_total_counts", length(zero_count_genes)),
    stringsAsFactors = FALSE
  ),
  file.path(result_tmp, "01_removed_zero_count_genes.tsv")
)
atomic_write_table(
  data.frame(
    gene = zero_variance_targets,
    reason = rep(
      "zero_variance_after_log_normalization",
      length(zero_variance_targets)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(result_tmp, "01_removed_zero_variance_targets.tsv")
)
atomic_write_table(donor_stage, file.path(result_tmp, "01_donor_stage_map.tsv"))
atomic_write_table(parameters, file.path(result_tmp, "01_run_parameters.tsv"))
atomic_write_table(packages, file.path(result_tmp, "01_package_versions.tsv"))
atomic_write_table(object_summary, file.path(result_tmp, "01_input_object_summary.tsv"))

session_lines <- capture.output(sessionInfo())
atomic_write_lines(session_lines, file.path(result_tmp, "01_sessionInfo.txt"))

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
result_manifest <- data.frame(
  file = basename(result_files),
  bytes = as.numeric(file.info(result_files)$size),
  sha256 = vapply(result_files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
atomic_write_table(result_manifest, file.path(result_tmp, "01_output_SHA256.tsv"))

significance_summary <- do.call(rbind, lapply(unique(contrast_stats$contrast), function(contrast_i) {
  x <- contrast_stats[contrast_stats$contrast == contrast_i, , drop = FALSE]
  sig <- x$pathway[is.finite(x$fdr_BH_primary) & x$fdr_BH_primary <= 0.05]
  data.frame(
    contrast = contrast_i,
    figure_role = unique(x$figure_role),
    pathways_tested = nrow(x),
    significant_pathways_FDR_le_0.05 = length(sig),
    significant_pathway_names = paste(sig, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

end_time <- Sys.time()
summary_lines <- c(
  "FLINC Figure 4 PROGENy execution summary",
  paste0("Started UTC: ", start_utc),
  paste0("Completed UTC: ", format(end_time, tz = "UTC", usetz = TRUE)),
  paste0("Elapsed seconds: ", signif(as.numeric(difftime(end_time, start_time, units = "secs")), 10)),
  paste0("Input object: ", input_path),
  paste0("Input SHA256: ", input_sha),
  paste0("Input spots: ", ncol(score_matrix)),
  paste0("Input Spatial features before zero-count filtering: ", object_summary$value[object_summary$metric == "spatial_features_before_filter"]),
  paste0("Removed zero-count genes: ", length(zero_count_genes)),
  paste0("Observed model targets: ", length(available_targets)),
  paste0("Removed zero-variance model targets: ", length(zero_variance_targets)),
  paste0("Scored pathways: ", nrow(score_matrix)),
  paste0("Selected decoupleR statistic: norm_wmean"),
  paste0("PROGENy model serialized SHA256: ", model_sha_serialized),
  paste0("Primary main contrast: F4_vs_F0"),
  paste0("Supplementary contrast: F1_vs_F0"),
  paste0("Command: ", paste(commandArgs(trailingOnly = FALSE), collapse = " "))
)
atomic_write_lines(summary_lines, file.path(log_tmp, "01_execution_summary.txt"))
atomic_write_table(parameters, file.path(log_tmp, "01_parameters.tsv"))
atomic_write_table(packages, file.path(log_tmp, "01_package_versions.tsv"))
atomic_write_table(significance_summary, file.path(log_tmp, "01_significance_summary.tsv"))
atomic_write_table(coverage, file.path(log_tmp, "01_model_coverage.tsv"))
atomic_write_lines(session_lines, file.path(log_tmp, "01_sessionInfo.txt"))

combined_hashes <- rbind(
  data.frame(
    role = "input_object", file = input_path,
    bytes = as.numeric(file.info(input_path)$size), sha256 = input_sha,
    stringsAsFactors = FALSE
  ),
  data.frame(
    role = "generated_output", file = file.path(result_final, result_manifest$file),
    bytes = result_manifest$bytes, sha256 = result_manifest$sha256,
    stringsAsFactors = FALSE
  )
)
atomic_write_table(combined_hashes, file.path(log_tmp, "01_input_output_SHA256.tsv"))
atomic_write_table(
  data.frame(
    resource = c("PROGENy_Human_top100_serialized_model", "DoRothEA_not_used_in_script_01"),
    package = c("progeny", "dorothea"),
    package_version = c(as.character(packageVersion("progeny")), as.character(packageVersion("dorothea"))),
    sha256 = c(model_sha_serialized, NA_character_),
    stringsAsFactors = FALSE
  ),
  file.path(log_tmp, "01_resource_versions.tsv")
)
atomic_write_lines("SUCCESS", file.path(result_tmp, "01_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "01_COMPLETED.txt"))

if (!file.rename(result_tmp, result_final)) fail("Could not atomically publish result directory")
if (!file.rename(log_tmp, log_final)) {
  rollback_ok <- file.rename(result_final, result_tmp)
  fail("Could not atomically publish log directory; result rollback status: ", rollback_ok)
}

completed <- TRUE
message("PROGENy analysis completed successfully.")
message("Results: ", result_final)
message("Compact logs: ", log_final)
