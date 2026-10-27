#!/usr/bin/env Rscript

# FLINC Figure 3 and Supplementary Figure 8 - final reproducibility audit.
#
# This script checks the completeness and internal consistency of operations
# 01 through 05. It performs no additional biological analysis. Generated
# figures and large result objects remain under ignored results/ directories;
# only compact execution records are copied into reproducibility_logs/.

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for SHA-256 verification.", call. = FALSE)
}

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_csv <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

safe_table <- function(path, reader) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, data = NULL, detail = paste("missing:", path)))
  }
  tryCatch(
    list(ok = TRUE, data = reader(path), detail = "read successfully"),
    error = function(e) list(
      ok = FALSE,
      data = NULL,
      detail = paste("could not read", path, "-", conditionMessage(e))
    )
  )
}

safe_rds <- function(path) safe_table(path, readRDS)

nonempty_files <- function(paths) {
  if (!length(paths)) return(logical(0))
  exists <- file.exists(paths) & !dir.exists(paths)
  sizes <- rep(NA_real_, length(paths))
  sizes[exists] <- file.info(paths[exists])$size
  exists & !is.na(sizes) & sizes > 0
}

safe_relative_paths <- function(paths) {
  paths <- as.character(paths)
  !is.na(paths) & nzchar(paths) &
    !grepl("^/", paths) &
    !grepl("(^|/)\\.\\.(/|$)", paths)
}

numeric_equal <- function(x, y, tolerance = 1e-12) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  length(x) == length(y) &&
    !anyNA(x) && !anyNA(y) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(y)))
}

key_from <- function(x, columns) {
  do.call(paste, c(x[columns], sep = "||"))
}

parameter_values <- function(x) {
  if (!all(c("parameter", "value") %in% colnames(x)) ||
      anyDuplicated(x$parameter)) {
    return(character(0))
  }
  stats::setNames(as.character(x$value), as.character(x$parameter))
}

hash_cache <- new.env(parent = emptyenv(), hash = TRUE)
sha256_file <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!exists(normalized, envir = hash_cache, inherits = FALSE)) {
    assign(
      normalized,
      digest::digest(
        object = normalized,
        algo = "sha256",
        serialize = FALSE,
        file = TRUE
      ),
      envir = hash_cache
    )
  }
  get(normalized, envir = hash_cache, inherits = FALSE)
}

atomic_write_table <- function(x, path, sep = "\t") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stopf("Refusing to replace an existing audit file: %s", path)
  }
  utils::write.table(
    x,
    partial,
    sep = sep,
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
  if (!file.rename(partial, path)) stopf("Could not finalize audit file: %s", path)
  invisible(path)
}

atomic_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stopf("Refusing to replace an existing audit file: %s", path)
  }
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stopf("Could not finalize audit file: %s", path)
  invisible(path)
}

atomic_copy <- function(source, target) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(target, ".partial")
  if (file.exists(target) || file.exists(partial)) {
    stopf("Refusing to replace an existing compact record: %s", target)
  }
  copied <- file.copy(
    source,
    partial,
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!isTRUE(copied) || !file.rename(partial, target)) {
    stopf("Could not copy compact record to: %s", target)
  }
  invisible(target)
}

repro_root <- env_or(
  "FLINC_REPRO_ROOT",
  "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)
repo_root <- normalizePath(file.path(repro_root, "..", ".."), mustWork = TRUE)
scripts_dir <- file.path(repro_root, "scripts", "figure3")
results_root <- file.path(repro_root, "results", "Fig3_CCC")
step01_dir <- file.path(results_root, "01_liana_inference")
step02_dir <- file.path(results_root, "02_main_and_supplementary_figures")
step04_dir <- file.path(results_root, "03_spatial_LR_all_samples")
step05_dir <- file.path(results_root, "04_main_figure3C_selected")
audit_parent <- file.path(repro_root, "reproducibility_logs", "Fig3_CCC")
audit_dir <- env_or(
  "FIG3_AUDIT_OUTPUT_DIR",
  file.path(audit_parent, "Fig3_final_audit")
)

if (dir.exists(audit_dir) && length(list.files(audit_dir, all.files = TRUE, no.. = TRUE))) {
  stopf(
    paste0(
      "Audit directory is not empty. Preserve it in a timestamped backup before ",
      "rerunning: %s"
    ),
    audit_dir
  )
}
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(audit_dir)) stopf("Could not create audit directory: %s", audit_dir)

status_rows <- list()
record_check <- function(component, pass, detail, required = TRUE) {
  pass <- isTRUE(pass)
  status_rows[[length(status_rows) + 1L]] <<- data.frame(
    component = component,
    required = required,
    pass = pass,
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  invisible(pass)
}

check_required_files <- function(component, paths, required = TRUE) {
  good <- nonempty_files(paths)
  missing <- paths[!good]
  detail <- if (!length(missing)) {
    sprintf("%d/%d required files present and non-empty", length(paths), length(paths))
  } else {
    paste("missing or empty:", paste(missing, collapse = "; "))
  }
  record_check(component, length(missing) == 0L, detail, required)
}

read_completion <- function(path) {
  result <- safe_table(path, read_tsv)
  if (!result$ok) return(result)
  x <- result$data
  if (!identical(colnames(x), c("field", "value")) || anyDuplicated(x$field)) {
    return(list(ok = FALSE, data = NULL, detail = "invalid completion schema"))
  }
  list(
    ok = TRUE,
    data = stats::setNames(as.character(x$value), as.character(x$field)),
    detail = sprintf("%d completion fields read", nrow(x))
  )
}

check_completion <- function(path, script_path, expected_numeric) {
  result <- read_completion(path)
  if (!result$ok) return(result)
  values <- result$data
  required <- c("status", "script_path", "script_sha256", names(expected_numeric))
  missing <- setdiff(required, names(values))
  if (length(missing)) {
    return(list(
      ok = FALSE,
      data = values,
      detail = paste("missing completion fields:", paste(missing, collapse = ", "))
    ))
  }
  numeric_ok <- all(vapply(
    names(expected_numeric),
    function(field) numeric_equal(values[[field]], expected_numeric[[field]], 1e-12),
    logical(1)
  ))
  script_ok <- file.exists(script_path) &&
    identical(basename(values[["script_path"]]), basename(script_path)) &&
    identical(tolower(values[["script_sha256"]]), tolower(sha256_file(script_path)))
  ok <- identical(values[["status"]], "complete") && script_ok && numeric_ok
  list(
    ok = ok,
    data = values,
    detail = sprintf(
      "status=%s; script checksum=%s; expected counts=%s",
      values[["status"]], script_ok, numeric_ok
    )
  )
}

verified_output_rows <- list()
verify_output_manifest <- function(
    step,
    directory,
    manifest_name,
    expected_rows,
    exact_directory = FALSE,
    required_prefix = "") {
  manifest_path <- file.path(directory, manifest_name)
  result <- safe_table(manifest_path, read_tsv)
  if (!result$ok) return(result)
  x <- result$data
  required <- c("file", "bytes", "sha256")
  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    return(list(ok = FALSE, data = x, detail = paste("missing columns:", paste(missing, collapse = ", "))))
  }
  paths_safe <- safe_relative_paths(x$file)
  prefix_ok <- !nzchar(required_prefix) || all(startsWith(basename(x$file), required_prefix))
  absolute <- file.path(directory, x$file)
  files_ok <- all(nonempty_files(absolute))
  bytes_ok <- files_ok && identical(
    as.numeric(file.info(absolute)$size),
    as.numeric(x$bytes)
  )
  actual_sha <- if (files_ok) vapply(absolute, sha256_file, character(1)) else rep(NA_character_, nrow(x))
  sha_ok <- files_ok && identical(
    tolower(unname(actual_sha)),
    tolower(as.character(x$sha256))
  )
  exact_ok <- TRUE
  if (exact_directory && dir.exists(directory)) {
    actual_files <- list.files(directory, recursive = TRUE, full.names = FALSE, all.files = FALSE)
    actual_files <- actual_files[
      file.info(file.path(directory, actual_files))$isdir %in% FALSE
    ]
    exact_ok <- setequal(
      setdiff(actual_files, c(manifest_name, paste0(step, "_COMPLETION.tsv"))),
      as.character(x$file)
    )
  }
  ok <- nrow(x) == expected_rows && !anyDuplicated(x$file) &&
    all(paths_safe) && prefix_ok && files_ok && bytes_ok && sha_ok && exact_ok
  verified_output_rows[[length(verified_output_rows) + 1L]] <<- data.frame(
    step = step,
    file = as.character(x$file),
    bytes = as.numeric(x$bytes),
    sha256 = tolower(as.character(x$sha256)),
    stringsAsFactors = FALSE
  )
  list(
    ok = ok,
    data = x,
    detail = sprintf(
      "%d/%d rows; files=%s; bytes=%s; SHA-256=%s; exact file set=%s",
      nrow(x), expected_rows, files_ok, bytes_ok, sha_ok, exact_ok
    )
  )
}

verified_input_rows <- list()
verify_input_manifest <- function(step, manifest_path) {
  result <- safe_table(manifest_path, read_tsv)
  if (!result$ok) return(result)
  x <- result$data
  required <- c("role", "bytes", "sha256")
  missing <- setdiff(required, colnames(x))
  path_column <- if ("relative_path" %in% colnames(x)) {
    "relative_path"
  } else if ("path" %in% colnames(x)) {
    "path"
  } else {
    NA_character_
  }
  if (length(missing) || is.na(path_column)) {
    return(list(
      ok = FALSE,
      data = x,
      detail = paste("missing columns:", paste(c(missing, if (is.na(path_column)) "path/relative_path"), collapse = ", "))
    ))
  }
  recorded_paths <- as.character(x[[path_column]])
  if (identical(path_column, "relative_path")) {
    paths_safe <- safe_relative_paths(recorded_paths)
    absolute <- file.path(repro_root, recorded_paths)
  } else {
    paths_safe <- !is.na(recorded_paths) & nzchar(recorded_paths)
    absolute <- ifelse(
      grepl("^/", recorded_paths),
      recorded_paths,
      file.path(repro_root, recorded_paths)
    )
  }
  files_ok <- all(paths_safe) && all(nonempty_files(absolute))
  bytes_ok <- files_ok && identical(
    as.numeric(file.info(absolute)$size),
    as.numeric(x$bytes)
  )
  actual_sha <- if (files_ok) vapply(absolute, sha256_file, character(1)) else rep(NA_character_, nrow(x))
  sha_ok <- files_ok && identical(
    tolower(unname(actual_sha)),
    tolower(as.character(x$sha256))
  )
  # A role can legitimately occur more than once (for example, the six
  # selected source LRlocal PDFs in Step 05). Input identity is defined by
  # the recorded path, which must remain unique.
  unique_paths_ok <- !anyDuplicated(recorded_paths)
  ok <- unique_paths_ok && files_ok && bytes_ok && sha_ok
  verified_input_rows[[length(verified_input_rows) + 1L]] <<- data.frame(
    step = step,
    role = as.character(x$role),
    recorded_path = recorded_paths,
    bytes = as.numeric(x$bytes),
    sha256 = tolower(as.character(x$sha256)),
    stringsAsFactors = FALSE
  )
  list(
    ok = ok,
    data = x,
    detail = sprintf(
      "%d inputs; unique paths=%s; files=%s; bytes=%s; SHA-256=%s",
      nrow(x), unique_paths_ok, files_ok, bytes_ok, sha_ok
    )
  )
}

pdf_magic_ok <- function(paths) {
  if (!all(nonempty_files(paths))) return(FALSE)
  all(vapply(paths, function(path) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    identical(rawToChar(readBin(connection, what = "raw", n = 4L)), "%PDF")
  }, logical(1)))
}

# -------------------------------------------------------------------------
# Public code and documentation
# -------------------------------------------------------------------------

canonical_script_names <- c(
  "01_run_liana_inference.R",
  "02_analyze_stage_activity.R",
  "03_plot_main_and_supplementary_ccc.R",
  "04_plot_spatial_lr_all_donors.R",
  "05_select_main_figure3C_panels.R",
  "06_validate_figure3_outputs.R"
)
canonical_scripts <- file.path(scripts_dir, canonical_script_names)
documentation <- c(
  file.path(repro_root, "README.md"),
  file.path(repro_root, "FIGURE3_SOURCE_WORKFLOW.md"),
  file.path(repro_root, "FIGURE3_SHA256.txt"),
  file.path(scripts_dir, "README.md")
)
check_required_files(
  "Public Figure 3 code and documentation",
  c(canonical_scripts, documentation)
)

public_expected_relative <- c(
  "README.md",
  "FIGURE3_SOURCE_WORKFLOW.md",
  file.path("scripts", "figure3", "README.md"),
  file.path("scripts", "figure3", canonical_script_names)
)
public_sha_path <- file.path(repro_root, "FIGURE3_SHA256.txt")
public_sha_ok <- FALSE
public_sha_detail <- paste("missing:", public_sha_path)
public_sha_rows <- data.frame()
if (file.exists(public_sha_path)) {
  lines <- trimws(readLines(public_sha_path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  valid_lines <- grepl("^[0-9A-Fa-f]{64}[[:space:]]+\\*?.+$", lines)
  if (length(lines) && all(valid_lines)) {
    hashes <- tolower(sub("^([0-9A-Fa-f]{64}).*$", "\\1", lines))
    paths <- sub("^[0-9A-Fa-f]{64}[[:space:]]+\\*?", "", lines)
    safe <- safe_relative_paths(paths)
    absolute <- file.path(repro_root, paths)
    files_ok <- all(safe) && all(nonempty_files(absolute))
    actual <- if (files_ok) vapply(absolute, sha256_file, character(1)) else rep(NA_character_, length(paths))
    public_sha_ok <- !anyDuplicated(paths) &&
      setequal(paths, public_expected_relative) &&
      files_ok && identical(hashes, tolower(unname(actual)))
    public_sha_detail <- sprintf(
      "%d public file checksums; complete coverage=%s; SHA-256=%s",
      length(paths), setequal(paths, public_expected_relative),
      files_ok && identical(hashes, tolower(unname(actual)))
    )
    public_sha_rows <- data.frame(
      relative_path = paths,
      bytes = if (files_ok) as.numeric(file.info(absolute)$size) else NA_real_,
      sha256 = hashes,
      stringsAsFactors = FALSE
    )
  } else {
    public_sha_detail <- "one or more FIGURE3_SHA256.txt lines are malformed"
  }
}
record_check("Public Figure 3 SHA-256 manifest", public_sha_ok, public_sha_detail)

# -------------------------------------------------------------------------
# Completion records, output manifests and input manifests
# -------------------------------------------------------------------------

completion_specs <- list(
  list(
    step = "01", directory = step01_dir, script = canonical_scripts[[1]],
    expected = c(
      aggregate_interactions = 609984,
      significant_interactions = 29205,
      pseudobulks = 288,
      retained_spatial_genes = 17775
    )
  ),
  list(
    step = "02", directory = step02_dir, script = canonical_scripts[[2]],
    expected = c(
      single_gene_candidates_tested = 3000,
      finite_spearman_tests = 2903,
      bh_significant_trends = 349,
      effect_size_pass = 315,
      selected_interactions = 60,
      chord_significant_interactions = 29205,
      cholangiocyte_participating_interactions = 15693,
      dotplot_input_tables = 4
    )
  ),
  list(
    step = "03", directory = step02_dir, script = canonical_scripts[[3]],
    expected = c(
      pdf_panels = 7,
      figure3_panels = 2,
      supplementary_figure8_panels = 5,
      selected_interactions = 60,
      chord_significant_interactions = 29205,
      cholangiocyte_partner_rows = 23
    )
  ),
  list(
    step = "04", directory = step04_dir, script = canonical_scripts[[4]],
    expected = c(
      donors = 24,
      spatial_lr_edges = 3,
      lrlocal_pdfs = 72,
      source_detect_pdfs = 0,
      target_detect_pdfs = 0,
      spot_level_score_rds_files = 0,
      spots = 26056,
      retained_spatial_genes = 17775,
      partial_files = 0
    )
  ),
  list(
    step = "05", directory = step05_dir, script = canonical_scripts[[5]],
    expected = c(
      selected_panels = 6,
      selected_interactions = 3,
      healthy_f0_panels = 3,
      f4_panels = 3,
      source_pdfs_verified = 6,
      byte_identical_copies = 6,
      partial_files = 0
    )
  )
)

completion_values <- list()
for (spec in completion_specs) {
  result <- check_completion(
    file.path(spec$directory, paste0(spec$step, "_COMPLETION.tsv")),
    spec$script,
    spec$expected
  )
  completion_values[[spec$step]] <- result$data
  record_check(
    paste0("Step ", spec$step, " completion record and script checksum"),
    result$ok,
    result$detail
  )
}

output_manifest_specs <- list(
  list("01", step01_dir, "01_OUTPUT_SHA256.tsv", 11L, TRUE, "01_"),
  list("02", step02_dir, "02_OUTPUT_SHA256.tsv", 17L, FALSE, "02_"),
  list("03", step02_dir, "03_OUTPUT_SHA256.tsv", 11L, FALSE, "03_"),
  list("04", step04_dir, "04_OUTPUT_SHA256.tsv", 79L, TRUE, ""),
  list("05", step05_dir, "05_OUTPUT_SHA256.tsv", 11L, TRUE, "")
)
for (spec in output_manifest_specs) {
  result <- verify_output_manifest(
    step = spec[[1]],
    directory = spec[[2]],
    manifest_name = spec[[3]],
    expected_rows = spec[[4]],
    exact_directory = spec[[5]],
    required_prefix = spec[[6]]
  )
  record_check(
    paste0("Step ", spec[[1]], " generated-output SHA-256 verification"),
    result$ok,
    result$detail
  )
}

input_manifest_specs <- list(
  list("01", file.path(step01_dir, "01_input_manifest.tsv")),
  list("02", file.path(step02_dir, "02_input_manifest.tsv")),
  list("03", file.path(step02_dir, "03_input_manifest.tsv")),
  list("04", file.path(step04_dir, "04_input_manifest.tsv")),
  list("05", file.path(step05_dir, "05_input_manifest.tsv"))
)
for (spec in input_manifest_specs) {
  say("Verifying Step", spec[[1]], "input files and SHA-256 checksums")
  result <- verify_input_manifest(spec[[1]], spec[[2]])
  record_check(
    paste0("Step ", spec[[1]], " input manifest and SHA-256 verification"),
    result$ok,
    result$detail
  )
}

source_sha_values <- vapply(
  completion_values[c("01", "04")],
  function(x) if (!is.null(x) && "source_object_sha256" %in% names(x)) {
    x[["source_object_sha256"]]
  } else if (!is.null(x) && "input_sha256" %in% names(x)) {
    x[["input_sha256"]]
  } else {
    NA_character_
  },
  character(1)
)
record_check(
  "Shared validated Figure 1 source-object identity",
  length(unique(na.omit(source_sha_values))) == 1L &&
    identical(
      unname(unique(na.omit(source_sha_values))),
      "478cbd318aa832381fd9a2cb1f0ecc7e21f7ce55c62b4a5d0de2d522df517fdd"
    ),
  paste("Step 01/04 source SHA-256:", paste(source_sha_values, collapse = "; "))
)

# -------------------------------------------------------------------------
# Step 01 LIANA inference integrity
# -------------------------------------------------------------------------

cell_types <- c(
  "Bcells", "CentralHep", "Cholangiocyte", "Erythroid", "InfMac",
  "InterHep", "NKTcell", "NonInfMac", "PortalEndo", "PortalHep",
  "Stellate", "cvLSECs"
)
stages <- c("F0", "F1", "F2", "F3", "F4")
expected_stage_counts <- c(F0 = 6L, F1 = 6L, F2 = 3L, F3 = 3L, F4 = 6L)

aggregate_result <- safe_rds(file.path(step01_dir, "01_liana_aggregate.rds"))
significant_result <- safe_rds(file.path(step01_dir, "01_liana_significant_rank_le_0.01.rds"))
liana_ok <- FALSE
liana_detail <- paste(aggregate_result$detail, significant_result$detail, sep = "; ")
if (aggregate_result$ok && significant_result$ok) {
  aggregate <- aggregate_result$data
  significant <- significant_result$data
  required_columns <- c(
    "source", "target", "ligand.complex", "receptor.complex",
    "aggregate_rank", "mean_rank", "natmi.edge_specificity", "natmi.rank",
    "sca.LRscore", "sca.rank", "logfc.logfc_comb", "logfc.rank",
    "connectome.weight_sc", "connectome.rank"
  )
  missing <- union(
    setdiff(required_columns, colnames(aggregate)),
    setdiff(required_columns, colnames(significant))
  )
  if (!length(missing)) {
    key_columns <- c("source", "target", "ligand.complex", "receptor.complex")
    aggregate_keys <- key_from(aggregate, key_columns)
    significant_keys <- key_from(significant, key_columns)
    expected_significant <- aggregate[aggregate$aggregate_rank <= 0.01, , drop = FALSE]
    liana_ok <- nrow(aggregate) == 609984L && ncol(aggregate) == 14L &&
      nrow(significant) == 29205L && ncol(significant) == 14L &&
      !anyDuplicated(aggregate_keys) && !anyDuplicated(significant_keys) &&
      setequal(unique(aggregate$source), cell_types) &&
      setequal(unique(aggregate$target), cell_types) &&
      all(is.finite(aggregate$aggregate_rank)) &&
      all(significant$aggregate_rank <= 0.01) &&
      isTRUE(all.equal(
        as.data.frame(significant),
        as.data.frame(expected_significant),
        tolerance = 0,
        check.attributes = FALSE
      ))
    liana_detail <- sprintf(
      "%d aggregate and %d rank<=0.01 interactions; exact significant subset=%s",
      nrow(aggregate), nrow(significant),
      isTRUE(all.equal(
        as.data.frame(significant), as.data.frame(expected_significant),
        tolerance = 0, check.attributes = FALSE
      ))
    )
  } else {
    liana_detail <- paste("missing LIANA columns:", paste(missing, collapse = ", "))
  }
}
record_check("LIANA aggregate and significant-interaction structure", liana_ok, liana_detail)

pseudobulk_result <- safe_table(
  file.path(step01_dir, "01_pseudobulk_metadata.csv"),
  read_csv
)
pseudobulk_ok <- FALSE
pseudobulk_detail <- pseudobulk_result$detail
if (pseudobulk_result$ok) {
  pb <- pseudobulk_result$data
  required <- c("pb_id", "celltype", "donor", "condition", "condition_stage")
  missing <- setdiff(required, colnames(pb))
  if (!length(missing)) {
    donor_stage <- unique(pb[c("donor", "condition_stage")])
    stage_counts <- table(donor_stage$condition_stage)
    expected_keys <- as.vector(outer(cell_types, donor_stage$donor, paste, sep = "||"))
    pseudobulk_ok <- nrow(pb) == 288L && !anyDuplicated(pb$pb_id) &&
      setequal(pb$pb_id, expected_keys) &&
      setequal(pb$celltype, cell_types) && nrow(donor_stage) == 24L &&
      !anyDuplicated(donor_stage$donor) &&
      setequal(names(stage_counts), stages) &&
      identical(as.integer(stage_counts[stages]), as.integer(expected_stage_counts)) &&
      identical(as.character(pb$condition), as.character(pb$condition_stage))
    pseudobulk_detail <- sprintf(
      "%d cell-type/donor pseudobulks; 12 cell types; 24 donors; stages=%s",
      nrow(pb), paste(as.integer(stage_counts[stages]), collapse = "/")
    )
  } else {
    pseudobulk_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("LIANA pseudobulk metadata structure", pseudobulk_ok, pseudobulk_detail)

step01_parameters_result <- safe_table(
  file.path(step01_dir, "01_analysis_parameters.tsv"), read_tsv
)
step01_parameters_ok <- FALSE
step01_parameters_detail <- step01_parameters_result$detail
if (step01_parameters_result$ok) {
  p <- parameter_values(step01_parameters_result$data)
  required <- c(
    "species", "resource", "methods", "liana_input_class",
    "pseudobulk_normalization_method", "pseudobulk_scale_factor",
    "significance_rank_threshold", "stage_levels"
  )
  methods_text <- if ("methods" %in% names(p)) p[["methods"]] else ""
  step01_parameters_ok <- all(required %in% names(p)) &&
    identical(tolower(p[["species"]]), "human") &&
    identical(tolower(p[["resource"]]), "omnipath") &&
    setequal(strsplit(methods_text, ";", fixed = TRUE)[[1]], c("cellchat", "natmi", "sca", "logfc", "connectome")) &&
    identical(p[["liana_input_class"]], "SingleCellExperiment") &&
    identical(p[["pseudobulk_normalization_method"]], "LogNormalize") &&
    numeric_equal(p[["pseudobulk_scale_factor"]], 1000000) &&
    numeric_equal(p[["significance_rank_threshold"]], 0.01) &&
    identical(p[["stage_levels"]], "F0;F1;F2;F3;F4")
  step01_parameters_detail <- sprintf(
    "%d recorded parameters; required LIANA specification present=%s",
    length(p), all(required %in% names(p))
  )
}
record_check("LIANA recorded analysis specification", step01_parameters_ok, step01_parameters_detail)

# -------------------------------------------------------------------------
# Step 02 stage-activity integrity
# -------------------------------------------------------------------------

stage_result <- safe_rds(file.path(step02_dir, "02_stage_activity.rds"))
selected_result <- safe_rds(file.path(step02_dir, "02_selected_interactions_top60.rds"))
row_z_result <- safe_table(
  file.path(step02_dir, "02_selected_interactions_row_z_scores.csv"), read_csv
)
stage_ok <- FALSE
stage_detail <- paste(stage_result$detail, selected_result$detail, row_z_result$detail, sep = "; ")
if (stage_result$ok && selected_result$ok && row_z_result$ok) {
  stage_activity <- stage_result$data
  selected <- selected_result$data
  row_z <- row_z_result$data
  key_columns <- c("source", "target", "ligand.complex", "receptor.complex")
  required_stage <- c(
    key_columns, "lig_gene", "rec_gene", stages,
    paste0("delta_", stages[-1], "_vs_F0"),
    "spearman_r", "spearman_p", "spearman_q"
  )
  required_selected <- c(required_stage, "rowlab")
  required_z <- c("source", "target", "ligand", "receptor", "rowlab", stages)
  missing <- union(
    union(setdiff(required_stage, colnames(stage_activity)), setdiff(required_selected, colnames(selected))),
    setdiff(required_z, colnames(row_z))
  )
  if (!length(missing)) {
    stage_keys <- key_from(stage_activity, key_columns)
    selected_keys <- key_from(selected, key_columns)
    z_keys <- key_from(
      transform(row_z, ligand.complex = ligand, receptor.complex = receptor),
      key_columns
    )
    selected_values <- as.matrix(selected[stages])
    recomputed_z <- t(apply(selected_values, 1, function(x) as.numeric(scale(x))))
    z_match <- match(z_keys, selected_keys)
    z_ok <- !anyNA(z_match) && numeric_equal(
      as.numeric(as.matrix(row_z[stages])),
      as.numeric(recomputed_z[z_match, , drop = FALSE]),
      tolerance = 1e-10
    )
    stage_ok <- nrow(stage_activity) == 3000L && ncol(stage_activity) == 18L &&
      nrow(selected) == 60L && ncol(selected) == 19L &&
      nrow(row_z) == 60L && ncol(row_z) == 10L &&
      !anyDuplicated(stage_keys) && !anyDuplicated(selected_keys) &&
      !anyDuplicated(z_keys) && all(selected_keys %in% stage_keys) &&
      setequal(z_keys, selected_keys) && z_ok &&
      sum(is.finite(stage_activity$spearman_p)) == 2903L &&
      sum(stage_activity$spearman_q <= 0.05, na.rm = TRUE) == 349L &&
      sum(
        stage_activity$spearman_q <= 0.05 &
          abs(stage_activity$delta_F4_vs_F0) >= 0.5,
        na.rm = TRUE
      ) == 315L
    stage_detail <- sprintf(
      "3000 tests; 2903 finite; 349 BH-significant; 315 effect-qualified; 60 selected; row-z consistency=%s",
      z_ok
    )
  } else {
    stage_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("Stage-activity tests, selected interactions and row-z scores", stage_ok, stage_detail)

summary_result <- safe_table(file.path(step02_dir, "02_statistical_summary.tsv"), read_tsv)
summary_ok <- FALSE
summary_detail <- summary_result$detail
if (summary_result$ok && identical(colnames(summary_result$data), c("statistic", "value"))) {
  summary_values <- stats::setNames(
    as.character(summary_result$data$value),
    as.character(summary_result$data$statistic)
  )
  expected <- c(
    top_ranked_candidates_requested = 3000,
    single_gene_candidates_tested = 3000,
    finite_spearman_tests = 2903,
    bh_significant_trends_q_le_0.05 = 349,
    significant_and_abs_f4_f0_delta_ge_0.5 = 315,
    selected_interactions = 60,
    chord_significant_interactions = 29205,
    cholangiocyte_participating_interactions = 15693
  )
  numeric_fields <- intersect(names(expected), names(summary_values))
  summary_ok <- nrow(summary_result$data) == 9L &&
    all(names(expected) %in% names(summary_values)) &&
    all(vapply(numeric_fields, function(field) {
      numeric_equal(summary_values[[field]], expected[[field]])
    }, logical(1))) &&
    identical(
      summary_values[["dotplot_interactions_per_stage"]],
      "F1:60;F2:60;F3:60;F4:60"
    )
  summary_detail <- sprintf("%d recorded statistical summaries; expected values=%s", nrow(summary_result$data), summary_ok)
}
record_check("Recorded statistical-test counts", summary_ok, summary_detail)

chord_result <- safe_table(file.path(step02_dir, "02_chord_significant_lr_counts.csv"), read_csv)
partner_result <- safe_table(file.path(step02_dir, "02_cholangiocyte_partner_counts.csv"), read_csv)
communication_ok <- FALSE
communication_detail <- paste(chord_result$detail, partner_result$detail, sep = "; ")
if (chord_result$ok && partner_result$ok) {
  chord <- chord_result$data
  partner <- partner_result$data
  chord_keys <- if (all(c("source", "target") %in% colnames(chord))) key_from(chord, c("source", "target")) else character(0)
  expected_chord_keys <- as.vector(outer(cell_types, cell_types, paste, sep = "||"))
  communication_ok <- identical(colnames(chord), c("source", "target", "n")) &&
    nrow(chord) == 144L && !anyDuplicated(chord_keys) &&
    setequal(chord_keys, expected_chord_keys) &&
    sum(chord$n) == 29205L &&
    identical(colnames(partner), c("direction", "partner", "n_lr_pairs")) &&
    nrow(partner) == 23L &&
    !anyDuplicated(paste(partner$direction, partner$partner, sep = "||")) &&
    sum(partner$n_lr_pairs) == 15693L
  communication_detail <- sprintf(
    "144 chord source-target counts sum to %d; 23 cholangiocyte partner counts sum to %d",
    sum(chord$n), sum(partner$n_lr_pairs)
  )
}
record_check("Chord and cholangiocyte-hub count consistency", communication_ok, communication_detail)

dotplot_ok <- TRUE
dotplot_details <- character(0)
if (selected_result$ok) {
  selected <- selected_result$data
  selected_keys <- key_from(selected, c("source", "target", "ligand.complex", "receptor.complex"))
  for (stage in stages[-1]) {
    path <- file.path(step02_dir, paste0("02_dotplot_input_", stage, "_vs_F0.csv"))
    result <- safe_table(path, read_csv)
    current_ok <- result$ok
    if (result$ok) {
      x <- result$data
      activity_column <- stage
      delta_column <- paste0("delta_", stage, "_vs_F0")
      required <- c(
        "source", "target", "ligand.complex", "receptor.complex",
        activity_column, delta_column, "comparison", "delta_pos"
      )
      current_ok <- all(required %in% colnames(x)) && nrow(x) == 60L
      if (current_ok) {
        keys <- key_from(x, c("source", "target", "ligand.complex", "receptor.complex"))
        idx <- match(keys, selected_keys)
        current_ok <- !anyDuplicated(keys) && !anyNA(idx) &&
          setequal(keys, selected_keys) &&
          numeric_equal(x[[activity_column]], selected[[activity_column]][idx]) &&
          numeric_equal(x[[delta_column]], selected[[delta_column]][idx]) &&
          numeric_equal(x$delta_pos, pmax(x[[delta_column]], 0)) &&
          all(x$comparison == paste(stage, "vs F0"))
      }
    }
    dotplot_ok <- dotplot_ok && current_ok
    dotplot_details <- c(dotplot_details, paste0(stage, "=", current_ok))
  }
} else {
  dotplot_ok <- FALSE
}
record_check(
  "Supplementary Figure 8 stage-versus-F0 dot-plot inputs",
  dotplot_ok,
  paste(dotplot_details, collapse = "; ")
)

# -------------------------------------------------------------------------
# Step 03 plot integrity
# -------------------------------------------------------------------------

plot_manifest_result <- safe_table(file.path(step02_dir, "03_plot_manifest.csv"), read_csv)
plot_manifest_ok <- FALSE
plot_manifest_detail <- plot_manifest_result$detail
if (plot_manifest_result$ok) {
  plots <- plot_manifest_result$data
  required <- c("figure", "panel", "content", "filename", "width_inches", "height_inches", "pages")
  missing <- setdiff(required, colnames(plots))
  if (!length(missing)) {
    expected_files <- c(
      "03_Figure3A_chord_LIANA_aggregate_rank_le_0.01.pdf",
      "03_Figure3B_stage_activity_heatmap_with_table.pdf",
      "03_SuppFigure8A_cholangiocyte_partner_counts.pdf",
      "03_SuppFigure8B_dotplot_F1_vs_F0.pdf",
      "03_SuppFigure8C_dotplot_F2_vs_F0.pdf",
      "03_SuppFigure8D_dotplot_F3_vs_F0.pdf",
      "03_SuppFigure8E_dotplot_F4_vs_F0.pdf"
    )
    pdf_paths <- file.path(step02_dir, plots$filename)
    plot_manifest_ok <- nrow(plots) == 7L && !anyDuplicated(plots$filename) &&
      identical(as.character(plots$filename), expected_files) &&
      sum(plots$figure == "Figure 3") == 2L &&
      sum(plots$figure == "Supplementary Figure 8") == 5L &&
      all(plots$pages == 1L) && pdf_magic_ok(pdf_paths)
    plot_manifest_detail <- sprintf(
      "%d one-page PDFs; Figure 3=%d; Supplementary Figure 8=%d; PDF signatures=%s",
      nrow(plots), sum(plots$figure == "Figure 3"),
      sum(plots$figure == "Supplementary Figure 8"), pdf_magic_ok(pdf_paths)
    )
  } else {
    plot_manifest_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("Figure 3A-B and Supplementary Figure 8A-E plot manifest", plot_manifest_ok, plot_manifest_detail)

# -------------------------------------------------------------------------
# Step 04 spatial-map integrity
# -------------------------------------------------------------------------

donor_map_result <- safe_table(file.path(step04_dir, "04_donor_stage_map.csv"), read_csv)
story_result <- safe_table(file.path(step04_dir, "04_story_interactions.csv"), read_csv)
limits_result <- safe_table(file.path(step04_dir, "04_spatial_lr_color_limits.csv"), read_csv)
spatial_manifest_result <- safe_table(file.path(step04_dir, "04_spatial_lr_manifest.csv"), read_csv)

expected_donor_map <- data.frame(
  stage = rep(stages, times = c(6, 6, 3, 3, 6)),
  donor = c(
    "6B", "5D", "4C", "3B", "2A", "1A",
    "6A", "5A", "4D", "3C", "2B", "1D",
    "6C", "5B", "1C",
    "4A", "3D", "2C",
    "6D", "5C", "4B", "3A", "2D", "1B"
  ),
  stringsAsFactors = FALSE
)
donor_map_ok <- FALSE
donor_map_detail <- donor_map_result$detail
if (donor_map_result$ok) {
  donor_map <- donor_map_result$data
  required <- c("stage_order", "stage", "stage_folder", "donor_order", "donor")
  missing <- setdiff(required, colnames(donor_map))
  if (!length(missing)) {
    donor_map_ok <- nrow(donor_map) == 24L && !anyDuplicated(donor_map$donor) &&
      identical(as.character(donor_map$stage), expected_donor_map$stage) &&
      identical(as.character(donor_map$donor), expected_donor_map$donor) &&
      identical(as.integer(table(donor_map$stage)[stages]), as.integer(expected_stage_counts)) &&
      identical(as.character(donor_map$stage_folder), ifelse(donor_map$stage == "F0", "Healthy_F0", donor_map$stage))
    donor_map_detail <- paste(
      "24 fixed donors; stage replication F0/F1/F2/F3/F4 =",
      paste(as.integer(table(donor_map$stage)[stages]), collapse = "/")
    )
  } else {
    donor_map_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("Fixed donor/stage map and ordering", donor_map_ok, donor_map_detail)

expected_story <- data.frame(
  source = c("Cholangiocyte", "Cholangiocyte", "Stellate"),
  target = c("Stellate", "InfMac", "PortalEndo"),
  ligand = c("NPNT", "MMP7", "CCL2"),
  receptor = c("ITGA8", "CD151", "ACKR1"),
  pair = c("NPNT_ITGA8", "MMP7_CD151", "CCL2_ACKR1"),
  direction = c(
    "Cholangiocyte_to_Stellate",
    "Cholangiocyte_to_InfMac",
    "Stellate_to_PortalEndo"
  ),
  stringsAsFactors = FALSE
)
story_ok <- FALSE
story_detail <- story_result$detail
if (story_result$ok) {
  story <- story_result$data
  required <- c("edge_order", names(expected_story), "note", "score_column")
  missing <- setdiff(required, colnames(story))
  if (!length(missing)) {
    story_ok <- nrow(story) == 3L && identical(story$edge_order, 1:3) &&
      all(vapply(names(expected_story), function(column) {
        identical(as.character(story[[column]]), expected_story[[column]])
      }, logical(1))) &&
      !anyDuplicated(story$score_column)
    story_detail <- paste("three fixed interactions:", paste(story$pair, collapse = "; "))
  } else {
    story_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("Fixed spatial LR interactions and source-target identities", story_ok, story_detail)

expected_q99 <- c(
  NPNT_ITGA8 = 1.56069258152837e-09,
  MMP7_CD151 = 3.02077740741239,
  CCL2_ACKR1 = 0.000991789688945745
)
limits_ok <- FALSE
limits_detail <- limits_result$detail
if (limits_result$ok) {
  limits <- limits_result$data
  required <- c(
    "edge_order", "source", "target", "ligand", "receptor", "pair",
    "score_column", "q01", "q99", "finite_spots"
  )
  missing <- setdiff(required, colnames(limits))
  if (!length(missing)) {
    idx <- match(names(expected_q99), limits$pair)
    limits_ok <- nrow(limits) == 3L && !anyDuplicated(limits$pair) &&
      !anyNA(idx) && all(limits$q01 == 0) &&
      numeric_equal(limits$q99[idx], expected_q99, tolerance = 1e-12) &&
      all(limits$finite_spots == 26056L)
    limits_detail <- paste(
      "pair-specific q01=0 and q99:",
      paste(sprintf("%s=%.15g", limits$pair, limits$q99), collapse = "; ")
    )
  } else {
    limits_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("Pair-specific shared spatial color limits", limits_ok, limits_detail)

spatial_ok <- FALSE
spatial_detail <- spatial_manifest_result$detail
if (spatial_manifest_result$ok && donor_map_result$ok && story_result$ok && limits_result$ok) {
  spatial <- spatial_manifest_result$data
  required <- c(
    "stage_order", "stage", "stage_folder", "donor_order", "donor",
    "edge_order", "source", "target", "ligand", "receptor", "pair",
    "direction", "relative_pdf_path", "spots", "image_key",
    "relative_aspect_ratio_difference", "spot_multiplier", "circle_vertices",
    "true_visium_spot_diameter_pixels", "true_visium_spot_radius_pixels",
    "maximum_center_overrun_pixels", "color_q01", "color_q99",
    "pdf_bytes", "pdf_sha256"
  )
  missing <- setdiff(required, colnames(spatial))
  if (!length(missing)) {
    spatial_keys <- paste(spatial$donor, spatial$pair, sep = "||")
    expected_keys <- as.vector(outer(donor_map$donor, story$pair, paste, sep = "||"))
    stage_rows <- table(spatial$stage)
    donor_spots <- spatial[!duplicated(spatial$donor), c("donor", "spots")]
    limit_idx <- match(spatial$pair, limits$pair)
    pdf_paths <- file.path(step04_dir, spatial$relative_pdf_path)
    geometry_ok <- all(spatial$relative_aspect_ratio_difference <= 0.002) &&
      all(spatial$maximum_center_overrun_pixels <= 1e-8) &&
      all(spatial$spot_multiplier == 1) &&
      all(spatial$circle_vertices == 48L) &&
      all(spatial$true_visium_spot_diameter_pixels > 15) &&
      all(spatial$true_visium_spot_diameter_pixels < 17) &&
      numeric_equal(
        spatial$true_visium_spot_radius_pixels,
        spatial$true_visium_spot_diameter_pixels / 2,
        tolerance = 1e-12
      )
    colors_ok <- !anyNA(limit_idx) &&
      numeric_equal(spatial$color_q01, limits$q01[limit_idx]) &&
      numeric_equal(spatial$color_q99, limits$q99[limit_idx], tolerance = 1e-12)
    pdfs_ok <- all(safe_relative_paths(spatial$relative_pdf_path)) &&
      all(nonempty_files(pdf_paths)) &&
      identical(as.numeric(file.info(pdf_paths)$size), as.numeric(spatial$pdf_bytes)) &&
      identical(
        tolower(unname(vapply(pdf_paths, sha256_file, character(1)))),
        tolower(as.character(spatial$pdf_sha256))
      ) && pdf_magic_ok(pdf_paths)
    spatial_ok <- nrow(spatial) == 72L && !anyDuplicated(spatial_keys) &&
      setequal(spatial_keys, expected_keys) &&
      identical(as.integer(stage_rows[stages]), 3L * as.integer(expected_stage_counts)) &&
      nrow(donor_spots) == 24L && sum(donor_spots$spots) == 26056L &&
      all(table(spatial$donor) == 3L) && all(table(spatial$pair) == 24L) &&
      geometry_ok && colors_ok && pdfs_ok
    spatial_detail <- sprintf(
      "72 PDFs; 24 donors x 3 pairs; 26056 spots; max aspect difference=%.6f%%; diameter=%.3f-%.3f px; colors=%s; checksums=%s",
      100 * max(spatial$relative_aspect_ratio_difference),
      min(spatial$true_visium_spot_diameter_pixels),
      max(spatial$true_visium_spot_diameter_pixels),
      colors_ok, pdfs_ok
    )
  } else {
    spatial_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("All-donor spatial LR manifest, geometry and PDF integrity", spatial_ok, spatial_detail)

step04_parameters_result <- safe_table(file.path(step04_dir, "04_analysis_parameters.tsv"), read_tsv)
step04_parameters_ok <- FALSE
step04_parameters_detail <- step04_parameters_result$detail
if (step04_parameters_result$ok) {
  p <- parameter_values(step04_parameters_result$data)
  required <- c(
    "input_assay", "input_layer", "spatial_normalization_method",
    "spatial_normalization_scale_factor", "nearest_neighbors_excluding_self",
    "include_self", "effective_neighbors_including_self", "gaussian_sigma",
    "score_transform", "fixed_color_limits", "spot_multiplier",
    "circle_vertices", "source_target_detection_outputs",
    "large_spot_level_score_object_saved"
  )
  step04_parameters_ok <- all(required %in% names(p)) &&
    identical(p[["input_assay"]], "Spatial") &&
    identical(p[["input_layer"]], "counts") &&
    identical(p[["spatial_normalization_method"]], "LogNormalize") &&
    numeric_equal(p[["spatial_normalization_scale_factor"]], 10000) &&
    numeric_equal(p[["nearest_neighbors_excluding_self"]], 30) &&
    tolower(p[["include_self"]]) == "true" &&
    numeric_equal(p[["effective_neighbors_including_self"]], 31) &&
    numeric_equal(p[["gaussian_sigma"]], 2) &&
    identical(p[["score_transform"]], "log1p(local_lr_score * 10000)") &&
    grepl("1st and 99th", p[["fixed_color_limits"]], fixed = TRUE) &&
    numeric_equal(p[["spot_multiplier"]], 1) &&
    numeric_equal(p[["circle_vertices"]], 48) &&
    identical(p[["source_target_detection_outputs"]], "not generated") &&
    tolower(p[["large_spot_level_score_object_saved"]]) == "false"
  step04_parameters_detail <- sprintf(
    "%d parameters; normalization/neighborhood/scoring/rendering specification=%s",
    length(p), step04_parameters_ok
  )
}
record_check("Spatial LR recorded analysis and rendering specification", step04_parameters_ok, step04_parameters_detail)

forbidden_pdfs <- if (dir.exists(step04_dir)) {
  list.files(
    step04_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "(sourceDetect|targetDetect|codetect).*\\.pdf$",
    ignore.case = TRUE
  )
} else character(0)
derived_rds <- if (dir.exists(step04_dir)) {
  list.files(step04_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.rds$", ignore.case = TRUE)
} else character(0)
record_check(
  "Absence of excluded spatial detection PDFs and derived RDS files",
  !length(forbidden_pdfs) && !length(derived_rds),
  sprintf("detection PDFs=%d; derived RDS files=%d", length(forbidden_pdfs), length(derived_rds))
)

# -------------------------------------------------------------------------
# Step 05 fixed Figure 3C panel selection integrity
# -------------------------------------------------------------------------

selected_manifest_result <- safe_table(
  file.path(step05_dir, "05_figure3C_selected_panel_manifest.csv"), read_csv
)
selected_limits_result <- safe_table(
  file.path(step05_dir, "05_selected_interactions_and_color_limits.csv"), read_csv
)
selection_parameters_result <- safe_table(
  file.path(step05_dir, "05_selection_parameters.tsv"), read_tsv
)

expected_panels <- data.frame(
  panel_order = 1:6,
  figure_row = c(rep("top", 3), rep("bottom", 3)),
  figure_column = rep(1:3, 2),
  stage = c(rep("F0", 3), rep("F4", 3)),
  donor = c("6B", "2A", "4C", "4B", "3A", "2D"),
  pair = rep(expected_story$pair, 2),
  stringsAsFactors = FALSE
)
expected_ranks <- c(
  NPNT_ITGA8 = 7.59661069840347e-06,
  MMP7_CD151 = 1.1479120956405e-05,
  CCL2_ACKR1 = 0.000101850349325487
)

selected_panels_ok <- FALSE
selected_panels_detail <- selected_manifest_result$detail
if (selected_manifest_result$ok && selected_limits_result$ok) {
  panels <- selected_manifest_result$data
  selected_limits <- selected_limits_result$data
  required_panels <- c(
    names(expected_panels), "source", "target", "ligand", "receptor",
    "direction", "liana_aggregate_rank", "color_q01", "color_q99",
    "source_relative_pdf_path", "selected_relative_pdf_path",
    "pdf_bytes", "pdf_sha256"
  )
  required_limits <- c(
    "edge_order", "source", "target", "ligand", "receptor", "pair",
    "direction", "liana_aggregate_rank", "aggregate_rank_threshold",
    "color_q01", "color_q99", "color_limits_calculated_across_donors",
    "finite_spots"
  )
  missing <- union(
    setdiff(required_panels, colnames(panels)),
    setdiff(required_limits, colnames(selected_limits))
  )
  if (!length(missing)) {
    fixed_composition_ok <- all(vapply(names(expected_panels), function(column) {
      identical(as.character(panels[[column]]), as.character(expected_panels[[column]]))
    }, logical(1)))
    source_paths <- file.path(repro_root, panels$source_relative_pdf_path)
    selected_paths <- file.path(repro_root, panels$selected_relative_pdf_path)
    files_ok <- all(safe_relative_paths(panels$source_relative_pdf_path)) &&
      all(safe_relative_paths(panels$selected_relative_pdf_path)) &&
      all(nonempty_files(c(source_paths, selected_paths)))
    source_sha <- if (files_ok) vapply(source_paths, sha256_file, character(1)) else rep(NA_character_, nrow(panels))
    selected_sha <- if (files_ok) vapply(selected_paths, sha256_file, character(1)) else rep(NA_character_, nrow(panels))
    byte_identical_ok <- files_ok &&
      identical(unname(source_sha), unname(selected_sha)) &&
      identical(
        tolower(unname(source_sha)),
        tolower(as.character(panels$pdf_sha256))
      ) &&
      identical(as.numeric(file.info(selected_paths)$size), as.numeric(panels$pdf_bytes))
    rank_idx <- match(panels$pair, names(expected_ranks))
    color_idx <- match(panels$pair, names(expected_q99))
    ranks_ok <- !anyNA(rank_idx) && numeric_equal(
      panels$liana_aggregate_rank,
      unname(expected_ranks[rank_idx]),
      tolerance = 1e-14
    )
    colors_ok <- !anyNA(color_idx) && all(panels$color_q01 == 0) &&
      numeric_equal(panels$color_q99, unname(expected_q99[color_idx]), tolerance = 1e-12)
    limit_idx <- match(selected_limits$pair, names(expected_ranks))
    selected_limits_ok <- nrow(selected_limits) == 3L && !anyNA(limit_idx) &&
      numeric_equal(selected_limits$liana_aggregate_rank, unname(expected_ranks[limit_idx]), tolerance = 1e-14) &&
      all(selected_limits$aggregate_rank_threshold == 0.01) &&
      all(selected_limits$color_limits_calculated_across_donors == 24L) &&
      all(selected_limits$finite_spots == 26056L)
    selected_panels_ok <- nrow(panels) == 6L && fixed_composition_ok &&
      byte_identical_ok && ranks_ok && colors_ok && selected_limits_ok &&
      pdf_magic_ok(selected_paths)
    selected_panels_detail <- sprintf(
      "six fixed F0/F4 panels; composition=%s; byte-identical copies=%s; ranks=%s; shared limits=%s",
      fixed_composition_ok, byte_identical_ok, ranks_ok, colors_ok
    )
  } else {
    selected_panels_detail <- paste("missing columns:", paste(missing, collapse = ", "))
  }
}
record_check("Fixed Figure 3C panel selection, ranks and byte-identical copies", selected_panels_ok, selected_panels_detail)

selection_record_ok <- FALSE
selection_record_detail <- selection_parameters_result$detail
if (selection_parameters_result$ok) {
  p <- parameter_values(selection_parameters_result$data)
  fixed_basis <- "selection_basis" %in% names(p) &&
    grepl("Fixed to reproduce the accepted manuscript Figure 3C", p[["selection_basis"]], fixed = TRUE)
  timing_parameter <- "historical_panel_selection_preceded_current_regenerated_output_visual_review"
  before_visual <- timing_parameter %in% names(p) &&
    identical(toupper(p[[timing_parameter]]), "TRUE")
  selection_record_ok <- length(p) == 14L && fixed_basis && before_visual
  selection_record_detail <- sprintf(
    "%d selection parameters; fixed manuscript composition=%s; selection before regenerated visual inspection=%s",
    length(p), fixed_basis, before_visual
  )
}
record_check("Pre-specified Figure 3C panel-selection record", selection_record_ok, selection_record_detail)

# -------------------------------------------------------------------------
# Session information and incomplete-output detection
# -------------------------------------------------------------------------

session_files <- c(
  file.path(step01_dir, "01_sessionInfo.txt"),
  file.path(step02_dir, "02_sessionInfo.txt"),
  file.path(step02_dir, "03_sessionInfo.txt"),
  file.path(step04_dir, "04_sessionInfo.txt"),
  file.path(step05_dir, "05_sessionInfo.txt")
)
session_ok <- all(nonempty_files(session_files)) && all(vapply(
  session_files,
  function(path) any(grepl("^R version ", readLines(path, warn = FALSE))),
  logical(1)
))
record_check(
  "Software session information for Steps 01-05",
  session_ok,
  sprintf("%d/%d non-empty sessionInfo files with R version", sum(nonempty_files(session_files)), length(session_files))
)

partial_roots <- c(results_root, scripts_dir, audit_parent)
partial_paths <- character(0)
for (root in partial_roots[dir.exists(partial_roots)]) {
  candidates <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  partial_paths <- c(
    partial_paths,
    candidates[grepl("\\.partial($|\\.)", basename(candidates))]
  )
}
partial_paths <- sort(unique(partial_paths))
record_check(
  "Absence of incomplete .partial outputs",
  length(partial_paths) == 0L,
  if (!length(partial_paths)) "no .partial files or directories detected" else paste(partial_paths, collapse = "; ")
)

# -------------------------------------------------------------------------
# Final decision, compact public records and audit outputs
# -------------------------------------------------------------------------

component_status <- do.call(rbind, status_rows)
rownames(component_status) <- NULL
audit_pass <- all(component_status$pass[component_status$required])

output_count_summary <- data.frame(
  item = c(
    "LIANA aggregate interactions",
    "LIANA significant interactions",
    "donor-by-cell-type pseudobulks",
    "stage-activity candidates tested",
    "finite Spearman tests",
    "BH-significant trends",
    "effect-size-qualified trends",
    "selected stage-progressive interactions",
    "Figure 3A-B and Supplementary Figure 8 PDFs",
    "all-donor LRlocal PDFs",
    "selected Figure 3C PDFs",
    "forbidden detection PDFs",
    "derived spatial-score RDS files",
    "incomplete .partial outputs"
  ),
  expected = c(609984, 29205, 288, 3000, 2903, 349, 315, 60, 7, 72, 6, 0, 0, 0),
  observed = c(
    if (aggregate_result$ok) nrow(aggregate_result$data) else NA,
    if (significant_result$ok) nrow(significant_result$data) else NA,
    if (pseudobulk_result$ok) nrow(pseudobulk_result$data) else NA,
    if (stage_result$ok) nrow(stage_result$data) else NA,
    if (stage_result$ok) sum(is.finite(stage_result$data$spearman_p)) else NA,
    if (stage_result$ok) sum(stage_result$data$spearman_q <= 0.05, na.rm = TRUE) else NA,
    if (stage_result$ok) sum(
      stage_result$data$spearman_q <= 0.05 &
        abs(stage_result$data$delta_F4_vs_F0) >= 0.5,
      na.rm = TRUE
    ) else NA,
    if (selected_result$ok) nrow(selected_result$data) else NA,
    if (plot_manifest_result$ok) nrow(plot_manifest_result$data) else NA,
    if (spatial_manifest_result$ok) nrow(spatial_manifest_result$data) else NA,
    if (selected_manifest_result$ok) nrow(selected_manifest_result$data) else NA,
    length(forbidden_pdfs), length(derived_rds), length(partial_paths)
  ),
  stringsAsFactors = FALSE
)
output_count_summary$pass <- output_count_summary$expected == output_count_summary$observed

atomic_write_table(
  component_status,
  file.path(audit_dir, "Fig3_component_status.tsv")
)
atomic_write_table(
  output_count_summary,
  file.path(audit_dir, "Fig3_expected_output_counts.tsv")
)
if (nrow(public_sha_rows)) {
  atomic_write_table(
    public_sha_rows,
    file.path(audit_dir, "Fig3_public_code_and_documentation_SHA256.tsv")
  )
}
if (length(verified_output_rows)) {
  atomic_write_table(
    do.call(rbind, verified_output_rows),
    file.path(audit_dir, "Fig3_verified_generated_output_SHA256.tsv")
  )
}
if (length(verified_input_rows)) {
  atomic_write_table(
    do.call(rbind, verified_input_rows),
    file.path(audit_dir, "Fig3_verified_input_SHA256.tsv")
  )
}

compact_records <- list(
  step01 = c(
    "01_COMPLETION.tsv", "01_OUTPUT_SHA256.tsv", "01_analysis_parameters.tsv",
    "01_donor_stage_map.csv", "01_input_manifest.tsv",
    "01_pseudobulk_metadata.csv", "01_sessionInfo.txt"
  ),
  step02 = c(
    "02_COMPLETION.tsv", "02_OUTPUT_SHA256.tsv", "02_analysis_parameters.tsv",
    "02_cholangiocyte_partner_counts.csv", "02_chord_significant_lr_counts.csv",
    paste0("02_dotplot_input_", stages[-1], "_vs_F0.csv"),
    "02_input_manifest.tsv", "02_selected_interactions_row_z_scores.csv",
    "02_selected_interactions_top60.csv", "02_sessionInfo.txt",
    "02_statistical_summary.tsv"
  ),
  step03 = c(
    "03_COMPLETION.tsv", "03_OUTPUT_SHA256.tsv", "03_input_manifest.tsv",
    "03_plot_manifest.csv", "03_plot_parameters.tsv", "03_sessionInfo.txt"
  ),
  step04 = c(
    "04_COMPLETION.tsv", "04_OUTPUT_SHA256.tsv", "04_analysis_parameters.tsv",
    "04_donor_stage_map.csv", "04_input_manifest.tsv", "04_sessionInfo.txt",
    "04_spatial_lr_color_limits.csv", "04_spatial_lr_manifest.csv",
    "04_story_interactions.csv"
  ),
  step05 = c(
    "05_COMPLETION.tsv", "05_OUTPUT_SHA256.tsv",
    "05_figure3C_selected_panel_manifest.csv", "05_input_manifest.tsv",
    "05_selected_interactions_and_color_limits.csv",
    "05_selection_parameters.tsv", "05_sessionInfo.txt"
  )
)
compact_sources <- list(
  step01 = step01_dir,
  step02 = step02_dir,
  step03 = step02_dir,
  step04 = step04_dir,
  step05 = step05_dir
)
compact_manifest_rows <- list()
if (audit_pass) {
  say("Copying verified compact execution records into the public audit directory")
  for (step in names(compact_records)) {
    for (filename in compact_records[[step]]) {
      source <- file.path(compact_sources[[step]], filename)
      target <- file.path(audit_dir, "records", step, filename)
      if (!isTRUE(nonempty_files(source))) stopf("Compact record is missing or empty: %s", source)
      atomic_copy(source, target)
      compact_manifest_rows[[length(compact_manifest_rows) + 1L]] <- data.frame(
        step = step,
        source_relative_path = substring(source, nchar(repro_root) + 2L),
        copied_relative_path = substring(target, nchar(audit_dir) + 2L),
        bytes = as.numeric(file.info(target)$size),
        sha256 = sha256_file(target),
        stringsAsFactors = FALSE
      )
    }
  }
  atomic_write_table(
    do.call(rbind, compact_manifest_rows),
    file.path(audit_dir, "Fig3_compact_record_manifest.tsv")
  )
}

summary_lines <- c(
  "FLINC Figure 3 and Supplementary Figure 8 final reproducibility audit",
  paste0("Timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Workflow root: ", repro_root),
  "",
  paste0(
    "Required checks passed: ",
    sum(component_status$pass & component_status$required),
    " / ",
    sum(component_status$required)
  ),
  "",
  sprintf(
    "%s: %s - %s",
    component_status$component,
    ifelse(component_status$pass, "PASS", "FAIL"),
    component_status$detail
  ),
  "",
  paste0("FINAL_FIGURE3_AUDIT_PASS = ", audit_pass),
  "",
  "Scope note: this audit checks required files, dimensions, statistical-test counts,",
  "internal consistency, completion records, SHA-256 checksums, software information,",
  "fixed Figure 3C panel identities and absence of incomplete or excluded outputs.",
  "It performs no additional biological analysis and contains no historical comparison."
)
atomic_write_lines(summary_lines, file.path(audit_dir, "audit_summary.txt"))
atomic_write_lines(capture.output(sessionInfo()), file.path(audit_dir, "sessionInfo.txt"))

completion_table <- data.frame(
  field = c(
    "status", "finished_at", "script_path", "script_sha256",
    "required_checks", "required_checks_passed", "audit_pass",
    "biological_analysis_performed", "compact_records_copied"
  ),
  value = c(
    if (audit_pass) "complete" else "failed",
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    canonical_scripts[[6]],
    if (file.exists(canonical_scripts[[6]])) sha256_file(canonical_scripts[[6]]) else NA_character_,
    sum(component_status$required),
    sum(component_status$pass & component_status$required),
    audit_pass,
    "none",
    length(compact_manifest_rows)
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(completion_table, file.path(audit_dir, "06_COMPLETION.tsv"))

audit_files_before_checksum <- list.files(
  audit_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)
audit_files_before_checksum <- audit_files_before_checksum[
  file.info(audit_files_before_checksum)$isdir %in% FALSE
]
audit_checksum_table <- data.frame(
  file = substring(audit_files_before_checksum, nchar(audit_dir) + 2L),
  bytes = as.numeric(file.info(audit_files_before_checksum)$size),
  sha256 = vapply(audit_files_before_checksum, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
audit_checksum_table <- audit_checksum_table[order(audit_checksum_table$file), , drop = FALSE]
atomic_write_table(
  audit_checksum_table,
  file.path(audit_dir, "06_OUTPUT_SHA256.tsv")
)

audit_partial <- list.files(
  audit_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  pattern = "\\.partial($|\\.)"
)
if (length(audit_partial)) {
  stopf("Audit created incomplete .partial files: %s", paste(audit_partial, collapse = "; "))
}

cat(paste(summary_lines, collapse = "\n"), "\n")

if (!audit_pass) {
  stopf("Figure 3 audit did not fully pass. Inspect: %s", audit_dir)
}

cat("\nFINAL FIGURE 3 AUDIT PASSED.\n")
