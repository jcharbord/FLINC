#!/usr/bin/env Rscript

# Figure 3C Step 05: select the six fixed manuscript panels from the complete
# 72-panel spatial LR output. This script performs no biological analysis and
# does not recalculate, rescale, or redraw any spatial LR score.

options(width = 180)

started_at <- Sys.time()

repo_root <- "/Users/jianwu/Documents/FLINC"
repro_root <- file.path(repo_root, "Jeremie", "Reproducible_code")

step01_dir <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "01_liana_inference"
)
step04_dir <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "03_spatial_LR_all_samples"
)
output_dir <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "04_main_figure3C_selected"
)
panel_dir <- file.path(output_dir, "panels")

step01_significant_csv <- file.path(
  step01_dir,
  "01_liana_significant_rank_le_0.01.csv"
)
step01_checksums_tsv <- file.path(step01_dir, "01_OUTPUT_SHA256.tsv")
step01_completion_tsv <- file.path(step01_dir, "01_COMPLETION.tsv")

step04_story_csv <- file.path(step04_dir, "04_story_interactions.csv")
step04_color_limits_csv <- file.path(
  step04_dir,
  "04_spatial_lr_color_limits.csv"
)
step04_donor_stage_csv <- file.path(step04_dir, "04_donor_stage_map.csv")
step04_spatial_manifest_csv <- file.path(
  step04_dir,
  "04_spatial_lr_manifest.csv"
)
step04_checksums_tsv <- file.path(step04_dir, "04_OUTPUT_SHA256.tsv")
step04_completion_tsv <- file.path(step04_dir, "04_COMPLETION.tsv")

aggregate_rank_threshold <- 0.01
expected_selected_panels <- 6L
expected_selected_interactions <- 3L

interaction_spec <- data.frame(
  edge_order = 1:3,
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
  expected_aggregate_rank = c(
    7.59661069840347e-06,
    1.1479120956405e-05,
    0.000101850349325487
  ),
  stringsAsFactors = FALSE
)

panel_spec <- data.frame(
  panel_order = 1:6,
  figure_row = c("top", "top", "top", "bottom", "bottom", "bottom"),
  figure_column = c(1L, 2L, 3L, 1L, 2L, 3L),
  stage = c("F0", "F0", "F0", "F4", "F4", "F4"),
  donor = c("6B", "2A", "4C", "4B", "3A", "2D"),
  pair = c(
    "NPNT_ITGA8", "MMP7_CD151", "CCL2_ACKR1",
    "NPNT_ITGA8", "MMP7_CD151", "CCL2_ACKR1"
  ),
  selected_pdf_filename = c(
    "05_Figure3C_01_top_F0_6B_NPNT_ITGA8.pdf",
    "05_Figure3C_02_top_F0_2A_MMP7_CD151.pdf",
    "05_Figure3C_03_top_F0_4C_CCL2_ACKR1.pdf",
    "05_Figure3C_04_bottom_F4_4B_NPNT_ITGA8.pdf",
    "05_Figure3C_05_bottom_F4_3A_MMP7_CD151.pdf",
    "05_Figure3C_06_bottom_F4_2D_CCL2_ACKR1.pdf"
  ),
  stringsAsFactors = FALSE
)

required_packages <- "digest"
if (!requireNamespace(required_packages, quietly = TRUE)) {
  stop("Missing required package: digest")
}

required_inputs <- c(
  step01_significant_csv,
  step01_checksums_tsv,
  step01_completion_tsv,
  step04_story_csv,
  step04_color_limits_csv,
  step04_donor_stage_csv,
  step04_spatial_manifest_csv,
  step04_checksums_tsv,
  step04_completion_tsv
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing Step 05 inputs: ", paste(missing_inputs, collapse = ", "))
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

read_completion <- function(path, step_label) {
  completion <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  status <- completion$value[completion$field == "status"]
  if (length(status) != 1L || status != "complete") {
    stop(step_label, " completion status is not 'complete'.")
  }
  invisible(completion)
}

verify_manifested_file <- function(path, checksum_table, step_label) {
  expected <- checksum_table$sha256[
    checksum_table$file == basename(path)
  ]
  if (length(expected) != 1L) {
    stop(
      "Expected one ", step_label, " checksum entry for: ", basename(path)
    )
  }
  observed <- sha256_file(path)
  if (!identical(tolower(observed), tolower(trimws(expected)))) {
    stop(step_label, " checksum mismatch for: ", basename(path))
  }
  invisible(observed)
}

atomic_write_csv <- function(object, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  utils::write.csv(object, partial, row.names = FALSE)
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

atomic_copy <- function(source, destination) {
  partial <- paste0(destination, ".partial")
  if (file.exists(destination) || file.exists(partial)) {
    stop("Refusing to replace existing selected panel: ", destination)
  }
  copied <- file.copy(
    source,
    partial,
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!isTRUE(copied) || !file.exists(partial)) {
    stop("Could not copy selected panel: ", source)
  }
  if (!file.rename(partial, destination)) {
    stop("Could not finalize selected panel: ", destination)
  }
}

relative_to_repro <- function(path) {
  normalized_path <- normalizePath(path, mustWork = TRUE)
  normalized_root <- normalizePath(repro_root, mustWork = TRUE)
  prefix <- paste0(normalized_root, .Platform$file.sep)
  if (!startsWith(normalized_path, prefix)) {
    stop("Path is not located under the reproducibility root: ", path)
  }
  substring(normalized_path, nchar(prefix) + 1L)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
} else {
  NA_character_
}
script_sha256 <- if (!is.na(script_path)) sha256_file(script_path) else NA_character_

selected_pdf_paths <- file.path(panel_dir, panel_spec$selected_pdf_filename)
metadata_paths <- c(
  panel_manifest_csv = file.path(
    output_dir,
    "05_figure3C_selected_panel_manifest.csv"
  ),
  selected_interactions_csv = file.path(
    output_dir,
    "05_selected_interactions_and_color_limits.csv"
  ),
  parameters_tsv = file.path(output_dir, "05_selection_parameters.tsv"),
  input_manifest_tsv = file.path(output_dir, "05_input_manifest.tsv"),
  session_info_txt = file.path(output_dir, "05_sessionInfo.txt"),
  output_checksums_tsv = file.path(output_dir, "05_OUTPUT_SHA256.tsv"),
  completion_tsv = file.path(output_dir, "05_COMPLETION.tsv")
)

planned_outputs <- c(selected_pdf_paths, unname(metadata_paths))
blocked_outputs <- planned_outputs[
  file.exists(planned_outputs) | file.exists(paste0(planned_outputs, ".partial"))
]
if (length(blocked_outputs)) {
  stop(
    "Refusing to replace existing Step 05 outputs: ",
    paste(blocked_outputs, collapse = ", ")
  )
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(panel_dir)) {
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
}

message("Figure 3C fixed-panel selection started at ", format(started_at))
message("Source: ", step04_dir)
message("Output: ", output_dir)
message("Verifying Step 01 and Step 04 completion records...")

read_completion(step01_completion_tsv, "Step 01")
read_completion(step04_completion_tsv, "Step 04")

step01_checksums <- utils::read.delim(
  step01_checksums_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
step04_checksums <- utils::read.delim(
  step04_checksums_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

verify_manifested_file(
  step01_significant_csv,
  step01_checksums,
  "Step 01"
)
for (path in c(
  step04_story_csv,
  step04_color_limits_csv,
  step04_donor_stage_csv,
  step04_spatial_manifest_csv
)) {
  verify_manifested_file(path, step04_checksums, "Step 04")
}

message("Resolving the fixed three interactions and six manuscript panels...")

story <- utils::read.csv(
  step04_story_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
color_limits <- utils::read.csv(
  step04_color_limits_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
donor_stage <- utils::read.csv(
  step04_donor_stage_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
spatial_manifest <- utils::read.csv(
  step04_spatial_manifest_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
liana_significant <- utils::read.csv(
  step01_significant_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

story_key <- paste(
  story$source,
  story$target,
  story$ligand,
  story$receptor,
  sep = "||"
)
spec_key <- paste(
  interaction_spec$source,
  interaction_spec$target,
  interaction_spec$ligand,
  interaction_spec$receptor,
  sep = "||"
)
story_hits <- match(spec_key, story_key)
if (anyNA(story_hits) ||
    anyDuplicated(story_key[story_hits]) ||
    sum(story_key %in% spec_key) != expected_selected_interactions) {
  stop("The Step 04 story table does not contain the unique fixed interactions.")
}
resolved_story <- story[story_hits, , drop = FALSE]
if (!identical(as.character(resolved_story$pair), interaction_spec$pair) ||
    !identical(as.character(resolved_story$direction), interaction_spec$direction)) {
  stop("The Step 04 pair or direction labels differ from the fixed specification.")
}

liana_key <- paste(
  liana_significant$source,
  liana_significant$target,
  liana_significant$ligand.complex,
  liana_significant$receptor.complex,
  sep = "||"
)
liana_hits <- match(spec_key, liana_key)
if (anyNA(liana_hits) ||
    anyDuplicated(liana_key[liana_hits]) ||
    sum(liana_key %in% spec_key) != expected_selected_interactions) {
  stop("The LIANA significant table does not contain the unique fixed interactions.")
}
aggregate_ranks <- as.numeric(liana_significant$aggregate_rank[liana_hits])
if (any(!is.finite(aggregate_ranks)) ||
    any(aggregate_ranks > aggregate_rank_threshold)) {
  stop("At least one fixed interaction does not satisfy aggregate_rank <= 0.01.")
}
rank_matches <- mapply(
  function(observed, expected) {
    isTRUE(all.equal(observed, expected, tolerance = 1e-12))
  },
  aggregate_ranks,
  interaction_spec$expected_aggregate_rank
)
if (!all(rank_matches)) {
  stop("At least one fixed interaction has an unexpected LIANA aggregate rank.")
}

color_hits <- match(interaction_spec$pair, color_limits$pair)
if (anyNA(color_hits) || anyDuplicated(color_limits$pair[color_hits])) {
  stop("The Step 04 color-limit table does not uniquely contain all three pairs.")
}
resolved_limits <- color_limits[color_hits, , drop = FALSE]
if (any(!is.finite(resolved_limits$q01)) ||
    any(!is.finite(resolved_limits$q99)) ||
    any(resolved_limits$q99 < resolved_limits$q01) ||
    any(resolved_limits$finite_spots != 26056L)) {
  stop("The pair-specific color limits are incomplete or invalid.")
}

donor_stage_key <- paste(donor_stage$donor, donor_stage$stage, sep = "||")
panel_stage_key <- paste(panel_spec$donor, panel_spec$stage, sep = "||")
donor_hits <- match(panel_stage_key, donor_stage_key)
if (anyNA(donor_hits) || anyDuplicated(donor_stage_key[donor_hits])) {
  stop("The fixed donor-stage assignments are absent or non-unique.")
}

manifest_key <- paste(spatial_manifest$donor, spatial_manifest$pair, sep = "||")
panel_key <- paste(panel_spec$donor, panel_spec$pair, sep = "||")
panel_hits <- match(panel_key, manifest_key)
if (anyNA(panel_hits) ||
    anyDuplicated(manifest_key[panel_hits]) ||
    sum(manifest_key %in% panel_key) != expected_selected_panels) {
  stop("The Step 04 manifest does not uniquely contain all six fixed panels.")
}
selected_source <- spatial_manifest[panel_hits, , drop = FALSE]
if (!identical(as.character(selected_source$stage), panel_spec$stage)) {
  stop("At least one fixed panel has an unexpected stage in the Step 04 manifest.")
}

pair_hits <- match(panel_spec$pair, interaction_spec$pair)
for (column in c("source", "target", "ligand", "receptor", "direction")) {
  expected_values <- interaction_spec[[column]][pair_hits]
  observed_values <- as.character(selected_source[[column]])
  if (!identical(observed_values, expected_values)) {
    stop("Selected-panel ", column, " values differ from the fixed specification.")
  }
}

selected_limit_hits <- match(selected_source$pair, resolved_limits$pair)
if (!isTRUE(all.equal(
  as.numeric(selected_source$color_q01),
  as.numeric(resolved_limits$q01[selected_limit_hits]),
  tolerance = 1e-15
)) || !isTRUE(all.equal(
  as.numeric(selected_source$color_q99),
  as.numeric(resolved_limits$q99[selected_limit_hits]),
  tolerance = 1e-15
))) {
  stop("Selected-panel color limits differ from the shared pair-specific limits.")
}

source_pdf_paths <- file.path(step04_dir, selected_source$relative_pdf_path)
if (any(!file.exists(source_pdf_paths)) ||
    any(file.info(source_pdf_paths)$size <= 0)) {
  stop("At least one fixed source PDF is missing or empty.")
}
source_pdf_sha256 <- vapply(source_pdf_paths, sha256_file, character(1))
if (!identical(
  tolower(unname(source_pdf_sha256)),
  tolower(trimws(as.character(selected_source$pdf_sha256)))
)) {
  stop("At least one fixed source PDF differs from its Step 04 checksum.")
}

message("Copying six byte-identical PDFs into the ordered main-panel folder...")
for (index in seq_len(expected_selected_panels)) {
  atomic_copy(source_pdf_paths[[index]], selected_pdf_paths[[index]])
}

selected_pdf_sha256 <- vapply(selected_pdf_paths, sha256_file, character(1))
if (!identical(unname(selected_pdf_sha256), unname(source_pdf_sha256))) {
  stop("At least one selected PDF is not byte-identical to its Step 04 source.")
}

selected_interactions <- data.frame(
  edge_order = interaction_spec$edge_order,
  source = interaction_spec$source,
  target = interaction_spec$target,
  ligand = interaction_spec$ligand,
  receptor = interaction_spec$receptor,
  pair = interaction_spec$pair,
  direction = interaction_spec$direction,
  liana_aggregate_rank = aggregate_ranks,
  aggregate_rank_threshold = aggregate_rank_threshold,
  color_q01 = resolved_limits$q01,
  color_q99 = resolved_limits$q99,
  color_limits_calculated_across_donors = 24L,
  finite_spots = resolved_limits$finite_spots,
  stringsAsFactors = FALSE
)

panel_pair_hits <- match(panel_spec$pair, selected_interactions$pair)
panel_manifest <- data.frame(
  panel_order = panel_spec$panel_order,
  figure_row = panel_spec$figure_row,
  figure_column = panel_spec$figure_column,
  stage = panel_spec$stage,
  donor = panel_spec$donor,
  source = selected_source$source,
  target = selected_source$target,
  ligand = selected_source$ligand,
  receptor = selected_source$receptor,
  pair = selected_source$pair,
  direction = selected_source$direction,
  liana_aggregate_rank = selected_interactions$liana_aggregate_rank[
    panel_pair_hits
  ],
  color_q01 = selected_source$color_q01,
  color_q99 = selected_source$color_q99,
  source_relative_pdf_path = vapply(
    source_pdf_paths,
    relative_to_repro,
    character(1)
  ),
  selected_relative_pdf_path = file.path(
    "results",
    "Fig3_CCC",
    "04_main_figure3C_selected",
    "panels",
    panel_spec$selected_pdf_filename
  ),
  pdf_bytes = unname(file.info(selected_pdf_paths)$size),
  pdf_sha256 = unname(selected_pdf_sha256),
  stringsAsFactors = FALSE
)

selection_parameters <- data.frame(
  parameter = c(
    "selection_basis",
    "historical_panel_selection_preceded_current_regenerated_output_visual_review",
    "regenerated_panel_substitution_based_on_signal_intensity",
    "selected_panels",
    "selected_interactions",
    "stages_displayed",
    "row_order",
    "column_order",
    "color_scale_scope",
    "color_limit_quantiles",
    "out_of_limit_values",
    "between_pair_absolute_score_comparison",
    "spatial_map_role",
    "scientific_analysis_performed_by_this_step"
  ),
  value = c(
    "Fixed to reproduce the accepted manuscript Figure 3C donor/pair composition",
    "TRUE",
    "FALSE",
    as.character(expected_selected_panels),
    as.character(expected_selected_interactions),
    "F0;F4",
    "top=F0;bottom=F4",
    "NPNT_ITGA8;MMP7_CD151;CCL2_ACKR1",
    "Identical limits across all 24 donors within each LR pair",
    "0.01;0.99",
    "clipped",
    "not valid",
    "descriptive; not used for statistical testing",
    "none; checksum-verified file selection only"
  ),
  stringsAsFactors = FALSE
)

input_manifest_paths <- c(
  step01_significant_csv,
  step01_checksums_tsv,
  step01_completion_tsv,
  step04_story_csv,
  step04_color_limits_csv,
  step04_donor_stage_csv,
  step04_spatial_manifest_csv,
  step04_checksums_tsv,
  step04_completion_tsv,
  source_pdf_paths
)
input_roles <- c(
  "LIANA significant interactions and aggregate ranks",
  "Step 01 output checksums",
  "Step 01 completion record",
  "fixed spatial LR interaction definitions",
  "pair-specific shared color limits",
  "donor-stage map",
  "72-panel spatial LR manifest",
  "Step 04 output checksums",
  "Step 04 completion record",
  rep("selected source LRlocal PDF", expected_selected_panels)
)
input_manifest <- data.frame(
  role = input_roles,
  relative_path = vapply(
    input_manifest_paths,
    relative_to_repro,
    character(1)
  ),
  bytes = unname(file.info(input_manifest_paths)$size),
  sha256 = vapply(input_manifest_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)

atomic_write_csv(panel_manifest, metadata_paths[["panel_manifest_csv"]])
atomic_write_csv(
  selected_interactions,
  metadata_paths[["selected_interactions_csv"]]
)
atomic_write_tsv(selection_parameters, metadata_paths[["parameters_tsv"]])
atomic_write_tsv(input_manifest, metadata_paths[["input_manifest_tsv"]])
atomic_write_lines(
  capture.output(utils::sessionInfo()),
  metadata_paths[["session_info_txt"]]
)

checksum_targets <- c(
  selected_pdf_paths,
  metadata_paths[["panel_manifest_csv"]],
  metadata_paths[["selected_interactions_csv"]],
  metadata_paths[["parameters_tsv"]],
  metadata_paths[["input_manifest_tsv"]],
  metadata_paths[["session_info_txt"]]
)
output_checksums <- data.frame(
  file = c(
    file.path("panels", basename(selected_pdf_paths)),
    basename(checksum_targets[-seq_along(selected_pdf_paths)])
  ),
  bytes = unname(file.info(checksum_targets)$size),
  sha256 = vapply(checksum_targets, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
atomic_write_tsv(
  output_checksums,
  metadata_paths[["output_checksums_tsv"]]
)

partial_files <- list.files(
  output_dir,
  pattern = "\\.partial$",
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE
)
if (length(partial_files)) {
  stop("Incomplete .partial outputs remain: ", paste(partial_files, collapse = ", "))
}

finished_at <- Sys.time()
completion <- data.frame(
  field = c(
    "status",
    "started_at",
    "finished_at",
    "elapsed_minutes",
    "script_path",
    "script_sha256",
    "selected_panels",
    "selected_interactions",
    "healthy_f0_panels",
    "f4_panels",
    "source_pdfs_verified",
    "byte_identical_copies",
    "partial_files"
  ),
  value = c(
    "complete",
    format(started_at, tz = "UTC", usetz = TRUE),
    format(finished_at, tz = "UTC", usetz = TRUE),
    sprintf("%.3f", as.numeric(difftime(finished_at, started_at, units = "mins"))),
    script_path,
    script_sha256,
    as.character(expected_selected_panels),
    as.character(expected_selected_interactions),
    as.character(sum(panel_spec$stage == "F0")),
    as.character(sum(panel_spec$stage == "F4")),
    as.character(length(source_pdf_paths)),
    as.character(sum(selected_pdf_sha256 == source_pdf_sha256)),
    as.character(length(partial_files))
  ),
  stringsAsFactors = FALSE
)
atomic_write_tsv(completion, metadata_paths[["completion_tsv"]])

message("Figure 3C fixed-panel selection completed at ", format(finished_at))
message("Selected PDFs copied: ", expected_selected_panels)
message("Scientific analysis performed: none")
message("Outputs: ", output_dir)
