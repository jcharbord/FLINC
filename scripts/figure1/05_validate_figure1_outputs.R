#!/usr/bin/env Rscript

repo_dir <- "/Users/jianwu/Documents/FLINC"

repro_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code"
)

scripts_dir <- file.path(
  repro_dir,
  "scripts"
)

results_dir <- file.path(
  repro_dir,
  "results"
)

logs_dir <- file.path(
  repro_dir,
  "reproducibility_logs"
)

audit_dir <- file.path(
  logs_dir,
  "Fig1_final_audit"
)

manifest_dir <- file.path(
  audit_dir,
  "small_manifests"
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  manifest_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


checksum_table <- function(paths, category) {

  paths <- unique(
    as.character(paths)
  )

  paths <- paths[
    !is.na(paths) &
      nzchar(paths)
  ]

  exists <- file.exists(paths)

  size_bytes <- rep(
    NA_real_,
    length(paths)
  )

  size_bytes[exists] <- file.info(
    paths[exists]
  )$size

  md5 <- rep(
    NA_character_,
    length(paths)
  )

  if (any(exists)) {

    cat(
      "MD5 ",
      category,
      ": ",
      sum(exists),
      " files\n",
      sep = ""
    )

    md5[exists] <- unname(
      tools::md5sum(
        paths[exists]
      )
    )
  }

  data.frame(
    category = category,
    path = paths,
    exists = exists,
    size_bytes = size_bytes,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# CANONICAL CODE
# ------------------------------------------------------------

canonical_scripts <- file.path(
  scripts_dir,
  c(
    "analysis/01_preprocess_integrate_and_reduce_spatial_data.R",
    "source_history/rebuild_clustering_source_from_historical_analysis.R",
    "analysis/02_cluster_annotate_and_characterize_spatial_domains.R",
    "source_history/clustering_annotation_characterization_source.Rmd",
    "source_history/clustering_annotation_characterization_extracted.R",
    "figure1/01_plot_spatial_clusters_on_histology.R",
    "figure1/02_plot_top_cluster_markers_spatially.R",
    "figure1/03_plot_final_annotated_umap.R",
    "figure1/05_validate_figure1_outputs.R"
  )
)

supporting_crop_scripts <- file.path(
  scripts_dir,
  "figure1/04_prepare_selected_histology_crop.R"
)

documentation <- c(
  file.path(repro_dir, "README.md"),
  file.path(repro_dir, "FIGURE1_SOURCE_WORKFLOW.md"),
  file.path(repro_dir, "FIGURE1_SHA256.txt"),
  file.path(scripts_dir, "analysis", "README.md"),
  file.path(scripts_dir, "figure1", "README.md"),
  file.path(scripts_dir, "source_history", "README.md")
)

code_checksums <- rbind(
  checksum_table(
    canonical_scripts,
    "canonical_script"
  ),
  checksum_table(
    supporting_crop_scripts,
    "supporting_crop_script"
  ),
  checksum_table(
    documentation,
    "documentation"
  )
)

write.csv(
  code_checksums,
  file.path(
    audit_dir,
    "code_and_documentation_checksums.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# CORE ANALYSIS FILES
# ------------------------------------------------------------

analysis_dir <- file.path(
  results_dir,
  "260701_approach1_cluster_go_marker_deconv"
)

final_object <- file.path(
  analysis_dir,
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)

top10_source <- file.path(
  analysis_dir,
  "Top10_markers_per_final_cluster_dotplot_data.csv"
)

color_map <- file.path(
  repo_dir,
  "data/st_analysis_pipeline/color_pseudotime.csv"
)

sample_order_file <- file.path(
  logs_dir,
  "sample_order_and_reference.csv"
)


# ------------------------------------------------------------
# SPACE RANGER H&E INPUTS
# ------------------------------------------------------------

hires_inputs <- character(0)

if (file.exists(sample_order_file)) {

  sample_table <- read.csv(
    sample_order_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if ("sample_id" %in% names(sample_table)) {

    hires_inputs <- file.path(
      "/Users/jianwu/Documents/DATA/spaceranger_output",
      sample_table$sample_id,
      "spatial",
      "tissue_hires_image.png"
    )
  }
}


critical_inputs <- c(
  final_object,
  top10_source,
  color_map,
  hires_inputs
)

critical_checksums <- checksum_table(
  critical_inputs,
  "critical_input"
)

write.csv(
  critical_checksums,
  file.path(
    audit_dir,
    "critical_input_checksums.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# FINAL UMAP
# ------------------------------------------------------------

umap_pdf <- file.path(
  analysis_dir,
  "UMAP_final_annotated_clusters_C4_on_top_C3green_C5red.pdf"
)

umap_png <- file.path(
  analysis_dir,
  "UMAP_final_annotated_clusters_C4_on_top_C3green_C5red.png"
)

umap_logs <- file.path(
  logs_dir,
  "UMAP_final_C5red",
  c(
    "run_parameters.txt",
    "sessionInfo.txt",
    "cluster_palette_C5red.csv",
    "UMAP_plot_order.csv",
    "final_cluster_counts.csv",
    "script_and_output_md5.csv"
  )
)


# ------------------------------------------------------------
# FINAL ALL-SAMPLE CLUSTER OVERLAYS
# ------------------------------------------------------------

cluster_dir <- file.path(
  results_dir,
  "all_samples_TRUE_Visium_clusters_C5red_circlelegend"
)

cluster_pdf <- list.files(
  cluster_dir,
  pattern = "\\.pdf$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

cluster_png <- list.files(
  cluster_dir,
  pattern = "\\.png$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

cluster_logs <- file.path(
  logs_dir,
  "all_samples_TRUE_Visium_clusters_C5red_circlelegend",
  c(
    "run_parameters.txt",
    "sessionInfo.txt",
    "cluster_palette_C5red.csv"
  )
)


# ------------------------------------------------------------
# FINAL DONOR-3D TOP10 SPATIAL MARKERS
# ------------------------------------------------------------

top10_dir <- file.path(
  results_dir,
  "3D_Top10_markers_by_cluster_TRUE_Visium_SCT"
)

top10_pdf <- list.files(
  top10_dir,
  pattern = "\\.pdf$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

top10_png <- list.files(
  top10_dir,
  pattern = "\\.png$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

top10_logs <- file.path(
  logs_dir,
  "3D_Top10_markers_by_cluster_TRUE_Visium_SCT",
  c(
    "run_parameters.txt",
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# CORE REPRODUCIBILITY LOGS
# ------------------------------------------------------------

core_logs <- file.path(
  logs_dir,
  c(
    "run_parameters_01_to_05.csv",
    "sample_order_and_reference.csv",
    "sessionInfo_start_01_to_05.txt",
    "ALL_main_object_md5_checksums_after_step05.csv",
    "ALL_outputs_06_to_10_as_Rscript_md5.csv",
    "sessionInfo_before_06_to_10_as_Rscript.txt",
    "sessionInfo_after_06_to_10_as_Rscript.txt"
  )
)

required_logs <- c(
  core_logs,
  cluster_logs,
  top10_logs,
  umap_logs
)

required_log_status <- data.frame(
  path = required_logs,
  exists = file.exists(required_logs),
  stringsAsFactors = FALSE
)

write.csv(
  required_log_status,
  file.path(
    audit_dir,
    "required_logs_status.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# SUPPORTING REPRESENTATIVE DONOR-3D FIGURE MATERIAL
# Not required for audit pass
# ------------------------------------------------------------

selected_crop_dir <- file.path(
  results_dir,
  "Fig1_Cropped_F3_sample3D"
)

selected_crop_marker_dir <- file.path(
  selected_crop_dir,
  "Top10_markers"
)

selected_crop_red_dir <- file.path(
  selected_crop_dir,
  "clusters_C5red"
)

selected_crop_yellow_dir <- file.path(
  selected_crop_dir,
  "clusters_C5yellow"
)

selected_crop_he_dir <- file.path(
  selected_crop_dir,
  "H_E"
)

selected_crop_records_dir <- file.path(
  selected_crop_dir,
  "reproducibility_records"
)

selected_crop_marker_pdf <- character(0)
selected_crop_marker_png <- character(0)

if (dir.exists(selected_crop_marker_dir)) {

  selected_crop_marker_pdf <- list.files(
    selected_crop_marker_dir,
    pattern = "\\.pdf$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  selected_crop_marker_png <- list.files(
    selected_crop_marker_dir,
    pattern = "\\.png$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

selected_crop_red_pdf <- file.path(
  selected_crop_red_dir,
  "Fig1_sample3D_clusters_C5red.pdf"
)

selected_crop_red_png <- file.path(
  selected_crop_red_dir,
  "Fig1_sample3D_clusters_C5red.png"
)

selected_crop_yellow_pdf <- file.path(
  selected_crop_yellow_dir,
  "Fig1_sample3D_clusters_C5yellow.pdf"
)

selected_crop_yellow_png <- file.path(
  selected_crop_yellow_dir,
  "Fig1_sample3D_clusters_C5yellow.png"
)

selected_crop_he_source <- file.path(
  selected_crop_he_dir,
  "3D_tissue_hires_image_SOURCE.png"
)

selected_crop_files <- character(0)

if (dir.exists(selected_crop_dir)) {

  selected_crop_files <- list.files(
    selected_crop_dir,
    pattern = "\\.(pdf|png|csv|rds|tsv|txt)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

selected_crop_checksums <- checksum_table(
  selected_crop_files,
  "supporting_representative_3D_output"
)

write.csv(
  selected_crop_checksums,
  file.path(
    audit_dir,
    "supporting_selected_crop_output_checksums.csv"
  ),
  row.names = FALSE
)

selected_crop_red_ok <- all(
  file.exists(
    c(
      selected_crop_red_pdf,
      selected_crop_red_png
    )
  )
)

selected_crop_yellow_ok <- all(
  file.exists(
    c(
      selected_crop_yellow_pdf,
      selected_crop_yellow_png
    )
  )
)

selected_crop_ok <- (
  dir.exists(selected_crop_dir) &&
    length(selected_crop_marker_pdf) == 70 &&
    length(selected_crop_marker_png) == 70 &&
    selected_crop_red_ok &&
    selected_crop_yellow_ok
)

selected_crop_status <- data.frame(
  component = "Fig1_Cropped_F3_sample3D",
  directory_exists = dir.exists(selected_crop_dir),
  marker_pdf_count = length(selected_crop_marker_pdf),
  marker_png_count = length(selected_crop_marker_png),
  C5red_pair_exists = selected_crop_red_ok,
  C5yellow_pair_exists = selected_crop_yellow_ok,
  registered_H_E_source_copy_exists =
    file.exists(selected_crop_he_source),
  required_for_audit = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  selected_crop_status,
  file.path(
    audit_dir,
    "supporting_selected_crop_status.csv"
  ),
  row.names = FALSE
)

selected_crop_logs <- file.path(
  selected_crop_records_dir,
  c(
    "run_parameters.txt",
    "sessionInfo.txt",
    "generation_validation.csv"
  )
)

write.csv(
  data.frame(
    path = selected_crop_logs,
    exists = file.exists(selected_crop_logs),
    required_for_audit = FALSE,
    stringsAsFactors = FALSE
  ),
  file.path(
    audit_dir,
    "supporting_selected_crop_logs_status.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# COPY SMALL RESULT CSV MANIFESTS INTO AUDIT AREA
# ------------------------------------------------------------

small_result_csvs <- c(
  list.files(
    cluster_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  ),
  list.files(
    top10_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
)

if (length(small_result_csvs) > 0) {

  for (src in small_result_csvs) {

    file.copy(
      src,
      file.path(
        manifest_dir,
        basename(src)
      ),
      overwrite = TRUE
    )
  }
}


# ------------------------------------------------------------
# CHECKSUM FINAL FIGURE SOURCE OUTPUTS
# ------------------------------------------------------------

final_source_outputs <- c(
  umap_pdf,
  umap_png,
  cluster_pdf,
  cluster_png,
  top10_pdf,
  top10_png
)

output_checksums <- checksum_table(
  final_source_outputs,
  "final_Fig1_source_output"
)

write.csv(
  output_checksums,
  file.path(
    audit_dir,
    "final_Fig1_source_output_checksums.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# PALETTE CONSISTENCY
# ------------------------------------------------------------

script14 <- paste(
  readLines(
    file.path(
      scripts_dir,
      "figure1/01_plot_spatial_clusters_on_histology.R"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)

script18 <- paste(
  readLines(
    file.path(
      scripts_dir,
      "figure1/03_plot_final_annotated_umap.R"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)

palette_check <- data.frame(
  component = c(
    "all_sample_TRUE_Visium_script",
    "final_UMAP_script"
  ),
  contains_C5_red_E15759 = c(
    grepl(
      "#E15759",
      script14,
      fixed = TRUE
    ),
    grepl(
      "#E15759",
      script18,
      fixed = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  palette_check,
  file.path(
    audit_dir,
    "C5red_palette_consistency.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# EXPECTED OUTPUT COUNTS
# ------------------------------------------------------------

cluster_output_ok <- (
  length(cluster_pdf) == 24 &&
    length(cluster_png) == 24
)

top10_output_ok <- (
  length(top10_pdf) == 70 &&
    length(top10_png) == 70
)

umap_output_ok <- (
  file.exists(umap_pdf) &&
    file.exists(umap_png)
)

canonical_code_ok <- all(
  file.exists(canonical_scripts)
)

critical_inputs_ok <- all(
  critical_checksums$exists
)

required_logs_ok <- all(
  required_log_status$exists
)

palette_ok <- all(
  palette_check$contains_C5_red_E15759
)


audit_pass <- all(
  canonical_code_ok,
  critical_inputs_ok,
  required_logs_ok,
  cluster_output_ok,
  top10_output_ok,
  umap_output_ok,
  palette_ok
)


component_status <- data.frame(
  component = c(
    "Analysis and figure source code",
    "Critical analytical inputs",
    "Core analysis and plotting run records",
    "Final C5-red UMAP",
    "24 all-sample cluster PDFs and PNGs",
    "70 top-marker PDFs and PNGs",
    "C5 red palette consistency",
    "Selected crop material"
  ),
  required = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE
  ),
  status = c(
    canonical_code_ok,
    critical_inputs_ok,
    required_logs_ok,
    umap_output_ok,
    cluster_output_ok,
    top10_output_ok,
    palette_ok,
    selected_crop_ok
  ),
  stringsAsFactors = FALSE
)

write.csv(
  component_status,
  file.path(
    audit_dir,
    "Fig1_component_status.csv"
  ),
  row.names = FALSE
)


summary_lines <- c(
  "FLINC Figure 1 source-workflow validation",
  paste0(
    "Timestamp: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S %Z"
    )
  ),
  "",
  paste0(
    "Canonical code OK: ",
    canonical_code_ok
  ),
  paste0(
    "Critical inputs OK: ",
    critical_inputs_ok
  ),
  paste0(
    "Required logs OK: ",
    required_logs_ok
  ),
  "",
  paste0(
    "Final C5-red UMAP PDF/PNG: ",
    umap_output_ok
  ),
  paste0(
    "Cluster PDFs: ",
    length(cluster_pdf),
    " / expected 24"
  ),
  paste0(
    "Cluster PNGs: ",
    length(cluster_png),
    " / expected 24"
  ),
  paste0(
    "Top10 marker PDFs: ",
    length(top10_pdf),
    " / expected 70"
  ),
  paste0(
    "Top10 marker PNGs: ",
    length(top10_png),
    " / expected 70"
  ),
  "",
  paste0(
    "C5 red palette consistency: ",
    palette_ok
  ),
  "",
  paste0(
    "Selected crop material present: ",
    selected_crop_ok
  ),
  "",
  "Final Figure 1 source scripts:",
  "figure1/01_plot_spatial_clusters_on_histology.R",
  "figure1/02_plot_top_cluster_markers_spatially.R",
  "figure1/03_plot_final_annotated_umap.R",
  "",
  "Supporting histology-crop script:",
  "figure1/04_prepare_selected_histology_crop.R",
  "",
  paste0(
    "FINAL_AUDIT_PASS = ",
    audit_pass
  )
)

writeLines(
  summary_lines,
  file.path(
    audit_dir,
    "audit_summary.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    audit_dir,
    "sessionInfo.txt"
  )
)

cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)

if (!audit_pass) {

  stop(
    "Figure 1 audit did not fully pass. ",
    "Inspect reproducibility_logs/Fig1_final_audit/"
  )
}

cat("\nFINAL AUDIT PASSED.\n")
