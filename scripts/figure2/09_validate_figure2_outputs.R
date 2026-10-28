#!/usr/bin/env Rscript

# FLINC Figure 2 - final reproducibility audit.
#
# This script checks the completeness and internal consistency of the public
# Figure 2 workflow after operations 01 through 08 have finished. It performs
# no additional biological analysis.

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

nonempty_files <- function(paths) {
  if (!length(paths)) return(logical(0))
  exists <- file.exists(paths) & !dir.exists(paths)
  sizes <- rep(NA_real_, length(paths))
  sizes[exists] <- file.info(paths[exists])$size
  exists & !is.na(sizes) & sizes > 0
}

checksum_table <- function(paths, category) {
  paths <- sort(unique(as.character(paths)))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  exists <- file.exists(paths) & !dir.exists(paths)
  size_bytes <- rep(NA_real_, length(paths))
  md5 <- rep(NA_character_, length(paths))
  if (any(exists)) {
    size_bytes[exists] <- file.info(paths[exists])$size
    md5[exists] <- unname(tools::md5sum(paths[exists]))
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

repro_root <- env_or(
  "FLINC_REPRO_ROOT",
  "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)

scripts_dir <- file.path(repro_root, "scripts", "figure2")
objects_dir <- file.path(repro_root, "objects", "Fig2_WGCNA")
results_dir <- file.path(repro_root, "results", "Fig2_WGCNA")
logs_dir <- file.path(repro_root, "reproducibility_logs", "Fig2_WGCNA")
audit_dir <- env_or(
  "FIG2_AUDIT_OUTPUT_DIR",
  file.path(logs_dir, "Fig2_final_audit")
)

enrichment_dir <- env_or(
  "FIG2_STEP05_OUTPUT_DIR",
  file.path(results_dir, "05_module_tables_and_enrichr")
)
enrichment_checkpoint_dir <- env_or(
  "FIG2_STEP05_CHECKPOINT_DIR",
  file.path(results_dir, "05_enrichr_checkpoints")
)
hub_network_dir <- env_or(
  "FIG2_STEP06_OUTPUT_DIR",
  file.path(results_dir, "06_colored_hub_networks")
)
module_expression_dir <- env_or(
  "FIG2_STEP07_OUTPUT_DIR",
  file.path(results_dir, "07_module_expression")
)
supplementary_umap_dir <- env_or(
  "FIG2_SUPPLEMENTARY_UMAP_OUTPUT_DIR",
  file.path(results_dir, "07_supplementary_wgcna_umaps")
)
gsea_results_root <- env_or(
  "FIG2_GSEA_RESULTS_DIR",
  file.path(repro_root, "results", "Fig2_GSEA", "08_stagewise_gsea")
)
gsea_logs_root <- env_or(
  "FIG2_GSEA_LOG_DIR",
  file.path(
    repro_root,
    "reproducibility_logs",
    "Fig2_GSEA",
    "08_stagewise_gsea"
  )
)

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

sha256_file <- function(path) {
  digest::digest(object = path, algo = "sha256", serialize = FALSE, file = TRUE)
}

resolve_gsea_run <- function(results_root, logs_root, requested_run = "") {
  run_pattern <- "^run_[0-9]{8}T[0-9]{6}Z$"
  result_runs <- if (dir.exists(results_root)) {
    basename(list.dirs(results_root, recursive = FALSE, full.names = TRUE))
  } else {
    character(0)
  }
  log_runs <- if (dir.exists(logs_root)) {
    basename(list.dirs(logs_root, recursive = FALSE, full.names = TRUE))
  } else {
    character(0)
  }
  result_runs <- result_runs[grepl(run_pattern, result_runs)]
  log_runs <- log_runs[grepl(run_pattern, log_runs)]
  paired <- sort(intersect(result_runs, log_runs))
  paired <- paired[
    file.exists(file.path(results_root, paired, "RESULTS_COMPLETE.txt")) &
      file.exists(file.path(logs_root, paired, "COMPLETED.txt"))
  ]

  if (nzchar(requested_run)) {
    if (!requested_run %in% paired) {
      return(list(
        ok = FALSE,
        run_name = requested_run,
        detail = paste("requested completed run pair not found:", requested_run)
      ))
    }
    selected <- requested_run
  } else if (length(paired)) {
    selected <- tail(paired, 1L)
  } else {
    return(list(
      ok = FALSE,
      run_name = NA_character_,
      detail = "no paired completed GSEA result/log run was found"
    ))
  }

  list(
    ok = TRUE,
    run_name = selected,
    detail = paste("selected paired completed GSEA run:", selected)
  )
}

gsea_run <- resolve_gsea_run(
  gsea_results_root,
  gsea_logs_root,
  Sys.getenv("FIG2_GSEA_RUN_NAME", unset = "")
)
gsea_results_dir <- if (gsea_run$ok) {
  file.path(gsea_results_root, gsea_run$run_name)
} else {
  file.path(gsea_results_root, "__NO_COMPLETED_RUN__")
}
gsea_logs_dir <- if (gsea_run$ok) {
  file.path(gsea_logs_root, gsea_run$run_name)
} else {
  file.path(gsea_logs_root, "__NO_COMPLETED_RUN__")
}
record_check(
  "Paired completed supplementary GSEA run",
  gsea_run$ok,
  gsea_run$detail
)

modules <- paste0("SM_", 1:6)
comparisons <- paste0(c("F1", "F2", "F3", "F4"), "_vs_F0")
cell_types <- c(
  "Bcells", "CentralHep", "Cholangiocyte", "Erythroid", "InfMac",
  "InterHep", "NKTcell", "NonInfMac", "PortalEndo", "PortalHep",
  "Stellate", "cvLSECs"
)
enrichr_databases <- c(
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "KEGG_2019_Human",
  "WikiPathways_2024_Human"
)

# -------------------------------------------------------------------------
# Public code and documentation
# -------------------------------------------------------------------------

canonical_scripts <- file.path(
  scripts_dir,
  c(
    "01_prepare_wgcna_input.R",
    "02_build_wgcna_network.R",
    "03_analyze_differential_module_eigengenes.R",
    "03_bis_analyze_differential_module_eigengenes_donor_level.R",
    "04_export_module_tables_and_run_enrichr.R",
    "05_generate_functionally_colored_hub_networks.R",
    "05_bis_build_ego_networks.R",
    "06_plot_harmonized_module_eigengenes_and_celltype_associations.R",
    "07_generate_supplementary_wgcna_umaps.R",
    "08_run_supplementary_stage_gsea.R",
    "09_validate_figure2_outputs.R"
  )
)
documentation <- c(
  file.path(repro_root, "README.md"),
  file.path(scripts_dir, "README.md"),
  file.path(repro_root, "FIGURE2_SOURCE_WORKFLOW.md"),
  file.path(repro_root, "FIGURE2_SHA256.txt")
)
reference_files <- file.path(
  scripts_dir,
  "reference",
  c("functional_subcategories.tsv", "ego_network_functional_context_annotations.tsv")
)

check_required_files(
  "Public Figure 2 code, documentation and reference table",
  c(canonical_scripts, documentation, reference_files)
)

say("Calculating public code, documentation and reference checksums")
code_checksums <- rbind(
  checksum_table(canonical_scripts, "canonical_script"),
  checksum_table(documentation, "documentation"),
  checksum_table(reference_files, "reference_table")
)
utils::write.table(
  code_checksums,
  file.path(audit_dir, "code_documentation_and_reference_checksums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

# -------------------------------------------------------------------------
# Required objects, tables, figures and execution records
# -------------------------------------------------------------------------

prepared_object <- file.path(objects_dir, "01_wgcna_input_prepared.rds")
network_object <- file.path(objects_dir, "02_wgcna_network_modules.rds")

network_output_files <- file.path(
  results_dir,
  "02_network",
  c(
    "02_soft_power_table.tsv",
    "02_soft_powers.pdf",
    "02_module_dendrogram.pdf",
    "02_module_assignments_and_kME.tsv",
    "02_module_sizes.tsv",
    "02_harmonized_module_eigengenes.tsv"
  )
)
network_log_files <- file.path(
  logs_dir,
  c(
    "01_wgcna_input_summary.tsv",
    "01_wgcna_input_sessionInfo.txt",
    "02_wgcna_parameters.txt",
    "02_wgcna_network_summary.tsv",
    "02_wgcna_network_sessionInfo.txt"
  )
)

dme_files <- c(
  file.path(
    results_dir,
    "03_DME_donor_level_2026-10-13",
    c(
      "03_DiffModuleEigengenes_donor_level.csv",
      "03_DiffModuleEigengenes_donor_level.tsv",
      "03_DiffModuleEigengenes_donor_level.xlsx",
      "03_donor_mean_module_eigengenes.csv",
      "03_QC_spots_per_donor_and_module.csv",
      "03_QC_donors_per_stage.csv",
      "03_module_assignments_used.csv",
      "03_DME_donor_level_lollipop_all_stages_vs_Healthy.pdf",
      "03_DME_donor_level_volcano_all_stages_vs_Healthy.pdf"
    )
  ),
  file.path(logs_dir, c("03_DME_donor_level_2026-10-13_summary.tsv", "03_DME_donor_level_2026-10-13_sessionInfo.txt"))
)

enrichment_root_files <- file.path(
  enrichment_dir,
  c(
    "05_input_manifest.tsv",
    "05_run_summary.tsv",
    "05_sessionInfo.txt",
    "05_completion_report.txt"
  )
)
enrichment_table_files <- file.path(
  enrichment_dir,
  "tables",
  c(
    "05_module_metadata.tsv",
    "05_all_genes_assignments_and_kME.tsv",
    "05_non_grey_genes_ranked_by_module_kME.tsv",
    "05_top70_genes_ranked_by_module_kME.tsv",
    "05_displayed_hubs_top50_SM1_top30_SM2_to_SM6.tsv",
    "05_enrichr_all_raw_results.tsv",
    "05_enrichr_significant_top10_for_plots.tsv",
    "05_enrichr_query_and_plot_index.tsv",
    "05_enrichr_database_catalog_snapshot.tsv"
  )
)
enrichment_checkpoint_files <- file.path(
  enrichment_checkpoint_dir,
  paste0("05_enrichr_", modules, ".rds")
)

panels <- setNames(LETTERS[3:8], modules)
hub_individual_bases <- vapply(
  modules,
  function(module) sprintf(
    "06_panel_%s_%s_colored_hub_network",
    panels[[module]],
    module
  ),
  character(1)
)
hub_plot_files <- c(
  file.path(hub_network_dir, "plots", paste0(hub_individual_bases, ".pdf")),
  file.path(hub_network_dir, "plots", paste0(hub_individual_bases, ".png")),
  file.path(
    hub_network_dir,
    "plots",
    c(
      "06_panels_C_to_H_colored_hub_networks.pdf",
      "06_panels_C_to_H_colored_hub_networks.png"
    )
  )
)
hub_table_files <- file.path(
  hub_network_dir,
  c(
    "06_hubs_with_functional_subcategories.tsv",
    "06_functional_subcategory_legend.tsv",
    "06_functional_subcategory_counts.tsv",
    "06_sessionInfo.txt",
    "06_completion_report.txt"
  )
)

spatial_map_bases <- paste0("07_", modules, "_hME_spatial_24_samples")
module_expression_files <- c(
  file.path(
    module_expression_dir,
    "dotplots",
    c("07_all_hMEs_by_cluster.pdf", "07_all_hMEs_by_cluster.png")
  ),
  file.path(
    module_expression_dir,
    "spatial_maps",
    c(
      paste0(spatial_map_bases, ".pdf"),
      paste0(spatial_map_bases, ".png"),
      "07_all_modules_hME_spatial_24_samples.pdf"
    )
  ),
  file.path(
    module_expression_dir,
    "correlations",
    c(
      "07_harmonized_module_eigengenes.tsv",
      "07_hME_celltype_spearman_long.csv",
      "07_hME_celltype_spearman_matrix.csv",
      "07_strongest_positive_driver_by_module.csv",
      "07_hME_celltype_association_heatmap.pdf"
    )
  ),
  file.path(
    module_expression_dir,
    "validation",
    c("07_summary.tsv", "07_sessionInfo.txt", "07_validation.tsv")
  )
)

supplementary_umap_files <- file.path(
  supplementary_umap_dir,
  c(
    "07_module_gene_umap_coordinates.tsv",
    "07_supplementary_wgcna_umap_parameters.tsv",
    "hubgene_umap_ggplot.pdf",
    "hubgene_umap_ggplot.png",
    "umap.MEs.pdf",
    "umap.MEs.png",
    "07_supplementary_wgcna_umaps_sessionInfo.txt",
    "07_supplementary_wgcna_umaps_completion_report.txt"
  )
)

gsea_collections <- c("H", "C2", "C3", "C5")
gsea_de_files <- file.path(
  gsea_results_dir,
  "de",
  paste0(comparisons, "_limma_voom_DE.tsv.gz")
)
gsea_ranking_files <- c(
  file.path(
    gsea_results_dir,
    "rankings",
    paste0(comparisons, "_moderated_t_ranking.tsv.gz")
  ),
  file.path(
    gsea_results_dir,
    "rankings",
    paste0(comparisons, "_gene_id_to_symbol_mapping.tsv.gz")
  )
)
gsea_keys <- expand.grid(
  comparison = comparisons,
  collection = gsea_collections,
  stringsAsFactors = FALSE
)
gsea_full_files <- c(
  file.path(
    gsea_results_dir,
    "gsea_full",
    paste0(
      gsea_keys$comparison,
      "_",
      gsea_keys$collection,
      "_GSEA_full.tsv.gz"
    )
  ),
  file.path(
    gsea_results_dir,
    "gsea_full",
    "all_comparisons_all_collections_GSEA_full.tsv.gz"
  )
)
gsea_significant_files <- c(
  file.path(
    gsea_results_dir,
    "gsea_significant",
    paste0(
      gsea_keys$comparison,
      "_",
      gsea_keys$collection,
      "_GSEA_FDR_lt_0.05.tsv.gz"
    )
  ),
  file.path(
    gsea_results_dir,
    "gsea_significant",
    "all_comparisons_all_collections_GSEA_FDR_lt_0.05.tsv.gz"
  )
)
gsea_figure_files <- file.path(
  gsea_results_dir,
  "figures",
  paste0(comparisons, "_H_dotplot.pdf")
)
gsea_qc_files <- file.path(
  gsea_results_dir,
  "qc",
  c("pseudobulk_MDS.pdf", "voom_mean_variance.pdf")
)
gsea_manifest_files <- file.path(
  gsea_results_dir,
  "manifests",
  c(
    "gene_filter_manifest.tsv.gz",
    "msigdb_term2gene_H_C2_C3_C5.tsv.gz",
    "msigdb_collection_manifest.tsv",
    "ranking_manifest.tsv"
  )
)
gsea_object_files <- file.path(
  gsea_results_dir,
  "objects",
  c("donor_pseudobulk_raw_counts.rds", "limma_voom_contrast_fit.rds")
)
gsea_root_files <- file.path(
  gsea_results_dir,
  c("GSEA_significance_counts.tsv", "RESULTS_COMPLETE.txt")
)
gsea_execution_files <- file.path(
  gsea_logs_dir,
  c(
    "analysis_scope.txt",
    "COMPLETED.txt",
    "execution_record_SHA256SUMS.tsv",
    "input_checksums.tsv",
    "input_object_summary.tsv",
    "package_versions.tsv",
    "result_SHA256SUMS.tsv",
    "run_parameters.txt",
    "sample_metadata.tsv",
    "sessionInfo.txt"
  )
)

check_required_files(
  "Prepared and completed WGCNA objects",
  c(prepared_object, network_object)
)
check_required_files(
  "Network outputs and execution records",
  c(network_output_files, network_log_files)
)
check_required_files("Differential hME outputs and execution records", dme_files)
check_required_files(
  "Module tables and EnrichR outputs",
  c(enrichment_root_files, enrichment_table_files, enrichment_checkpoint_files)
)
check_required_files(
  "Functionally colored hub-network outputs",
  c(hub_plot_files, hub_table_files)
)
check_required_files(
  "Harmonized ME plots and cell-type association outputs",
  module_expression_files
)
check_required_files(
  "Supplementary module-gene and harmonized-ME UMAP outputs",
  supplementary_umap_files
)
check_required_files(
  "Supplementary donor-pseudobulk GSEA outputs",
  c(
    gsea_de_files,
    gsea_ranking_files,
    gsea_full_files,
    gsea_significant_files,
    gsea_figure_files,
    gsea_qc_files,
    gsea_manifest_files,
    gsea_object_files,
    gsea_root_files
  )
)
check_required_files(
  "Supplementary donor-pseudobulk GSEA execution records",
  gsea_execution_files
)

# -------------------------------------------------------------------------
# Tabular integrity checks
# -------------------------------------------------------------------------

module_table_file <- file.path(
  results_dir,
  "02_network",
  "02_module_assignments_and_kME.tsv"
)
module_table_result <- safe_table(module_table_file, read_tsv)
module_table_ok <- FALSE
module_table_detail <- module_table_result$detail
if (module_table_result$ok) {
  module_table <- module_table_result$data
  required_columns <- c(
    "gene_name", "module", "color", paste0("kME_SM_", 1:6)
  )
  missing_columns <- setdiff(required_columns, colnames(module_table))
  if (length(missing_columns)) {
    module_table_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    expected_sizes <- c(
      SM_1 = 1652L, SM_2 = 50L, SM_3 = 69L,
      SM_4 = 90L, SM_5 = 53L, SM_6 = 292L
    )
    observed_sizes <- table(as.character(module_table$module))
    size_ok <- all(
      names(expected_sizes) %in% names(observed_sizes)
    ) && identical(
      as.integer(observed_sizes[names(expected_sizes)]),
      as.integer(expected_sizes)
    )
    module_table_ok <- nrow(module_table) == 11515L &&
      setequal(unique(as.character(module_table$module)), c(modules, "grey")) &&
      !anyDuplicated(module_table$gene_name) &&
      size_ok
    module_table_detail <- sprintf(
      "%d genes; modules: %s",
      nrow(module_table),
      paste(sort(unique(as.character(module_table$module))), collapse = ", ")
    )
  }
}
record_check("Network module assignments and harmonized kME structure", module_table_ok, module_table_detail)

hme_table_file <- file.path(
  results_dir,
  "02_network",
  "02_harmonized_module_eigengenes.tsv"
)
hme_table_result <- safe_table(hme_table_file, read_tsv)
hme_table_ok <- FALSE
hme_table_detail <- hme_table_result$detail
if (hme_table_result$ok) {
  hme_table <- hme_table_result$data
  required_columns <- c("barcode", modules)
  missing_columns <- setdiff(required_columns, colnames(hme_table))
  if (length(missing_columns)) {
    hme_table_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    hme_values <- suppressWarnings(
      apply(hme_table[, modules, drop = FALSE], 2, as.numeric)
    )
    hme_table_ok <- nrow(hme_table) == 26056L &&
      !anyDuplicated(hme_table$barcode) &&
      all(is.finite(hme_values))
    hme_table_detail <- sprintf(
      "%d spots; expected hME columns present: TRUE",
      nrow(hme_table)
    )
  }
}
record_check("Harmonized module eigengene export", hme_table_ok, hme_table_detail)

dme_table_file <- file.path(
  results_dir,
  "03_DME_donor_level_2026-10-13",
  "03_DiffModuleEigengenes_donor_level.csv"
)
dme_result <- safe_table(dme_table_file, read_csv)
dme_ok <- FALSE
dme_detail <- dme_result$detail
if (dme_result$ok) {
  dme <- dme_result$data
  required_columns <- c(
    "module", "comparison", "delta_mean", "p_wilcox", "p_adj",
    "significant_FDR_0.05", "paired", "n_baseline_used", "n_stage_used",
    "mean_baseline", "mean_stage", "me_source", "donor_column"
  )
  missing_columns <- setdiff(required_columns, colnames(dme))
  if (length(missing_columns)) {
    dme_detail <- paste("missing columns:", paste(missing_columns, collapse = ", "))
  } else {
    keys <- paste(dme$module, dme$comparison, sep = "::")
    expected_keys <- as.vector(outer(modules, comparisons, paste, sep = "::"))
    dme_ok <- nrow(dme) == 24L && !anyDuplicated(keys) &&
      setequal(keys, expected_keys) &&
      all(dme$paired == FALSE) && all(dme$me_source == "harmonized") &&
      all(dme$n_baseline_used == 6L) &&
      all(dme$n_stage_used == ifelse(dme$comparison %in% c("F2_vs_F0", "F3_vs_F0"), 3L, 6L)) &&
      isTRUE(all.equal(dme$delta_mean, dme$mean_stage - dme$mean_baseline)) &&
      isTRUE(all.equal(dme$p_adj, ave(dme$p_wilcox, dme$comparison,
        FUN = function(x) stats::p.adjust(x, method = "BH"))))
    dme_detail <- sprintf("%d module-by-stage tests", nrow(dme))
  }
}
record_check("Donor-level differential harmonized ME table", dme_ok, dme_detail)

enrichment_index_file <- file.path(
  enrichment_dir,
  "tables",
  "05_enrichr_query_and_plot_index.tsv"
)
enrichment_index_result <- safe_table(enrichment_index_file, read_tsv)
enrichment_ok <- FALSE
enrichment_detail <- enrichment_index_result$detail
if (enrichment_index_result$ok) {
  index <- enrichment_index_result$data
  required_columns <- c("module", "database", "n_query_genes")
  missing_columns <- setdiff(required_columns, colnames(index))
  if (length(missing_columns)) {
    enrichment_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    keys <- paste(index$module, index$database, sep = "::")
    expected_keys <- as.vector(outer(modules, enrichr_databases, paste, sep = "::"))
    raw_files <- list.files(
      file.path(enrichment_dir, "raw_enrichr"),
      pattern = "__raw\\.tsv$",
      full.names = TRUE
    )
    gene_lists <- list.files(
      file.path(enrichment_dir, "gene_lists"),
      pattern = "\\.txt$",
      full.names = TRUE
    )
    enrichment_ok <- nrow(index) == 30L && !anyDuplicated(keys) &&
      setequal(keys, expected_keys) && length(raw_files) == 30L &&
      length(gene_lists) == 18L && all(nonempty_files(c(raw_files, gene_lists)))
    enrichment_detail <- sprintf(
      "%d module/database queries; %d raw tables; %d gene lists",
      nrow(index), length(raw_files), length(gene_lists)
    )
  }
}
record_check("EnrichR query and export structure", enrichment_ok, enrichment_detail)

hub_table_file <- file.path(
  hub_network_dir,
  "06_hubs_with_functional_subcategories.tsv"
)
hub_result <- safe_table(hub_table_file, read_tsv)
hub_ok <- FALSE
hub_detail <- hub_result$detail
if (hub_result$ok) {
  hubs <- hub_result$data
  required_columns <- c("module", "gene_name")
  missing_columns <- setdiff(required_columns, colnames(hubs))
  if (length(missing_columns)) {
    hub_detail <- paste("missing columns:", paste(missing_columns, collapse = ", "))
  } else {
    expected_hub_counts <- c(
      SM_1 = 50L, SM_2 = 30L, SM_3 = 30L,
      SM_4 = 30L, SM_5 = 30L, SM_6 = 30L
    )
    observed_hub_counts <- table(as.character(hubs$module))
    hub_keys <- paste(hubs$module, hubs$gene_name, sep = "::")
    hub_ok <- nrow(hubs) == 200L && !anyDuplicated(hub_keys) &&
      setequal(names(observed_hub_counts), modules) &&
      identical(
        as.integer(observed_hub_counts[modules]),
        as.integer(expected_hub_counts[modules])
      )
    hub_detail <- sprintf("%d unique displayed module/gene pairs", nrow(hubs))
  }
}
record_check("Displayed hub-gene output structure", hub_ok, hub_detail)

correlation_file <- file.path(
  module_expression_dir,
  "correlations",
  "07_hME_celltype_spearman_long.csv"
)
correlation_result <- safe_table(correlation_file, read_csv)
correlation_ok <- FALSE
correlation_detail <- correlation_result$detail
if (correlation_result$ok) {
  correlations <- correlation_result$data
  required_columns <- c(
    "module", "cell_type", "rho", "p_value", "p_adj_global",
    "p_adj_within_module"
  )
  missing_columns <- setdiff(required_columns, colnames(correlations))
  if (length(missing_columns)) {
    correlation_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    keys <- paste(correlations$module, correlations$cell_type, sep = "::")
    expected_keys <- as.vector(outer(modules, cell_types, paste, sep = "::"))
    correlation_ok <- nrow(correlations) == 72L && !anyDuplicated(keys) &&
      setequal(keys, expected_keys)
    correlation_detail <- sprintf(
      "%d module-by-cell-type correlations",
      nrow(correlations)
    )
  }
}
record_check("Harmonized ME-cell-type correlation table", correlation_ok, correlation_detail)

validation_file <- file.path(
  module_expression_dir,
  "validation",
  "07_validation.tsv"
)
validation_result <- safe_table(validation_file, read_tsv)
validation_ok <- FALSE
validation_detail <- validation_result$detail
if (validation_result$ok) {
  validation <- validation_result$data
  required_columns <- c("check", "pass")
  missing_columns <- setdiff(required_columns, colnames(validation))
  if (length(missing_columns)) {
    validation_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    pass_values <- as.logical(validation$pass)
    validation_ok <- nrow(validation) >= 8L &&
      !anyNA(pass_values) && all(pass_values)
    validation_detail <- sprintf(
      "%d/%d operation-level checks passed",
      sum(pass_values, na.rm = TRUE),
      length(pass_values)
    )
  }
}
record_check("Operation 06 internal validation record", validation_ok, validation_detail)

umap_coordinates_file <- file.path(
  supplementary_umap_dir,
  "07_module_gene_umap_coordinates.tsv"
)
umap_parameters_file <- file.path(
  supplementary_umap_dir,
  "07_supplementary_wgcna_umap_parameters.tsv"
)
umap_coordinates_result <- safe_table(umap_coordinates_file, read_tsv)
umap_parameters_result <- safe_table(umap_parameters_file, read_tsv)
umap_ok <- FALSE
umap_detail <- paste(
  umap_coordinates_result$detail,
  umap_parameters_result$detail,
  sep = "; "
)
if (umap_coordinates_result$ok && umap_parameters_result$ok) {
  coordinates <- umap_coordinates_result$data
  umap_parameters <- umap_parameters_result$data
  coordinate_columns <- c("gene", "module", "color", "kME", "UMAP1", "UMAP2")
  parameter_columns <- c("parameter", "value")
  missing_columns <- c(
    setdiff(coordinate_columns, colnames(coordinates)),
    setdiff(parameter_columns, colnames(umap_parameters))
  )
  if (length(missing_columns)) {
    umap_detail <- paste("missing columns:", paste(unique(missing_columns), collapse = ", "))
  } else {
    expected_umap_sizes <- c(
      SM_1 = 1652L, SM_2 = 50L, SM_3 = 69L,
      SM_4 = 90L, SM_5 = 53L, SM_6 = 292L
    )
    observed_umap_sizes <- table(as.character(coordinates$module))
    coordinate_values <- suppressWarnings(as.matrix(data.frame(
      kME = as.numeric(coordinates$kME),
      UMAP1 = as.numeric(coordinates$UMAP1),
      UMAP2 = as.numeric(coordinates$UMAP2)
    )))
    parameter_values <- setNames(
      as.character(umap_parameters$value),
      as.character(umap_parameters$parameter)
    )
    required_parameters <- c(
      "n_embedded_genes", "gene_labels", "hub_circles",
      "ME_feature_type", "ME_reduction"
    )
    umap_ok <- nrow(coordinates) == 2206L &&
      !anyDuplicated(coordinates$gene) &&
      all(is.finite(coordinate_values)) &&
      setequal(names(observed_umap_sizes), modules) &&
      identical(
        as.integer(observed_umap_sizes[modules]),
        as.integer(expected_umap_sizes[modules])
      ) &&
      !anyDuplicated(umap_parameters$parameter) &&
      all(required_parameters %in% names(parameter_values)) &&
      identical(parameter_values[["n_embedded_genes"]], "2206") &&
      identical(parameter_values[["gene_labels"]], "none") &&
      identical(parameter_values[["hub_circles"]], "none") &&
      identical(
        parameter_values[["ME_feature_type"]],
        "harmonized module eigengenes (hMEs)"
      ) &&
      identical(parameter_values[["ME_reduction"]], "umap_SCT")
    umap_detail <- sprintf(
      "%d embedded genes; no gene labels or hub circles; hMEs on umap_SCT",
      nrow(coordinates)
    )
  }
}
record_check(
  "Supplementary module-gene and harmonized-ME UMAP structure",
  umap_ok,
  umap_detail
)

# -------------------------------------------------------------------------
# Supplementary donor-pseudobulk GSEA integrity checks
# -------------------------------------------------------------------------

gsea_sample_file <- file.path(gsea_logs_dir, "sample_metadata.tsv")
gsea_sample_result <- safe_table(gsea_sample_file, read_tsv)
gsea_sample_ok <- FALSE
gsea_sample_detail <- gsea_sample_result$detail
if (gsea_sample_result$ok) {
  gsea_samples <- gsea_sample_result$data
  required_columns <- c("donor", "stage", "n_spots", "library_size")
  missing_columns <- setdiff(required_columns, colnames(gsea_samples))
  if (length(missing_columns)) {
    gsea_sample_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    expected_stage_counts <- c(F0 = 6L, F1 = 6L, F2 = 3L, F3 = 3L, F4 = 6L)
    observed_stage_counts <- table(as.character(gsea_samples$stage))
    n_spots <- suppressWarnings(as.numeric(gsea_samples$n_spots))
    library_size <- suppressWarnings(as.numeric(gsea_samples$library_size))
    gsea_sample_ok <- nrow(gsea_samples) == 24L &&
      !anyDuplicated(gsea_samples$donor) &&
      setequal(names(observed_stage_counts), names(expected_stage_counts)) &&
      identical(
        as.integer(observed_stage_counts[names(expected_stage_counts)]),
        as.integer(expected_stage_counts)
      ) &&
      all(is.finite(n_spots) & n_spots > 0) &&
      sum(n_spots) == 26056L &&
      all(is.finite(library_size) & library_size > 0)
    gsea_sample_detail <- sprintf(
      "%d donors and %d spots; stage replication F0/F1/F2/F3/F4 = %s",
      nrow(gsea_samples),
      sum(n_spots),
      paste(as.integer(observed_stage_counts[names(expected_stage_counts)]), collapse = "/")
    )
  }
}
record_check(
  "GSEA donor-level sample structure",
  gsea_sample_ok,
  gsea_sample_detail
)

ranking_manifest_file <- file.path(
  gsea_results_dir,
  "manifests",
  "ranking_manifest.tsv"
)
ranking_manifest_result <- safe_table(ranking_manifest_file, read_tsv)
ranking_manifest_ok <- FALSE
ranking_manifest_detail <- ranking_manifest_result$detail
if (ranking_manifest_result$ok) {
  ranking_manifest <- ranking_manifest_result$data
  required_columns <- c(
    "comparison",
    "expression_filtered_genes",
    "ranked_unique_symbols",
    "ranking_statistic",
    "DEG_significance_filter_applied"
  )
  missing_columns <- setdiff(required_columns, colnames(ranking_manifest))
  if (length(missing_columns)) {
    ranking_manifest_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    filtered_n <- suppressWarnings(
      as.numeric(ranking_manifest$expression_filtered_genes)
    )
    ranked_n <- suppressWarnings(
      as.numeric(ranking_manifest$ranked_unique_symbols)
    )
    deg_filter <- as.logical(ranking_manifest$DEG_significance_filter_applied)
    ranking_manifest_ok <- nrow(ranking_manifest) == 4L &&
      !anyDuplicated(ranking_manifest$comparison) &&
      setequal(ranking_manifest$comparison, comparisons) &&
      all(is.finite(filtered_n) & filtered_n > 1000) &&
      all(is.finite(ranked_n) & ranked_n > 1000 & ranked_n <= filtered_n) &&
      all(ranking_manifest$ranking_statistic == "limma moderated t") &&
      !anyNA(deg_filter) && !any(deg_filter)
    ranking_manifest_detail <- paste0(
      "four complete moderated-t rankings; DEG significance filtering before GSEA: ",
      if (all(!deg_filter, na.rm = TRUE)) "none" else "detected"
    )
  }
}
record_check(
  "GSEA whole-list ranking specification",
  ranking_manifest_ok,
  ranking_manifest_detail
)

gsea_package_file <- file.path(gsea_logs_dir, "package_versions.tsv")
gsea_package_result <- safe_table(gsea_package_file, read_tsv)
gsea_package_ok <- FALSE
gsea_package_detail <- gsea_package_result$detail
if (gsea_package_result$ok) {
  gsea_packages <- gsea_package_result$data
  required_columns <- c("package", "version")
  required_components <- c(
    "R", "R_platform", "Seurat", "SeuratObject", "Matrix", "edgeR",
    "limma", "statmod", "msigdbr", "clusterProfiler", "fgsea",
    "AnnotationDbi", "org.Hs.eg.db", "ggplot2", "stringr", "digest"
  )
  missing_columns <- setdiff(required_columns, colnames(gsea_packages))
  if (length(missing_columns)) {
    gsea_package_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    gsea_package_ok <- !anyDuplicated(gsea_packages$package) &&
      all(required_components %in% gsea_packages$package) &&
      all(!is.na(gsea_packages$version) & nzchar(gsea_packages$version))
    gsea_package_detail <- sprintf(
      "%d R/platform/package versions recorded",
      nrow(gsea_packages)
    )
  }
}
record_check(
  "GSEA package and platform version record",
  gsea_package_ok,
  gsea_package_detail
)

msigdb_manifest_file <- file.path(
  gsea_results_dir,
  "manifests",
  "msigdb_collection_manifest.tsv"
)
msigdb_manifest_result <- safe_table(msigdb_manifest_file, read_tsv)
msigdb_manifest_ok <- FALSE
msigdb_manifest_detail <- msigdb_manifest_result$detail
if (msigdb_manifest_result$ok) {
  msigdb_manifest <- msigdb_manifest_result$data
  required_columns <- c(
    "collection", "n_gene_sets", "n_term_gene_pairs",
    "msigdbr_package_version", "msigdb_version"
  )
  missing_columns <- setdiff(required_columns, colnames(msigdb_manifest))
  if (length(missing_columns)) {
    msigdb_manifest_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    n_gene_sets <- suppressWarnings(as.numeric(msigdb_manifest$n_gene_sets))
    n_pairs <- suppressWarnings(as.numeric(msigdb_manifest$n_term_gene_pairs))
    msigdb_manifest_ok <- nrow(msigdb_manifest) == 4L &&
      !anyDuplicated(msigdb_manifest$collection) &&
      setequal(msigdb_manifest$collection, gsea_collections) &&
      all(is.finite(n_gene_sets) & n_gene_sets > 0) &&
      all(is.finite(n_pairs) & n_pairs > 0) &&
      all(n_pairs >= n_gene_sets) &&
      all(nzchar(msigdb_manifest$msigdbr_package_version)) &&
      all(nzchar(msigdb_manifest$msigdb_version))
    msigdb_manifest_detail <- sprintf(
      "MSigDB H/C2/C3/C5 recorded; database version: %s",
      paste(unique(msigdb_manifest$msigdb_version), collapse = ",")
    )
  }
}
record_check(
  "GSEA MSigDB collection and database version record",
  msigdb_manifest_ok,
  msigdb_manifest_detail
)

gsea_counts_file <- file.path(gsea_results_dir, "GSEA_significance_counts.tsv")
gsea_counts_result <- safe_table(gsea_counts_file, read_tsv)
gsea_counts_ok <- FALSE
gsea_counts_detail <- gsea_counts_result$detail
if (gsea_counts_result$ok) {
  gsea_counts <- gsea_counts_result$data
  required_columns <- c(
    "comparison", "collection", "tested_gene_sets",
    "significant_FDR_lt_0.05", "activated", "suppressed"
  )
  missing_columns <- setdiff(required_columns, colnames(gsea_counts))
  if (length(missing_columns)) {
    gsea_counts_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    keys <- paste(gsea_counts$comparison, gsea_counts$collection, sep = "::")
    expected_keys <- as.vector(
      outer(comparisons, gsea_collections, paste, sep = "::")
    )
    tested <- suppressWarnings(as.numeric(gsea_counts$tested_gene_sets))
    significant <- suppressWarnings(
      as.numeric(gsea_counts$significant_FDR_lt_0.05)
    )
    activated <- suppressWarnings(as.numeric(gsea_counts$activated))
    suppressed <- suppressWarnings(as.numeric(gsea_counts$suppressed))
    gsea_counts_ok <- nrow(gsea_counts) == 16L &&
      !anyDuplicated(keys) && setequal(keys, expected_keys) &&
      all(is.finite(tested) & tested > 0) &&
      all(is.finite(significant) & significant >= 0 & significant <= tested) &&
      all(is.finite(activated) & activated >= 0) &&
      all(is.finite(suppressed) & suppressed >= 0) &&
      all(significant == activated + suppressed)
    gsea_counts_detail <- sprintf(
      "%d comparison/collection tests; %d significant gene sets",
      nrow(gsea_counts),
      sum(significant)
    )
  }
}
record_check(
  "GSEA comparison and collection result structure",
  gsea_counts_ok,
  gsea_counts_detail
)

hallmark_significant_file <- file.path(
  gsea_results_dir,
  "gsea_significant",
  "all_comparisons_all_collections_GSEA_FDR_lt_0.05.tsv.gz"
)
hallmark_result <- safe_table(hallmark_significant_file, read_tsv)
hallmark_ok <- FALSE
hallmark_detail <- hallmark_result$detail
if (hallmark_result$ok && gsea_counts_result$ok) {
  all_significant <- hallmark_result$data
  required_columns <- c(
    "comparison", "collection", "Description", "NES", "p.adjust",
    "setSize", "direction"
  )
  missing_columns <- setdiff(required_columns, colnames(all_significant))
  if (length(missing_columns)) {
    hallmark_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else {
    hallmarks <- all_significant[all_significant$collection == "H", , drop = FALSE]
    hallmark_expected <- sum(
      gsea_counts$significant_FDR_lt_0.05[gsea_counts$collection == "H"]
    )
    nes <- suppressWarnings(as.numeric(hallmarks$NES))
    fdr <- suppressWarnings(as.numeric(hallmarks$p.adjust))
    expected_direction <- ifelse(nes > 0, "Activated", "Suppressed")
    hallmark_ok <- nrow(hallmarks) == hallmark_expected &&
      all(hallmarks$comparison %in% comparisons) &&
      all(is.finite(nes) & nes != 0) &&
      all(is.finite(fdr) & fdr >= 0 & fdr < 0.05) &&
      identical(as.character(hallmarks$direction), expected_direction)
    hallmark_detail <- sprintf(
      "%d significant Hallmark results are valid and represented without a top-N filter",
      nrow(hallmarks)
    )
  }
}
record_check(
  "Significant Hallmark table and plot source",
  hallmark_ok,
  hallmark_detail
)

run_parameters_file <- file.path(gsea_logs_dir, "run_parameters.txt")
run_parameter_lines <- if (file.exists(run_parameters_file)) {
  readLines(run_parameters_file, warn = FALSE)
} else {
  character(0)
}
required_parameter_lines <- c(
  "seed_base=260814",
  "computational_threads=1",
  "independent_biological_unit=donor",
  "de_method=donor_pseudobulk_limma_voom",
  "ranking_statistic=limma_moderated_t",
  "ranking_gene_universe=all_expression_filtered_genes",
  "DEG_significance_filter_before_GSEA=none",
  "raw_count_assay=Spatial",
  "SCT_values_or_models_used=FALSE",
  "filter=edgeR_filterByExpr",
  "normalization=edgeR_TMM",
  "clusterProfiler_GSEA_nproc_argument=not_overridden",
  "Hallmark_plot_selection=all_significant_Hallmark_gene_sets"
)
run_parameters_ok <- length(run_parameter_lines) &&
  all(required_parameter_lines %in% run_parameter_lines) &&
  paste0("run_name=", gsea_run$run_name) %in% run_parameter_lines
record_check(
  "GSEA recorded analysis specification",
  run_parameters_ok,
  if (run_parameters_ok) {
    paste(
      "donor pseudobulk, raw counts, moderated-t whole-list ranking",
      "and all-significant Hallmark plotting recorded"
    )
  } else {
    "one or more required GSEA parameter records are absent"
  }
)

result_checksum_file <- file.path(gsea_logs_dir, "result_SHA256SUMS.tsv")
result_checksum_result <- safe_table(result_checksum_file, read_tsv)
result_checksum_ok <- FALSE
result_checksum_detail <- result_checksum_result$detail
if (result_checksum_result$ok) {
  checksum_manifest <- result_checksum_result$data
  required_columns <- c("sha256", "path")
  missing_columns <- setdiff(required_columns, colnames(checksum_manifest))
  unsafe_paths <- if ("path" %in% colnames(checksum_manifest)) {
    grepl("^/|(^|/)\\.\\.(/|$)", checksum_manifest$path)
  } else {
    logical(0)
  }
  if (length(missing_columns)) {
    result_checksum_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else if (any(unsafe_paths)) {
    result_checksum_detail <- "checksum manifest contains an unsafe relative path"
  } else {
    listed_files <- file.path(gsea_results_dir, checksum_manifest$path)
    actual_result_files <- list.files(
      gsea_results_dir,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )
    actual_result_files <- actual_result_files[
      file.info(actual_result_files)$isdir %in% FALSE
    ]
    actual_relative_paths <- substring(
      actual_result_files,
      nchar(gsea_results_dir) + 2L
    )
    files_ok <- all(nonempty_files(listed_files))
    actual_sha256 <- if (files_ok) {
      vapply(listed_files, sha256_file, character(1))
    } else {
      rep(NA_character_, length(listed_files))
    }
    result_checksum_ok <- !anyDuplicated(checksum_manifest$path) &&
      setequal(checksum_manifest$path, actual_relative_paths) &&
      files_ok &&
      identical(
        tolower(as.character(checksum_manifest$sha256)),
        tolower(unname(actual_sha256))
      )
    result_checksum_detail <- sprintf(
      "%d generated GSEA result checksums verified",
      nrow(checksum_manifest)
    )
  }
}
record_check(
  "GSEA generated-result SHA-256 verification",
  result_checksum_ok,
  result_checksum_detail
)

execution_checksum_file <- file.path(
  gsea_logs_dir,
  "execution_record_SHA256SUMS.tsv"
)
execution_checksum_result <- safe_table(execution_checksum_file, read_tsv)
execution_checksum_ok <- FALSE
execution_checksum_detail <- execution_checksum_result$detail
if (execution_checksum_result$ok) {
  execution_manifest <- execution_checksum_result$data
  required_columns <- c("sha256", "path")
  missing_columns <- setdiff(required_columns, colnames(execution_manifest))
  unsafe_paths <- if ("path" %in% colnames(execution_manifest)) {
    grepl("^/|(^|/)\\.\\.(/|$)", execution_manifest$path)
  } else {
    logical(0)
  }
  if (length(missing_columns)) {
    execution_checksum_detail <- paste(
      "missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  } else if (any(unsafe_paths)) {
    execution_checksum_detail <- "execution manifest contains an unsafe relative path"
  } else {
    listed_files <- file.path(gsea_logs_dir, execution_manifest$path)
    actual_log_files <- list.files(
      gsea_logs_dir,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )
    actual_log_files <- actual_log_files[
      file.info(actual_log_files)$isdir %in% FALSE
    ]
    actual_relative_paths <- substring(
      actual_log_files,
      nchar(gsea_logs_dir) + 2L
    )
    expected_manifest_paths <- setdiff(
      actual_relative_paths,
      "execution_record_SHA256SUMS.tsv"
    )
    files_ok <- all(nonempty_files(listed_files))
    actual_sha256 <- if (files_ok) {
      vapply(listed_files, sha256_file, character(1))
    } else {
      rep(NA_character_, length(listed_files))
    }
    execution_checksum_ok <- !anyDuplicated(execution_manifest$path) &&
      setequal(execution_manifest$path, expected_manifest_paths) &&
      files_ok &&
      identical(
        tolower(as.character(execution_manifest$sha256)),
        tolower(unname(actual_sha256))
      )
    execution_checksum_detail <- sprintf(
      "%d execution-record checksums verified",
      nrow(execution_manifest)
    )
  }
}
record_check(
  "GSEA execution-record SHA-256 verification",
  execution_checksum_ok,
  execution_checksum_detail
)

# -------------------------------------------------------------------------
# Incomplete-output detection
# -------------------------------------------------------------------------

partial_roots <- c(
  objects_dir,
  results_dir,
  logs_dir,
  gsea_results_root,
  gsea_logs_root
)
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
  if (!length(partial_paths)) {
    "no .partial files or directories detected"
  } else {
    paste(partial_paths, collapse = "; ")
  }
)

# -------------------------------------------------------------------------
# Checksums, summary and final decision
# -------------------------------------------------------------------------

recursive_output_files <- unlist(lapply(
  c(
    file.path(results_dir, "02_network"),
    file.path(results_dir, "03_DME_donor_level_2026-10-13"),
    enrichment_dir,
    hub_network_dir,
    module_expression_dir,
    supplementary_umap_dir
  ),
  function(path) {
    if (!dir.exists(path)) return(character(0))
    list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  }
), use.names = FALSE)

gsea_selected_run_files <- unlist(lapply(
  c(gsea_results_dir, gsea_logs_dir),
  function(path) {
    if (!dir.exists(path)) return(character(0))
    files <- list.files(
      path,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE
    )
    files[file.info(files)$isdir %in% FALSE]
  }
), use.names = FALSE)

final_output_files <- unique(c(
  prepared_object,
  network_object,
  network_log_files,
  file.path(logs_dir, c("03_DME_donor_level_2026-10-13_summary.tsv", "03_DME_donor_level_2026-10-13_sessionInfo.txt")),
  enrichment_checkpoint_files,
  recursive_output_files,
  gsea_selected_run_files
))

say("Calculating final Figure 2 output checksums; large RDS files may take time")
output_checksums <- checksum_table(final_output_files, "final_Fig2_workflow_output")
utils::write.table(
  output_checksums,
  file.path(audit_dir, "final_output_checksums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

component_status <- do.call(rbind, status_rows)
rownames(component_status) <- NULL
audit_pass <- all(component_status$pass[component_status$required])

utils::write.table(
  component_status,
  file.path(audit_dir, "Fig2_component_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

summary_lines <- c(
  "FLINC Figure 2 final reproducibility audit",
  paste0("Timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Workflow root: ", repro_root),
  paste0(
    "Supplementary GSEA run: ",
    if (gsea_run$ok) gsea_run$run_name else "not resolved"
  ),
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
  paste0("FINAL_FIGURE2_AUDIT_PASS = ", audit_pass),
  "",
  paste0(
    "Code/reference checksum table: ",
    file.path(audit_dir, "code_documentation_and_reference_checksums.tsv")
  ),
  paste0("Output checksum table: ", file.path(audit_dir, "final_output_checksums.tsv")),
  "",
  paste0(
    "Scope note: this audit checks output completeness and internal consistency; ",
    "it performs no additional biological analysis."
  )
)

writeLines(summary_lines, file.path(audit_dir, "audit_summary.txt"), useBytes = TRUE)
writeLines(
  capture.output(sessionInfo()),
  file.path(audit_dir, "sessionInfo.txt"),
  useBytes = TRUE
)

cat(paste(summary_lines, collapse = "\n"), "\n")

if (!audit_pass) {
  stopf("Figure 2 audit did not fully pass. Inspect: %s", audit_dir)
}

cat("\nFINAL FIGURE 2 AUDIT PASSED.\n")

