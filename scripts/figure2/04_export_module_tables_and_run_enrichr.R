#!/usr/bin/env Rscript

# FLINC Figure 2 — export deterministic module and hub-gene tables and run
# EnrichR on every gene assigned to each non-grey module.
#
# Restart safety:
#   * EnrichR responses are checkpointed atomically, one module at a time.
#   * A retry validates and reuses compatible checkpoints.
#   * Final outputs are built in a sibling .partial directory and renamed into
#     place only after all tables and plots have been generated successfully.
#   * Existing final or partial output directories are never overwritten.

options(stringsAsFactors = FALSE, warn = 1)

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

move_into_place <- function(partial, path) {
  if (file.exists(path) || dir.exists(path)) {
    stopf("Refusing to overwrite existing output: %s", path)
  }
  if (!file.rename(partial, path)) {
    stopf("Could not move partial output into place: %s", path)
  }
}

save_rds_atomic <- function(object, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path)) {
    stopf("Refusing to overwrite existing checkpoint: %s", path)
  }
  if (file.exists(partial)) {
    stopf("Partial checkpoint exists; inspect before retrying: %s", partial)
  }
  saveRDS(object, partial, compress = "gzip", version = 3)
  move_into_place(partial, path)
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
}

write_text <- function(x, path) {
  writeLines(x, path, useBytes = TRUE)
}

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

clean_term <- function(x) {
  x <- gsub("\\s*\\((GO:?\\s*[0-9]+|GO:\\s*[0-9]+)\\)\\s*$", "", x, perl = TRUE)
  x <- gsub("\\s*WP[0-9]+\\s*$", "", x, perl = TRUE)
  x <- gsub("\\s*\\([^\\)]+\\)\\s*$", "", x, perl = TRUE)
  trimws(gsub("[[:space:]]+", " ", x, perl = TRUE))
}

rbind_fill <- function(x) {
  if (length(x) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(x, colnames), use.names = FALSE))
  x <- lapply(x, function(one) {
    missing_names <- setdiff(all_names, colnames(one))
    for (name in missing_names) {
      one[[name]] <- rep(NA, nrow(one))
    }
    one[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

save_pdf <- function(plot, path, width = 8, height = 6) {
  grDevices::pdf(
    file = path,
    width = width,
    height = height,
    useDingbats = FALSE,
    onefile = TRUE
  )
  tryCatch(
    print(plot),
    finally = grDevices::dev.off()
  )
}

make_module_barplot <- function(
  table,
  module,
  database,
  fill_color,
  n_terms = 10L,
  p_cutoff = 0.05
) {
  required <- c("Term", "Adjusted.P.value", "Combined.Score")
  if (!all(required %in% colnames(table))) {
    return(list(plot = NULL, selected = table[0, , drop = FALSE]))
  }

  adjusted_p <- suppressWarnings(as.numeric(table[["Adjusted.P.value"]]))
  combined_score <- suppressWarnings(as.numeric(table[["Combined.Score"]]))
  keep <- !is.na(adjusted_p) & adjusted_p <= p_cutoff &
    !is.na(combined_score) & is.finite(combined_score) & combined_score > 0
  selected <- table[keep, , drop = FALSE]

  if (nrow(selected) == 0L) {
    return(list(plot = NULL, selected = selected))
  }

  selected <- selected[
    order(
      -as.numeric(selected[["Combined.Score"]]),
      as.numeric(selected[["Adjusted.P.value"]]),
      as.character(selected[["Term"]])
    ),
    ,
    drop = FALSE
  ]
  selected <- head(selected, n_terms)
  selected$Term_clean <- clean_term(as.character(selected$Term))
  selected$plot_score <- log(as.numeric(selected[["Combined.Score"]]))
  selected$plot_rank <- seq_len(nrow(selected))

  plot_data <- selected
  plot_data$plot_label <- factor(
    plot_data$Term_clean,
    levels = rev(unique(plot_data$Term_clean))
  )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = plot_label, y = plot_score)
  ) +
    ggplot2::geom_col(
      fill = fill_color,
      color = "black",
      linewidth = 0.35,
      width = 0.85
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0(module, " - ", database),
      x = NULL,
      y = "Enrichment log(Combined.Score)"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(color = "black", face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black")
    )

  list(plot = plot, selected = selected)
}

required_packages <- c("enrichR", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stopf(
    paste0(
      "Missing required R package(s): %s. Install them in ",
      "hdWGCNA_pipeline-env before running Step 05."
    ),
    paste(missing_packages, collapse = ", ")
  )
}

# enrichR 3.4 initializes these options in .onAttach(), not .onLoad().
# This script deliberately uses namespace-qualified calls rather than attaching
# packages, so set the human Enrichr endpoint explicitly and reproducibly.
enrichr_base_address <- "https://maayanlab.cloud/Enrichr/"
options(
  enrichR.sites.base.address = "https://maayanlab.cloud/",
  enrichR.base.address = enrichr_base_address,
  speedrichr.base.address = "https://maayanlab.cloud/speedrichr/api/",
  enrichR.live = TRUE,
  enrichR.quiet = FALSE
)
if (!identical(getOption("enrichR.base.address"), enrichr_base_address)) {
  stopf("Could not initialize the human Enrichr endpoint")
}

observed_r_version <- paste(R.version$major, R.version$minor, sep = ".")
if (!identical(observed_r_version, "4.4.1")) {
  stopf("R version mismatch: expected 4.4.1, observed %s", R.version.string)
}

repro_root <- env_or(
  "FLINC_REPRO_ROOT",
  "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)
module_table_file <- env_or(
  "FIG2_MODULE_TABLE_INPUT",
  file.path(
    repro_root,
    "results",
    "Fig2_WGCNA",
    "02_network",
    "02_module_assignments_and_kME.tsv"
  )
)
output_dir <- env_or(
  "FIG2_STEP05_OUTPUT_DIR",
  file.path(
    repro_root,
    "results",
    "Fig2_WGCNA",
    "05_module_tables_and_enrichr"
  )
)
checkpoint_dir <- env_or(
  "FIG2_STEP05_CHECKPOINT_DIR",
  file.path(
    repro_root,
    "results",
    "Fig2_WGCNA",
    "05_enrichr_checkpoints"
  )
)
partial_output_dir <- paste0(output_dir, ".partial")
sleep_seconds <- suppressWarnings(as.numeric(env_or("FIG2_ENRICHR_SLEEP_SECONDS", "1")))
if (length(sleep_seconds) != 1L || is.na(sleep_seconds) || sleep_seconds < 0) {
  stopf("FIG2_ENRICHR_SLEEP_SECONDS must be one non-negative number")
}

if (!file.exists(module_table_file)) {
  stopf("Step-02 module table does not exist: %s", module_table_file)
}
if (dir.exists(output_dir) || file.exists(output_dir)) {
  stopf("Final Step-05 output already exists; do not rerun: %s", output_dir)
}
if (dir.exists(partial_output_dir) || file.exists(partial_output_dir)) {
  stopf(
    "Partial Step-05 output exists; inspect before retrying: %s",
    partial_output_dir
  )
}

dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(checkpoint_dir)) {
  stopf("Could not create checkpoint directory: %s", checkpoint_dir)
}

module_metadata <- data.frame(
  module = paste0("SM_", 1:6),
  biological_name = c(
    "Fibro-immune ECM remodeling",
    "Interferon antiviral",
    "Cholesterol metabolism and lipogenesis",
    "Hepatocyte homeostatic function",
    "Stress and protein homeostasis",
    "Hepatic plasma"
  ),
  module_color = c(
    "#CAB2D6", "#6A3D9A", "#E69F00",
    "#7EC0EE", "#993F00", "#023FA5"
  ),
  expected_n_genes = c(1652L, 50L, 69L, 90L, 53L, 292L),
  displayed_hub_n = c(50L, rep(30L, 5L)),
  stringsAsFactors = FALSE
)
module_metadata$top70_n <- ceiling(module_metadata$expected_n_genes * 0.7)

say("Reading module assignments and harmonized kME values")
modules <- read_tsv(module_table_file)
required_columns <- c(
  "gene_name", "module", "color", paste0("kME_SM_", 1:6)
)
missing_columns <- setdiff(required_columns, colnames(modules))
if (length(missing_columns) > 0L) {
  stopf(
    "Step-02 module table is missing: %s",
    paste(missing_columns, collapse = ", ")
  )
}
if (nrow(modules) != 11515L) {
  stopf("Expected 11,515 WGCNA genes; observed %d", nrow(modules))
}
if (anyNA(modules$gene_name) || any(!nzchar(as.character(modules$gene_name)))) {
  stopf("Module table contains missing or empty gene names")
}
if (anyDuplicated(modules$gene_name)) {
  stopf("Module table contains duplicated gene names")
}
modules$gene_name <- as.character(modules$gene_name)
modules$module <- as.character(modules$module)
modules$color <- as.character(modules$color)

expected_modules <- c(module_metadata$module, "grey")
unexpected_modules <- setdiff(unique(modules$module), expected_modules)
missing_modules <- setdiff(expected_modules, unique(modules$module))
if (length(unexpected_modules) > 0L || length(missing_modules) > 0L) {
  stopf(
    "Unexpected module set. Unexpected: [%s]; missing: [%s]",
    paste(unexpected_modules, collapse = ", "),
    paste(missing_modules, collapse = ", ")
  )
}

observed_sizes <- table(modules$module)
for (i in seq_len(nrow(module_metadata))) {
  module_name <- module_metadata$module[i]
  observed_size <- unname(observed_sizes[module_name])
  expected_size <- module_metadata$expected_n_genes[i]
  if (!identical(as.integer(observed_size), as.integer(expected_size))) {
    stopf(
      "Module-size mismatch for %s: expected %d, observed %d",
      module_name,
      expected_size,
      observed_size
    )
  }
  observed_colors <- unique(modules$color[modules$module == module_name])
  if (!identical(observed_colors, module_metadata$module_color[i])) {
    stopf(
      "Module-colour mismatch for %s: expected %s, observed [%s]",
      module_name,
      module_metadata$module_color[i],
      paste(observed_colors, collapse = ", ")
    )
  }
}

ranked_parts <- vector("list", nrow(module_metadata))
genes_by_module <- setNames(vector("list", nrow(module_metadata)), module_metadata$module)
for (i in seq_len(nrow(module_metadata))) {
  module_name <- module_metadata$module[i]
  kme_column <- paste0("kME_", module_name)
  part <- modules[modules$module == module_name, , drop = FALSE]
  kme <- suppressWarnings(as.numeric(part[[kme_column]]))
  if (anyNA(kme) || any(!is.finite(kme))) {
    stopf("Missing or non-finite %s values within %s", kme_column, module_name)
  }
  part$module_kME <- kme
  part <- part[order(-part$module_kME, part$gene_name), , drop = FALSE]
  part$module_rank <- seq_len(nrow(part))
  part$module_size <- nrow(part)
  part$biological_name <- module_metadata$biological_name[i]
  part$selected_top70 <- part$module_rank <= module_metadata$top70_n[i]
  part$selected_displayed_hub <-
    part$module_rank <= module_metadata$displayed_hub_n[i]
  ranked_parts[[i]] <- part
  genes_by_module[[module_name]] <- part$gene_name
}
ranked_non_grey <- do.call(rbind, ranked_parts)
rownames(ranked_non_grey) <- NULL
top70 <- ranked_non_grey[ranked_non_grey$selected_top70, , drop = FALSE]
displayed_hubs <- ranked_non_grey[
  ranked_non_grey$selected_displayed_hub,
  ,
  drop = FALSE
]

if (nrow(top70) != sum(module_metadata$top70_n)) {
  stopf("Top-70%% export has an unexpected number of rows")
}
if (nrow(displayed_hubs) != sum(module_metadata$displayed_hub_n)) {
  stopf("Displayed-hub export has an unexpected number of rows")
}

biological_names <- setNames(
  module_metadata$biological_name,
  module_metadata$module
)
modules$biological_name <- unname(biological_names[modules$module])
modules$biological_name[modules$module == "grey"] <- NA_character_

enrichr_databases <- c(
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023",
  "KEGG_2019_Human",
  "WikiPathways_2024_Human"
)

say("Checking availability of the five EnrichR databases")
database_catalog <- tryCatch(
  enrichR::listEnrichrDbs(),
  error = function(e) {
    stopf("Could not retrieve the EnrichR database catalog: %s", conditionMessage(e))
  }
)
if (!is.data.frame(database_catalog) || !("libraryName" %in% colnames(database_catalog))) {
  stopf("EnrichR database catalog does not contain libraryName")
}
unavailable_databases <- setdiff(
  enrichr_databases,
  as.character(database_catalog$libraryName)
)
if (length(unavailable_databases) > 0L) {
  stopf(
    "Required EnrichR database(s) are unavailable: %s",
    paste(unavailable_databases, collapse = ", ")
  )
}

say("Running or resuming EnrichR queries")
checkpoints <- setNames(vector("list", nrow(module_metadata)), module_metadata$module)
for (module_name in module_metadata$module) {
  genes <- genes_by_module[[module_name]]
  checkpoint_file <- file.path(
    checkpoint_dir,
    paste0("05_enrichr_", module_name, ".rds")
  )
  checkpoint_partial <- paste0(checkpoint_file, ".partial")
  if (file.exists(checkpoint_partial)) {
    stopf(
      "Partial EnrichR checkpoint exists; inspect before retrying: %s",
      checkpoint_partial
    )
  }

  if (file.exists(checkpoint_file)) {
    say("Validating and reusing checkpoint for", module_name)
    checkpoint <- readRDS(checkpoint_file)
    required_checkpoint_fields <- c(
      "module", "genes", "databases", "enrichr_base_address",
      "queried_at_utc", "results"
    )
    if (!is.list(checkpoint) ||
        !all(required_checkpoint_fields %in% names(checkpoint)) ||
        !identical(checkpoint$module, module_name) ||
        !identical(checkpoint$genes, genes) ||
        !identical(checkpoint$databases, enrichr_databases) ||
        !identical(checkpoint$enrichr_base_address, enrichr_base_address)) {
      stopf("Existing checkpoint is incompatible: %s", checkpoint_file)
    }
    result_names <- names(checkpoint$results)
    if (!setequal(result_names, enrichr_databases)) {
      stopf("Checkpoint has an unexpected database set: %s", checkpoint_file)
    }
    checkpoint$results <- checkpoint$results[enrichr_databases]
  } else {
    say("Querying", module_name, "with", length(genes), "genes")
    results <- tryCatch(
      enrichR::enrichr(genes, enrichr_databases),
      error = function(e) {
        stopf("EnrichR query failed for %s: %s", module_name, conditionMessage(e))
      }
    )
    if (!is.list(results) || !setequal(names(results), enrichr_databases)) {
      stopf("EnrichR returned an unexpected database set for %s", module_name)
    }
    results <- results[enrichr_databases]
    checkpoint <- list(
      module = module_name,
      genes = genes,
      databases = enrichr_databases,
      enrichr_base_address = enrichr_base_address,
      queried_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      R_version = R.version.string,
      enrichR_version = as.character(utils::packageVersion("enrichR")),
      results = results
    )
    save_rds_atomic(checkpoint, checkpoint_file)
    say("Checkpoint saved:", checkpoint_file)
    if (sleep_seconds > 0) {
      Sys.sleep(sleep_seconds)
    }
  }
  checkpoints[[module_name]] <- checkpoint
}

say("All EnrichR responses are available; finalizing Step-05 outputs")
dir.create(partial_output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(partial_output_dir)) {
  stopf("Could not create partial output directory: %s", partial_output_dir)
}

table_dir <- file.path(partial_output_dir, "tables")
gene_list_dir <- file.path(partial_output_dir, "gene_lists")
raw_dir <- file.path(partial_output_dir, "raw_enrichr")
plot_dir <- file.path(partial_output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_list_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(
  module_metadata,
  file.path(table_dir, "05_module_metadata.tsv")
)
write_tsv(
  modules,
  file.path(table_dir, "05_all_genes_assignments_and_kME.tsv")
)
write_tsv(
  ranked_non_grey,
  file.path(table_dir, "05_non_grey_genes_ranked_by_module_kME.tsv")
)
write_tsv(
  top70,
  file.path(table_dir, "05_top70_genes_ranked_by_module_kME.tsv")
)
write_tsv(
  displayed_hubs,
  file.path(
    table_dir,
    "05_displayed_hubs_top50_SM1_top30_SM2_to_SM6.tsv"
  )
)

for (module_name in module_metadata$module) {
  module_part <- ranked_non_grey[
    ranked_non_grey$module == module_name,
    ,
    drop = FALSE
  ]
  write_text(
    module_part$gene_name,
    file.path(gene_list_dir, paste0(module_name, "_all_genes.txt"))
  )
  write_text(
    module_part$gene_name[module_part$selected_top70],
    file.path(gene_list_dir, paste0(module_name, "_top70_genes.txt"))
  )
  write_text(
    module_part$gene_name[module_part$selected_displayed_hub],
    file.path(gene_list_dir, paste0(module_name, "_displayed_hubs.txt"))
  )
}

raw_tables <- list()
selected_tables <- list()
plot_index_rows <- list()
raw_index <- 0L
selected_index <- 0L
plot_index <- 0L

for (i in seq_len(nrow(module_metadata))) {
  module_name <- module_metadata$module[i]
  biological_name <- module_metadata$biological_name[i]
  module_color <- module_metadata$module_color[i]
  checkpoint <- checkpoints[[module_name]]
  module_plot_dir <- file.path(plot_dir, module_name)
  dir.create(module_plot_dir, recursive = TRUE, showWarnings = FALSE)

  for (database in enrichr_databases) {
    result <- checkpoint$results[[database]]
    if (is.null(result)) {
      result <- data.frame()
    }
    result <- as.data.frame(result, check.names = FALSE, stringsAsFactors = FALSE)
    result$module <- rep(module_name, nrow(result))
    result$biological_name <- rep(biological_name, nrow(result))
    result$database <- rep(database, nrow(result))
    result$n_query_genes <- rep(length(checkpoint$genes), nrow(result))
    result$queried_at_utc <- rep(checkpoint$queried_at_utc, nrow(result))

    raw_file <- file.path(
      raw_dir,
      paste0(module_name, "__", database, "__raw.tsv")
    )
    write_tsv(result, raw_file)

    raw_index <- raw_index + 1L
    raw_tables[[raw_index]] <- result

    plotted <- make_module_barplot(
      table = result,
      module = module_name,
      database = database,
      fill_color = module_color,
      n_terms = 10L,
      p_cutoff = 0.05
    )

    selected <- plotted$selected
    if (nrow(selected) > 0L) {
      selected$module <- module_name
      selected$biological_name <- biological_name
      selected$database <- database
      selected_index <- selected_index + 1L
      selected_tables[[selected_index]] <- selected
    }

    plot_relative <- NA_character_
    if (!is.null(plotted$plot)) {
      plot_file <- file.path(module_plot_dir, paste0(database, ".pdf"))
      save_pdf(plotted$plot, plot_file, width = 8, height = 6)
      plot_relative <- file.path("plots", module_name, paste0(database, ".pdf"))
    }

    adjusted_p <- if ("Adjusted.P.value" %in% colnames(result)) {
      suppressWarnings(as.numeric(result[["Adjusted.P.value"]]))
    } else {
      rep(NA_real_, nrow(result))
    }
    n_significant <- sum(!is.na(adjusted_p) & adjusted_p <= 0.05)

    plot_index <- plot_index + 1L
    plot_index_rows[[plot_index]] <- data.frame(
      module = module_name,
      biological_name = biological_name,
      database = database,
      n_query_genes = length(checkpoint$genes),
      n_raw_terms = nrow(result),
      n_adjusted_p_le_0.05 = n_significant,
      n_terms_plotted = nrow(selected),
      plot_file = plot_relative,
      queried_at_utc = checkpoint$queried_at_utc,
      stringsAsFactors = FALSE
    )
  }
}

combined_raw <- rbind_fill(raw_tables)
combined_selected <- rbind_fill(selected_tables)
if (ncol(combined_selected) == 0L) {
  combined_selected <- data.frame(
    module = character(0),
    biological_name = character(0),
    database = character(0),
    Term = character(0),
    Adjusted.P.value = numeric(0),
    Combined.Score = numeric(0),
    Term_clean = character(0),
    plot_score = numeric(0),
    plot_rank = integer(0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
plot_index_table <- do.call(rbind, plot_index_rows)
rownames(plot_index_table) <- NULL

write_tsv(
  combined_raw,
  file.path(table_dir, "05_enrichr_all_raw_results.tsv")
)
write_tsv(
  combined_selected,
  file.path(table_dir, "05_enrichr_significant_top10_for_plots.tsv")
)
write_tsv(
  plot_index_table,
  file.path(table_dir, "05_enrichr_query_and_plot_index.tsv")
)
write_tsv(
  database_catalog,
  file.path(table_dir, "05_enrichr_database_catalog_snapshot.tsv")
)

input_manifest <- data.frame(
  input = "module_table",
  path = module_table_file,
  md5 = unname(tools::md5sum(module_table_file)),
  stringsAsFactors = FALSE
)
write_tsv(
  input_manifest,
  file.path(partial_output_dir, "05_input_manifest.tsv")
)

query_times <- vapply(
  checkpoints[module_metadata$module],
  function(x) as.character(x$queried_at_utc),
  character(1)
)
summary_table <- rbind(
  data.frame(metric = "module_table_input", value = module_table_file),
  data.frame(metric = "R_version", value = observed_r_version),
  data.frame(metric = "enrichR_version", value = as.character(utils::packageVersion("enrichR"))),
  data.frame(metric = "ggplot2_version", value = as.character(utils::packageVersion("ggplot2"))),
  data.frame(metric = "n_wgcna_genes", value = nrow(modules)),
  data.frame(metric = "n_non_grey_genes", value = nrow(ranked_non_grey)),
  data.frame(metric = "n_top70_genes", value = nrow(top70)),
  data.frame(metric = "n_displayed_hub_genes", value = nrow(displayed_hubs)),
  data.frame(metric = "gene_selection_for_enrichr", value = "all non-grey module genes; no cap"),
  data.frame(metric = "gene_ranking", value = "descending signed module-specific kME; gene name tie-break"),
  data.frame(metric = "enrichr_base_address", value = enrichr_base_address),
  data.frame(metric = "enrichr_databases", value = paste(enrichr_databases, collapse = ",")),
  data.frame(metric = "adjusted_p_cutoff", value = "0.05"),
  data.frame(metric = "terms_per_plot", value = "10"),
  data.frame(metric = "plot_ranking", value = "descending Combined.Score"),
  data.frame(metric = "n_module_database_queries", value = nrow(plot_index_table)),
  data.frame(metric = "n_raw_enrichment_rows", value = nrow(combined_raw)),
  data.frame(metric = "n_top_plot_rows", value = nrow(combined_selected)),
  data.frame(
    metric = paste0("queried_at_utc_", names(query_times)),
    value = unname(query_times)
  )
)
write_tsv(
  summary_table,
  file.path(partial_output_dir, "05_run_summary.tsv")
)
write_text(
  capture.output(sessionInfo()),
  file.path(partial_output_dir, "05_sessionInfo.txt")
)

report_lines <- c(
  "FLINC Fig. 2 Step 05: deterministic module tables and EnrichR",
  paste0("Finalized: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  paste0("WGCNA genes: ", nrow(modules)),
  paste0("Non-grey genes queried: ", nrow(ranked_non_grey)),
  paste0("Top-70% table rows: ", nrow(top70)),
  paste0("Displayed hub-gene rows: ", nrow(displayed_hubs)),
  paste0("Module/database queries: ", nrow(plot_index_table)),
  paste0("Raw enrichment rows saved: ", nrow(combined_raw)),
  paste0("Top plot rows saved: ", nrow(combined_selected)),
  "",
  "EnrichR parameters:",
  "  Species/databases: human database set",
  paste0("  Endpoint: ", enrichr_base_address),
  paste0("  Databases: ", paste(enrichr_databases, collapse = ", ")),
  "  Gene list: every gene in each non-grey module (no cap)",
  "  Plot filter: Adjusted.P.value <= 0.05",
  "  Plot terms: top 10 by descending Combined.Score",
  "",
  "STEP05_MODULE_TABLES_PASS: TRUE",
  "STEP05_ENRICHR_PASS: TRUE"
)
write_text(
  report_lines,
  file.path(partial_output_dir, "05_completion_report.txt")
)

move_into_place(partial_output_dir, output_dir)

cat("\n=== STEP 05 SUMMARY ===\n")
print(module_metadata, row.names = FALSE)
cat("\nModule/database query index:\n")
print(plot_index_table, row.names = FALSE)
cat("\nSTEP05_MODULE_TABLES_PASS: TRUE\n")
cat("STEP05_ENRICHR_PASS: TRUE\n")
say("Final output:", output_dir)
say("EnrichR checkpoints:", checkpoint_dir)
say("FLINC Fig. 2 Step 05 completed successfully")
