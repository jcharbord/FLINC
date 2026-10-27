options(stringsAsFactors = FALSE, digits = 10, scipen = 999)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
set.seed(42L)

fail <- function(...) stop(paste0(...), call. = FALSE)

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table",
  "ggplot2", "cowplot", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  fail("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) fail("Arguments must be supplied as --key value pairs")
arg_list <- as.list(args[seq(2L, length(args), by = 2L)])
names(arg_list) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])

required_args <- c(
  "input", "module-table", "uniprot-snapshot",
  "historical-counts", "historical-topgenes",
  "output-dir", "log-dir",
  "expected-input-sha256", "expected-uniprot-sha256"
)
missing_args <- setdiff(required_args, names(arg_list))
if (length(missing_args) > 0L) {
  fail("Missing arguments: ", paste(missing_args, collapse = ", "))
}

input_path <- normalizePath(arg_list[["input"]], mustWork = TRUE)
module_path <- normalizePath(arg_list[["module-table"]], mustWork = TRUE)
uniprot_path <- normalizePath(arg_list[["uniprot-snapshot"]], mustWork = TRUE)
historical_counts_path <- normalizePath(arg_list[["historical-counts"]], mustWork = TRUE)
historical_topgenes_path <- normalizePath(arg_list[["historical-topgenes"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
expected_input_sha <- tolower(arg_list[["expected-input-sha256"]])
expected_uniprot_sha <- tolower(arg_list[["expected-uniprot-sha256"]])

if (dir.exists(output_dir) || file.exists(output_dir)) {
  fail("Final output path already exists: ", output_dir)
}
if (dir.exists(log_dir) || file.exists(log_dir)) {
  fail("Final log path already exists: ", log_dir)
}

dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_dir), recursive = TRUE, showWarnings = FALSE)

run_stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
result_tmp <- paste0(output_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
log_tmp <- paste0(log_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
dir.create(result_tmp, recursive = FALSE, showWarnings = FALSE)
dir.create(log_tmp, recursive = FALSE, showWarnings = FALSE)

completed <- FALSE
on.exit(
  if (!completed) {
    message("Incomplete run retained: ", result_tmp, "\n", log_tmp)
  },
  add = TRUE
)

atomic_write_table <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(
    x, temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

started_at <- Sys.time()

input_sha <- sha256(input_path)
uniprot_sha <- sha256(uniprot_path)
module_sha <- sha256(module_path)

if (!identical(tolower(input_sha), expected_input_sha)) {
  fail("Input SHA256 mismatch. Expected ", expected_input_sha, "; observed ", input_sha)
}
if (!identical(tolower(uniprot_sha), expected_uniprot_sha)) {
  fail("UniProt snapshot SHA256 mismatch. Expected ", expected_uniprot_sha, "; observed ", uniprot_sha)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

message("Reading regenerated WGCNA module assignments...")
module_table <- data.table::fread(module_path, data.table = FALSE)
if (!all(c("gene_name", "module") %in% names(module_table))) {
  fail("Module table lacks gene_name/module columns")
}

module_table$gene <- toupper(trimws(as.character(module_table$gene_name)))
module_table$module <- gsub("-", "_", trimws(as.character(module_table$module)))
module_table <- module_table[
  grepl("^SM_[0-9]+$", module_table$module) &
    !is.na(module_table$gene) &
    nzchar(module_table$gene),
  c("gene", "module"),
  drop = FALSE
]

if (anyDuplicated(module_table$gene)) {
  fail("Module table has duplicate gene assignments after normalization")
}

module_levels <- paste0(
  "SM_",
  sort(unique(as.integer(sub("^SM_", "", module_table$module))))
)
if (!identical(module_levels, paste0("SM_", 1:6))) {
  fail("Expected regenerated SM_1 through SM_6 modules; observed ",
       paste(module_levels, collapse = ", "))
}

message("Reading frozen historical UniProt snapshot...")
uniprot <- data.table::fread(uniprot_path, data.table = FALSE)
column_names_lower <- tolower(names(uniprot))

find_column <- function(exact_name, pattern = NULL) {
  hit <- names(uniprot)[match(exact_name, column_names_lower)]
  if (!is.na(hit)) return(hit)
  if (!is.null(pattern)) {
    hits <- names(uniprot)[grepl(pattern, column_names_lower)]
    if (length(hits) == 1L) return(hits)
  }
  NA_character_
}

gene_column <- find_column("gene names (primary)", "gene.*primary")
signal_column <- find_column("signal peptide", "^signal")
tm_column <- find_column("transmembrane", "^transmem")
location_column <- find_column("subcellular location [cc]", "subcellular")

if (anyNA(c(gene_column, signal_column, tm_column, location_column))) {
  fail("Could not identify required columns in frozen UniProt snapshot")
}

uniprot$gene <- toupper(trimws(as.character(uniprot[[gene_column]])))
uniprot$has_signal <- !is.na(uniprot[[signal_column]]) &
  nzchar(as.character(uniprot[[signal_column]]))
uniprot$has_transmembrane <- !is.na(uniprot[[tm_column]]) &
  nzchar(as.character(uniprot[[tm_column]]))
uniprot$membrane_location <- grepl(
  "membrane",
  tolower(ifelse(is.na(uniprot[[location_column]]), "", uniprot[[location_column]]))
)
uniprot$secreted_real <- uniprot$has_signal &
  !uniprot$has_transmembrane &
  !uniprot$membrane_location

secreted_real_genes <- sort(unique(
  uniprot$gene[
    !is.na(uniprot$gene) &
      nzchar(uniprot$gene) &
      uniprot$secreted_real
  ]
))

if (length(secreted_real_genes) == 0L) {
  fail("Frozen UniProt snapshot produced zero strict secreted genes")
}

message("Loading regenerated Seurat object read-only...")
object <- readRDS(input_path)
if (!("Spatial" %in% names(object@assays))) {
  fail("Regenerated object lacks the Spatial assay")
}

counts <- SeuratObject::LayerData(object[["Spatial"]], layer = "counts")
if (!(inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix"))) {
  fail("Spatial counts are not a sparse matrix")
}
if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  fail("Spatial counts lack feature or spot names")
}

feature_symbols <- toupper(rownames(counts))
if (anyDuplicated(feature_symbols)) {
  fail("Upper-cased Spatial feature names are duplicated")
}
rownames(counts) <- feature_symbols

library_sizes <- Matrix::colSums(counts)
if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
  fail("Invalid spot library sizes")
}

module_sizes <- as.data.frame(table(module_table$module), stringsAsFactors = FALSE)
names(module_sizes) <- c("module", "n_genes_module")
module_sizes$module <- as.character(module_sizes$module)
module_sizes$module_number <- as.integer(sub("^SM_", "", module_sizes$module))
module_sizes <- module_sizes[order(module_sizes$module_number), , drop = FALSE]

secreted_by_module <- lapply(module_sizes$module, function(module_name) {
  module_genes <- module_table$gene[module_table$module == module_name]
  sort(intersect(module_genes, secreted_real_genes))
})
names(secreted_by_module) <- module_sizes$module

summary_table <- data.frame(
  module = module_sizes$module,
  n_genes_module = module_sizes$n_genes_module,
  n_secreted_real = vapply(secreted_by_module, length, integer(1)),
  stringsAsFactors = FALSE
)
summary_table$pct_secreted_real <- 100 *
  summary_table$n_secreted_real / summary_table$n_genes_module

all_selected_genes <- sort(unique(unlist(secreted_by_module, use.names = FALSE)))
available_selected_genes <- intersect(all_selected_genes, rownames(counts))
missing_selected_genes <- setdiff(all_selected_genes, available_selected_genes)

message(
  "Calculating historical-equivalent log-normalized mean expression for ",
  length(available_selected_genes), " strict secreted factors..."
)

selected_counts <- counts[available_selected_genes, , drop = FALSE]
normalization_multiplier <- 10000 / library_sizes
normalized_selected <- selected_counts %*%
  Matrix::Diagonal(
    n = ncol(selected_counts),
    x = normalization_multiplier
  )
normalized_selected <- log1p(normalized_selected)

mean_expression <- Matrix::rowMeans(normalized_selected)
expression_table <- data.frame(
  gene = names(mean_expression),
  mean_normalized_expression = as.numeric(mean_expression),
  stringsAsFactors = FALSE
)

top30_by_module <- do.call(
  rbind,
  lapply(summary_table$module, function(module_name) {
    genes <- secreted_by_module[[module_name]]
    current <- expression_table[
      expression_table$gene %in% genes,
      ,
      drop = FALSE
    ]
    current <- current[
      order(-current$mean_normalized_expression, current$gene),
      ,
      drop = FALSE
    ]
    current <- head(current, 30L)
    if (nrow(current) == 0L) {
      return(data.frame(
        gene = character(),
        mean_normalized_expression = numeric(),
        module = character(),
        rank = integer(),
        stringsAsFactors = FALSE
      ))
    }
    current$module <- module_name
    current$rank <- seq_len(nrow(current))
    current
  })
)

if (is.null(top30_by_module) || nrow(top30_by_module) == 0L) {
  fail("No strict secreted factors were available for expression ranking")
}

historical_counts <- data.table::fread(historical_counts_path, data.table = FALSE)
historical_topgenes <- data.table::fread(historical_topgenes_path, data.table = FALSE)

historical_count_comparison <- merge(
  summary_table,
  historical_counts[, c("module", "n_genes_module", "n_secreted_real"), drop = FALSE],
  by = "module",
  all = TRUE,
  suffixes = c("_regenerated", "_historical")
)

historical_count_comparison$n_genes_module_difference <-
  historical_count_comparison$n_genes_module_regenerated -
  historical_count_comparison$n_genes_module_historical
historical_count_comparison$n_secreted_real_difference <-
  historical_count_comparison$n_secreted_real_regenerated -
  historical_count_comparison$n_secreted_real_historical

historical_top_overlap <- do.call(
  rbind,
  lapply(c("SM_1", "SM_6"), function(module_name) {
    historical_current <- historical_topgenes[
      historical_topgenes$module == module_name,
      ,
      drop = FALSE
    ]
    current <- top30_by_module[top30_by_module$module == module_name, , drop = FALSE]
    data.frame(
      module = module_name,
      regenerated_top_n = nrow(current),
      historical_top_n = nrow(historical_current),
      shared_gene_count = length(intersect(current$gene, historical_current$gene)),
      stringsAsFactors = FALSE
    )
  })
)

module_colours <- c(
  SM_1 = "#CAB2D6",
  SM_2 = "#6A3D9A",
  SM_3 = "#E69F00",
  SM_4 = "#7EC0EE",
  SM_5 = "#993F00",
  SM_6 = "#023FA5"
)

plot_table <- summary_table
plot_table$module <- factor(plot_table$module, levels = summary_table$module)

count_plot <- ggplot(
  plot_table,
  aes(x = module, y = n_secreted_real, fill = module)
) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n_secreted_real), hjust = -0.12, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = module_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Strict secreted factors per WGCNA module",
    subtitle = "Frozen reviewed-human UniProt snapshot; signal peptide, no transmembrane region, no membrane localization",
    x = NULL,
    y = "Number of secreted factors"
  ) +
  cowplot::theme_cowplot()

percentage_plot <- ggplot(
  plot_table,
  aes(x = module, y = pct_secreted_real, fill = module)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_secreted_real)),
    hjust = -0.12,
    size = 3.5
  ) +
  coord_flip() +
  scale_fill_manual(values = module_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Proportion of strict secreted factors per WGCNA module",
    subtitle = "Frozen reviewed-human UniProt snapshot",
    x = NULL,
    y = "Secreted factors (% of module genes)"
  ) +
  cowplot::theme_cowplot()

write_top_gene_plot <- function(module_name, panel_letter) {
  current <- top30_by_module[top30_by_module$module == module_name, , drop = FALSE]
  if (nrow(current) == 0L) fail("No ranked factors available for ", module_name)

  current <- current[
    order(current$mean_normalized_expression, current$gene),
    ,
    drop = FALSE
  ]
  current$gene <- factor(current$gene, levels = current$gene)

  plot <- ggplot(
    current,
    aes(x = gene, y = mean_normalized_expression)
  ) +
    geom_col(fill = module_colours[[module_name]], width = 0.78) +
    coord_flip() +
    labs(
      title = paste0(
        panel_letter, ". Top ", nrow(current),
        " strict secreted factors in ", module_name
      ),
      subtitle = "Ranked by mean log-normalized Spatial expression across all spots",
      x = NULL,
      y = "Mean log1p(count / spot library size × 10,000)"
    ) +
    cowplot::theme_cowplot()

  stem <- paste0(
    "Fig4", panel_letter, "_Top",
    nrow(current), "_", module_name,
    "_strict_secreted_factors_mean_normalized_expression"
  )

  ggplot2::ggsave(
    file.path(result_tmp, paste0(stem, ".pdf")),
    plot,
    width = 8.5,
    height = 6.8,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    file.path(result_tmp, paste0(stem, ".png")),
    plot,
    width = 8.5,
    height = 6.8,
    dpi = 300,
    bg = "white"
  )
}

message("Rendering secreted-factor count and manuscript candidate panels...")
ggplot2::ggsave(
  file.path(result_tmp, "SupFig_secreted_factor_counts_by_module.pdf"),
  count_plot,
  width = 8.6,
  height = 4.8,
  device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  file.path(result_tmp, "SupFig_secreted_factor_percent_by_module.pdf"),
  percentage_plot,
  width = 8.6,
  height = 4.8,
  device = grDevices::cairo_pdf
)
write_top_gene_plot("SM_1", "H")
write_top_gene_plot("SM_6", "I")

definition_table <- data.frame(
  parameter = c(
    "secreted_factor_definition",
    "uniprot_resource",
    "uniprot_snapshot_sha256",
    "organism",
    "reviewed_records_only",
    "expression_input",
    "expression_normalization",
    "expression_summary",
    "module_source",
    "neighbour_or_hub_selection"
  ),
  value = c(
    "signal peptide AND no transmembrane region AND no membrane localization",
    uniprot_path,
    uniprot_sha,
    "Homo sapiens",
    "TRUE; frozen historical snapshot",
    "Spatial counts layer of regenerated Seurat object",
    "log1p(count / spot library size * 10000), calculated without modifying input object",
    "mean across all spots",
    module_path,
    "none; all genes assigned to each regenerated WGCNA module are screened"
  ),
  stringsAsFactors = FALSE
)

uniprot_summary <- data.frame(
  snapshot_rows = nrow(uniprot),
  unique_primary_gene_symbols = length(unique(uniprot$gene)),
  strict_secreted_real_gene_symbols = length(secreted_real_genes),
  genes_available_for_expression_ranking = length(available_selected_genes),
  selected_secreted_genes_absent_from_object = length(missing_selected_genes),
  stringsAsFactors = FALSE
)

input_hashes <- data.frame(
  input = c(
    "regenerated_seurat_object",
    "regenerated_module_assignment_table",
    "frozen_historical_uniprot_snapshot",
    "historical_count_table",
    "historical_topgene_table"
  ),
  path = c(
    input_path,
    module_path,
    uniprot_path,
    historical_counts_path,
    historical_topgenes_path
  ),
  sha256 = c(
    input_sha,
    module_sha,
    uniprot_sha,
    sha256(historical_counts_path),
    sha256(historical_topgenes_path)
  ),
  stringsAsFactors = FALSE
)

atomic_write_table(summary_table, file.path(result_tmp, "12_secreted_factor_counts_and_percentages.tsv"))
atomic_write_table(top30_by_module, file.path(result_tmp, "12_top30_strict_secreted_factors_by_module.tsv"))
atomic_write_table(
  data.frame(gene = missing_selected_genes, stringsAsFactors = FALSE),
  file.path(result_tmp, "12_selected_secreted_genes_absent_from_expression_object.tsv")
)
atomic_write_table(
  historical_count_comparison,
  file.path(result_tmp, "12_historical_vs_regenerated_module_count_comparison.tsv")
)
atomic_write_table(
  historical_top_overlap,
  file.path(result_tmp, "12_historical_vs_regenerated_SM1_SM6_top30_overlap.tsv")
)
atomic_write_table(definition_table, file.path(result_tmp, "12_secreted_factor_definition.tsv"))
atomic_write_table(uniprot_summary, file.path(result_tmp, "12_uniprot_snapshot_summary.tsv"))
atomic_write_table(input_hashes, file.path(log_tmp, "12_input_SHA256.tsv"))
atomic_write_table(
  data.frame(
    package = c("R", required_packages),
    version = c(
      R.version.string,
      vapply(
        required_packages,
        function(package_name) as.character(utils::packageVersion(package_name)),
        character(1)
      )
    ),
    stringsAsFactors = FALSE
  ),
  file.path(log_tmp, "12_package_versions.tsv")
)

finished_at <- Sys.time()

atomic_write_lines(
  c(
    "FLINC Figure 4 strict secreted-factor execution summary",
    paste("Started UTC:", format(started_at, tz = "UTC", usetz = TRUE)),
    paste("Completed UTC:", format(finished_at, tz = "UTC", usetz = TRUE)),
    paste("Elapsed seconds:", as.numeric(difftime(finished_at, started_at, units = "secs"))),
    paste("Regenerated object SHA256:", input_sha),
    paste("Regenerated module-table SHA256:", module_sha),
    paste("Frozen UniProt snapshot SHA256:", uniprot_sha),
    paste("Secreted definition: signal peptide AND no transmembrane AND no membrane localization"),
    paste("Expression normalization: log1p(count / library size * 10000)"),
    paste("SM_1 historical/regenerated top30 overlap:",
          historical_top_overlap$shared_gene_count[
            historical_top_overlap$module == "SM_1"
          ]),
    paste("SM_6 historical/regenerated top30 overlap:",
          historical_top_overlap$shared_gene_count[
            historical_top_overlap$module == "SM_6"
          ]),
    paste("Command:", paste(commandArgs(), collapse = " "))
  ),
  file.path(log_tmp, "12_execution_summary.txt")
)
atomic_write_lines(capture.output(sessionInfo()), file.path(log_tmp, "12_sessionInfo.txt"))

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
manifest <- data.frame(
  file = basename(result_files),
  bytes = file.info(result_files)$size,
  sha256 = vapply(result_files, sha256, character(1)),
  stringsAsFactors = FALSE
)
atomic_write_table(manifest, file.path(result_tmp, "12_output_SHA256.tsv"))

atomic_write_lines("SUCCESS", file.path(result_tmp, "12_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "12_COMPLETED.txt"))

if (!file.rename(result_tmp, output_dir)) {
  fail("Could not atomically publish results")
}
if (!file.rename(log_tmp, log_dir)) {
  fail("Could not atomically publish logs")
}

completed <- TRUE
message("Strict secreted-factor analysis completed successfully.")
