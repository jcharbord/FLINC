#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 10, scipen = 999)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
set.seed(42L)

fail <- function(...) stop(paste0(...), call. = FALSE)
required_packages <- c("hdWGCNA", "digest", "ggplot2", "ggrepel", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) fail("Missing packages: ", paste(missing_packages, collapse = ", "))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) fail("Arguments must be supplied as --key value pairs")
arg_list <- as.list(args[seq(2L, length(args), by = 2L)])
names(arg_list) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])
required_args <- c("network", "tom", "module-table", "annotation", "output-dir", "log-dir")
if (length(setdiff(required_args, names(arg_list))) > 0L) fail("Missing arguments: ", paste(setdiff(required_args, names(arg_list)), collapse = ", "))

network_path <- normalizePath(arg_list[["network"]], mustWork = TRUE)
tom_path <- normalizePath(arg_list[["tom"]], mustWork = TRUE)
module_path <- normalizePath(arg_list[["module-table"]], mustWork = TRUE)
annotation_path <- normalizePath(arg_list[["annotation"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(log_dir) || file.exists(log_dir)) fail("Final output or log path already exists")
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_dir), recursive = TRUE, showWarnings = FALSE)
run_stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
result_tmp <- paste0(output_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
log_tmp <- paste0(log_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
dir.create(result_tmp, recursive = FALSE, showWarnings = FALSE)
dir.create(log_tmp, recursive = FALSE, showWarnings = FALSE)
completed <- FALSE
on.exit(if (!completed) message("Incomplete run retained: ", result_tmp, "\n", log_tmp), add = TRUE)

atomic_write_table <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(x, temporary, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}
atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}
sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
script_path <- {
  found <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(found) == 1L) normalizePath(sub("^--file=", "", found), mustWork = TRUE) else NA_character_
}
started_at <- Sys.time()

message("Hashing and loading regenerated Figure 2 inputs read-only...")
input_hashes <- data.frame(
  input = c("network_rds", "tom_rda", "module_table", "annotation_table"),
  path = c(network_path, tom_path, module_path, annotation_path),
  sha256 = c(sha256(network_path), sha256(tom_path), sha256(module_path), sha256(annotation_path)),
  stringsAsFactors = FALSE
)
module_table <- utils::read.delim(module_path, check.names = FALSE, stringsAsFactors = FALSE)
annotation <- utils::read.delim(annotation_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene_name", "module") %in% names(module_table))) fail("Module table lacks gene_name/module")
if (!all(c("module", "gene", "annotation", "curation_basis") %in% names(annotation))) fail("Annotation table lacks required columns")
if (anyDuplicated(module_table$gene_name)) fail("Module table has duplicate gene names")
if (any(!nzchar(annotation$annotation))) fail("Annotation labels must be non-empty")
if (anyDuplicated(annotation[c("module", "gene")])) fail("Extended annotation table has duplicate module/gene entries")

suppressPackageStartupMessages(library(hdWGCNA))
network_object <- readRDS(network_path)
tom_dissimilarity <- GetTOM(network_object, wgcna_name = "vis")
if (!is.matrix(tom_dissimilarity) || !identical(rownames(tom_dissimilarity), colnames(tom_dissimilarity))) fail("GetTOM did not return named square matrix")
if (any(!is.finite(tom_dissimilarity)) || min(tom_dissimilarity) < -1e-12 || max(tom_dissimilarity) > 1 + 1e-12) fail("Invalid TOM dissimilarity values")
if (any(abs(diag(tom_dissimilarity)) > 1e-12)) fail("Expected GetTOM diagonal to be zero for a dissimilarity matrix")

tom_holder <- new.env(parent = emptyenv())
load(tom_path, envir = tom_holder)
if (!exists("consTomDS", envir = tom_holder, inherits = FALSE)) fail("TOM RDA lacks consTomDS")
stored_dissimilarity <- get("consTomDS", envir = tom_holder)
if (!inherits(stored_dissimilarity, "dist") || attr(stored_dissimilarity, "Size") != nrow(tom_dissimilarity)) fail("Stored TOM is not matching consTomDS distance representation")

dist_entry <- function(x, i, j) {
  if (i == j) return(0)
  if (i > j) { temporary <- i; i <- j; j <- temporary }
  n <- attr(x, "Size")
  x[n * (i - 1L) - i * (i - 1L) / 2L + j - i]
}

role_for_gene <- function(module, gene, index_gene) {
  if (identical(gene, index_gene)) return("Index")
  labels <- unique(annotation$annotation[annotation$module == module & annotation$gene == gene])
  if (length(labels) == 0L) return("Other")
  if (length(labels) != 1L) fail("Expected one extended functional-context label for ", module, "/", gene)
  labels[[1]]
}

module_border <- c(SM_1 = "#4B1D68", SM_3 = "#D69200")
role_order <- c(
  "Fibrosis / ECM",
  "Immune / inflammation",
  "Stress / signalling",
  "Structure / adhesion",
  "Transport / homeostasis",
  "Lipid metabolism",
  "Metabolism / detoxification",
  "Other"
)
role_colours <- c(
  Index = "#F05A28",
  "Fibrosis / ECM" = "#1976D2",
  "Immune / inflammation" = "#43A047",
  "Stress / signalling" = "#8E24AA",
  "Structure / adhesion" = "#00897B",
  "Transport / homeostasis" = "#F9A825",
  "Lipid metabolism" = "#D81B60",
  "Metabolism / detoxification" = "#6A1B9A",
  Other = "#BDBDBD"
)
# Legacy F/G values in exported panel columns are retained for traceability;
# output filenames use gene names without assigning new manuscript panel letters.
specifications <- data.frame(index_gene = c("PCSK7", "PNPLA3"), expected_module = c("SM_1", "SM_3"), panel = c("F", "G"), stringsAsFactors = FALSE)
top_k <- 20L
all_ranked <- list(); all_nodes <- list(); all_edges <- list(); plots <- list()

for (i in seq_len(nrow(specifications))) {
  index_gene <- specifications$index_gene[[i]]
  expected_module <- specifications$expected_module[[i]]
  panel <- specifications$panel[[i]]
  actual_module <- module_table$module[match(index_gene, module_table$gene_name)]
  if (is.na(actual_module) || actual_module != expected_module) fail(index_gene, " is not assigned to expected ", expected_module)
  module_genes <- intersect(module_table$gene_name[module_table$module == actual_module], rownames(tom_dissimilarity))
  if (!(index_gene %in% module_genes)) fail(index_gene, " is absent from its TOM module")
  similarity <- 1 - tom_dissimilarity[index_gene, module_genes]
  similarity <- similarity[names(similarity) != index_gene]
  similarity <- similarity[order(-similarity, names(similarity))]
  if (length(similarity) < top_k) fail("Fewer than ", top_k, " available neighbours for ", index_gene)
  selected <- head(similarity, top_k)
  selected_genes <- names(selected)
  positions <- match(c(index_gene, selected_genes), rownames(tom_dissimilarity))
  stored_selected <- vapply(positions[-1], function(j) dist_entry(stored_dissimilarity, positions[[1]], j), numeric(1))
  if (!identical(unname(tom_dissimilarity[index_gene, selected_genes]), unname(stored_selected))) fail("GetTOM and consTomDS disagree for selected ", index_gene, " edges")
  ranked <- data.frame(panel = panel, index_gene = index_gene, module = actual_module, rank = seq_len(top_k), neighbor = selected_genes, TOM_dissimilarity = as.numeric(tom_dissimilarity[index_gene, selected_genes]), TOM_similarity = as.numeric(selected), stringsAsFactors = FALSE)
  ranked$role <- vapply(ranked$neighbor, function(gene) role_for_gene(actual_module, gene, index_gene), character(1))
  all_ranked[[index_gene]] <- ranked
  nodes <- data.frame(panel = panel, index_gene = index_gene, module = actual_module, gene = c(index_gene, selected_genes), role = c("Index", ranked$role), stringsAsFactors = FALSE)
  angles <- seq(pi / 2, pi / 2 - 2 * pi + 2 * pi / top_k, length.out = top_k)
  nodes$x <- c(0, cos(angles)); nodes$y <- c(0, sin(angles))
  nodes$TOM_similarity <- c(NA_real_, ranked$TOM_similarity)
  all_nodes[[index_gene]] <- nodes
  edges <- data.frame(panel = panel, index_gene = index_gene, module = actual_module, from = index_gene, to = selected_genes, TOM_similarity = ranked$TOM_similarity, stringsAsFactors = FALSE)
  edge_xy <- merge(edges, nodes[, c("gene", "x", "y")], by.x = "from", by.y = "gene", all.x = TRUE)
  edge_xy <- merge(edge_xy, nodes[, c("gene", "x", "y")], by.x = "to", by.y = "gene", all.x = TRUE, suffixes = c("_from", "_to"))
  all_edges[[index_gene]] <- edges
  present_role_breaks <- c("Index", role_order[role_order %in% unique(ranked$role)])
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_xy, ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = TOM_similarity), colour = "#636363", alpha = 0.78) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y, colour = role), size = 5.5) +
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x = x, y = y, label = gene), seed = 42, size = 3.6, max.overlaps = Inf, box.padding = 0.25, point.padding = 0.1, min.segment.length = 0, segment.colour = NA) +
    ggplot2::scale_colour_manual(values = role_colours, breaks = present_role_breaks, name = "functional context") +
    ggplot2::scale_linewidth_continuous(range = c(0.45, 2.7), name = "TOM similarity") +
    ggplot2::coord_equal(xlim = c(-1.38, 1.38), ylim = c(-1.25, 1.25), clip = "off") +
    ggplot2::labs(title = paste0(index_gene, " ego (k = ", top_k, ") in module ", actual_module), subtitle = "Strict top-20 TOM neighbours; colours show curated functional context") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8), legend.position = "right", plot.margin = ggplot2::margin(10, 24, 10, 24), panel.border = ggplot2::element_rect(colour = module_border[[actual_module]], fill = NA, linewidth = 1.2))
  plots[[index_gene]] <- p
  stem <- paste0("Fig2_", index_gene, "_ego_top20_extended_functional_context_", actual_module)
  ggplot2::ggsave(file.path(result_tmp, paste0(stem, ".pdf")), p, width = 8.5, height = 6.6, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(result_tmp, paste0(stem, ".png")), p, width = 8.5, height = 6.6, dpi = 300, bg = "white")
}

ranked_table <- do.call(rbind, all_ranked)
node_table <- do.call(rbind, all_nodes)
edge_table <- do.call(rbind, all_edges)
atomic_write_table(ranked_table, file.path(result_tmp, "11_ego_top20_neighbors.tsv"))
atomic_write_table(node_table, file.path(result_tmp, "11_ego_network_nodes.tsv"))
atomic_write_table(edge_table, file.path(result_tmp, "11_ego_network_edges.tsv"))
atomic_write_table(annotation, file.path(result_tmp, "11_historical_curated_functional_annotations.tsv"))
atomic_write_table(input_hashes, file.path(log_tmp, "11_input_SHA256.tsv"))

parameters <- data.frame(parameter = c("network_source", "tom_storage", "GetTOM_interpretation", "edge_weight", "selection_scope", "selection_rule", "top_k", "edge_set", "role_source", "role_rule", "layout", "layout_seed", "blas_threads"), value = c(network_path, "consTomDS dist object", "TOM dissimilarity", "1 - GetTOM()", "same regenerated module as index gene", "descending TOM similarity; gene-name ascending tie break", as.character(top_k), "index-to-top20-neighbour edges only", "explicit extended functional-context annotation table", "Index overrides other role; all strict neighbours receive one predefined functional-context label", "deterministic radial star", "42", "1"), stringsAsFactors = FALSE)
atomic_write_table(parameters, file.path(result_tmp, "11_run_parameters.tsv"))
atomic_write_table(data.frame(package = c("R", required_packages), version = c(R.version.string, vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))), stringsAsFactors = FALSE), file.path(log_tmp, "11_package_versions.tsv"))

finished_at <- Sys.time()
atomic_write_lines(c("FLINC Figure 2 extended-functional-context ego-network candidate execution summary", paste("Started UTC:", format(started_at, tz = "UTC", usetz = TRUE)), paste("Completed UTC:", format(finished_at, tz = "UTC", usetz = TRUE)), paste("Elapsed seconds:", as.numeric(difftime(finished_at, started_at, units = "secs"))), paste("Network RDS SHA256:", input_hashes$sha256[[1]]), paste("TOM RDA SHA256:", input_hashes$sha256[[2]]), paste("GetTOM interpretation: dissimilarity"), paste("Edge weight: TOM similarity = 1 - GetTOM()"), paste("PCSK7-PNPLA3 similarity audit:", 1 - tom_dissimilarity["PCSK7", "PNPLA3"]), paste("Command:", paste(commandArgs(), collapse = " "))), file.path(log_tmp, "11_execution_summary.txt"))
atomic_write_lines(capture.output(sessionInfo()), file.path(log_tmp, "11_sessionInfo.txt"))

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
manifest <- data.frame(file = basename(result_files), bytes = file.info(result_files)$size, sha256 = vapply(result_files, sha256, character(1)), stringsAsFactors = FALSE)
atomic_write_table(manifest, file.path(result_tmp, "11_output_SHA256.tsv"))
atomic_write_lines("SUCCESS", file.path(result_tmp, "11_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "11_COMPLETED.txt"))
if (!file.rename(result_tmp, output_dir)) fail("Could not atomically publish results")
if (!file.rename(log_tmp, log_dir)) fail("Could not atomically publish logs")
completed <- TRUE
message("Figure 2 ego networks completed successfully.")
