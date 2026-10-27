#!/usr/bin/env Rscript

# FLINC Figure 4: render the frozen historical secreted-factor enrichment panel.
# The historical GO:CC/KEGG result is treated as a frozen resource.  This script
# never queries KEGG and does not perform a new enrichment calculation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L || any(!startsWith(args[seq(1L, length(args), 2L)], "--"))) {
  stop("Arguments must be supplied as --key value pairs", call. = FALSE)
}
arg_list <- setNames(args[seq(2L, length(args), 2L)], sub("^--", "", args[seq(1L, length(args), 2L)]))
required <- c("frozen-top2", "input-audit", "output-dir", "log-dir")
missing <- setdiff(required, names(arg_list))
if (length(missing) > 0L) stop("Missing argument(s): ", paste(missing, collapse = ", "), call. = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("cowplot", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
  stop("Packages ggplot2, cowplot, and digest are required", call. = FALSE)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256")
need_file <- function(path) normalizePath(path, mustWork = TRUE)
frozen_top2_path <- need_file(arg_list[["frozen-top2"]])
input_audit_path <- need_file(arg_list[["input-audit"]])
output_dir <- normalizePath(arg_list[["output-dir"]], mustWork = FALSE)
log_dir <- normalizePath(arg_list[["log-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

atomic_write_table <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.table(x, tmp, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
}
atomic_write_lines <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not publish ", path, call. = FALSE)
}

audit <- utils::read.delim(input_audit_path, check.names = FALSE, stringsAsFactors = FALSE)
required_audit <- c("module", "module_genes_identical", "secreted_genes_identical", "module_colour")
if (!all(required_audit %in% names(audit))) stop("Input audit lacks required columns", call. = FALSE)
audit <- audit[audit$module %in% c("SM_1", "SM_6"), , drop = FALSE]
if (nrow(audit) != 2L || !all(audit$module_genes_identical) || !all(audit$secreted_genes_identical)) {
  stop("Cannot render frozen enrichment: exact input-identity audit did not pass", call. = FALSE)
}
palette <- stats::setNames(audit$module_colour, audit$module)

top <- utils::read.csv(frozen_top2_path, check.names = FALSE, stringsAsFactors = FALSE)
name_for <- function(options) {
  present <- options[options %in% names(top)]
  if (length(present) == 0L) NULL else present[[1L]]
}
term_col <- name_for(c("term", "Description"))
p_col <- name_for(c("p_adj", "p.adjust"))
score_col <- name_for(c("score"))
if (is.null(term_col) || is.null(p_col) || !all(c("module", "db") %in% names(top))) {
  stop("Frozen top-2 table lacks module, db, term/Description, or p_adj/p.adjust", call. = FALSE)
}
top$term <- as.character(top[[term_col]])
top$p_adj <- as.numeric(top[[p_col]])
top$score <- if (is.null(score_col)) -log10(top$p_adj) else as.numeric(top[[score_col]])
top <- top[top$module %in% c("SM_1", "SM_6") & top$db %in% c("GO:CC", "KEGG"), c("module", "db", "term", "p_adj", "score"), drop = FALSE]
if (nrow(top) != 8L || any(table(top$module, top$db) != 2L) || any(!is.finite(top$score))) {
  stop("Frozen table must contain exactly two finite-ranked terms for each SM_1/SM_6 × GO:CC/KEGG panel", call. = FALSE)
}

make_panel <- function(module, database) {
  d <- top[top$module == module & top$db == database, , drop = FALSE]
  d <- d[order(d$score, d$term), , drop = FALSE]
  d$term_wrapped <- vapply(d$term, function(x) paste(strwrap(x, width = 34L), collapse = "\n"), character(1))
  d$term_wrapped <- factor(d$term_wrapped, levels = d$term_wrapped)
  ggplot2::ggplot(d, ggplot2::aes(x = score, y = term_wrapped)) +
    ggplot2::geom_col(fill = palette[[module]], width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("FDR %.2g", p_adj)), hjust = -0.08, size = 3.1) +
    ggplot2::coord_cartesian(clip = "off", xlim = c(0, max(d$score) * 1.25)) +
    ggplot2::labs(title = paste0(module, " — ", database), x = expression(-log[10](BH~adjusted~P)), y = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.92),
      plot.margin = ggplot2::margin(5.5, 35, 5.5, 110)
    )
}

panels <- list(
  make_panel("SM_1", "GO:CC"), make_panel("SM_6", "GO:CC"),
  make_panel("SM_1", "KEGG"), make_panel("SM_6", "KEGG")
)
figure <- cowplot::plot_grid(
  plotlist = panels,
  ncol = 1,
  rel_heights = c(1, 1, 1, 1)
)
figure <- cowplot::ggdraw() +
  cowplot::draw_label("Enrichment of strict secreted factors in SM_1 and SM_6", x = 0.5, y = 0.99, fontface = "bold", size = 15) +
  cowplot::draw_label("GO Cellular Component and KEGG; top 2 BH-ranked terms per module", x = 0.5, y = 0.962, size = 10) +
  cowplot::draw_plot(figure, x = 0, y = 0, width = 1, height = 0.94)

stem <- "SupFig_secreted_factor_enrichment_SM1_SM6_GOCC_KEGG_top2"
pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
png_path <- file.path(output_dir, paste0(stem, ".png"))
ggplot2::ggsave(pdf_path, figure, width = 11, height = 14, device = grDevices::cairo_pdf)
ggplot2::ggsave(png_path, figure, width = 11, height = 14, dpi = 300)

atomic_write_table(top, file.path(output_dir, "14_frozen_top2_enrichment_terms.tsv"))
atomic_write_table(
  data.frame(
    input = c("frozen_historical_top2_enrichment", "passed_input_identity_audit"),
    path = c(frozen_top2_path, input_audit_path),
    sha256 = c(sha256(frozen_top2_path), sha256(input_audit_path)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "14_input_checksums.tsv")
)
manifest_paths <- c(pdf_path, png_path, file.path(output_dir, "14_frozen_top2_enrichment_terms.tsv"), file.path(output_dir, "14_input_checksums.tsv"))
atomic_write_table(
  data.frame(
    file = basename(manifest_paths),
    bytes = file.info(manifest_paths)$size,
    sha256 = vapply(manifest_paths, sha256, character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "14_output_SHA256.tsv")
)
atomic_write_lines(c(
  "FLINC Figure 4 frozen secreted-factor enrichment rendering summary",
  "Method: visual rendering only; no live GO or KEGG query was made.",
  paste("Frozen top-2 enrichment SHA256:", sha256(frozen_top2_path)),
  paste("Verified input-audit SHA256:", sha256(input_audit_path)),
  "SM_1 and SM_6 strict-secreted gene lists and WGCNA universes were exactly identical to the historical inputs.",
  paste("R:", R.version.string)
), file.path(log_dir, "14_execution_summary.txt"))
atomic_write_lines("SUCCESS", file.path(output_dir, "14_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_dir, "14_COMPLETED.txt"))
cat("FROZEN_SECRETED_ENRICHMENT_RENDERING_OK\n")
