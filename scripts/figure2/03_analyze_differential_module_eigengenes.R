# Historical spot-level DME; retained for traceability.
# For final manuscript inference use 03_bis_analyze_differential_module_eigengenes_donor_level.R.
#!/usr/bin/env Rscript

# FLINC Figure 2 — differential harmonized module eigengene analysis.
#
# Input:  step-02 hdWGCNA object with harmonized MEs and module assignments.
# Output: stage-vs-F0 Wilcoxon results and summary plots.
# Statistical calculations use the six harmonized module eigengenes (hMEs).

options(stringsAsFactors = FALSE, warn = 1)
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(ggplot2)
  library(ggrepel)
})

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

move_into_place <- function(partial, path) {
  if (file.exists(path)) stopf("Refusing to overwrite existing output: %s", path)
  if (!file.rename(partial, path)) stopf("Could not move partial file into place: %s", path)
}

write_csv_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) stopf("Partial output already exists: %s", partial)
  utils::write.csv(x, partial, row.names = FALSE, na = "NA")
  move_into_place(partial, path)
}

write_tsv_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) stopf("Partial output already exists: %s", partial)
  utils::write.table(x, partial, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
  move_into_place(partial, path)
}

write_lines_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) stopf("Partial output already exists: %s", partial)
  writeLines(x, partial, useBytes = TRUE)
  move_into_place(partial, path)
}

save_plot_atomic <- function(plot, path, width, height) {
  partial <- paste0(path, ".partial.pdf")
  if (file.exists(partial)) stopf("Partial plot already exists: %s", partial)
  ggplot2::ggsave(partial, plot, device = grDevices::cairo_pdf,
                  width = width, height = height, units = "in")
  move_into_place(partial, path)
}

repro_root <- Sys.getenv(
  "FLINC_REPRO_ROOT",
  unset = "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)

input_file <- Sys.getenv(
  "FIG2_NETWORK_INPUT",
  unset = file.path(repro_root, "objects", "Fig2_WGCNA", "02_wgcna_network_modules.rds")
)

results_dir <- file.path(repro_root, "results", "Fig2_WGCNA", "03_DME")
log_dir <- file.path(repro_root, "reproducibility_logs", "Fig2_WGCNA")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

csv_file <- file.path(results_dir, "03_DiffModuleEigengenes_harmonized.csv")
tsv_file <- file.path(results_dir, "03_DiffModuleEigengenes_harmonized.tsv")
xlsx_file <- file.path(results_dir, "03_DiffModuleEigengenes_harmonized.xlsx")
lollipop_file <- file.path(results_dir, "03_DME_lollipop_all_stages_vs_F0.pdf")
volcano_file <- file.path(results_dir, "03_DME_volcano_all_stages_vs_F0.pdf")
summary_file <- file.path(log_dir, "03_DME_summary.tsv")
session_file <- file.path(log_dir, "03_DME_sessionInfo.txt")

all_outputs <- c(csv_file, tsv_file, xlsx_file, lollipop_file, volcano_file,
                 summary_file, session_file)
existing <- all_outputs[file.exists(all_outputs)]
if (length(existing) > 0L) {
  stopf("Step 03 output already exists; do not rerun. Existing: %s",
        paste(existing, collapse = ", "))
}

if (!file.exists(input_file)) stopf("Step-02 input does not exist: %s", input_file)

say("Loading completed step-02 object:", input_file)
object <- readRDS(input_file)
if (!inherits(object, "Seurat")) stopf("Input is not a Seurat object")
if (!("F_score" %in% colnames(object[[]]))) stopf("Input lacks F_score metadata")
if (!("vis" %in% names(object@misc))) stopf("Input lacks hdWGCNA experiment 'vis'")

MEs <- GetMEs(object, harmonized = TRUE, wgcna_name = "vis")
if (is.null(MEs)) stopf("No harmonized module eigengenes found")
MEs <- as.data.frame(MEs)
if (nrow(MEs) != ncol(object)) {
  stopf("hME rows (%d) do not match object spots (%d)", nrow(MEs), ncol(object))
}

modules <- grep("^SM_[0-9]+$", colnames(MEs), value = TRUE)
modules <- modules[order(as.integer(sub("^SM_", "", modules)))]
expected_modules <- paste0("SM_", 1:6)
if (!identical(modules, expected_modules)) {
  stopf("Expected modules [%s], observed [%s]",
        paste(expected_modules, collapse = ", "), paste(modules, collapse = ", "))
}

stage <- as.character(object$F_score)
expected_stages <- c("F0", "F1", "F2", "F3", "F4")
observed_stages <- sort(unique(stage))
if (!identical(observed_stages, expected_stages)) {
  stopf("Expected stages [%s], observed [%s]",
        paste(expected_stages, collapse = ", "), paste(observed_stages, collapse = ", "))
}

stage_counts <- table(factor(stage, levels = expected_stages))
if (any(stage_counts == 0L)) stopf("At least one expected stage contains no spots")
say("Spot counts:", paste(names(stage_counts), stage_counts, sep = "=", collapse = ", "))

comparisons <- paste0(c("F1", "F2", "F3", "F4"), "_vs_F0")
rows <- vector("list", length(comparisons) * length(modules))
k <- 0L
for (comparison in comparisons) {
  current_stage <- sub("_vs_F0$", "", comparison)
  i0 <- which(stage == "F0")
  i1 <- which(stage == current_stage)
  for (module in modules) {
    k <- k + 1L
    x0 <- MEs[i0, module]
    x1 <- MEs[i1, module]
    test <- suppressWarnings(stats::wilcox.test(x1, x0, alternative = "two.sided"))
    rows[[k]] <- data.frame(
      module = module,
      comparison = comparison,
      n_F0_spots = sum(!is.na(x0)),
      n_stage_spots = sum(!is.na(x1)),
      mean_F0 = mean(x0, na.rm = TRUE),
      mean_stage = mean(x1, na.rm = TRUE),
      delta_mean = mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE),
      p_wilcox = test$p.value,
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, rows)
results$p_adj <- ave(results$p_wilcox, results$comparison,
                     FUN = function(x) stats::p.adjust(x, method = "BH"))
results$significant_FDR_0.05 <- results$p_adj < 0.05
results <- results[order(match(results$comparison, comparisons),
                         match(results$module, modules)), ]
rownames(results) <- NULL

say("Writing DME result tables")
write_csv_atomic(results, csv_file)
write_tsv_atomic(results, tsv_file)
if (!requireNamespace("writexl", quietly = TRUE)) {
  stopf("Package 'writexl' is required for the XLSX export")
}
xlsx_partial <- paste0(xlsx_file, ".partial.xlsx")
if (file.exists(xlsx_partial)) stopf("Partial output already exists: %s", xlsx_partial)
writexl::write_xlsx(results, xlsx_partial)
move_into_place(xlsx_partial, xlsx_file)

plot_df <- results
plot_df$comparison <- factor(plot_df$comparison, levels = comparisons)
plot_df$module <- factor(plot_df$module, levels = rev(modules))
offsets <- c(F1_vs_F0 = -0.24, F2_vs_F0 = -0.08,
             F3_vs_F0 = 0.08, F4_vs_F0 = 0.24)
plot_df$module_index <- as.numeric(plot_df$module)
plot_df$y <- plot_df$module_index + offsets[as.character(plot_df$comparison)]

stage_palette <- c(F1_vs_F0 = "#fcbba1", F2_vs_F0 = "#fc9272",
                   F3_vs_F0 = "#fb6a4a", F4_vs_F0 = "#cb181d")

p_lollipop <- ggplot(plot_df) +
  geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
  geom_segment(aes(x = 0, xend = delta_mean, y = y, yend = y,
                   colour = comparison), linewidth = 1) +
  geom_point(aes(x = delta_mean, y = y, fill = comparison),
             shape = 21, size = 3.2, stroke = 0.2) +
  scale_colour_manual(values = stage_palette, guide = "none") +
  scale_fill_manual(values = stage_palette, name = "Contrast",
                    labels = c("F1 vs F0", "F2 vs F0", "F3 vs F0", "F4 vs F0")) +
  scale_y_continuous(breaks = seq_along(modules), labels = rev(modules),
                     expand = expansion(add = 0.7)) +
  labs(x = expression(Delta*" mean hME (stage - F0)"), y = NULL,
       title = "Differential harmonized module eigengenes") +
  theme_minimal(base_family = "Helvetica") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
save_plot_atomic(p_lollipop, lollipop_file, 10, 7)

module_table <- GetModules(object, wgcna_name = "vis")
module_colors <- unique(module_table[c("module", "color")])
module_palette <- setNames(module_colors$color, module_colors$module)
plot_df$neglog10_padj <- -log10(pmax(plot_df$p_adj, .Machine$double.xmin))
plot_df$point_color <- ifelse(plot_df$significant_FDR_0.05,
                             unname(module_palette[as.character(plot_df$module)]),
                             "grey70")
plot_df$point_color[is.na(plot_df$point_color)] <- "black"

p_volcano <- ggplot(plot_df, aes(x = delta_mean, y = neglog10_padj)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = point_color), size = 2.8) +
  ggrepel::geom_text_repel(aes(label = module, color = point_color), size = 3,
                           min.segment.length = 0, max.overlaps = Inf,
                           show.legend = FALSE) +
  scale_color_identity() +
  facet_wrap(~ comparison, scales = "free") +
  labs(x = expression(Delta*" mean hME (stage - F0)"),
       y = expression(-log[10](p[adj])),
       title = "Differential harmonized module eigengenes") +
  theme_bw(base_size = 12)
save_plot_atomic(p_volcano, volcano_file, 12, 10)

summary <- data.frame(
  input_file = input_file,
  input_size_bytes = file.info(input_file)$size,
  n_spots = ncol(object),
  n_modules = length(modules),
  n_tests = nrow(results),
  n_significant_FDR_0.05 = sum(results$significant_FDR_0.05),
  random_seed = 12345,
  stringsAsFactors = FALSE
)
write_tsv_atomic(summary, summary_file)
write_lines_atomic(capture.output(sessionInfo()), session_file)

say("Significant DME tests at BH FDR < 0.05:", sum(results$significant_FDR_0.05))
say("Fig. 2 differential module eigengene analysis completed successfully")
