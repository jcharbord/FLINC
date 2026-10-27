#!/usr/bin/env Rscript

# FLINC Figure 2 -- donor-level differential module eigengene analysis.
# Based on 03_analyze_differential_module_eigengenes.R and DME_stat_donor_level_.Rmd.
#
# Input: completed step-02 Seurat/hdWGCNA object. No network reconstruction.
# Analysis: mean spot score per donor and module; independent donor Wilcoxon
# comparisons, F1-F4 versus Healthy; BH across SIX modules within each contrast.
# Outputs: CSV/TSV/XLSX, donor means and QC, volcano/lollipop PDF+PNG, run logs.
#
# Run in RStudio: source("scripts/figure2/03_bis_analyze_differential_module_eigengenes_donor_level.R")
# Or: Rscript --vanilla scripts/figure2/03_bis_analyze_differential_module_eigengenes_donor_level.R
#
# Optional overrides, set BEFORE running (defaults match the GitHub workflow):
# Sys.setenv(FLINC_REPRO_ROOT = "/path/to/FLINC/Jeremie/Reproducible_code")
# Sys.setenv(FIG2_NETWORK_INPUT = "/path/to/02_wgcna_network_modules.rds")
# Sys.setenv(FIG2_DONOR_COLUMN = "donor", FIG2_STAGE_COLUMN = "F_score")
# Sys.setenv(FIG2_DME_RESULTS_DIR = "/path/to/a/new/output/folder")
# FIG2_ME_SOURCE: "harmonized" (default), "raw", or "metadata".
# Use metadata ONLY if the six metadata SM_* columns are the intended scores.
# There is deliberately no automatic fallback between different ME sources.
# To reproduce the OLD Rmd numerically, select its SAME object and ME source.
# The old CSV alone does not establish which source that Rmd selected.
#
# Results folder: results/Fig2_WGCNA/03_DME_donor_level_2026-10-13
# F0 is retained in machine-readable stage identifiers and displayed as Healthy.
# Donor IDs in exported tables must be study IDs, not participant names.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(12345)

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
  ext <- tolower(tools::file_ext(path))
  partial <- paste0(path, ".partial.", ext)
  if (file.exists(partial)) stopf("Partial plot already exists: %s", partial)
  device <- if (ext == "pdf") grDevices::cairo_pdf else "png"
  ggplot2::ggsave(partial, plot, device = device, width = width, height = height,
                  units = "in", dpi = 300, bg = "white")
  move_into_place(partial, path)
}

# Pure calculation function: scores must already be aligned to metadata by barcode.
calculate_donor_dme <- function(MEs, stage, donor, modules) {
  expected_stages <- paste0("F", 0:4)
  if (nrow(MEs) != length(stage) || length(stage) != length(donor)) {
    stopf("Score and metadata lengths differ")
  }
  if (anyNA(stage) || anyNA(donor) || any(!nzchar(stage)) || any(!nzchar(donor))) {
    stopf("Stage and donor IDs must be present for every spot")
  }
  if (!identical(sort(unique(stage)), expected_stages)) {
    stopf("Expected stages F0-F4, observed: %s", paste(sort(unique(stage)), collapse = ", "))
  }
  if (!all(vapply(MEs[modules], is.numeric, logical(1)))) stopf("ME columns must be numeric")
  n_stages <- vapply(split(stage, donor), function(x) length(unique(x)), integer(1))
  if (any(n_stages != 1L)) {
    stopf("Donor IDs occur in multiple stages: %s. Check identifiers; this is an unpaired cross-sectional analysis.",
          paste(names(n_stages)[n_stages != 1L], collapse = ", "))
  }
  # All available spots per donor are averaged, including multiple sections if present.
  # This is a tissue-wide donor mean; no equal weighting of sections is introduced.
  donors <- unique(data.frame(stage = stage, donor = donor))
  donors <- donors[order(match(donors$stage, expected_stages), donors$donor), , drop = FALSE]
  rownames(donors) <- NULL
  donor_means <- donors
  donor_qc <- donors
  donor_qc$n_spots <- vapply(donors$donor, function(id) sum(donor == id), integer(1))
  for (m in modules) {
    groups <- lapply(donors$donor, function(id) MEs[[m]][donor == id])
    donor_qc[[paste0(m, "_n_finite_spots")]] <- vapply(groups, function(x) sum(is.finite(x)), integer(1))
    donor_means[[m]] <- vapply(groups, function(x) {
      x <- x[is.finite(x)]
      if (length(x)) mean(x) else NA_real_
    }, numeric(1))
  }
  if (any(!is.finite(as.matrix(MEs[modules])))) {
    warning("Non-finite spot scores excluded within donor; see donor spot-count QC.", call. = FALSE)
  }
  stage_qc <- data.frame(stage = expected_stages)
  stage_qc$n_donor <- vapply(expected_stages, function(s) sum(donors$stage == s), integer(1))
  stage_qc$n_spots <- vapply(expected_stages, function(s) sum(stage == s), integer(1))
  for (m in modules) {
    stage_qc[[paste0(m, "_n_finite_donors")]] <- vapply(expected_stages, function(s) {
      sum(is.finite(donor_means[[m]][donor_means$stage == s]))
    }, integer(1))
  }
  rows <- list()
  for (s in expected_stages[-1L]) {
    for (m in modules) {
      x0 <- donor_means[[m]][donor_means$stage == "F0"]
      x1 <- donor_means[[m]][donor_means$stage == s]
      x0 <- x0[is.finite(x0)]
      x1 <- x1[is.finite(x1)]
      if (length(x0) < 2L || length(x1) < 2L) {
        stopf("%s versus Healthy / %s has fewer than two usable donors in a group", s, m)
      }
      # R's default exact-test selection and continuity correction match the old
      # donor Rmd (exact when eligible, normal approximation with ties).
      test <- stats::wilcox.test(x1, x0, alternative = "two.sided", paired = FALSE,
                                exact = NULL, correct = TRUE)
      if (!is.finite(test$p.value)) stopf("Undefined Wilcoxon P value for %s / %s", s, m)
      rows[[length(rows) + 1L]] <- data.frame(
        comparison = paste0(s, "_vs_F0"), stage = s, baseline = "F0", module = m,
        paired = FALSE, n_baseline_used = length(x0), n_stage_used = length(x1),
        mean_baseline = mean(x0), mean_stage = mean(x1),
        delta_mean = mean(x1) - mean(x0), p_wilcox = test$p.value,
        test_method = test$method)
    }
  }
  results <- do.call(rbind, rows)
  results$p_adj <- ave(results$p_wilcox, results$comparison,
                       FUN = function(x) stats::p.adjust(x, method = "BH"))
  results$significant_FDR_0.05 <- results$p_adj < 0.05
  results$stars <- ifelse(results$p_adj < 0.001, "***",
                         ifelse(results$p_adj < 0.01, "**",
                                ifelse(results$p_adj < 0.05, "*", "")))
  rownames(results) <- NULL
  list(results = results, donor_means = donor_means, donor_qc = donor_qc, stage_qc = stage_qc)
}

make_dme_plots <- function(results, modules) {
  # Final manuscript labels and asterisks, formerly applied by a separate replot.
  module_names <- c(SM_1 = "Fibro\u2013immune ECM\nremodeling",
    SM_2 = "Interferon", SM_3 = "Cholesterol metabolism\nand lipogenesis",
    SM_4 = "Hepatocyte homeostatic\nfunction",
    SM_5 = "Stress and protein\nhomeostasis", SM_6 = "Hepatic plasma")
  if (!all(modules %in% names(module_names))) stop("Unknown module label.")
  module_axis_labels <- paste(rev(modules), unname(module_names[rev(modules)]), sep = "\n")
  comparisons <- paste0(paste0("F", 1:4), "_vs_F0")
  contrast_labels <- setNames(paste0("F", 1:4, " vs Healthy"), comparisons)
  plot_df <- results
  plot_df$comparison <- factor(plot_df$comparison, levels = comparisons)
  plot_df$module <- factor(plot_df$module, levels = rev(modules))
  plot_df$module_label <- gsub("_", "-", as.character(plot_df$module))
  offsets <- c(F1_vs_F0 = -0.24, F2_vs_F0 = -0.08, F3_vs_F0 = 0.08, F4_vs_F0 = 0.24)
  plot_df$y <- as.numeric(plot_df$module) + offsets[as.character(plot_df$comparison)]
  stage_palette <- c(F1_vs_F0 = "#fcbba1", F2_vs_F0 = "#fc9272",
                     F3_vs_F0 = "#fb6a4a", F4_vs_F0 = "#cb181d")
  star_dx <- 0.025 * max(abs(plot_df$delta_mean), 1)
  plot_df$star_x <- plot_df$delta_mean + ifelse(plot_df$delta_mean >= 0, star_dx, -star_dx)
  plot_df$star_hjust <- ifelse(plot_df$delta_mean >= 0, 0, 1)
  p_lollipop <- ggplot(plot_df) +
    geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
    geom_segment(aes(x = 0, xend = delta_mean, y = y, yend = y,
                     colour = comparison), linewidth = 1) +
    geom_point(aes(x = delta_mean, y = y, fill = comparison),
               shape = 21, size = 3.2, stroke = 0.2) +
    geom_text(aes(x = star_x, y = y, label = stars, hjust = star_hjust), size = 7, vjust = 0.78) +
    scale_colour_manual(values = stage_palette, guide = "none") +
    scale_fill_manual(values = stage_palette, name = "Contrast", labels = contrast_labels) +
    scale_y_continuous(breaks = seq_along(modules), labels = module_axis_labels,
                       expand = expansion(add = 0.7)) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(x = expression(Delta*" mean ME (stage - Healthy) at donor level"), y = NULL,
         title = "Donor-level differential module eigengenes",
         caption = "BH-adjusted P: * < 0.05; ** < 0.01; *** < 0.001; unmarked, not significant") +
    theme_minimal(base_family = "sans") +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "plain",
            hjust = 1, vjust = 0.5, lineheight = 0.95, margin = margin(r = 6)),
          axis.title.x = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  # Established FLINC module colours, independent of colours stored in the object.
  module_palette <- c(SM_1 = "#CAB2D6", SM_2 = "#6A3D9A", SM_3 = "#E69F00",
                      SM_4 = "#7EC0EE", SM_5 = "#993F00", SM_6 = "#023FA5")
  plot_df$neglog10_padj <- -log10(pmax(plot_df$p_adj, .Machine$double.xmin))
  plot_df$point_color <- ifelse(plot_df$significant_FDR_0.05,
                               unname(module_palette[as.character(plot_df$module)]), "#A6A6A6")
  # SM_3 and SM_4 almost overlap at F1 in the recovered results. An open SM_4
  # ring distinguishes them without changing their statistical coordinates.
  ring <- plot_df$stage == "F1" & as.character(plot_df$module) == "SM_4"
  p_volcano <- ggplot(plot_df, aes(x = delta_mean, y = neglog10_padj)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(data = plot_df[!ring, , drop = FALSE], aes(colour = point_color), size = 2.8) +
    geom_point(data = plot_df[ring, , drop = FALSE], aes(colour = point_color),
               shape = 1, size = 3.5, stroke = 1) +
    ggrepel::geom_text_repel(aes(label = module_label, colour = point_color),
                             size = 3.5, seed = 12345, min.segment.length = 0,
                             max.overlaps = Inf, box.padding = 0.5,
                             point.padding = 0.3, show.legend = FALSE) +
    scale_colour_identity() +
    facet_wrap(~ comparison, ncol = 2, scales = "free",
               labeller = as_labeller(contrast_labels)) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.16))) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.22))) +
    labs(x = expression(Delta*" mean ME (stage - Healthy) at donor level"),
         y = expression(-log[10](P[adj])),
         title = "Donor-level differential module eigengenes",
         subtitle = "Module colors: adjusted P < 0.05; grey: not significant") +
    theme_classic(base_size = 12, base_family = "sans") +
    theme(strip.background = element_rect(fill = "#CECECE", colour = NA),
          strip.text = element_text(size = 12), panel.spacing = grid::unit(1, "lines"),
          plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 11), axis.title = element_text(size = 14),
          legend.position = "none")
  list(lollipop = p_lollipop, volcano = p_volcano)
}

# ---- Paths and packages: same convention as the original GitHub step 03 ----
repro_root <- Sys.getenv("FLINC_REPRO_ROOT",
  unset = "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code")
input_file <- Sys.getenv("FIG2_NETWORK_INPUT",
  unset = file.path(repro_root, "objects", "Fig2_WGCNA", "02_wgcna_network_modules.rds"))
run_tag <- "03_DME_donor_level_2026-10-13"
results_dir <- Sys.getenv("FIG2_DME_RESULTS_DIR",
  unset = file.path(repro_root, "results", "Fig2_WGCNA", run_tag))
log_dir <- file.path(repro_root, "reproducibility_logs", "Fig2_WGCNA")
stage_col <- Sys.getenv("FIG2_STAGE_COLUMN", unset = "F_score")
donor_col <- Sys.getenv("FIG2_DONOR_COLUMN", unset = "donor")
me_source <- Sys.getenv("FIG2_ME_SOURCE", unset = "harmonized")
wgcna_name <- "vis"
modules <- paste0("SM_", 1:6)
if (!(me_source %in% c("harmonized", "raw", "metadata"))) {
  stopf("FIG2_ME_SOURCE must be harmonized, raw, or metadata")
}
packages <- c("Seurat", "hdWGCNA", "ggplot2", "ggrepel", "writexl")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stopf("Required packages missing: %s", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(ggplot2)
  library(ggrepel)
})
if (!capabilities("cairo")) stopf("R requires Cairo support for the PDF/PNG exports")
if (!file.exists(input_file)) stopf("Step-02 input does not exist: %s", input_file)

csv_file <- file.path(results_dir, "03_DiffModuleEigengenes_donor_level.csv")
tsv_file <- file.path(results_dir, "03_DiffModuleEigengenes_donor_level.tsv")
xlsx_file <- file.path(results_dir, "03_DiffModuleEigengenes_donor_level.xlsx")
donor_file <- file.path(results_dir, "03_donor_mean_module_eigengenes.csv")
donor_qc_file <- file.path(results_dir, "03_QC_spots_per_donor_and_module.csv")
stage_qc_file <- file.path(results_dir, "03_QC_donors_per_stage.csv")
modules_file <- file.path(results_dir, "03_module_assignments_used.csv")
lollipop_file <- file.path(results_dir, "03_DME_donor_level_lollipop_all_stages_vs_Healthy.pdf")
volcano_file <- file.path(results_dir, "03_DME_donor_level_volcano_all_stages_vs_Healthy.pdf")
summary_file <- file.path(log_dir, paste0(run_tag, "_summary.tsv"))
session_file <- file.path(log_dir, paste0(run_tag, "_sessionInfo.txt"))
manifest_file <- file.path(results_dir, "03_output_manifest.tsv")
all_outputs <- c(csv_file, tsv_file, xlsx_file, donor_file, donor_qc_file, stage_qc_file,
                 modules_file, lollipop_file, volcano_file,
                 sub("\\.pdf$", ".png", lollipop_file), sub("\\.pdf$", ".png", volcano_file),
                 summary_file, session_file, manifest_file)
existing <- all_outputs[file.exists(all_outputs)]
if (length(existing)) stopf("Outputs already exist; keep this run or archive it before rerunning: %s",
                            paste(existing, collapse = ", "))

# ---- Read and align scores explicitly by spot barcode ----
say("Loading completed step-02 object:", input_file)
object <- readRDS(input_file)
if (!inherits(object, "Seurat")) stopf("Input is not a Seurat object")
md <- object[[]]
if (!all(c(stage_col, donor_col) %in% colnames(md))) {
  stopf("Input must contain metadata columns '%s' and '%s'. Set FIG2_DONOR_COLUMN only to a verified biological donor ID.",
        stage_col, donor_col)
}
if (!(wgcna_name %in% names(object@misc))) stopf("Input lacks hdWGCNA experiment 'vis'")
say("Using ME source:", me_source, "| donor column:", donor_col)
if (me_source == "metadata") {
  if (!all(modules %in% colnames(md))) stopf("Metadata lacks SM_1-SM_6 score columns")
  MEs <- md[, modules, drop = FALSE]
} else {
  MEs <- GetMEs(object, harmonized = (me_source == "harmonized"), wgcna_name = wgcna_name)
}
if (is.null(MEs)) stopf("Selected ME source is empty; no fallback was applied")
MEs <- as.data.frame(MEs)
found_modules <- grep("^SM_[0-9]+$", colnames(MEs), value = TRUE)
if (!setequal(found_modules, modules)) stopf("Expected exactly SM_1-SM_6 in selected ME source")
spots <- colnames(object)
if (is.null(rownames(MEs)) || anyDuplicated(rownames(MEs)) ||
    !setequal(rownames(MEs), spots) || anyDuplicated(rownames(md)) ||
    !setequal(rownames(md), spots)) stopf("ME/metadata rows must match all object spot barcodes uniquely")
MEs <- MEs[spots, modules, drop = FALSE]
md <- md[spots, , drop = FALSE]
stage <- trimws(as.character(md[[stage_col]]))
stage[!is.na(stage) & stage == "Healthy"] <- "F0"
donor <- trimws(as.character(md[[donor_col]]))
module_table <- as.data.frame(GetModules(object, wgcna_name = wgcna_name))
if (!all(c("gene_name", "module") %in% colnames(module_table))) {
  stopf("Module assignments must contain gene_name and module columns")
}
if (!all(modules %in% as.character(module_table$module))) stopf("Module assignments lack one or more SM_1-SM_6 modules")
module_table <- module_table[as.character(module_table$module) %in% modules, , drop = FALSE]

# ---- Average within donor, then test donors ----
say("Calculating donor means and stage-versus-Healthy tests")
analysis <- calculate_donor_dme(MEs, stage, donor, modules)
results <- analysis$results
results$me_source <- me_source
results$donor_column <- donor_col
say("Donor counts:", paste(analysis$stage_qc$stage, analysis$stage_qc$n_donor,
                          sep = "=", collapse = ", "))
plots <- make_dme_plots(results, modules)

# ---- Export tables, figures and provenance ----
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
say("Writing donor-level results to:", results_dir)
write_csv_atomic(results, csv_file)
write_tsv_atomic(results, tsv_file)
write_csv_atomic(analysis$donor_means, donor_file)
write_csv_atomic(analysis$donor_qc, donor_qc_file)
write_csv_atomic(analysis$stage_qc, stage_qc_file)
write_csv_atomic(module_table, modules_file)
xlsx_partial <- paste0(xlsx_file, ".partial.xlsx")
if (file.exists(xlsx_partial)) stopf("Partial output already exists: %s", xlsx_partial)
writexl::write_xlsx(list(DME = results, donor_means = analysis$donor_means,
                       donor_QC = analysis$donor_qc, stage_QC = analysis$stage_qc), xlsx_partial)
move_into_place(xlsx_partial, xlsx_file)
save_plot_atomic(plots$lollipop, lollipop_file, 11.5, 7)
save_plot_atomic(plots$lollipop, sub("\\.pdf$", ".png", lollipop_file), 11.5, 7)
save_plot_atomic(plots$volcano, volcano_file, 12, 10)
save_plot_atomic(plots$volcano, sub("\\.pdf$", ".png", volcano_file), 12, 10)
summary <- data.frame(
  input_file = normalizePath(input_file), input_size_bytes = file.info(input_file)$size,
  input_md5 = unname(tools::md5sum(input_file)), run_time = as.character(Sys.time()),
  me_source = me_source, stage_column = stage_col, donor_column = donor_col,
  wgcna_name = wgcna_name, n_spots = ncol(object), n_donors = nrow(analysis$donor_means),
  n_modules = length(modules), n_tests = nrow(results),
  n_significant_FDR_0.05 = sum(results$significant_FDR_0.05),
  aggregation = "mean of finite spot MEs within each donor; equal donor weight",
  test = "two-sided unpaired Wilcoxon rank-sum; exact=NULL; correct=TRUE",
  adjustment = "BH across 6 modules separately within each stage vs Healthy",
  random_seed = 12345)
write_tsv_atomic(summary, summary_file)
write_lines_atomic(capture.output(sessionInfo()), session_file)
manifest_paths <- setdiff(all_outputs, manifest_file)
write_tsv_atomic(data.frame(file = normalizePath(manifest_paths),
                            size_bytes = file.info(manifest_paths)$size,
                            md5 = unname(tools::md5sum(manifest_paths))), manifest_file)
say("Significant DME tests at BH FDR < 0.05:", sum(results$significant_FDR_0.05))
say("Fig. 2 donor-level differential module eigengene analysis completed successfully")
