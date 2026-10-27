#!/usr/bin/env Rscript

# FLINC Figure 2 — plot harmonized module eigengenes and cell-type associations:
#   1) harmonized-ME dot plot by the final seven spatial clusters;
#   2) harmonized-ME tissue maps for SM_1 through SM_6;
#   3) harmonized module eigengene correlations with 12 cell-type fractions.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(viridis)
  library(ggrastr)
})

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

write_tsv <- function(x, path) utils::write.table(
  x, path, sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = TRUE, na = "NA"
)
write_csv <- function(x, path, row.names = FALSE) utils::write.csv(
  x, path, row.names = row.names, na = "NA"
)

make_hex_spots <- function(cdata, fill) {
  r <- 1 / 2
  R <- (2 / sqrt(3)) * r
  spot_positions <- select_spot_positions(cdata, fill = fill)
  spot_positions <- adjust_hex_centers(spot_positions)
  vertex_offsets <- data.frame(
    x.offset = c(0, r, r, 0, -r, -r),
    y.offset = c(-R, -R / 2, R / 2, R, R / 2, -R / 2)
  )
  spot_vertices <- make_spot_vertices(spot_positions, vertex_offsets)
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  spot_vertices
}

select_spot_positions <- function(cdata, x = "col", y = "row", fill) {
  if (length(fill) != nrow(cdata)) stopf("fill length does not match coordinate rows")
  spot_positions <- cdata[, c(x, y), drop = FALSE]
  colnames(spot_positions) <- c("x.pos", "y.pos")
  spot_positions$fill <- as.numeric(fill)
  spot_positions$spot <- rownames(spot_positions)
  spot_positions
}

adjust_hex_centers <- function(spot_positions) {
  r <- 1 / 2
  R <- (2 / sqrt(3)) * r
  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos * R * (3 / 2)
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
  spot_positions
}

make_spot_vertices <- function(spot_positions, vertex_offsets) {
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  as.data.frame(spot_vertices)
}

sample_feature_plot <- function(object, feature, sample_col = "F_score_donor",
                                ncol = 4L, plot_max = "q99", plot_min = 0,
                                dpi = 600) {
  metadata <- object[[]]
  needed <- c("row", "col", sample_col, feature)
  missing <- setdiff(needed, colnames(metadata))
  if (length(missing)) stopf("Spatial plot metadata missing: %s", paste(missing, collapse = ", "))
  if (!is.numeric(metadata[[feature]])) stopf("Spatial feature is not numeric: %s", feature)

  # Preserve the sample order stored in the object metadata.
  sample_order <- unique(as.character(metadata[[sample_col]]))

  values <- metadata[[feature]]
  upper <- if (is.character(plot_max)) {
    as.numeric(stats::quantile(values, as.numeric(sub("^q", "", plot_max)) / 100,
                               na.rm = TRUE, names = FALSE))
  } else as.numeric(plot_max)
  lower <- if (is.character(plot_min)) {
    as.numeric(stats::quantile(values, as.numeric(sub("^q", "", plot_min)) / 100,
                               na.rm = TRUE, names = FALSE))
  } else as.numeric(plot_min)
  clipped <- pmin(values, upper)
  metadata$.plot_value <- clipped

  colors <- rev(viridis::inferno(256))
  plots <- lapply(sample_order, function(sample_id) {
    current <- metadata[as.character(metadata[[sample_col]]) == sample_id, , drop = FALSE]
    vertices <- make_hex_spots(current, current$.plot_value)
    positive_vertices <- vertices[vertices$fill > lower, , drop = FALSE]
    background_vertices <- vertices[vertices$fill <= lower, , drop = FALSE]
    ggplot(positive_vertices,
           aes(x = x.vertex, y = y.vertex, group = spot, fill = fill)) +
      ggrastr::rasterise(
        geom_polygon(
          data = background_vertices,
          aes(x = x.vertex, y = y.vertex, group = spot),
          inherit.aes = FALSE, fill = "lightgrey", linewidth = 0
        ),
        dpi = dpi / 2
      ) +
      ggrastr::rasterise(geom_polygon(linewidth = 0), dpi = dpi) +
      scale_fill_gradientn(colors = colors, limits = c(lower, upper), oob = scales::squish) +
      coord_equal() +
      ggtitle(sample_id) +
      labs(fill = feature) +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 8))
  })
  names(plots) <- sample_order
  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = feature) &
    theme(plot.title = element_text(hjust = 0.5))
}

repro_root <- env_or(
  "FLINC_REPRO_ROOT",
  "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)
network_file <- env_or(
  "FIG2_NETWORK_OBJECT_INPUT",
  file.path(repro_root, "objects", "Fig2_WGCNA", "02_wgcna_network_modules.rds")
)
deconvolution_file <- env_or(
  "FIG2_DECONVOLUTION_OBJECT_INPUT",
  paste0(
    "/Users/jianwu/Documents/FLINC/Jeremie/Figures_new/WGCNA/",
    "FLINC6_all_stages/data_harmonized/hdWGCNA.anno_norm.rds"
  )
)
output_dir <- env_or(
  "FIG2_STEP07_OUTPUT_DIR",
  file.path(repro_root, "results", "Fig2_WGCNA", "07_module_expression")
)
partial_dir <- paste0(output_dir, ".partial")
log_dir <- file.path(repro_root, "reproducibility_logs", "Fig2_WGCNA")
wgcna_name <- env_or("FIG2_WGCNA_NAME", "vis")

if (!file.exists(network_file)) stopf("Network object does not exist: %s", network_file)
if (dir.exists(output_dir) || file.exists(output_dir)) stopf("Final output already exists: %s", output_dir)
if (dir.exists(partial_dir) || file.exists(partial_dir)) stopf("Partial output already exists: %s", partial_dir)
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(partial_dir, "dotplots"), recursive = TRUE)
dir.create(file.path(partial_dir, "spatial_maps"), recursive = TRUE)
dir.create(file.path(partial_dir, "correlations"), recursive = TRUE)
dir.create(file.path(partial_dir, "validation"), recursive = TRUE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

modules <- paste0("SM_", 1:6)
cell_types <- c(
  "Bcells", "CentralHep", "Cholangiocyte", "Erythroid", "InfMac",
  "InterHep", "NKTcell", "NonInfMac", "PortalEndo", "PortalHep",
  "Stellate", "cvLSECs"
)
module_labels <- c(
  SM_1 = "Fibro-immune ECM remodeling",
  SM_2 = "Interferon antiviral",
  SM_3 = "Cholesterol metabolism and lipogenesis",
  SM_4 = "Hepatocyte homeostatic function",
  SM_5 = "Stress and protein homeostasis",
  SM_6 = "Hepatic plasma"
)

say("Loading WGCNA object")
object <- readRDS(network_file)
if (!inherits(object, "Seurat")) stopf("Network input is not a Seurat object")
if (!(wgcna_name %in% names(object@misc))) stopf("Missing hdWGCNA experiment: %s", wgcna_name)
needed_metadata <- c("cluster_named", "F_score", "F_score_donor", "row", "col", "barcode")
missing_metadata <- setdiff(needed_metadata, colnames(object[[]]))
if (length(missing_metadata)) stopf("Network metadata missing: %s", paste(missing_metadata, collapse = ", "))
if (ncol(object) != 26056L) stopf("Expected 26056 spots; observed %d", ncol(object))
if (length(unique(object$F_score_donor)) != 24L) stopf("Expected 24 samples")

say("Retrieving harmonized module eigengenes")
hmes <- as.data.frame(GetMEs(object, harmonized = TRUE, wgcna_name = wgcna_name), check.names = FALSE)
if (!setequal(rownames(hmes), colnames(object))) stopf("hME rows do not match object spots")
hmes <- hmes[colnames(object), , drop = FALSE]
if (!all(modules %in% colnames(hmes))) stopf("hMEs are missing SM_1 through SM_6")
hmes <- hmes[, modules, drop = FALSE]
if (any(!is.finite(as.matrix(hmes)))) stopf("hMEs contain non-finite values")
for (module in modules) object[[module]] <- hmes[[module]]

cluster_levels <- unique(as.character(object$cluster_named))
expected_clusters <- c(
  "C0 Portal tracts", "C1 Periportal hepatocytes",
  "C2 Centrilobular hepatocytes", "C3 Immune cells", "C4 Myofibroblasts",
  "C5 Cholangiocytes", "C6 Matrix-rich fibroblast"
)
if (!setequal(cluster_levels, expected_clusters)) {
  stopf("Unexpected cluster_named values: %s", paste(sort(cluster_levels), collapse = "; "))
}
object$cluster_named <- factor(as.character(object$cluster_named), levels = expected_clusters)
Idents(object) <- "cluster_named"

say("Writing all-module dot plot by final spatial cluster")
dotplot <- DotPlot(object, features = rev(modules), group.by = "cluster_named", dot.min = 0.1) +
  coord_flip(clip = "off") +
  scale_color_gradient2(low = "blue", mid = "grey95", high = "red", midpoint = 0) +
  labs(x = NULL, y = NULL, color = "Average expression", size = "Percent expressed") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        plot.margin = margin(t = 8, r = 8, b = 95, l = 8))
ggsave(file.path(partial_dir, "dotplots", "07_all_hMEs_by_cluster.pdf"), dotplot,
       device = cairo_pdf, width = 11, height = 8.2, units = "in")
ggsave(file.path(partial_dir, "dotplots", "07_all_hMEs_by_cluster.png"), dotplot,
       width = 11, height = 8.2, units = "in", dpi = 300)

say("Writing SM_1 through SM_6 tissue maps")
spatial_plots <- list()
for (module in modules) {
  say("Spatial map:", module)
  plot <- sample_feature_plot(object, module, ncol = 4L, plot_max = "q99", plot_min = 0)
  spatial_plots[[module]] <- plot
  ggsave(file.path(partial_dir, "spatial_maps", paste0("07_", module, "_hME_spatial_24_samples.pdf")),
         plot, device = cairo_pdf, width = 20, height = 30, units = "in")
  ggsave(file.path(partial_dir, "spatial_maps", paste0("07_", module, "_hME_spatial_24_samples.png")),
         plot, width = 20, height = 30, units = "in", dpi = 300)
}
pdf(file.path(partial_dir, "spatial_maps", "07_all_modules_hME_spatial_24_samples.pdf"),
    width = 20, height = 30, onefile = TRUE)
for (module in modules) print(spatial_plots[[module]])
dev.off()

say("Preparing the 12 deconvolution fractions")
missing_cell_types <- setdiff(cell_types, colnames(object[[]]))
deconv_source <- "network object metadata"
if (length(missing_cell_types)) {
  if (!file.exists(deconvolution_file)) {
    stopf("Cell-type fractions are missing and deconvolution object does not exist: %s", deconvolution_file)
  }
  deconv <- readRDS(deconvolution_file)
  if (!inherits(deconv, "Seurat")) stopf("Deconvolution input is not a Seurat object")
  deconv_metadata <- deconv[[]]
  deconv_needed <- c("barcode", cell_types)
  deconv_missing <- setdiff(deconv_needed, colnames(deconv_metadata))
  if (length(deconv_missing)) stopf("Deconvolution metadata missing: %s", paste(deconv_missing, collapse = ", "))
  if (anyDuplicated(as.character(object$barcode))) stopf("Network object has duplicated barcode values")
  if (anyDuplicated(as.character(deconv_metadata$barcode))) stopf("Deconvolution object has duplicated barcode values")
  match_index <- match(as.character(object$barcode), as.character(deconv_metadata$barcode))
  if (anyNA(match_index)) stopf("Deconvolution join left %d unmatched spots", sum(is.na(match_index)))
  fractions <- deconv_metadata[match_index, cell_types, drop = FALSE]
  rownames(fractions) <- colnames(object)
  object <- AddMetaData(object, metadata = fractions)
  deconv_source <- deconvolution_file
  rm(deconv, deconv_metadata, fractions)
  invisible(gc())
}
fraction_matrix <- as.matrix(object[[]][, cell_types, drop = FALSE])
storage.mode(fraction_matrix) <- "double"
if (any(!is.finite(fraction_matrix))) stopf("Cell-type fractions contain non-finite values")
if (any(fraction_matrix < 0)) stopf("Cell-type fractions contain negative values")

say("Preparing harmonized module eigengenes for cell-type correlations")
hME_output <- data.frame(barcode = object$barcode, hmes,
                         row.names = NULL, check.names = FALSE)
write_tsv(hME_output,
          file.path(partial_dir, "correlations", "07_harmonized_module_eigengenes.tsv"))

say("Computing 6 x 12 spot-level Spearman correlations")
cor_parts <- lapply(modules, function(module) {
  do.call(rbind, lapply(cell_types, function(cell_type) {
    x <- hmes[[module]]
    y <- fraction_matrix[, cell_type]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
      return(data.frame(module = module, cell_type = cell_type, rho = NA_real_,
                        p_value = NA_real_, n = sum(ok)))
    }
    test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.frame(module = module, cell_type = cell_type,
               rho = unname(test$estimate), p_value = test$p.value, n = sum(ok))
  }))
})
correlations <- do.call(rbind, cor_parts)
correlations$p_adj_global <- p.adjust(correlations$p_value, method = "BH")
correlations$p_adj_within_module <- ave(
  correlations$p_value, correlations$module,
  FUN = function(x) p.adjust(x, method = "BH")
)
write_csv(correlations,
          file.path(partial_dir, "correlations", "07_hME_celltype_spearman_long.csv"))

cor_matrix <- matrix(NA_real_, nrow = length(modules), ncol = length(cell_types),
                     dimnames = list(modules, cell_types))
for (i in seq_len(nrow(correlations))) {
  cor_matrix[correlations$module[i], correlations$cell_type[i]] <- correlations$rho[i]
}
write_csv(data.frame(module = rownames(cor_matrix), cor_matrix, check.names = FALSE),
          file.path(partial_dir, "correlations", "07_hME_celltype_spearman_matrix.csv"))

drivers <- do.call(rbind, lapply(modules, function(module) {
  part <- correlations[correlations$module == module & is.finite(correlations$rho) & correlations$rho > 0, ]
  if (!nrow(part)) return(data.frame(module = module, driver = NA_character_, rho = NA_real_, p_adj = NA_real_))
  best <- part[which.max(part$rho), , drop = FALSE]
  data.frame(module = module, driver = best$cell_type, rho = best$rho,
             p_adj = best$p_adj_within_module)
}))
write_csv(drivers,
          file.path(partial_dir, "correlations", "07_strongest_positive_driver_by_module.csv"))

row_annotation <- data.frame(Driver = drivers$driver, Rho = drivers$rho,
                             row.names = drivers$module)
row_labels <- paste0(modules, ":\n", unname(module_labels[modules]))
pdf(file.path(partial_dir, "correlations", "07_hME_celltype_association_heatmap.pdf"),
    width = 10, height = 6)
pheatmap(
  cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  annotation_row = row_annotation,
  labels_row = row_labels,
  main = "Harmonized module eigengene-cell-type associations\nSpearman rho across spatial spots"
)
dev.off()

summary <- data.frame(
  metric = c(
    "network_input", "deconvolution_source", "n_spots", "n_samples",
    "n_clusters", "n_modules", "n_cell_types", "module_activity_for_plots",
    "module_activity_for_correlations",
    "correlation_method", "n_correlations"
  ),
  value = c(
    network_file, deconv_source, ncol(object), length(unique(object$F_score_donor)),
    length(levels(object$cluster_named)), length(modules), length(cell_types),
    "harmonized module eigengenes", "harmonized module eigengenes",
    "Spearman spot-level", nrow(correlations)
  )
)
write_tsv(summary, file.path(partial_dir, "validation", "07_summary.tsv"))
writeLines(capture.output(sessionInfo()),
           file.path(partial_dir, "validation", "07_sessionInfo.txt"))

checks <- c(
  STEP07_HME_MODULES_PASS = identical(colnames(hmes), modules),
  STEP07_SEVEN_CLUSTERS_PASS = identical(levels(object$cluster_named), expected_clusters),
  STEP07_24_SAMPLES_PASS = length(unique(object$F_score_donor)) == 24L,
  STEP07_ALL_SPATIAL_MAPS_PASS = length(spatial_plots) == 6L,
  STEP07_12_CELLTYPES_PASS = ncol(fraction_matrix) == 12L,
  STEP07_HME_CORRELATION_INPUT_PASS = identical(colnames(hmes), modules),
  STEP07_72_CORRELATIONS_PASS = nrow(correlations) == 72L,
  STEP07_PLOTS_PASS = all(file.exists(c(
    file.path(partial_dir, "dotplots", "07_all_hMEs_by_cluster.pdf"),
    file.path(partial_dir, "spatial_maps", "07_all_modules_hME_spatial_24_samples.pdf"),
    file.path(partial_dir, "correlations", "07_hME_celltype_association_heatmap.pdf")
  )))
)
write_tsv(data.frame(check = names(checks), pass = unname(checks)),
          file.path(partial_dir, "validation", "07_validation.tsv"))
if (!all(checks)) stopf("At least one Step 07 validation failed")

if (!file.rename(partial_dir, output_dir)) stopf("Could not finalize output directory: %s", output_dir)
for (name in names(checks)) cat(name, ": ", checks[[name]], "\n", sep = "")
say("Step 07 completed successfully")
say("Output:", output_dir)
