#!/usr/bin/env Rscript

# FLINC Figure 2 — regenerate two supplementary hdWGCNA UMAP outputs.
#
# Output 1: hubgene_umap_ggplot.pdf
#   Each point is one non-grey WGCNA gene, colored by module and sized by its
#   module-specific kME. No individual hub genes are circled or labelled.
#
# Output 2: umap.MEs.pdf
#   Harmonized module eigengenes (hMEs) for SM_1 through SM_6 plotted on the
#   retained annotated-spot UMAP. The standard, unharmonized MEs are not used.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
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

check_exact_version <- function(package, expected) {
  observed <- as.character(packageVersion(package))
  if (!identical(observed, expected)) {
    stopf(
      "Package version mismatch for %s: expected %s, observed %s",
      package, expected, observed
    )
  }
}

write_table_file <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

repro_root <- env_or(
  "FLINC_REPRO_ROOT",
  "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)
input_file <- env_or(
  "FIG2_NETWORK_OBJECT_INPUT",
  file.path(repro_root, "objects", "Fig2_WGCNA", "02_wgcna_network_modules.rds")
)
output_dir <- env_or(
  "FIG2_SUPPLEMENTARY_UMAP_OUTPUT_DIR",
  file.path(repro_root, "results", "Fig2_WGCNA", "07_supplementary_wgcna_umaps")
)
partial_dir <- paste0(output_dir, ".partial")

if (!file.exists(input_file) || dir.exists(input_file)) {
  stopf("Completed Figure 2 network object is missing: %s", input_file)
}
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stopf("Refusing to overwrite existing supplementary UMAP output: %s", output_dir)
}
if (file.exists(partial_dir) || dir.exists(partial_dir)) {
  stopf("Partial supplementary UMAP output exists; inspect it first: %s", partial_dir)
}
if (!dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)) {
  stopf("Could not create partial output directory: %s", partial_dir)
}

if (!identical(paste(R.version$major, R.version$minor, sep = "."), "4.4.1")) {
  stopf("R version mismatch: expected 4.4.1, observed %s", R.version.string)
}
check_exact_version("Seurat", "5.3.0")
check_exact_version("SeuratObject", "5.2.0")
check_exact_version("hdWGCNA", "0.4.6")

wgcna_name <- "vis"
modules <- paste0("SM_", 1:6)
expected_module_sizes <- c(
  SM_1 = 1652L,
  SM_2 = 50L,
  SM_3 = 69L,
  SM_4 = 90L,
  SM_5 = 53L,
  SM_6 = 292L
)
expected_gene_count <- sum(expected_module_sizes)
umap_seed <- 3210L
spot_reduction <- "umap_SCT"

say("Loading completed WGCNA object:", input_file)
object <- readRDS(input_file)
if (!inherits(object, "Seurat")) {
  stopf("Input is not a Seurat object")
}
if (!(wgcna_name %in% names(object@misc))) {
  stopf("Input does not contain the hdWGCNA experiment '%s'", wgcna_name)
}
if (!(spot_reduction %in% Reductions(object))) {
  stopf(
    "The spot-level reduction required for umap.MEs.pdf is missing: %s",
    spot_reduction
  )
}

module_table <- GetModules(object, wgcna_name = wgcna_name)
if (!all(c("gene_name", "module") %in% colnames(module_table))) {
  stopf("GetModules output is missing gene_name or module")
}
observed_module_sizes <- table(as.character(module_table$module))
if (!all(modules %in% names(observed_module_sizes)) ||
    !identical(
      as.integer(observed_module_sizes[modules]),
      as.integer(expected_module_sizes[modules])
    )) {
  stopf("The six retained module sizes do not match the validated network")
}

hmes <- as.data.frame(
  GetMEs(object, harmonized = TRUE, wgcna_name = wgcna_name),
  check.names = FALSE
)
if (!all(modules %in% colnames(hmes)) ||
    nrow(hmes) != ncol(object) ||
    any(!is.finite(as.matrix(hmes[, modules, drop = FALSE])))) {
  stopf("The harmonized module eigengene matrix failed validation")
}

network_data <- GetNetworkData(object, wgcna_name = wgcna_name)
tom_files <- unique(as.character(unlist(network_data$TOMFiles)))
if (!length(tom_files) || any(!file.exists(tom_files))) {
  stopf(
    "The network TOM required by RunModuleUMAP is missing. Recorded paths: %s",
    paste(tom_files, collapse = ", ")
  )
}
say("Using TOM:", paste(tom_files, collapse = ", "))

say("Running module-gene UMAP with historical parameters and seed 3210")
set.seed(umap_seed)
object_with_umap <- RunModuleUMAP(
  object,
  n_hubs = 5,
  exclude_grey = TRUE,
  wgcna_name = wgcna_name,
  n_neighbors = 15,
  min_dist = 0.3,
  spread = 2,
  seed = umap_seed
)

umap_df <- as.data.frame(
  GetModuleUMAP(object_with_umap, wgcna_name = wgcna_name),
  stringsAsFactors = FALSE
)
required_umap_columns <- c("gene", "module", "color", "kME", "UMAP1", "UMAP2")
missing_umap_columns <- setdiff(required_umap_columns, colnames(umap_df))
if (length(missing_umap_columns)) {
  stopf(
    "GetModuleUMAP output is missing columns: %s",
    paste(missing_umap_columns, collapse = ", ")
  )
}

umap_df$module <- as.character(umap_df$module)
umap_df$gene <- as.character(umap_df$gene)
if (nrow(umap_df) != expected_gene_count ||
    !setequal(unique(umap_df$module), modules) ||
    anyDuplicated(umap_df$gene) ||
    any(!is.finite(umap_df$UMAP1)) ||
    any(!is.finite(umap_df$UMAP2)) ||
    any(!is.finite(umap_df$kME))) {
  stopf("The module-gene UMAP table failed structural validation")
}

# Reproduce hubgene_umap_ggplot.pdf without individual hub circles or labels.
point_sizes <- pmax(umap_df$kME, 0) * 3
module_gene_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  ggrastr::rasterise(
    geom_point(
      color = umap_df$color,
      size = point_sizes,
      alpha = 0.95
    ),
    dpi = 800,
    dpi_scale = 0.5
  ) +
  coord_equal() +
  theme_void(base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

say("Constructing harmonized-ME feature maps on", spot_reduction)
hme_plot_list <- ModuleFeaturePlot(
  object,
  module_names = modules,
  wgcna_name = wgcna_name,
  reduction = spot_reduction,
  features = "hMEs",
  order = TRUE,
  point_size = 0.2
)
if (!is.list(hme_plot_list) || length(hme_plot_list) != 6L) {
  stopf("ModuleFeaturePlot did not return six harmonized-ME panels")
}
hme_plot <- patchwork::wrap_plots(hme_plot_list, ncol = 3)

coordinates_file <- file.path(partial_dir, "07_module_gene_umap_coordinates.tsv")
parameters_file <- file.path(partial_dir, "07_supplementary_wgcna_umap_parameters.tsv")
module_gene_pdf <- file.path(partial_dir, "hubgene_umap_ggplot.pdf")
module_gene_png <- file.path(partial_dir, "hubgene_umap_ggplot.png")
hme_pdf <- file.path(partial_dir, "umap.MEs.pdf")
hme_png <- file.path(partial_dir, "umap.MEs.png")

write_table_file(umap_df, coordinates_file)
write_table_file(
  data.frame(
    parameter = c(
      "input_file", "wgcna_name", "n_hubs_for_embedding",
      "n_neighbors", "min_dist", "spread", "seed", "exclude_grey",
      "n_embedded_genes", "gene_labels", "hub_circles",
      "ME_feature_type", "ME_reduction", "ME_point_size",
      "R_version", "Seurat_version", "SeuratObject_version",
      "hdWGCNA_version", "ggplot2_version", "ggrastr_version",
      "patchwork_version"
    ),
    value = c(
      input_file, wgcna_name, "5", "15", "0.3", "2", as.character(umap_seed),
      "TRUE", as.character(nrow(umap_df)), "none", "none",
      "harmonized module eigengenes (hMEs)", spot_reduction, "0.2",
      R.version.string,
      as.character(packageVersion("Seurat")),
      as.character(packageVersion("SeuratObject")),
      as.character(packageVersion("hdWGCNA")),
      as.character(packageVersion("ggplot2")),
      as.character(packageVersion("ggrastr")),
      as.character(packageVersion("patchwork"))
    ),
    stringsAsFactors = FALSE
  ),
  parameters_file
)

ggsave(
  filename = module_gene_pdf,
  plot = module_gene_plot,
  device = grDevices::cairo_pdf,
  width = 8,
  height = 8,
  units = "in",
  bg = "white"
)
ggsave(
  filename = module_gene_png,
  plot = module_gene_plot,
  device = "png",
  width = 8,
  height = 8,
  units = "in",
  dpi = 400,
  bg = "white"
)
ggsave(
  filename = hme_pdf,
  plot = hme_plot,
  device = grDevices::cairo_pdf,
  width = 18,
  height = 8,
  units = "in",
  bg = "white"
)
ggsave(
  filename = hme_png,
  plot = hme_plot,
  device = "png",
  width = 18,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(partial_dir, "07_supplementary_wgcna_umaps_sessionInfo.txt"),
  useBytes = TRUE
)
writeLines(
  c(
    "FLINC Figure 2 operation 07: supplementary hdWGCNA UMAP outputs",
    paste0("Completed: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("Input network object: ", input_file),
    paste0("TOM: ", paste(tom_files, collapse = ", ")),
    "",
    "hubgene_umap_ggplot.pdf:",
    "  RunModuleUMAP parameters: n_hubs=5; n_neighbors=15; min_dist=0.3; spread=2; seed=3210",
    paste0("  Embedded non-grey genes: ", nrow(umap_df)),
    "  Individual gene labels: none",
    "  Hub circles: none",
    "",
    "umap.MEs.pdf:",
    "  ModuleFeaturePlot features: hMEs",
    "  Eigengene source: GetMEs(..., harmonized=TRUE)",
    paste0("  Reduction: ", spot_reduction),
    "  Layout: 3 columns",
    "",
    "The input RDS was not overwritten; the transient UMAP-bearing object was not saved."
  ),
  file.path(partial_dir, "07_supplementary_wgcna_umaps_completion_report.txt"),
  useBytes = TRUE
)

required_outputs <- list.files(partial_dir, full.names = TRUE)
if (length(required_outputs) != 8L ||
    any(file.info(required_outputs)$size <= 0)) {
  stopf("The supplementary WGCNA UMAP output set is incomplete")
}

if (!file.rename(partial_dir, output_dir)) {
  stopf("Could not move completed supplementary UMAP output into place: %s", output_dir)
}

say("Supplementary WGCNA UMAP outputs completed successfully")
say("Output:", output_dir)
say("Module-gene UMAP genes:", nrow(umap_df))
say("Module-gene UMAP labels: none")
say("Harmonized-ME panels:", length(hme_plot_list))
