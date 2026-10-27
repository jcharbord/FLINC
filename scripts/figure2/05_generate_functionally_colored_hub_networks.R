#!/usr/bin/env Rscript

# FLINC Figure 2 — generate functionally colored hub-gene networks.
# Select hub genes with hdWGCNA::GetHubGenes, attach the accepted
# module-specific functional classification, and draw the six module networks.

options(stringsAsFactors = FALSE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
say <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
read_tsv <- function(path) utils::read.delim(
  path, sep = "\t", header = TRUE, quote = "", comment.char = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
write_tsv <- function(x, path) utils::write.table(
  x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
key_for <- function(module, gene) paste(module, gene, sep = "::")

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(igraph)
})

repro_root <- env_or(
  "FLINC_REPRO_ROOT",
  "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)
network_object_file <- env_or(
  "FIG2_NETWORK_OBJECT_INPUT",
  file.path(repro_root, "objects", "Fig2_WGCNA", "02_wgcna_network_modules.rds")
)
annotation_file <- env_or(
  "FIG2_FUNCTIONAL_ANNOTATION_INPUT",
  file.path(
    repro_root, "scripts", "figure2", "reference",
    "functional_subcategories.tsv"
  )
)
output_dir <- env_or(
  "FIG2_STEP06_OUTPUT_DIR",
  file.path(repro_root, "results", "Fig2_WGCNA", "06_colored_hub_networks")
)
partial_dir <- paste0(output_dir, ".partial")
wgcna_name <- env_or("FIG2_WGCNA_NAME", "vis")

for (path in c(network_object_file, annotation_file)) {
  if (!file.exists(path)) stopf("Required input does not exist: %s", path)
}
if (dir.exists(output_dir) || file.exists(output_dir)) {
  stopf("Final output already exists; do not overwrite: %s", output_dir)
}
if (dir.exists(partial_dir) || file.exists(partial_dir)) {
  stopf("Partial output exists; inspect it before retrying: %s", partial_dir)
}

modules <- paste0("SM_", 1:6)
expected_n <- c(SM_1 = 50L, SM_2 = 30L, SM_3 = 30L,
                SM_4 = 30L, SM_5 = 30L, SM_6 = 30L)
inner_n <- c(SM_1 = 20L, SM_2 = 10L, SM_3 = 10L,
             SM_4 = 10L, SM_5 = 10L, SM_6 = 10L)
panels <- c(SM_1 = "C", SM_2 = "D", SM_3 = "E",
            SM_4 = "F", SM_5 = "G", SM_6 = "H")
module_titles <- c(
  SM_1 = "Fibro-immune ECM remodeling",
  SM_2 = "Interferon antiviral",
  SM_3 = "Cholesterol metabolism and lipogenesis",
  SM_4 = "Hepatocyte homeostatic function",
  SM_5 = "Stress and protein homeostasis",
  SM_6 = "Hepatic plasma"
)
palette <- c(
  `1` = "#1432F3", `2` = "#4FAC5B", `3` = "#F4C243",
  `4` = "#FEFE55", `5` = "#E93323", `6` = "#69349B"
)
allowed_categories <- list(
  SM_1 = c("ECM", "Cytoskeleton", "Immune process", "Oxidative stress",
           "Secreted lipid mediator", "Regulation"),
  SM_2 = c("ISGs", "OAS", "Restriction", "ISGylation", "IFN signaling", "Stress"),
  SM_3 = c("Cholesterol synthesis", "Sterol processing", "Sterol regulation",
           "Mevalonate pathway", "Lipogenesis", "Hepatocyte genes"),
  SM_4 = c("Secreted", "Transporters", "Lipid", "Amino acid/urea",
           "Small-molecule", "Regulatory"),
  SM_5 = c("Early stress", "Membrane signaling", "Vesicle/lysosome",
           "Ubiquitin/proteostasis", "Cytoskeletal remodeling",
           "Transcription/metabolism"),
  SM_6 = c("Coagulation factors", "Complement system", "Acute-phase response",
           "Plasma carriers", "Antioxidant/heme binding", "Metabolism/transport")
)

say("Loading WGCNA object")
object <- readRDS(network_object_file)

say("Selecting hubs with GetHubGenes")
hub50 <- GetHubGenes(object, n_hubs = 50, wgcna_name = wgcna_name)
hub30 <- GetHubGenes(object, n_hubs = 30, wgcna_name = wgcna_name)
hub_parts <- lapply(modules, function(module) {
  source <- if (module == "SM_1") hub50 else hub30
  part <- source[as.character(source$module) == module, , drop = FALSE]
  if (nrow(part) != expected_n[[module]]) {
    stopf("%s: expected %d hubs, observed %d", module, expected_n[[module]], nrow(part))
  }
  part$module <- as.character(part$module)
  part$hub_rank <- seq_len(nrow(part))
  part
})
hubs <- do.call(rbind, hub_parts)
rownames(hubs) <- NULL
if (!all(c("module", "gene_name") %in% colnames(hubs))) {
  stopf("GetHubGenes output must contain module and gene_name")
}
hub_keys <- key_for(hubs$module, hubs$gene_name)
if (length(hub_keys) != 200L || anyDuplicated(hub_keys)) {
  stopf("The selected hub table must contain 200 unique module/gene pairs")
}

say("Validating the functional classification")
annotation <- read_tsv(annotation_file)
required_cols <- c(
  "figure_panel", "module", "gene_name", "functional_subcategory",
  "category_order", "category_color_hex"
)
missing_cols <- setdiff(required_cols, colnames(annotation))
if (length(missing_cols)) {
  stopf("Annotation is missing columns: %s", paste(missing_cols, collapse = ", "))
}
annotation$module <- as.character(annotation$module)
annotation$gene_name <- as.character(annotation$gene_name)
annotation$figure_panel <- as.character(annotation$figure_panel)
annotation$functional_subcategory <- as.character(annotation$functional_subcategory)
annotation$category_order <- as.integer(annotation$category_order)
annotation$category_color_hex <- toupper(as.character(annotation$category_color_hex))
annotation_keys <- key_for(annotation$module, annotation$gene_name)

if (nrow(annotation) != 200L || anyDuplicated(annotation_keys)) {
  stopf("Annotation must contain 200 unique module/gene pairs")
}
if (!setequal(annotation_keys, hub_keys)) {
  missing_annotation <- setdiff(hub_keys, annotation_keys)
  extra_annotation <- setdiff(annotation_keys, hub_keys)
  cat("Hubs missing annotation:\n", paste(missing_annotation, collapse = "\n"), "\n")
  cat("Annotations not selected as hubs:\n", paste(extra_annotation, collapse = "\n"), "\n")
  stopf("Functional annotation does not exactly match the current GetHubGenes selection")
}

for (module in modules) {
  part <- annotation[annotation$module == module, , drop = FALSE]
  if (nrow(part) != expected_n[[module]]) stopf("Wrong annotation count for %s", module)
  if (!all(part$figure_panel == panels[[module]])) stopf("Wrong panel for %s", module)
  expected_category <- allowed_categories[[module]][part$category_order]
  if (anyNA(expected_category) || !identical(part$functional_subcategory, expected_category)) {
    stopf("Category name/order mismatch in %s", module)
  }
}
expected_colors <- unname(palette[as.character(annotation$category_order)])
if (anyNA(expected_colors) || !identical(annotation$category_color_hex, expected_colors)) {
  stopf("Category colours do not match the six-colour palette")
}

annotation <- annotation[match(hub_keys, annotation_keys), , drop = FALSE]
annotated_hubs <- cbind(
  hubs,
  annotation[, c("figure_panel", "functional_subcategory",
                 "category_order", "category_color_hex"), drop = FALSE]
)

say("Reading the topological overlap matrix")
tom <- GetTOM(object, wgcna_name = wgcna_name)
if (is.null(colnames(tom)) || is.null(rownames(tom))) stopf("TOM has no gene names")
if (!all(hubs$gene_name %in% colnames(tom))) {
  stopf("Not all selected hubs occur in the TOM")
}

dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(partial_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

draw_module <- function(module, include_legend = TRUE, compact = FALSE) {
  info <- annotated_hubs[annotated_hubs$module == module, , drop = FALSE]
  genes <- info$gene_name
  reduced <- as.matrix(tom[genes, genes, drop = FALSE])
  reduced[!is.finite(reduced)] <- 0
  diag(reduced) <- 0
  graph <- igraph::graph_from_adjacency_matrix(
    reduced, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  weights <- igraph::E(graph)$weight
  scaled <- if (!length(weights) || diff(range(weights)) == 0) {
    rep(0.5, length(weights))
  } else {
    (weights - min(weights)) / diff(range(weights))
  }
  # SM_1 contains 50 hubs (versus 30 in SM_2-SM_6), producing a much denser
  # graph. Use slightly lighter edges for SM_1 so its visual density is
  # comparable while preserving all TOM weights, edges, and geometry.
  edge_alpha <- if (module == "SM_1") {
    0.05 + 0.20 * scaled
  } else {
    0.08 + 0.32 * scaled
  }

  n_inner <- inner_n[[module]]
  n_total <- nrow(info)
  inner_graph <- igraph::make_ring(n_inner)
  outer_graph <- igraph::make_ring(n_total - n_inner)
  layout <- rbind(
    igraph::layout_in_circle(inner_graph) * 0.48,
    igraph::layout_in_circle(outer_graph)
  )
  rownames(layout) <- genes
  layout <- layout[igraph::V(graph)$name, , drop = FALSE]

  vertex_index <- match(igraph::V(graph)$name, info$gene_name)
  node_colors <- info$category_color_hex[vertex_index]
  node_labels <- info$gene_name[vertex_index]

  old_xpd <- par(xpd = NA)
  on.exit(par(xpd = old_xpd), add = TRUE)
  plot(
    graph,
    layout = layout,
    rescale = TRUE,
    edge.color = vapply(edge_alpha, function(alpha) grDevices::adjustcolor("black", alpha.f = alpha), character(1)),
    edge.width = if (compact) 0.45 else 0.65,
    edge.curved = 0,
    vertex.color = node_colors,
    vertex.frame.color = NA,
    vertex.size = if (compact) 4.8 else 5.6,
    vertex.label = node_labels,
    vertex.label.color = "black",
    vertex.label.family = "Helvetica",
    vertex.label.font = 3,
    vertex.label.cex = if (compact) 0.53 else 0.72,
    vertex.label.dist = if (compact) 0.72 else 0.82,
    vertex.label.degree = -pi / 4,
    main = sprintf("%s. %s\n%s", panels[[module]], module, module_titles[[module]])
  )

  if (include_legend) {
    legend(
      "bottom", inset = c(0, if (compact) -0.075 else -0.085),
      horiz = FALSE, ncol = 2, bty = "n",
      legend = allowed_categories[[module]],
      col = unname(palette), pch = 16, pt.cex = if (compact) 1.25 else 1.4,
      cex = if (compact) 0.58 else 0.76, xpd = NA
    )
  }
}

say("Writing individual PDF and PNG panels")
for (module in modules) {
  pdf_file <- file.path(plot_dir, sprintf("06_panel_%s_%s_colored_hub_network.pdf",
                                         panels[[module]], module))
  grDevices::pdf(pdf_file, width = 10, height = 10, useDingbats = FALSE)
  par(mar = c(8.4, 1, 3.5, 1))
  draw_module(module, include_legend = TRUE, compact = FALSE)
  grDevices::dev.off()

  png_file <- sub("\\.pdf$", ".png", pdf_file)
  grDevices::png(png_file, width = 3000, height = 3000, res = 300)
  par(mar = c(8.4, 1, 3.5, 1))
  draw_module(module, include_legend = TRUE, compact = FALSE)
  grDevices::dev.off()
}

say("Writing the combined six-panel figure")
combined_pdf <- file.path(plot_dir, "06_panels_C_to_H_colored_hub_networks.pdf")
grDevices::pdf(combined_pdf, width = 18, height = 17, useDingbats = FALSE)
par(mfrow = c(2, 3), mar = c(8.2, 0.8, 3.0, 0.8), oma = c(0.5, 0.5, 1, 0.5))
for (module in modules) draw_module(module, include_legend = TRUE, compact = TRUE)
grDevices::dev.off()

combined_png <- sub("\\.pdf$", ".png", combined_pdf)
grDevices::png(combined_png, width = 5400, height = 5100, res = 300)
par(mfrow = c(2, 3), mar = c(8.2, 0.8, 3.0, 0.8), oma = c(0.5, 0.5, 1, 0.5))
for (module in modules) draw_module(module, include_legend = TRUE, compact = TRUE)
grDevices::dev.off()

legend_table <- do.call(rbind, lapply(modules, function(module) data.frame(
  figure_panel = panels[[module]],
  module = module,
  biological_name = module_titles[[module]],
  category_order = 1:6,
  functional_subcategory = allowed_categories[[module]],
  category_color_hex = unname(palette),
  stringsAsFactors = FALSE
)))
category_counts <- as.data.frame(table(
  module = factor(annotated_hubs$module, levels = modules),
  functional_subcategory = annotated_hubs$functional_subcategory
), stringsAsFactors = FALSE)
category_counts <- category_counts[category_counts$Freq > 0, , drop = FALSE]
names(category_counts)[3] <- "n_genes"

write_tsv(annotated_hubs, file.path(partial_dir, "06_hubs_with_functional_subcategories.tsv"))
write_tsv(legend_table, file.path(partial_dir, "06_functional_subcategory_legend.tsv"))
write_tsv(category_counts, file.path(partial_dir, "06_functional_subcategory_counts.tsv"))
writeLines(capture.output(sessionInfo()), file.path(partial_dir, "06_sessionInfo.txt"))
writeLines(c(
  "FLINC Fig. 2 Step 06: functionally coloured hub-gene networks",
  paste0("Completed: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "Hub selection: hdWGCNA::GetHubGenes",
  "SM_1: 50 genes (20 inner, 30 outer)",
  "SM_2-SM_6: 30 genes each (10 inner, 20 outer)",
  "Functional classification: module-specific categories",
  "Network weights: topological overlap matrix",
  "",
  "STEP06_GETHUBGENES_COUNTS_PASS: TRUE",
  "STEP06_EXACT_200_GENE_ANNOTATION_PASS: TRUE",
  "STEP06_FUNCTIONAL_CATEGORY_PASS: TRUE",
  "STEP06_CATEGORY_PALETTE_PASS: TRUE",
  "STEP06_NETWORK_PLOTS_PASS: TRUE"
), file.path(partial_dir, "06_completion_report.txt"))

if (!file.rename(partial_dir, output_dir)) {
  stopf("Could not finalize output directory: %s", output_dir)
}

cat("\nSTEP06_GETHUBGENES_COUNTS_PASS: TRUE\n")
cat("STEP06_EXACT_200_GENE_ANNOTATION_PASS: TRUE\n")
cat("STEP06_FUNCTIONAL_CATEGORY_PASS: TRUE\n")
cat("STEP06_CATEGORY_PALETTE_PASS: TRUE\n")
cat("STEP06_NETWORK_PLOTS_PASS: TRUE\n")
say("Step 06 completed successfully")
say("Output:", output_dir)
