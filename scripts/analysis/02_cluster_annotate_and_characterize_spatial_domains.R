# ============================================================
# 02_cluster_annotate_and_characterize_spatial_domains.R
#
# Purpose:
#   Perform the principal biological analysis of the integrated
#   FLINC spatial transcriptomics dataset.
#
# Main operations:
#   - graph-based clustering
#   - hierarchical subclustering
#   - final biological annotation
#   - cluster marker identification
#   - Gene Ontology enrichment
#   - cell-type deconvolution summaries
#   - cluster-cell-type chord relationships
#   - stage-associated cluster proportion analyses
#   - generation of final annotated analysis objects
#
# Implementation:
#   All clustering and characterization code is included below as R code.
#   Run this script directly; no notebook extraction or rendering is required.
#
# Reference execution on the validated Mac environment:
#
#   cd "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
#
#   FLINC_THREADS=4 caffeinate -dimsu Rscript \
#     scripts/analysis/02_cluster_annotate_and_characterize_spatial_domains.R \
#     > logs/clustering_characterization_$(date +%Y%m%d_%H%M%S).log 2>&1
# ============================================================

rm(list = ls())
gc()

N_THREADS <- as.integer(Sys.getenv("FLINC_THREADS", unset = "4"))
if (is.na(N_THREADS) || N_THREADS < 1) N_THREADS <- 4

Sys.setenv(
  OMP_NUM_THREADS = as.character(N_THREADS),
  OPENBLAS_NUM_THREADS = as.character(N_THREADS),
  MKL_NUM_THREADS = as.character(N_THREADS),
  VECLIB_MAXIMUM_THREADS = as.character(N_THREADS),
  NUMEXPR_NUM_THREADS = as.character(N_THREADS)
)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 30 * 1024^3)

REPRO_ROOT <- "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
SCRIPT_DIR <- file.path(REPRO_ROOT, "scripts")
OBJECT_DIR <- file.path(REPRO_ROOT, "objects")
RESULTS_DIR <- file.path(REPRO_ROOT, "results")
LOG_DIR <- file.path(REPRO_ROOT, "logs")
REPRO_LOG_DIR <- file.path(REPRO_ROOT, "reproducibility_logs")

dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPRO_LOG_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(REPRO_ROOT)

analysis_script <- file.path(
  SCRIPT_DIR, "analysis", "02_cluster_annotate_and_characterize_spatial_domains.R"
)

input_rds <- file.path(OBJECT_DIR, "object.ffpe.rds")
if (!file.exists(input_rds)) {
  stop("Missing integrated spatial analysis object:\n", input_rds)
}

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

write_checksum_csv <- function(files, out_file) {
  files <- files[!is.na(files) & file.exists(files)]
  if (length(files) == 0) {
    write.csv(data.frame(file = character(), md5 = character()), out_file, row.names = FALSE)
    return(invisible(NULL))
  }
  files_norm <- normalizePath(files, mustWork = TRUE)
  md5 <- tools::md5sum(files_norm)
  write.csv(
    data.frame(file = names(md5), md5 = unname(md5), stringsAsFactors = FALSE),
    out_file,
    row.names = FALSE
  )
  invisible(md5)
}

log_msg("Running spatial clustering, annotation and biological characterization")
log_msg("Analysis script:", analysis_script)
log_msg("Input object:", input_rds)
log_msg("Threads:", N_THREADS)

sink(file.path(REPRO_LOG_DIR, "sessionInfo_before_06_to_10_as_Rscript.txt"))
cat("Timestamp:\n")
print(Sys.time())
cat("\nR version:\n")
print(R.version.string)
cat("\n.libPaths():\n")
print(.libPaths())
cat("\nSys.info():\n")
print(Sys.info())
cat("\nThread settings:\n")
print(Sys.getenv(c(
  "OMP_NUM_THREADS",
  "OPENBLAS_NUM_THREADS",
  "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS",
  "NUMEXPR_NUM_THREADS"
)))
cat("\nsessionInfo():\n")
print(sessionInfo())
sink()

write.csv(
  installed.packages()[, c("Package", "Version", "LibPath", "Built")],
  file.path(REPRO_LOG_DIR, "installed_packages_before_06_to_10_as_Rscript.csv"),
  row.names = FALSE
)

write_checksum_csv(
  c(analysis_script, input_rds),
  file.path(REPRO_LOG_DIR, "input_md5_before_06_to_10_as_Rscript.csv")
)

log_msg("Running the clustering and characterization analysis")
log_msg("This will run clustering, annotation, markers, GO, deconvolution, chord diagrams, proportions, and final checksums.")

# Keep analysis objects and helpers separate from the execution-record wrapper.
run_env <- new.env(parent = globalenv())

evalq({
# BEGIN SPATIAL DOMAIN ANALYSIS
# ============================================================
# 260628_approach1_cluster_go_marker_deconv
# Reviewer-proof continuation after:
# 01 load data
# 02 SCTransform
# 03 integration
# 04 deconvolution assay
# 05 PCA / Harmony / UMAP
#
# Approach 1:
#   Main clustering resolution = 0.3
#   Main seed = 39152
#   Neighbor graph from harmony_SCT
#   Subcluster main C2 at resolution = 0.4
#   Sub-subcluster C2/S2 at resolution = 0.3
# ============================================================


# ============================================================
# 06. Setup and load post-reduction object
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(grid)
})

setwd("/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code")

input_rds <- "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code/objects/object.ffpe.rds"

out_dir <- "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code/results/260701_approach1_cluster_go_marker_deconv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

repro_dir <- file.path(out_dir, "reproducibility_logs")
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

write_checksum_csv <- function(files, out_file) {
  files <- files[!is.na(files) & file.exists(files)]
  if (length(files) == 0) {
    write.csv(data.frame(file = character(), md5 = character()), out_file, row.names = FALSE)
    return(invisible(NULL))
  }
  files_norm <- normalizePath(files, mustWork = TRUE)
  md5 <- tools::md5sum(files_norm)
  write.csv(data.frame(file = names(md5), md5 = unname(md5), stringsAsFactors = FALSE), out_file, row.names = FALSE)
  invisible(md5)
}

sink(file.path(repro_dir, "sessionInfo_before_06_to_10.txt"))
cat("Timestamp:\n"); print(Sys.time())
cat("\nR version:\n"); print(R.version.string)
cat("\n.libPaths():\n"); print(.libPaths())
cat("\nSys.info():\n"); print(Sys.info())
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

write.csv(installed.packages()[, c("Package", "Version", "LibPath", "Built")], file.path(repro_dir, "installed_packages_before_06_to_10.csv"), row.names = FALSE)

write_checksum_csv(input_rds, file.path(repro_dir, "input_object_ffpe_md5_before_06_to_10.csv"))


approach_name <- "Approach 1: main res 0.3 seed 39152"

cat("============================================================\n")
cat("Approach 1 clustering / markers / GO / deconvolution\n")
cat("============================================================\n")

cat("\nR version:\n")
print(R.version.string)

cat("\nLibrary paths:\n")
print(.libPaths())

cat("\nInput object:\n")
cat(input_rds, "\n")

cat("\nOutput directory:\n")
cat(out_dir, "\n")

object <- readRDS(input_rds)

cat("\nObject loaded successfully.\n")

cat("\nAssays in object:\n")
print(Assays(object))

cat("\nReductions in object:\n")
print(Reductions(object))

stopifnot("SCT" %in% Assays(object))
stopifnot("harmony_SCT" %in% Reductions(object))
stopifnot("umap_SCT" %in% Reductions(object))

cat("\nRequired assays/reductions found. Proceeding.\n")


# ============================================================
# 07. Approach 1 main clustering
# ============================================================
# Important:
# FindNeighbors uses harmony_SCT to build graph_SCT.
# FindClusters then clusters graph_SCT.
# This is the original Approach 1 logic.
# ============================================================

set.seed(39152)

object <- FindNeighbors(
  object,
  assay = "SCT",
  dims = 1:15,
  k.param = 23,
  graph.name = "graph_SCT",
  reduction = "harmony_SCT"
)

object <- FindClusters(
  object,
  resolution = 0.3,
  save.SNN = TRUE,
  do.sparse = TRUE,
  graph.name = "graph_SCT",
  random.seed = 39152,
  group.singletons = TRUE
)

# Seurat 4.2.3-compatible replacement for cluster.name = "cluster"
object$cluster <- as.character(Idents(object))
Idents(object) <- "cluster"

cat("\nMain cluster table:\n")
print(table(object$cluster))

write.csv(
  as.data.frame(table(object$cluster)),
  file = file.path(out_dir, "main_cluster_counts.csv"),
  row.names = FALSE
)

p_main <- DimPlot(
  object,
  reduction = "umap_SCT",
  group.by = "cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(paste0(approach_name, " - main clusters")) +
  theme_classic()

print(p_main)

ggsave(
  filename = file.path(out_dir, "UMAP_main_clusters.pdf"),
  plot = p_main,
  width = 8,
  height = 6
)


# ============================================================
# 08. Subcluster main cluster 2
# ============================================================
# Original Approach 1:
#   subset cluster == 2
#   FindNeighbors on harmony_SCT dims 1:40
#   k.param = 30
#   graph.name = sub_clustering
#   FindClusters resolution = 0.4
# ============================================================

if (!"2" %in% as.character(unique(object$cluster))) {
  stop("Cluster 2 was not found in object$cluster.")
}

object.cluster2 <- subset(object, subset = cluster == 2)

object.cluster2 <- FindNeighbors(
  object.cluster2,
  reduction = "harmony_SCT",
  dims = 1:40,
  assay = "SCT",
  k.param = 30,
  graph.name = "sub_clustering"
)

object.cluster2 <- FindClusters(
  object.cluster2,
  graph.name = "sub_clustering",
  resolution = 0.4
)

# Seurat 4.2.3-compatible replacement for cluster.name = "sub_cluster"
object.cluster2$sub_cluster <- paste0("S", as.character(Idents(object.cluster2)))
Idents(object.cluster2) <- "sub_cluster"

cat("\nSubcluster table for main cluster 2:\n")
print(table(object.cluster2$sub_cluster))

write.csv(
  as.data.frame(table(object.cluster2$sub_cluster)),
  file = file.path(out_dir, "cluster2_subcluster_counts.csv"),
  row.names = FALSE
)

p_sub <- DimPlot(
  object.cluster2,
  reduction = "umap_SCT",
  group.by = "sub_cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(paste0(approach_name, " - subclusters inside C2")) +
  theme_classic()

print(p_sub)

ggsave(
  filename = file.path(out_dir, "UMAP_cluster2_subclusters.pdf"),
  plot = p_sub,
  width = 8,
  height = 6
)


# ============================================================
# 09. Sub-subcluster C2/S2
# ============================================================
# Original Approach 1:
#   subset sub_cluster == S2
#   FindNeighbors on harmony_SCT dims 1:40
#   k.param = 30
#   graph.name = sub_sub_clustering
#   FindClusters resolution = 0.3
# ============================================================

if (!"S2" %in% as.character(unique(object.cluster2$sub_cluster))) {
  stop("Subcluster S2 was not found in object.cluster2$sub_cluster.")
}

object.cluster3 <- subset(object.cluster2, subset = sub_cluster == "S2")

object.cluster3 <- FindNeighbors(
  object.cluster3,
  reduction = "harmony_SCT",
  dims = 1:40,
  assay = "SCT",
  k.param = 30,
  graph.name = "sub_sub_clustering"
)

object.cluster3 <- FindClusters(
  object.cluster3,
  graph.name = "sub_sub_clustering",
  resolution = 0.3
)

# Seurat 4.2.3-compatible replacement for cluster.name = "sub_sub_cluster"
object.cluster3$sub_sub_cluster <- paste0("K", as.character(Idents(object.cluster3)))
Idents(object.cluster3) <- "sub_sub_cluster"

cat("\nSub-subcluster table for C2/S2:\n")
print(table(object.cluster3$sub_sub_cluster))

write.csv(
  as.data.frame(table(object.cluster3$sub_sub_cluster)),
  file = file.path(out_dir, "cluster2_S2_subsubcluster_counts.csv"),
  row.names = FALSE
)

p_subsub <- DimPlot(
  object.cluster3,
  reduction = "umap_SCT",
  group.by = "sub_sub_cluster",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(paste0(approach_name, " - sub-subclusters inside C2/S2")) +
  theme_classic()

print(p_subsub)

ggsave(
  filename = file.path(out_dir, "UMAP_cluster2_S2_subsubclusters.pdf"),
  plot = p_subsub,
  width = 8,
  height = 6
)

# ============================================================
# FINAL ANNOTATION + MARKERS + GO + DECONV + CHORD DIAGRAMS
# 10 Approach 1 reviewer-facing output
#
# Changes included:
#   - C3 Immune cells = green
#   - C5 Cholangiocytes = yellow
#   - C6 is NOT plotted in front on UMAP
#   - C4 Myofibroblasts plotted last/on top
#   - UMAP uses manual ggplot/ggrepel style
#   - GO colour scale made darker
#   - deconvolution names fixed for current assay
#   - immune/stromal and parenchymal chord diagrams generated
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(grid)
  library(ggrepel)
  library(circlize)
})

# ============================================================
# 10.0. Safety checks
# ============================================================

stopifnot(exists("object"))
stopifnot(exists("object.cluster2"))
stopifnot(exists("object.cluster3"))
stopifnot(exists("out_dir"))
stopifnot(exists("approach_name"))

stopifnot("SCT" %in% Assays(object))
stopifnot("umap_SCT" %in% Reductions(object))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("Starting final annotation + markers + GO + deconv + chord\n")
cat("Output folder:\n")
cat(out_dir, "\n")
cat("============================================================\n")


# ============================================================
# 10.1. Create raw Approach-1 labels
# ============================================================

object$cluster_final_raw <- paste0("C", object$cluster)

object$sub_cluster <- NA_character_
object$sub_cluster[colnames(object.cluster2)] <- as.character(object.cluster2$sub_cluster)

object$sub_sub_cluster <- NA_character_
object$sub_sub_cluster[colnames(object.cluster3)] <- as.character(object.cluster3$sub_sub_cluster)

object$cluster_sub_all_raw <- object$cluster_final_raw

object$cluster_sub_all_raw[!is.na(object$sub_cluster)] <-
  paste0("C2_", object$sub_cluster[!is.na(object$sub_cluster)])

object$cluster_sub_all_raw[!is.na(object$sub_sub_cluster)] <-
  paste0(
    "C2_S2_",
    object$sub_sub_cluster[!is.na(object$sub_sub_cluster)]
  )

cat("\nRaw Approach-1 combined labels:\n")
print(table(object$cluster_sub_all_raw))


# ============================================================
# 10.2. Rename to final reviewer-facing clusters and annotations
# ============================================================

cluster_rename_map <- c(
  "C0"       = "C1",  # Periportal hepatocytes
  "C1"       = "C2",  # Centrilobular hepatocytes
  "C2_S0"    = "C0",  # Portal tracts
  "C2_S1"    = "C6",  # Matrix-rich fibroblast
  "C2_S2_K0" = "C5",  # Cholangiocytes
  "C2_S2_K1" = "C4",  # Myofibroblasts
  "C2_S3"    = "C3"   # Immune cells
)

annotation_map <- c(
  "C0" = "Portal tracts",
  "C1" = "Periportal hepatocytes",
  "C2" = "Centrilobular hepatocytes",
  "C3" = "Immune cells",
  "C4" = "Myofibroblasts",
  "C5" = "Cholangiocytes",
  "C6" = "Matrix-rich fibroblast"
)

final_order <- c(
  "C0 Portal tracts",
  "C1 Periportal hepatocytes",
  "C2 Centrilobular hepatocytes",
  "C3 Immune cells",
  "C4 Myofibroblasts",
  "C5 Cholangiocytes",
  "C6 Matrix-rich fibroblast"
)

# C3 immune = green
# C5 cholangiocytes = yellow
final_palette <- c(
  C0 = "#56B4E9",
  C1 = "#FF9DA7",
  C2 = "#B189A8",
  C3 = "#59A14F",
  C4 = "#F28E2B",
  C5 = "#EDC948",
  C6 = "#E5C494"
)

palette_named <- setNames(unname(final_palette), final_order)

if (any(is.na(cluster_rename_map[object$cluster_sub_all_raw]))) {
  missing_raw <- unique(object$cluster_sub_all_raw[is.na(cluster_rename_map[object$cluster_sub_all_raw])])
  stop(
    "Some raw cluster labels are missing from cluster_rename_map:\n",
    paste(missing_raw, collapse = ", ")
  )
}

object$cluster <- unname(cluster_rename_map[object$cluster_sub_all_raw])
object$annotation <- unname(annotation_map[object$cluster])
object$cluster_named <- paste(object$cluster, object$annotation)

object$cluster <- factor(object$cluster, levels = names(annotation_map))
object$annotation <- factor(object$annotation, levels = unname(annotation_map))
object$cluster_named <- factor(object$cluster_named, levels = final_order)

Idents(object) <- "cluster"

cat("\nFinal reviewer-facing cluster table:\n")
print(table(object$cluster))

cat("\nFinal annotation table:\n")
print(table(object$annotation))

cat("\nFinal named cluster table:\n")
print(table(object$cluster_named))

write.csv(
  as.data.frame(table(object$cluster_sub_all_raw)),
  file = file.path(out_dir, "raw_approach1_cluster_sub_all_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(object$cluster)),
  file = file.path(out_dir, "final_cluster_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(object$cluster_named)),
  file = file.path(out_dir, "final_cluster_named_counts.csv"),
  row.names = FALSE
)

cluster_mapping_table <- data.frame(
  raw_approach1_label = names(cluster_rename_map),
  final_cluster = unname(cluster_rename_map),
  annotation = unname(annotation_map[unname(cluster_rename_map)]),
  final_cluster_named = paste(
    unname(cluster_rename_map),
    unname(annotation_map[unname(cluster_rename_map)])
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_mapping_table,
  file = file.path(out_dir, "approach1_raw_to_final_cluster_mapping.csv"),
  row.names = FALSE
)

cat("\nCluster mapping table:\n")
print(cluster_mapping_table)


# ============================================================
# 10.3. UMAP with small niche clusters plotted on top
# ============================================================

emb <- Embeddings(object, "umap_SCT")

umap_df <- data.frame(
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  cluster = as.character(object$cluster),
  cluster_named = factor(as.character(object$cluster_named), levels = final_order),
  stringsAsFactors = FALSE
)

# Plot order:
# Big/background clusters first.
# C6 is drawn before the small niche clusters, so it is NOT in front.
# C4 is drawn last, so myofibroblasts stay visible.
plot_order <- c(
  "C2 Centrilobular hepatocytes",
  "C1 Periportal hepatocytes",
  "C6 Matrix-rich fibroblast",
  "C0 Portal tracts",
  "C3 Immune cells",
  "C5 Cholangiocytes",
  "C4 Myofibroblasts"
)

umap_df$plot_order <- factor(
  as.character(umap_df$cluster_named),
  levels = plot_order
)

umap_df <- umap_df[order(umap_df$plot_order), ]

centers <- aggregate(
  umap_df[, c("UMAP_1", "UMAP_2")],
  by = list(cluster = umap_df$cluster),
  FUN = median,
  na.rm = TRUE
)

p_umap_front <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = cluster_named)
) +
  geom_point(size = 0.55, alpha = 0.9) +
  scale_color_manual(
    values = palette_named,
    name = "Annotation",
    drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = cluster),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.alpha = 0.6,
    fontface = "bold"
  ) +
  ggtitle("Approach 1 final clustering") +
  theme_minimal(base_size = 14) +
  coord_cartesian(clip = "off") +
  guides(
    color = guide_legend(
      override.aes = list(size = 5, alpha = 1)
    )
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    legend.key.width  = unit(8, "pt"),
    legend.key.height = unit(8, "pt"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

print(p_umap_front)

ggsave(
  filename = file.path(out_dir, "UMAP_final_annotated_clusters_C4_on_top_C3green_C5yellow.pdf"),
  plot = p_umap_front,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "UMAP_final_annotated_clusters_C4_on_top_C3green_C5yellow.png"),
  plot = p_umap_front,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

saveRDS(
  object,
  file = file.path(out_dir, "object_full_with_final_cluster_annotation_labels.rds")
)


# ============================================================
# 10.4. Marker analysis using final C0-C6 clusters
# ============================================================

DefaultAssay(object) <- "SCT"
Idents(object) <- "cluster"

cat("\nFinal clusters used for marker analysis:\n")
print(table(Idents(object)))

write.csv(
  as.data.frame(table(Idents(object))),
  file = file.path(out_dir, "marker_input_final_cluster_counts.csv"),
  row.names = FALSE
)

object <- PrepSCTFindMarkers(
  object = object,
  assay = "SCT",
  verbose = TRUE
)

cluster_marker <- FindAllMarkers(
  object = object,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  cluster_marker,
  file = file.path(out_dir, "cluster_marker_final_clusters.csv"),
  row.names = FALSE
)

cluster_marker_top20 <- cluster_marker %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  cluster_marker_top20,
  file = file.path(out_dir, "cluster_marker_final_clusters_top20_per_cluster.csv"),
  row.names = FALSE
)

marker_summary <- cluster_marker %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_markers_padj_005 = sum(p_val_adj < 0.05, na.rm = TRUE),
    max_avg_log2FC = max(avg_log2FC, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  marker_summary,
  file = file.path(out_dir, "cluster_marker_summary.csv"),
  row.names = FALSE
)

cat("\nMarker summary:\n")
print(marker_summary)

sink(file.path(out_dir, "sessionInfo_markers.txt"))
print(sessionInfo())
sink()


# ============================================================
# 10.5. Top 10 marker dotplot
# ============================================================

dat <- cluster_marker %>%
  filter(
    !is.na(cluster),
    !is.na(gene),
    gene != "",
    !is.na(avg_log2FC),
    !is.na(p_val_adj),
    avg_log2FC > 0.5
  )

if (nrow(dat) == 0) {
  stop("No markers passed avg_log2FC > 0.5 for top marker dotplot.")
}

eps <- suppressWarnings(min(dat$p_val_adj[dat$p_val_adj > 0], na.rm = TRUE))
if (!is.finite(eps)) eps <- 1e-300

plot_df <- dat %>%
  mutate(
    cluster = as.character(cluster),
    p_adj_clipped = if_else(
      p_val_adj <= 0 | is.na(p_val_adj),
      eps,
      p_val_adj
    )
  ) %>%
  arrange(cluster, p_adj_clipped, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    negLog10Padj = -log10(p_adj_clipped),
    cluster = factor(cluster, levels = names(annotation_map)),
    gene = factor(gene, levels = rev(unique(gene)))
  )

write.csv(
  plot_df,
  file = file.path(out_dir, "Top10_markers_per_final_cluster_dotplot_data.csv"),
  row.names = FALSE
)

p_top10 <- ggplot(plot_df, aes(x = cluster, y = gene)) +
  geom_point(aes(size = negLog10Padj, color = avg_log2FC)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "avg_log2FC"
  ) +
  scale_size_area(max_size = 8, name = "-log10 adj.p") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(expand = expansion(mult = 0, add = 0.60)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Final cluster",
    y = "Top markers",
    title = "Top 10 markers per final cluster",
    subtitle = approach_name
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 14, b = 10, l = 14)
  )

print(p_top10)

ggsave(
  filename = file.path(out_dir, "Top10_markers_per_final_cluster_dotplot.pdf"),
  plot = p_top10,
  width = 7,
  height = 9,
  units = "in"
)


# ============================================================
# 10.6. GO:BP enrichment with darker colour scale
# ============================================================

markers_all <- cluster_marker

markers <- markers_all %>%
  mutate(
    cluster = as.character(cluster),
    gene = trimws(as.character(gene))
  ) %>%
  filter(
    !is.na(cluster),
    !is.na(gene),
    gene != "",
    !is.na(avg_log2FC),
    avg_log2FC > 0.5,
    !is.na(p_val_adj),
    p_val_adj < 0.05
  )

if (nrow(markers) == 0) {
  stop("No marker genes passed avg_log2FC > 0.5 and p_val_adj < 0.05 for GO.")
}

cat("\nFiltered marker counts per cluster for GO:\n")
print(table(markers$cluster))

bg_symbols <- unique(trimws(as.character(markers_all$gene)))
bg_symbols <- bg_symbols[!is.na(bg_symbols) & nchar(bg_symbols) > 0]

bg_map <- suppressMessages(
  bitr(
    bg_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

bg_entrez <- unique(bg_map$ENTREZID)

sym2ent <- suppressMessages(
  bitr(
    unique(markers$gene),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

get_cluster_entrez <- function(df, cl, map_tbl) {
  syms <- df %>%
    filter(cluster == cl) %>%
    distinct(gene) %>%
    pull(gene)
  
  map_tbl %>%
    semi_join(tibble::tibble(SYMBOL = syms), by = "SYMBOL") %>%
    pull(ENTREZID) %>%
    unique()
}

cluster_levels <- names(annotation_map)
cluster_levels <- cluster_levels[cluster_levels %in% unique(markers$cluster)]

ego_list <- list()

for (cl in cluster_levels) {
  
  g <- get_cluster_entrez(markers, cl, sym2ent)
  
  cat("\nCluster:", cl, "- Entrez genes:", length(g), "\n")
  
  if (length(g) < 5) {
    cat("Skipping", cl, "because fewer than 5 Entrez genes.\n")
    next
  }
  
  ego <- suppressMessages(
    enrichGO(
      gene          = g,
      universe      = bg_entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pvalueCutoff  = 0.1,
      qvalueCutoff  = 0.2,
      pAdjustMethod = "BH",
      readable      = TRUE,
      minGSSize     = 10,
      maxGSSize     = 500
    )
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    
    ego <- simplify(
      ego,
      cutoff = 0.6,
      by = "p.adjust",
      select_fun = min
    )
    
    ego_list[[cl]] <- as.data.frame(ego) %>%
      mutate(cluster = cl)
  }
}

ego_all <- bind_rows(ego_list)

if (!nrow(ego_all)) {
  stop("No enriched GO terms found at current cutoffs.")
}

write.csv(
  ego_all,
  file = file.path(out_dir, "GO_BP_all_terms_final_clusters.csv"),
  row.names = FALSE
)

shorten_terms <- function(x) {
  x <- str_replace(x, regex("^peptide antigen assembly with MHC protein complex$", TRUE), "MHC peptide loading")
  x <- str_replace(x, regex("peptide antigen assembly with MHC class I protein complex", TRUE), "MHC I peptide loading")
  x <- str_replace(x, regex("peptide antigen assembly with MHC class II protein complex", TRUE), "MHC II peptide loading")
  x <- str_replace(x, regex("^enzyme[-\\s]?linked receptor (protein )?signall?ing pathway$", TRUE), "enzyme-linked receptor signaling")
  x <- str_replace(x, regex("^positive regulation of immune system process$", TRUE), "immune system activation")
  x <- str_replace(x, regex("^positive regulation of immune response$", TRUE), "immune activation")
  x
}

ego_top <- ego_all %>%
  mutate(
    Description = str_trunc(Description, 70),
    Description = shorten_terms(Description),
    Cluster = as.character(cluster)
  ) %>%
  group_by(cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

write.csv(
  ego_top,
  file = file.path(out_dir, "GO_BP_top_terms_final_clusters.csv"),
  row.names = FALSE
)

cluster_levels_plot <- cluster_levels[cluster_levels %in% unique(ego_top$Cluster)]

term_order <- ego_top %>%
  mutate(Cluster = factor(Cluster, levels = cluster_levels_plot)) %>%
  group_by(Description) %>%
  summarise(
    best_cluster_index = min(as.numeric(Cluster), na.rm = TRUE),
    best_p = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(best_cluster_index, best_p) %>%
  pull(Description)

ego_top <- ego_top %>%
  mutate(
    Description = factor(Description, levels = rev(unique(term_order))),
    Cluster = factor(Cluster, levels = cluster_levels_plot),
    negLog10Padj = -log10(p.adjust)
  )

p_go <- ggplot(
  ego_top,
  aes(x = Cluster, y = Description)
) +
  geom_vline(
    xintercept = seq_along(cluster_levels_plot),
    linewidth = 0.4,
    colour = "grey88"
  ) +
  geom_point(
    aes(size = Count, colour = negLog10Padj),
    alpha = 0.98
  ) +
  scale_colour_gradientn(
    colours = c("#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"),
    name = expression(-log[10]("adj. p"))
  ) +
  scale_size_area(
    max_size = 7,
    name = "Count"
  ) +
  scale_x_discrete(
    limits = names(annotation_map),
    drop = FALSE,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "GO:BP top 5 terms per final cluster",
    subtitle = paste0(approach_name, " | colour = -log10(adj. p)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 10,
      angle = 20,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    plot.margin = margin(t = 10, r = 14, b = 70, l = 32)
  ) +
  coord_cartesian(clip = "off")

print(p_go)

ggsave(
  filename = file.path(out_dir, "GO_matrix_terms_by_final_cluster_darker_blue.pdf"),
  plot = p_go,
  width = 8.5,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "GO_matrix_terms_by_final_cluster_darker_blue.png"),
  plot = p_go,
  width = 8.5,
  height = 10,
  units = "in",
  dpi = 300
)

sink(file.path(out_dir, "sessionInfo_GO.txt"))
print(sessionInfo())
sink()


# ============================================================
# 10.7. Deconvolution stacked barplot using final clusters
# ============================================================

stopifnot("singlecell" %in% Assays(object))

DefaultAssay(object) <- "singlecell"

deconv_features <- rownames(object[["singlecell"]])

cat("\nDeconvolution features:\n")
print(deconv_features)

# Correct colours for current deconvolution feature names
cols_celltype_all <- c(
  CentralHep    = "#76B7B2",
  Cholangiocyte = "#EDC948",
  Erythroid     = "#B07AA1",
  InterHep      = "#BD9BB5",
  PortalEndo    = "#4E79A7",
  PortalHep     = "#FFB4BC",
  Stellate      = "#F8AB61",
  cvLSECs       = "#BAB0AC",
  Bcell         = "#59A14F",
  Bcells        = "#59A14F",
  InfMac        = "#59A14F",
  NKTcell       = "#59A14F",
  NonInfMac     = "#59A14F",
  Tcell         = "#59A14F",
  Macrophage    = "#59A14F",
  Monocyte      = "#59A14F",
  Immune        = "#59A14F"
)

scores <- FetchData(object, vars = deconv_features)
scores <- as.data.frame(scores)
scores[] <- lapply(scores, as.numeric)

# Normalize each spot to sum to 1
rs <- rowSums(scores, na.rm = TRUE)
rs[rs == 0 | is.na(rs)] <- 1
scores <- scores / rs

# Merge immune classes if present, but keep all structural/parenchymal classes
immune_names <- intersect(
  c("Bcell", "Bcells", "InfMac", "NKTcell", "NonInfMac", "Tcell", "Macrophage", "Monocyte"),
  colnames(scores)
)

cat("\nImmune classes merged for deconv stacked plot:\n")
print(immune_names)

if (length(immune_names) > 0) {
  scores$Immune <- rowSums(scores[, immune_names, drop = FALSE], na.rm = TRUE)
  keep_features <- setdiff(colnames(scores), c(immune_names, "Immune"))
  keep_features <- unique(c(keep_features, "Immune"))
} else {
  keep_features <- unique(colnames(scores))
}

scores <- scores[, keep_features, drop = FALSE]

scores$cluster <- as.character(object$cluster_named)
scores$cluster <- factor(scores$cluster, levels = final_order)

cluster_comp <- scores %>%
  group_by(cluster) %>%
  summarise(
    across(
      all_of(keep_features),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

write.csv(
  cluster_comp,
  file = file.path(out_dir, "deconvolution_mean_fraction_per_final_cluster_wide_FIXED.csv"),
  row.names = FALSE
)

plot_data <- cluster_comp %>%
  pivot_longer(
    cols = all_of(keep_features),
    names_to = "celltype",
    values_to = "fraction"
  ) %>%
  mutate(
    celltype = factor(celltype, levels = keep_features),
    cluster = factor(cluster, levels = final_order)
  )

write.csv(
  plot_data,
  file = file.path(out_dir, "deconvolution_mean_fraction_per_final_cluster_long_FIXED.csv"),
  row.names = FALSE
)

dominant <- plot_data %>%
  group_by(cluster) %>%
  arrange(desc(fraction), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

write.csv(
  dominant,
  file = file.path(out_dir, "deconvolution_top3_celltypes_per_final_cluster_FIXED.csv"),
  row.names = FALSE
)

cat("\nTop 3 deconvolution classes per final cluster:\n")
print(dominant)

missing_cols <- setdiff(levels(plot_data$celltype), names(cols_celltype_all))

if (length(missing_cols) > 0) {
  auto_cols <- setNames(
    grDevices::hcl.colors(length(missing_cols), palette = "Dark 3"),
    missing_cols
  )
  cols_celltype_plot <- c(cols_celltype_all, auto_cols)
} else {
  cols_celltype_plot <- cols_celltype_all
}

cols_celltype_plot <- cols_celltype_plot[levels(plot_data$celltype)]

p_deconv_fixed <- ggplot(
  plot_data,
  aes(x = cluster, y = fraction, fill = celltype)
) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(
    values = cols_celltype_plot,
    drop = FALSE,
    guide = guide_legend(nrow = 2, byrow = TRUE, title = NULL)
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    y = "Mean deconvolution fraction",
    x = NULL,
    title = "Deconvolution composition per final cluster",
    subtitle = approach_name
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.55, "lines"),
    legend.key.width = unit(0.9, "lines"),
    legend.box.margin = margin(b = 6),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

print(p_deconv_fixed)

ggsave(
  filename = file.path(out_dir, "deconvolution_stacked_barplot_per_final_cluster_FIXED.pdf"),
  plot = p_deconv_fixed,
  width = 9,
  height = 5.8,
  units = "in",
  device = cairo_pdf,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "deconvolution_stacked_barplot_per_final_cluster_FIXED.png"),
  plot = p_deconv_fixed,
  width = 9,
  height = 5.8,
  units = "in",
  dpi = 300
)


# ============================================================
# 10.8. Chord diagram helper
# ============================================================

make_chord_diagram <- function(
    object,
    celltype_keep,
    clusters_named_keep,
    output_pdf,
    cols_celltype_all,
    final_palette,
    celltype_label_override = NULL
) {
  
  DefaultAssay(object) <- "singlecell"
  
  metadata <- FetchData(object, vars = rownames(object[["singlecell"]]))
  metadata$cluster <- as.character(object$cluster)
  metadata$cluster_named <- as.character(object$cluster_named)
  
  celltype <- intersect(celltype_keep, colnames(metadata))
  
  if (length(celltype) == 0) {
    warning("No requested cell types found for chord diagram: ", output_pdf)
    return(NULL)
  }
  
  metadata <- metadata[
    metadata$cluster_named %in% clusters_named_keep,
    ,
    drop = FALSE
  ]
  
  if (nrow(metadata) == 0) {
    warning("No requested clusters found for chord diagram: ", output_pdf)
    return(NULL)
  }
  
  cluster_named_vec <- final_order[final_order %in% unique(metadata$cluster_named)]
  cluster_ids_vec <- sub(" .*", "", cluster_named_vec)
  
  cols_cluster <- setNames(final_palette[cluster_ids_vec], cluster_named_vec)
  
  missing_cols <- setdiff(celltype, names(cols_celltype_all))
  if (length(missing_cols) > 0) {
    auto_cols <- setNames(
      grDevices::hcl.colors(length(missing_cols), palette = "Dark 3"),
      missing_cols
    )
    cols_celltype_local <- c(cols_celltype_all, auto_cols)
  } else {
    cols_celltype_local <- cols_celltype_all
  }
  
  cols_celltype <- cols_celltype_local[celltype]
  
  temp <- data.frame()
  
  for (i in seq_along(celltype)) {
    for (j in seq_along(cluster_named_vec)) {
      flag <- which(metadata$cluster_named == cluster_named_vec[j])
      temp[j, celltype[i]] <- sum(metadata[flag, celltype[i]], na.rm = TRUE)
    }
  }
  
  row_sums <- rowSums(temp, na.rm = TRUE)
  row_sums[row_sums == 0 | is.na(row_sums)] <- 1
  temp <- temp / row_sums
  temp$cluster <- cluster_named_vec
  
  chord_data <- pivot_longer(
    temp,
    cols = all_of(celltype),
    names_to = "celltype",
    values_to = "cell_frac"
  ) %>%
    filter(!is.na(cell_frac), cell_frac > 0)
  
  write.csv(
    chord_data,
    file = sub("\\.pdf$", "_data.csv", output_pdf),
    row.names = FALSE
  )
  
  pdf(output_pdf, width = 7, height = 7)
  
  grid.col <- c(cols_celltype, cols_cluster)
  df <- as.data.frame(chord_data)
  link.col <- cols_celltype[df$celltype]
  
  circos.clear()
  circos.par(
    track.height = 0.16,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.06)
  )
  
  chordDiagram(
    df,
    grid.col              = grid.col,
    col                   = link.col,
    transparency          = 0.25,
    annotationTrack       = "grid",
    annotationTrackHeight = 0.12,
    preAllocateTracks     = list(track.height = 0.08)
  )
  
  highlight.sector(
    sector.index = celltype,
    track.index  = 1,
    col          = "white",
    border       = "black",
    text         = "cell-type",
    text.col     = "black",
    cex          = 1.2,
    font         = 2,
    niceFacing   = TRUE,
    facing       = "bending.inside"
  )
  
  highlight.sector(
    sector.index = cluster_named_vec,
    track.index  = 1,
    col          = "white",
    border       = "black",
    text         = "spatial cluster",
    text.col     = "black",
    cex          = 1.2,
    font         = 2,
    niceFacing   = TRUE,
    facing       = "bending.inside"
  )
  
  circos.trackPlotRegion(
    track.index = 2,
    panel.fun = function(x, y) {
      
      sector <- CELL_META$sector.index
      ylim <- CELL_META$ylim
      y_min <- min(ylim)
      y_max <- max(ylim)
      
      y_inner_pos <- y_min + 0.18 * (y_max - y_min)
      y_outer_pos <- y_min + 0.58 * (y_max - y_min)
      y_mid_pos   <- y_min + 0.40 * (y_max - y_min)
      
      if (sector %in% cluster_named_vec) {
        
        lab <- sub(" .*", "", sector)
        
        circos.text(
          CELL_META$xcenter, y_mid_pos, lab,
          facing     = "bending.inside",
          niceFacing = TRUE,
          adj        = c(0.5, 0),
          cex        = 1,
          font       = 2
        )
        
      } else {
        
        inner_lab <- sector
        outer_lab <- NULL
        
        if (!is.null(celltype_label_override) && sector %in% names(celltype_label_override)) {
          inner_lab <- celltype_label_override[[sector]]
        }
        
        if (sector == "PortalEndo") {
          inner_lab <- "E"
          outer_lab <- "P"
        }
        
        if (sector == "Cholangiocyte") {
          inner_lab <- "Chol"
        }
        
        if (is.null(outer_lab)) {
          circos.text(
            CELL_META$xcenter, y_mid_pos, inner_lab,
            facing     = "bending.inside",
            niceFacing = TRUE,
            adj        = c(0.5, 0),
            cex        = 1,
            font       = 2
          )
        } else {
          circos.text(
            CELL_META$xcenter, y_inner_pos, inner_lab,
            facing     = "bending.inside",
            niceFacing = TRUE,
            adj        = c(0.5, 0),
            cex        = 1,
            font       = 2
          )
          circos.text(
            CELL_META$xcenter, y_outer_pos, outer_lab,
            facing     = "bending.inside",
            niceFacing = TRUE,
            adj        = c(0.5, 0),
            cex        = 1,
            font       = 2
          )
        }
      }
    },
    bg.border = NA
  )
  
  circos.clear()
  dev.off()
  
  cat("\nSaved chord diagram:\n")
  cat(output_pdf, "\n")
  
  return(chord_data)
}


# ============================================================
# 10.9. Chord diagrams
# ============================================================

# Immune / stromal chord
chord_immune_stromal <- make_chord_diagram(
  object = object,
  celltype_keep = c("NonInfMac", "NKTcell", "InfMac", "Bcell", "Bcells", "Stellate"),
  clusters_named_keep = c(
    "C3 Immune cells",
    "C4 Myofibroblasts",
    "C6 Matrix-rich fibroblast"
  ),
  output_pdf = file.path(out_dir, "chordDiagram_immune_stromal_current_names_C3green_C5yellow.pdf"),
  cols_celltype_all = cols_celltype_all,
  final_palette = final_palette
)

# Parenchymal / structural chord
chord_parenchymal <- make_chord_diagram(
  object = object,
  celltype_keep = c("PortalHep", "CentralHep", "InterHep", "Cholangiocyte", "PortalEndo"),
  clusters_named_keep = c(
    "C0 Portal tracts",
    "C1 Periportal hepatocytes",
    "C2 Centrilobular hepatocytes",
    "C5 Cholangiocytes"
  ),
  output_pdf = file.path(out_dir, "chordDiagram_parenchymal_current_names_C3green_C5yellow.pdf"),
  cols_celltype_all = cols_celltype_all,
  final_palette = final_palette,
  celltype_label_override = c(
    PortalHep = "PortHep",
    CentralHep = "CentrHep",
    InterHep = "InterHep",
    Cholangiocyte = "Chol"
  )
)

# Optional third chord: all major cell types versus all final clusters
chord_all_major <- make_chord_diagram(
  object = object,
  celltype_keep = c(
    "PortalHep", "CentralHep", "InterHep",
    "Cholangiocyte", "PortalEndo", "cvLSECs",
    "Stellate", "NonInfMac", "NKTcell", "InfMac", "Bcell", "Bcells"
  ),
  clusters_named_keep = final_order,
  output_pdf = file.path(out_dir, "chordDiagram_all_major_celltypes_all_final_clusters_C3green_C5yellow.pdf"),
  cols_celltype_all = cols_celltype_all,
  final_palette = final_palette,
  celltype_label_override = c(
    PortalHep = "PortHep",
    CentralHep = "CentrHep",
    InterHep = "InterHep",
    Cholangiocyte = "Chol"
  )
)




# ============================================================
# 10.9b. Add fibrosis stage F_score from donor ID
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

stopifnot("donor" %in% colnames(object@meta.data))

donor_to_F_score <- c(
  # F0
  "6B" = "F0",
  "5D" = "F0",
  "4C" = "F0",
  "3B" = "F0",
  "2A" = "F0",
  "1A" = "F0",
  
  # F1
  "6A" = "F1",
  "5A" = "F1",
  "4D" = "F1",
  "3C" = "F1",
  "2B" = "F1",
  "1D" = "F1",
  
  # F2
  "6C" = "F2",
  "5B" = "F2",
  "1C" = "F2",
  
  # F3
  "4A" = "F3",
  "3D" = "F3",
  "2C" = "F3",
  
  # F4
  "6D" = "F4",
  "5C" = "F4",
  "4B" = "F4",
  "3A" = "F4",
  "2D" = "F4",
  "1B" = "F4"
)

object$donor <- as.character(object$donor)

missing_donors <- setdiff(unique(object$donor), names(donor_to_F_score))

if (length(missing_donors) > 0) {
  stop(
    "Some donors in object$donor are missing from donor_to_F_score:\n",
    paste(missing_donors, collapse = ", ")
  )
}

object$F_score <- unname(donor_to_F_score[object$donor])

object$F_score <- factor(
  object$F_score,
  levels = c("F0", "F1", "F2", "F3", "F4")
)

cat("\nF_score added from donor metadata.\n")

cat("\nDonor to F_score table:\n")
print(
  object@meta.data %>%
    dplyr::distinct(donor, F_score) %>%
    dplyr::arrange(F_score, donor)
)

cat("\nSpot counts per F_score:\n")
print(table(object$F_score))

cat("\nDonor counts per F_score:\n")
print(
  object@meta.data %>%
    dplyr::distinct(donor, F_score) %>%
    dplyr::count(F_score)
)

write.csv(
  object@meta.data %>%
    dplyr::distinct(donor, F_score) %>%
    dplyr::arrange(F_score, donor),
  file = file.path(out_dir, "donor_to_F_score_mapping_added_to_object.csv"),
  row.names = FALSE
)

saveRDS(
  object,
  file = file.path(out_dir, "object_full_with_final_annotation_and_F_score.rds")
)


# ============================================================
# 10.10. Cluster proportions by fibrosis stage
# Donor-level proportions + Wilcoxon tests + saved plots/tables
# ============================================================

prop_out_dir <- file.path(out_dir, "cluster_proportion_by_F_score")
dir.create(prop_out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot("cluster" %in% colnames(object@meta.data))
stopifnot("cluster_named" %in% colnames(object@meta.data))
stopifnot("F_score" %in% colnames(object@meta.data))
stopifnot("donor" %in% colnames(object@meta.data))

stage_palette <- c(
  "F0" = "#BDBDBD",
  "F1" = "#F8C8D8",
  "F2" = "#F06292",
  "F3" = "#FF7043",
  "F4" = "#C62828"
)

p_to_stars <- function(p) {
  ifelse(
    is.na(p), "ns",
    ifelse(
      p < 0.001, "***",
      ifelse(
        p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns")
      )
    )
  )
}

make_cluster_prop_analysis <- function(
    object,
    cluster_ids,
    label,
    comparison = c("F0", "F4"),
    stage_palette,
    out_dir,
    y_percent_accuracy = 1
) {
  
  meta <- object@meta.data
  
  meta <- meta %>%
    dplyr::filter(
      !is.na(.data$F_score),
      !is.na(.data$donor),
      !is.na(.data$cluster)
    ) %>%
    dplyr::mutate(
      F_score = as.character(.data$F_score),
      donor = as.character(.data$donor),
      cluster = as.character(.data$cluster)
    )
  
  stage_levels <- intersect(
    c("F0", "F1", "F2", "F3", "F4"),
    sort(unique(meta$F_score))
  )
  
  if (length(stage_levels) < 2) {
    warning("Skipping ", label, ": fewer than two F_score levels.")
    return(NULL)
  }
  
  meta$F_score <- factor(meta$F_score, levels = stage_levels)
  
  cluster_ids <- as.character(cluster_ids)
  
  prop_by_donor <- meta %>%
    dplyr::group_by(F_score, donor) %>%
    dplyr::summarise(
      total_spots = dplyr::n(),
      selected_spots = sum(cluster %in% cluster_ids),
      prop = ifelse(total_spots > 0, selected_spots / total_spots, NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(prop))
  
  if (nrow(prop_by_donor) == 0) {
    warning("Skipping ", label, ": no donor-level proportions.")
    return(NULL)
  }
  
  write.csv(
    prop_by_donor,
    file = file.path(out_dir, paste0(label, "_donor_level_proportions.csv")),
    row.names = FALSE
  )
  
  stage_summary <- prop_by_donor %>%
    dplyr::group_by(F_score) %>%
    dplyr::summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      sd_prop = sd(prop, na.rm = TRUE),
      n = dplyr::n(),
      se_prop = sd_prop / sqrt(n),
      .groups = "drop"
    )
  
  write.csv(
    stage_summary,
    file = file.path(out_dir, paste0(label, "_stage_summary_mean_SEM.csv")),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # Selected comparison, e.g. F0/F4 or F1/F4
  # ------------------------------------------------------------
  
  comparison <- as.character(comparison)
  comparison <- comparison[comparison %in% stage_levels]
  
  pval_main <- NA_real_
  lab_text <- "p=NA (ns)"
  
  if (length(comparison) == 2) {
    sub_comp <- prop_by_donor %>%
      dplyr::filter(F_score %in% comparison) %>%
      droplevels()
    
    if (length(unique(sub_comp$F_score)) == 2) {
      wt <- suppressWarnings(
        wilcox.test(
          prop ~ F_score,
          data = sub_comp,
          exact = FALSE,
          correct = FALSE
        )
      )
      pval_main <- wt$p.value
      lab_text <- paste0("p=", signif(pval_main, 2), " (", p_to_stars(pval_main), ")")
    }
  }
  
  main_test_df <- data.frame(
    label = label,
    clusters = paste(cluster_ids, collapse = "+"),
    group1 = ifelse(length(comparison) >= 1, comparison[1], NA),
    group2 = ifelse(length(comparison) >= 2, comparison[2], NA),
    p_raw = pval_main,
    p_label = lab_text,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    main_test_df,
    file = file.path(out_dir, paste0(label, "_main_comparison_wilcox.csv")),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # All pairwise Wilcoxon tests
  # ------------------------------------------------------------
  
  pair_tests <- lapply(combn(stage_levels, 2, simplify = FALSE), function(g) {
    
    g1 <- g[1]
    g2 <- g[2]
    
    sub <- prop_by_donor %>%
      dplyr::filter(F_score %in% c(g1, g2)) %>%
      droplevels()
    
    if (length(unique(sub$F_score)) < 2) {
      return(data.frame(
        group1 = g1,
        group2 = g2,
        p_raw = NA_real_,
        n1 = sum(sub$F_score == g1),
        n2 = sum(sub$F_score == g2)
      ))
    }
    
    wt <- suppressWarnings(
      wilcox.test(
        prop ~ F_score,
        data = sub,
        exact = FALSE,
        correct = FALSE
      )
    )
    
    data.frame(
      group1 = g1,
      group2 = g2,
      p_raw = wt$p.value,
      n1 = sum(sub$F_score == g1),
      n2 = sum(sub$F_score == g2)
    )
  })
  
  pair_df <- dplyr::bind_rows(pair_tests) %>%
    dplyr::mutate(
      p_adj = p.adjust(p_raw, method = "holm"),
      label = paste0("p=", signif(p_adj, 2), " (", p_to_stars(p_adj), ")")
    )
  
  write.csv(
    pair_df,
    file = file.path(out_dir, paste0(label, "_all_pairwise_Wilcoxon_recap.csv")),
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # Plot 1: selected comparison bracket
  # ------------------------------------------------------------
  
  tops <- stage_summary %>%
    dplyr::mutate(top = mean_prop + ifelse(is.na(se_prop), 0, se_prop))
  
  y_top <- max(tops$top, na.rm = TRUE)
  if (!is.finite(y_top) || y_top <= 0) y_top <- max(stage_summary$mean_prop, na.rm = TRUE)
  if (!is.finite(y_top) || y_top <= 0) y_top <- 0.01
  
  y_bracket <- min(1, y_top * 1.12)
  y_label <- min(1, y_bracket * 1.04)
  y_max <- min(1, y_label * 1.08)
  if (y_max <= y_label) y_max <- y_label + 0.01
  
  p_main <- ggplot(
    stage_summary,
    aes(x = F_score, y = mean_prop, fill = F_score)
  ) +
    geom_col(width = 0.65) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_prop - ifelse(is.na(se_prop), 0, se_prop), 0),
        ymax = pmin(mean_prop + ifelse(is.na(se_prop), 0, se_prop), 1)
      ),
      width = 0.18
    ) +
    scale_fill_manual(values = stage_palette, guide = "none", drop = FALSE) +
    scale_y_continuous(
      labels = percent_format(accuracy = y_percent_accuracy),
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0(label, " proportion across fibrosis stages"),
      x = NULL,
      y = paste0("n(", label, " spots) / n(total spots)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 11)
    )
  
  if (length(comparison) == 2 && all(comparison %in% stage_levels)) {
    
    x1 <- which(stage_levels == comparison[1])
    x2 <- which(stage_levels == comparison[2])
    h_leg <- y_top * 0.03
    
    p_main <- p_main +
      annotate(
        "segment",
        x = x1,
        xend = x1,
        y = y_bracket - h_leg,
        yend = y_bracket
      ) +
      annotate(
        "segment",
        x = x2,
        xend = x2,
        y = y_bracket - h_leg,
        yend = y_bracket
      ) +
      annotate(
        "segment",
        x = x1,
        xend = x2,
        y = y_bracket,
        yend = y_bracket
      ) +
      annotate(
        "text",
        x = mean(c(x1, x2)),
        y = y_label,
        label = lab_text,
        size = 4.2
      )
  }
  
  print(p_main)
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(label, "_", paste(comparison, collapse = "_vs_"), "_wilcox.pdf")
    ),
    plot = p_main,
    width = 4.6,
    height = 4.2,
    units = "in",
    device = cairo_pdf,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(label, "_", paste(comparison, collapse = "_vs_"), "_wilcox.png")
    ),
    plot = p_main,
    width = 4.6,
    height = 4.2,
    units = "in",
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # Plot 2: adjacent comparisons only
  # ------------------------------------------------------------
  
  adj_df <- pair_df %>%
    dplyr::mutate(
      idx1 = match(group1, stage_levels),
      idx2 = match(group2, stage_levels)
    ) %>%
    dplyr::filter(idx2 - idx1 == 1) %>%
    dplyr::arrange(idx1)
  
  if (nrow(adj_df) > 0) {
    
    tops2 <- stage_summary %>%
      dplyr::mutate(top = mean_prop + ifelse(is.na(se_prop), 0, se_prop)) %>%
      dplyr::select(F_score, top)
    
    y_base <- vapply(2:length(stage_levels), function(i) {
      max(tops2$top[tops2$F_score %in% c(stage_levels[i - 1], stage_levels[i])], na.rm = TRUE)
    }, numeric(1))
    
    y_pos <- y_base * 1.08 + (seq_along(y_base) - 1) * 0.035
    adj_df$y.position <- pmin(y_pos, 1.0)
    
    y_max_adj <- max(adj_df$y.position, na.rm = TRUE) * 1.08
    if (!is.finite(y_max_adj) || y_max_adj <= 0) y_max_adj <- y_max
    
    p_adjacent <- ggplot(
      stage_summary,
      aes(x = F_score, y = mean_prop, fill = F_score)
    ) +
      geom_col(width = 0.65) +
      geom_errorbar(
        aes(
          ymin = pmax(mean_prop - ifelse(is.na(se_prop), 0, se_prop), 0),
          ymax = pmin(mean_prop + ifelse(is.na(se_prop), 0, se_prop), 1)
        ),
        width = 0.18
      ) +
      scale_fill_manual(values = stage_palette, guide = "none", drop = FALSE) +
      scale_y_continuous(
        labels = percent_format(accuracy = y_percent_accuracy),
        limits = c(0, min(1, y_max_adj)),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(
        title = paste0(label, " proportion across fibrosis stages"),
        x = "F score",
        y = paste0(label, " spots (proportion per donor)")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )
    
    # Manual adjacent brackets, avoiding ggpubr aesthetic inheritance problems

for (kk in seq_len(nrow(adj_df))) {
  
  x1 <- match(adj_df$group1[kk], stage_levels)
  x2 <- match(adj_df$group2[kk], stage_levels)
  y  <- adj_df$y.position[kk]
  h  <- y * 0.025
  
  p_adjacent <- p_adjacent +
    annotate(
      "segment",
      x = x1,
      xend = x1,
      y = y - h,
      yend = y
    ) +
    annotate(
      "segment",
      x = x2,
      xend = x2,
      y = y - h,
      yend = y
    ) +
    annotate(
      "segment",
      x = x1,
      xend = x2,
      y = y,
      yend = y
    ) +
    annotate(
      "text",
      x = mean(c(x1, x2)),
      y = y + h * 1.5,
      label = adj_df$label[kk],
      size = 3.4
    )
}
    
    print(p_adjacent)
    
    ggsave(
      filename = file.path(out_dir, paste0(label, "_adjacent_stage_wilcox.pdf")),
      plot = p_adjacent,
      width = 5.5,
      height = 4.2,
      units = "in",
      device = cairo_pdf,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(out_dir, paste0(label, "_adjacent_stage_wilcox.png")),
      plot = p_adjacent,
      width = 5.5,
      height = 4.2,
      units = "in",
      dpi = 300
    )
  }
  
  return(
    list(
      donor_level = prop_by_donor,
      stage_summary = stage_summary,
      main_test = main_test_df,
      pairwise_tests = pair_df,
      plot_main = p_main
    )
  )
}

# ============================================================
# Run proportion analyses
# ============================================================

cluster_prop_specs <- list(
  list(cluster_ids = "C0", label = "C0_Portal_tracts", comparison = c("F1", "F4")),
  list(cluster_ids = "C1", label = "C1_Periportal_hepatocytes", comparison = c("F0", "F4")),
  list(cluster_ids = "C2", label = "C2_Centrilobular_hepatocytes", comparison = c("F0", "F4")),
  list(cluster_ids = "C3", label = "C3_Immune_cells", comparison = c("F0", "F3")),
  list(cluster_ids = "C4", label = "C4_Myofibroblasts", comparison = c("F0", "F4")),
  list(cluster_ids = "C5", label = "C5_Cholangiocytes", comparison = c("F0", "F4")),
  list(cluster_ids = "C6", label = "C6_Matrix_rich_fibroblast", comparison = c("F0", "F4")),
  list(cluster_ids = c("C4", "C6"), label = "C4C6_Myofibroblast_MatrixRich", comparison = c("F0", "F4"))
)

cluster_prop_results <- list()

for (spec in cluster_prop_specs) {
  
  cat("\n============================================================\n")
  cat("Running proportion analysis:", spec$label, "\n")
  cat("Clusters:", paste(spec$cluster_ids, collapse = ", "), "\n")
  cat("Comparison:", paste(spec$comparison, collapse = " vs "), "\n")
  cat("============================================================\n")
  
  cluster_prop_results[[spec$label]] <- make_cluster_prop_analysis(
    object = object,
    cluster_ids = spec$cluster_ids,
    label = spec$label,
    comparison = spec$comparison,
    stage_palette = stage_palette,
    out_dir = prop_out_dir,
    y_percent_accuracy = 1
  )
}

# ============================================================
# Combined tables across all proportion analyses
# ============================================================

all_main_tests <- dplyr::bind_rows(
  lapply(names(cluster_prop_results), function(nm) {
    x <- cluster_prop_results[[nm]]
    if (is.null(x)) return(NULL)
    x$main_test
  })
)

all_pairwise_tests <- dplyr::bind_rows(
  lapply(names(cluster_prop_results), function(nm) {
    x <- cluster_prop_results[[nm]]
    if (is.null(x)) return(NULL)
    x$pairwise_tests %>%
      dplyr::mutate(label_cluster = nm, .before = 1)
  })
)

all_stage_summaries <- dplyr::bind_rows(
  lapply(names(cluster_prop_results), function(nm) {
    x <- cluster_prop_results[[nm]]
    if (is.null(x)) return(NULL)
    x$stage_summary %>%
      dplyr::mutate(label_cluster = nm, .before = 1)
  })
)

all_donor_level_props <- dplyr::bind_rows(
  lapply(names(cluster_prop_results), function(nm) {
    x <- cluster_prop_results[[nm]]
    if (is.null(x)) return(NULL)
    x$donor_level %>%
      dplyr::mutate(label_cluster = nm, .before = 1)
  })
)

write.csv(
  all_main_tests,
  file = file.path(prop_out_dir, "ALL_selected_comparison_Wilcoxon_tests.csv"),
  row.names = FALSE
)

write.csv(
  all_pairwise_tests,
  file = file.path(prop_out_dir, "ALL_pairwise_Wilcoxon_tests_Holm_corrected.csv"),
  row.names = FALSE
)

write.csv(
  all_stage_summaries,
  file = file.path(prop_out_dir, "ALL_stage_summary_mean_SEM.csv"),
  row.names = FALSE
)

write.csv(
  all_donor_level_props,
  file = file.path(prop_out_dir, "ALL_donor_level_cluster_proportions.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("DONE: cluster proportion by F_score analyses\n")
cat("Saved in:\n")
cat(prop_out_dir, "\n")
cat("============================================================\n")

# ============================================================
# 10.11. Final save and file check
# ============================================================

saveRDS(
  object,
  file = file.path(out_dir, "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds")
)

sink(file.path(out_dir, "sessionInfo_final_full_chunk_C3green_C5yellow.txt"))
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("DONE: final annotation + markers + GO + deconv + chord diagrams + proportions\n")
cat("All outputs saved in:\n")
cat(out_dir, "\n")
cat("============================================================\n")

cat("\nGenerated files in main output directory:\n")
print(list.files(out_dir))

if (exists("prop_out_dir")) {
  cat("\nGenerated files in proportion output directory:\n")
  print(list.files(prop_out_dir))
}

# ============================================================
# 10.12. Reproducibility checksum of all final outputs
# ============================================================

all_output_files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
all_output_files <- all_output_files[file.exists(all_output_files)]
all_output_files <- all_output_files[file.info(all_output_files)$isdir == FALSE]

write_checksum_csv(
  all_output_files,
  file.path(repro_dir, "ALL_output_file_md5_checksums_after_06_to_10.csv")
)

sink(file.path(repro_dir, "sessionInfo_after_06_to_10.txt"))
cat("Timestamp:\n"); print(Sys.time())
cat("\nR version:\n"); print(R.version.string)
cat("\n.libPaths():\n"); print(.libPaths())
cat("\nSys.info():\n"); print(Sys.info())
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

write.csv(
  installed.packages()[, c("Package", "Version", "LibPath", "Built")],
  file.path(repro_dir, "installed_packages_after_06_to_10.csv"),
  row.names = FALSE
)

cat("\nReproducibility checksums and sessionInfo saved in:\n")
cat(repro_dir, "\n")
# END SPATIAL DOMAIN ANALYSIS
}, envir = run_env)

log_msg("Completed spatial clustering, annotation and biological characterization")

out_dir <- file.path(RESULTS_DIR, "260701_approach1_cluster_go_marker_deconv")
if (dir.exists(out_dir)) {
  all_outputs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  all_outputs <- all_outputs[file.exists(all_outputs)]
  all_outputs <- all_outputs[file.info(all_outputs)$isdir == FALSE]

  write_checksum_csv(
    all_outputs,
    file.path(REPRO_LOG_DIR, "ALL_outputs_06_to_10_as_Rscript_md5.csv")
  )
}

sink(file.path(REPRO_LOG_DIR, "sessionInfo_after_06_to_10_as_Rscript.txt"))
cat("Timestamp:\n")
print(Sys.time())
cat("\nR version:\n")
print(R.version.string)
cat("\n.libPaths():\n")
print(.libPaths())
cat("\nSys.info():\n")
print(Sys.info())
cat("\nsessionInfo():\n")
print(sessionInfo())
sink()

write.csv(
  installed.packages()[, c("Package", "Version", "LibPath", "Built")],
  file.path(REPRO_LOG_DIR, "installed_packages_after_06_to_10_as_Rscript.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("DONE: spatial clustering, annotation and biological characterization completed\n")
cat("Results folder:\n")
cat(out_dir, "\n")
cat("Reproducibility logs:\n")
cat(REPRO_LOG_DIR, "\n")
cat("============================================================\n")
