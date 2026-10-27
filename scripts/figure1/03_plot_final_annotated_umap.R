#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

repo_dir <- "/Users/jianwu/Documents/FLINC"

analysis_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/results/260701_approach1_cluster_go_marker_deconv"
)

object_path <- file.path(
  analysis_dir,
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)

script_file <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/scripts/figure1/03_plot_final_annotated_umap.R"
)

log_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/reproducibility_logs/UMAP_final_C5red"
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(object_path)) {
  stop("Final annotated object not found: ", object_path)
}

cat("Loading final annotated object...\n")
object <- readRDS(object_path)

stopifnot(
  "umap_SCT" %in% Reductions(object)
)

required_meta <- c(
  "cluster",
  "cluster_named"
)

missing_meta <- setdiff(
  required_meta,
  colnames(object@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

final_order <- c(
  "C0 Portal tracts",
  "C1 Periportal hepatocytes",
  "C2 Centrilobular hepatocytes",
  "C3 Immune cells",
  "C4 Myofibroblasts",
  "C5 Cholangiocytes",
  "C6 Matrix-rich fibroblast"
)

# Same final palette as the TRUE Visium cluster overlays.
# Only C5 differs from the historical yellow version.
final_palette <- c(
  C0 = "#56B4E9",
  C1 = "#FF9DA7",
  C2 = "#B189A8",
  C3 = "#59A14F",
  C4 = "#F28E2B",
  C5 = "#E15759",
  C6 = "#E5C494"
)

palette_named <- setNames(
  unname(final_palette),
  final_order
)

emb <- Embeddings(
  object,
  "umap_SCT"
)

umap_df <- data.frame(
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  cluster = as.character(object$cluster),
  cluster_named = factor(
    as.character(object$cluster_named),
    levels = final_order
  ),
  stringsAsFactors = FALSE
)

# Preserve EXACT plotting order from the original final UMAP:
# large/background populations first,
# C6 before niche populations,
# C4 last/on top.
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

umap_df <- umap_df[
  order(umap_df$plot_order),
  ,
  drop = FALSE
]

centers <- aggregate(
  umap_df[, c("UMAP_1", "UMAP_2")],
  by = list(
    cluster = umap_df$cluster
  ),
  FUN = median,
  na.rm = TRUE
)

p_umap <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = cluster_named
  )
) +
  geom_point(
    size = 0.55,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = palette_named,
    name = "Annotation",
    drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = cluster
    ),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.alpha = 0.6,
    fontface = "bold"
  ) +
  ggtitle(
    "Final annotated spatial clusters"
  ) +
  theme_minimal(
    base_size = 14
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5,
        alpha = 1
      )
    )
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    legend.key.width = unit(
      8,
      "pt"
    ),
    legend.key.height = unit(
      8,
      "pt"
    ),
    legend.text = element_text(
      size = 11
    ),
    legend.title = element_text(
      size = 12,
      face = "bold"
    )
  )

pdf_file <- file.path(
  analysis_dir,
  "UMAP_final_annotated_clusters_C4_on_top_C3green_C5red.pdf"
)

png_file <- file.path(
  analysis_dir,
  "UMAP_final_annotated_clusters_C4_on_top_C3green_C5red.png"
)

ggsave(
  filename = pdf_file,
  plot = p_umap,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = p_umap,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

write.csv(
  data.frame(
    cluster = names(final_palette),
    color = unname(final_palette),
    stringsAsFactors = FALSE
  ),
  file.path(
    log_dir,
    "cluster_palette_C5red.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    plotting_position = seq_along(plot_order),
    cluster_named = plot_order,
    stringsAsFactors = FALSE
  ),
  file.path(
    log_dir,
    "UMAP_plot_order.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    table(object$cluster)
  ),
  file.path(
    log_dir,
    "final_cluster_counts.csv"
  ),
  row.names = FALSE
)

output_md5 <- tools::md5sum(
  c(
    script_file,
    pdf_file,
    png_file
  )
)

write.csv(
  data.frame(
    file = names(output_md5),
    md5 = unname(output_md5),
    stringsAsFactors = FALSE
  ),
  file.path(
    log_dir,
    "script_and_output_md5.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    paste0(
      "timestamp: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S %Z"
      )
    ),
    paste0(
      "input_object: ",
      object_path
    ),
    "reduction: umap_SCT",
    "UMAP coordinates recalculated: NO",
    "cluster assignments modified: NO",
    "plotting-only change: YES",
    "C5 historical color: #EDC948",
    "C5 final color: #E15759",
    "point_size: 0.55",
    "point_alpha: 0.9",
    "C4 plotted last/on top: YES",
    paste0(
      "pdf_file: ",
      pdf_file
    ),
    paste0(
      "png_file: ",
      png_file
    )
  ),
  file.path(
    log_dir,
    "run_parameters.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    log_dir,
    "sessionInfo.txt"
  )
)

cat("\nDONE: FINAL C5-RED UMAP\n")
cat("PDF:\n", pdf_file, "\n")
cat("PNG:\n", png_file, "\n")
cat("No Seurat object was modified or saved.\n")
