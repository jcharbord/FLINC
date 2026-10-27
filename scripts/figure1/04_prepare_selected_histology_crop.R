#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(png)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- match(name, args)
  if (is.na(hit)) {
    return(default)
  }
  if (hit == length(args)) {
    stop("Missing value after ", name)
  }
  args[[hit + 1]]
}

has_flag <- function(name) {
  name %in% args
}

script_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(script_arg) != 1) {
  stop("Could not determine the executing script path")
}

script_path <- normalizePath(
  sub("^--file=", "", script_arg),
  mustWork = TRUE
)

default_repro_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = FALSE
)

repro_root <- normalizePath(
  arg_value("--repro-root", default_repro_root),
  mustWork = TRUE
)

flinc_root <- normalizePath(
  file.path(repro_root, "..", ".."),
  mustWork = FALSE
)

if (has_flag("--help")) {
  cat(paste(
    "Usage:",
    "  Rscript 04_prepare_selected_histology_crop.R [options]",
    "",
    "Options:",
    "  --repro-root PATH  Reproducible_code directory",
    "  --object PATH      Final annotated Seurat object",
    "  --hires-png PATH   Sample 3D tissue_hires_image.png",
    "  --fullres-he-png PATH  Optional H&E source override; default = registered Space Ranger hires PNG",
    "  --cmap PATH        color_pseudotime.csv",
    "  --output-dir PATH  Stable output directory",
    "  --overwrite        Allow replacement of expected outputs",
    "  --help             Show this message",
    sep = "\n"
  ))
  quit(save = "no", status = 0)
}

# ============================================================
# SOURCE REPRODUCIBLE OBJECT
# ============================================================

object_path <- file.path(
  repro_root,
  "results/260701_approach1_cluster_go_marker_deconv",
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)

object_path <- normalizePath(
  arg_value("--object", object_path),
  mustWork = TRUE
)

donor_use <- "3D"
image_key <- "X3D"

crop_id <- "F3_sample3D_final"

crop_coordinates <- data.frame(
  crop_id = crop_id,
  xmin = 520,
  xmax = 1040,
  ymin = 1020,
  ymax = 1440,
  stringsAsFactors = FALSE
)

cluster_colors_red <- c(
  C0 = "#56B4E9",
  C1 = "#FF9DA7",
  C2 = "#B189A8",
  C3 = "#59A14F",
  C4 = "#F28E2B",
  C5 = "#E15759",
  C6 = "#E5C494"
)

cluster_colors_yellow <- cluster_colors_red
cluster_colors_yellow[["C5"]] <- "#EDC948"

cluster_labels_short <- setNames(
  names(cluster_colors_red),
  names(cluster_colors_red)
)

# ============================================================
# EXACT TOP 10 MARKERS PER FINAL CLUSTER
# from Top10_markers_per_final_cluster_dotplot_data.csv
# ============================================================

markers_by_cluster <- list(

  C0 = c(
    "AQP1",
    "KRT7",
    "C7",
    "IGFBP7",
    "SFRP5",
    "CHST4",
    "LTBP4",
    "FXYD2",
    "DCDC2",
    "ANXA4"
  ),

  C1 = c(
    "SDS",
    "HAL",
    "HAMP",
    "GLS2",
    "NNMT",
    "AMDHD1",
    "CXCL2",
    "SDSL",
    "SLPI",
    "CRP"
  ),

  C2 = c(
    "SLCO1B3",
    "CYP1A2",
    "CYP3A4",
    "GLUL",
    "SLCO1B7",
    "NOTUM",
    "OAT",
    "TBX3",
    "PLIN1",
    "CYP2E1"
  ),

  C3 = c(
    "TMSB4X",
    "CD74",
    "HLA-DRA",
    "LSP1",
    "ACTB",
    "LYZ",
    "HLA-DPB1",
    "HLA-DPA1",
    "PTPRC",
    "VIM"
  ),

  C4 = c(
    "LMOD1",
    "SEMA3G",
    "CCL19",
    "SMOC2",
    "MAMDC2",
    "ACTA2",
    "ID4",
    "ITGB4",
    "PARM1",
    "NOTCH3"
  ),

  C5 = c(
    "SOD3",
    "S100A6",
    "BICC1",
    "KRT7",
    "SFRP5",
    "AQP1",
    "DPT",
    "CHST4",
    "C7",
    "IGFBP7"
  ),

  C6 = c(
    "PTGDS",
    "AEBP1",
    "IGFBP7",
    "C7",
    "COL1A2",
    "CD74",
    "DCN",
    "DPT",
    "COL3A1",
    "HLA-DRA"
  )
)

# ============================================================
# HIRES H&E + EXPRESSION COLOR MAP
# ============================================================

hires_png <- file.path(
  dirname(flinc_root),
  "DATA/spaceranger_output",
  donor_use,
  "spatial",
  "tissue_hires_image.png"
)

hires_png <- normalizePath(
  arg_value("--hires-png", hires_png),
  mustWork = TRUE
)

# ============================================================
# CANONICAL H&E SOURCE
#
# Use the registered Space Ranger hires H&E for sample 3D.
# This is the same coordinate system used for the cluster and
# marker overlays.
#
# A copy of this source image is also placed in H_E/ so the
# final manuscript crop can be prepared manually in Photoshop
# without relying on an external duplicate image.
# ============================================================

fullres_he_png <- normalizePath(
  arg_value("--fullres-he-png", hires_png),
  mustWork = TRUE
)

cmap_file <- file.path(
  flinc_root,
  "data/st_analysis_pipeline/color_pseudotime.csv"
)

cmap_file <- normalizePath(
  arg_value("--cmap", cmap_file),
  mustWork = TRUE
)

# ============================================================
# OUTPUT FOLDERS
# ============================================================

output_root <- normalizePath(
  arg_value(
    "--output-dir",
    file.path(repro_root, "results/Fig1_Cropped_F3_sample3D")
  ),
  mustWork = FALSE
)

overwrite <- has_flag("--overwrite")

if (
  dir.exists(output_root) &&
    length(list.files(output_root, all.files = FALSE)) > 0 &&
    !overwrite
) {
  stop(
    "Output directory is not empty: ",
    output_root,
    "\nUse --overwrite only after verifying the target."
  )
}

markers_dir <- file.path(output_root, "Top10_markers")
he_dir <- file.path(output_root, "H_E")
clusters_red_dir <- file.path(output_root, "clusters_C5red")
clusters_yellow_dir <- file.path(output_root, "clusters_C5yellow")
records_dir <- file.path(output_root, "reproducibility_records")

for (path_use in c(
  output_root,
  markers_dir,
  he_dir,
  clusters_red_dir,
  clusters_yellow_dir,
  records_dir
)) {
  dir.create(path_use, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# COPY REGISTERED H&E SOURCE FOR MANUAL FIGURE ASSEMBLY
# ============================================================

he_source_copy <- file.path(
  he_dir,
  paste0(donor_use, "_tissue_hires_image_SOURCE.png")
)

if (!file.copy(
  hires_png,
  he_source_copy,
  overwrite = TRUE
)) {
  stop(
    "Could not copy registered Space Ranger H&E source to: ",
    he_source_copy
  )
}

out_dir <- markers_dir
log_dir <- records_dir

# Make C0-C6 subfolders
for (cluster_use in names(markers_by_cluster)) {

  dir.create(
    file.path(out_dir, cluster_use),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

write.csv(
  crop_coordinates,
  file.path(output_root, "crop_coordinates.csv"),
  row.names = FALSE
)

saveRDS(
  crop_coordinates,
  file.path(records_dir, "crop_coordinates.rds")
)

write.csv(
  data.frame(
    cluster = names(cluster_colors_red),
    color = unname(cluster_colors_red),
    stringsAsFactors = FALSE
  ),
  file.path(records_dir, "cluster_palette_C5red.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    cluster = names(cluster_colors_yellow),
    color = unname(cluster_colors_yellow),
    stringsAsFactors = FALSE
  ),
  file.path(records_dir, "cluster_palette_C5yellow.csv"),
  row.names = FALSE
)

# ============================================================
# LOAD OBJECT
# ============================================================

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for SHA256 validation")
}

object_sha256_before <- digest::digest(
  file = object_path,
  algo = "sha256"
)

hires_png_sha256 <- digest::digest(
  file = hires_png,
  algo = "sha256"
)

he_source_copy_sha256 <- digest::digest(
  file = he_source_copy,
  algo = "sha256"
)

if (!identical(hires_png_sha256, he_source_copy_sha256)) {
  stop("Copied H&E source does not match the Space Ranger source image")
}

fullres_he_png_sha256 <- digest::digest(
  file = fullres_he_png,
  algo = "sha256"
)

script_sha256 <- digest::digest(
  file = script_path,
  algo = "sha256"
)

cat("Loading final annotated Seurat object...\n")
cat(object_path, "\n")

object <- readRDS(object_path)

stopifnot("SCT" %in% Assays(object))
stopifnot(image_key %in% names(object@images))
stopifnot(file.exists(hires_png))
stopifnot(file.exists(fullres_he_png))
stopifnot(file.exists(cmap_file))

cat("\nSCT layers:\n")
print(Layers(object[["SCT"]]))

# ============================================================
# SCT NORMALIZED EXPRESSION
# ============================================================

expr_mat <- GetAssayData(
  object,
  assay = "SCT",
  layer = "data"
)

cat("\nSCT expression matrix dimensions:\n")
print(dim(expr_mat))

# ============================================================
# SAVE EXACT TOP-10 TABLE USED BY THIS SCRIPT
# ============================================================

marker_table <- do.call(
  rbind,
  lapply(
    names(markers_by_cluster),
    function(cl) {

      data.frame(
        cluster = cl,
        rank = seq_along(markers_by_cluster[[cl]]),
        gene = markers_by_cluster[[cl]],
        stringsAsFactors = FALSE
      )
    }
  )
)

marker_table$found_in_SCT <-
  marker_table$gene %in% rownames(expr_mat)

write.csv(
  marker_table,
  file.path(
    out_dir,
    "Top10_markers_per_cluster_used_for_spatial_plots.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("TOP MARKER CHECK\n")
cat("============================================\n")

print(
  table(
    marker_table$cluster,
    marker_table$found_in_SCT
  )
)

missing_table <- marker_table[
  !marker_table$found_in_SCT,
  ,
  drop = FALSE
]

if (nrow(missing_table) > 0) {

  cat("\nMarkers missing from SCT assay:\n")
  print(missing_table)

} else {

  cat("\nAll 70 cluster-marker entries were found in SCT.\n")
}

# ============================================================
# DONOR 3D COORDINATES
# ============================================================

img <- object@images[[image_key]]

coords0 <- img@coordinates
coords0$cell <- rownames(coords0)

shared_cells <- intersect(
  coords0$cell,
  colnames(expr_mat)
)

if (length(shared_cells) == 0) {
  stop(
    "No shared donor 3D spots between image coordinates and SCT matrix."
  )
}

coords0 <- coords0[
  shared_cells,
  ,
  drop = FALSE
]

cat(
  "\nNumber of donor 3D spots: ",
  nrow(coords0),
  "\n",
  sep = ""
)

# ============================================================
# LOAD REGISTERED SPACE RANGER H&E
# ============================================================

he <- readPNG(hires_png, native = TRUE)

img_h <- dim(he)[1]
img_w <- dim(he)[2]

sf_hires <- img@scale.factors$hires

coords0$x_hires <-
  coords0$imagecol * sf_hires

coords0$y_hires <-
  coords0$imagerow * sf_hires

# Image coordinates start top-left.
# ggplot coordinates start bottom-left.
coords0$y_plot <-
  img_h - coords0$y_hires

# ============================================================
# TRUE VISIUM SPOT GEOMETRY
# ============================================================

spot_diameter_fullres <-
  img@scale.factors$spot

spot_diameter_hires <-
  spot_diameter_fullres * sf_hires

spot_radius_hires <-
  spot_diameter_hires / 2

cat("\n============================================\n")
cat("TRUE VISIUM SPOT GEOMETRY\n")
cat("============================================\n")

cat(
  "spot_diameter_fullres: ",
  spot_diameter_fullres,
  " pixels\n",
  sep = ""
)

cat(
  "hires scale factor: ",
  sf_hires,
  "\n",
  sep = ""
)

cat(
  "TRUE Visium diameter on hires image: ",
  spot_diameter_hires,
  " pixels\n",
  sep = ""
)

cat(
  "spot multiplier: 1.00\n"
)

# ============================================================
# EXPRESSION COLOR MAP
# ============================================================

cmap <- read.csv(
  cmap_file,
  stringsAsFactors = FALSE
)

if (!("x" %in% colnames(cmap))) {
  stop(
    "Expected column 'x' in color_pseudotime.csv"
  )
}

hex_cols <- cmap$x

# Same behavior as previous marker-expression plots
if (length(hex_cols) > 2) {

  hex_cols <- hex_cols[
    2:(length(hex_cols) - 1)
  ]
}

# ============================================================
# TRUE CIRCULAR SPOT GEOMETRY
# ============================================================

n_circle_points <- 48

theta <- seq(
  0,
  2 * pi,
  length.out = n_circle_points + 1
)[-(n_circle_points + 1)]

# ============================================================
# FINAL CROP, H&E, AND CLUSTER PANELS
# ============================================================

cc <- crop_coordinates[1, , drop = FALSE]

crop_width <- cc$xmax - cc$xmin
crop_height <- cc$ymax - cc$ymin

panel_width_in <- 6
panel_height_in <- panel_width_in * crop_height / crop_width
dpi_out <- 300

# The registered Space Ranger hires H&E is the source for the H&E-only
# panel and the coordinate reference for all spatial overlays.
fullres_he <- readPNG(fullres_he_png, native = TRUE)

fullres_img_h <- dim(fullres_he)[1]
fullres_img_w <- dim(fullres_he)[2]

hires_to_fullres_x <- fullres_img_w / img_w
hires_to_fullres_y <- fullres_img_h / img_h

fullres_col_start <- max(
  1L,
  floor(cc$xmin * hires_to_fullres_x) + 1L
)

fullres_col_end <- min(
  fullres_img_w,
  ceiling(cc$xmax * hires_to_fullres_x)
)

# Convert bottom-left ggplot y coordinates to top-left raster row indices.
fullres_row_start <- max(
  1L,
  floor((img_h - cc$ymax) * hires_to_fullres_y) + 1L
)

fullres_row_end <- min(
  fullres_img_h,
  ceiling((img_h - cc$ymin) * hires_to_fullres_y)
)

fullres_he_crop <- fullres_he[
  fullres_row_start:fullres_row_end,
  fullres_col_start:fullres_col_end,
  drop = FALSE
]

# Matrix subsetting can drop the nativeRaster class; restore it explicitly so
# grid renders the packed RGBA values directly and without resampling first.
class(fullres_he_crop) <- "nativeRaster"

rm(fullres_he)
invisible(gc())

cat("\n============================================\n")
cat("REGISTERED SPACE RANGER H&E CROP\n")
cat("============================================\n")
cat("Space Ranger hires dimensions: ", img_w, " x ", img_h, "\n", sep = "")
cat(
  "H&E source dimensions: ",
  fullres_img_w,
  " x ",
  fullres_img_h,
  "\n",
  sep = ""
)
cat(
  "Coordinate scale x/y: ",
  format(hires_to_fullres_x, digits = 10),
  " / ",
  format(hires_to_fullres_y, digits = 10),
  "\n",
  sep = ""
)
cat(
  "H&E source crop columns: ",
  fullres_col_start,
  "-",
  fullres_col_end,
  "; rows: ",
  fullres_row_start,
  "-",
  fullres_row_end,
  "\n",
  sep = ""
)

save_plot_pair <- function(plot_use, directory, stem) {
  pdf_file <- file.path(directory, paste0(stem, ".pdf"))
  png_file <- file.path(directory, paste0(stem, ".png"))

  ggsave(
    filename = pdf_file,
    plot = plot_use,
    width = panel_width_in,
    height = panel_height_in,
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = png_file,
    plot = plot_use,
    width = panel_width_in,
    height = panel_height_in,
    dpi = dpi_out,
    device = "png",
    bg = "white"
  )

  c(pdf_file, png_file)
}

# A Visium spot is 55 micrometres in diameter. This converts the
# Space Ranger spot diameter into a 200-micrometre H&E scale bar.
visium_spot_diameter_um <- 55
scale_bar_um <- 200
scale_bar_px <- spot_diameter_hires * (
  scale_bar_um / visium_spot_diameter_um
)

scale_bar_xend <- cc$xmax - 20
scale_bar_xstart <- scale_bar_xend - scale_bar_px
scale_bar_y <- cc$ymin + 25
scale_bar_label_y <- scale_bar_y + 14

p_he <- ggplot() +
  annotation_raster(
    fullres_he_crop,
    xmin = cc$xmin,
    xmax = cc$xmax,
    ymin = cc$ymin,
    ymax = cc$ymax,
    interpolate = TRUE
  ) +
  annotate(
    "segment",
    x = scale_bar_xstart,
    xend = scale_bar_xend,
    y = scale_bar_y,
    yend = scale_bar_y,
    linewidth = 1.2,
    lineend = "butt",
    colour = "black"
  ) +
  annotate(
    "text",
    x = (scale_bar_xstart + scale_bar_xend) / 2,
    y = scale_bar_label_y,
    label = paste0(scale_bar_um, " µm"),
    size = 3.2,
    fontface = "bold",
    vjust = 0
  ) +
  coord_fixed(
    xlim = c(cc$xmin, cc$xmax),
    ylim = c(cc$ymin, cc$ymax),
    expand = FALSE
  ) +
  theme_void(base_size = 14) +
  theme(plot.margin = margin(0, 0, 0, 0))

he_files <- save_plot_pair(
  p_he,
  he_dir,
  "Fig1_sample3D_HE"
)

cluster_coords <- coords0
cluster_coords$cluster <- factor(
  as.character(
    object@meta.data[cluster_coords$cell, "cluster"]
  ),
  levels = names(cluster_colors_red)
)

if (anyNA(cluster_coords$cluster)) {
  stop("Missing or unexpected cluster assignments for sample 3D")
}

n_cluster_spots <- nrow(cluster_coords)

cluster_circle_df <- data.frame(
  cell = rep(cluster_coords$cell, each = n_circle_points),
  cluster = factor(
    rep(as.character(cluster_coords$cluster), each = n_circle_points),
    levels = names(cluster_colors_red)
  ),
  x = rep(cluster_coords$x_hires, each = n_circle_points) +
    rep(cos(theta), times = n_cluster_spots) * spot_radius_hires,
  y = rep(cluster_coords$y_plot, each = n_circle_points) +
    rep(sin(theta), times = n_cluster_spots) * spot_radius_hires,
  stringsAsFactors = FALSE
)

cluster_legend_df <- data.frame(
  x = rep(cc$xmin, length(cluster_colors_red)),
  y = rep(cc$ymin, length(cluster_colors_red)),
  cluster = factor(
    names(cluster_colors_red),
    levels = names(cluster_colors_red)
  )
)

make_cluster_plot <- function(palette) {
  ggplot() +
    annotation_raster(
      he,
      xmin = 0,
      xmax = img_w,
      ymin = 0,
      ymax = img_h
    ) +
    geom_polygon(
      data = cluster_circle_df,
      aes(x = x, y = y, group = cell, fill = cluster),
      colour = NA,
      show.legend = FALSE
    ) +
    geom_point(
      data = cluster_legend_df,
      aes(x = x, y = y, fill = cluster),
      shape = 21,
      size = 4.5,
      colour = NA,
      alpha = 0,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = palette,
      breaks = names(palette),
      labels = cluster_labels_short,
      drop = FALSE,
      name = "Cluster"
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          size = 5,
          colour = NA,
          alpha = 1
        ),
        ncol = 4,
        byrow = TRUE
      )
    ) +
    coord_fixed(
      xlim = c(cc$xmin, cc$xmax),
      ylim = c(cc$ymin, cc$ymax),
      expand = FALSE
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "inside",
      legend.position.inside = c(0.76, 0.10),
      legend.justification = c(0.5, 0.5),
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(
        fill = "white",
        colour = "grey70",
        linewidth = 0.4
      ),
      legend.margin = margin(3, 5, 3, 5)
    )
}

p_clusters_red <- make_cluster_plot(cluster_colors_red)
p_clusters_yellow <- make_cluster_plot(cluster_colors_yellow)

cluster_red_files <- save_plot_pair(
  p_clusters_red,
  clusters_red_dir,
  "Fig1_sample3D_clusters_C5red"
)

cluster_yellow_files <- save_plot_pair(
  p_clusters_yellow,
  clusters_yellow_dir,
  "Fig1_sample3D_clusters_C5yellow"
)

# ============================================================
# EXPORT SETTINGS
# ============================================================

summary_list <- list()

plot_counter <- 0

# ============================================================
# LOOP:
# CLUSTER -> TOP 10 MARKERS
# ============================================================

for (cluster_use in names(markers_by_cluster)) {

  cat("\n\n####################################################\n")
  cat("CLUSTER: ", cluster_use, "\n", sep = "")
  cat("####################################################\n")

  cluster_dir <- file.path(
    out_dir,
    cluster_use
  )

  genes_use <-
    markers_by_cluster[[cluster_use]]

  for (gene_rank in seq_along(genes_use)) {

    gene_use <-
      genes_use[gene_rank]

    cat("\n====================================================\n")
    cat(
      cluster_use,
      " | rank ",
      gene_rank,
      " | ",
      gene_use,
      "\n",
      sep = ""
    )
    cat("====================================================\n")

    # ========================================================
    # Check gene
    # ========================================================

    if (!(gene_use %in% rownames(expr_mat))) {

      warning(
        cluster_use,
        " marker ",
        gene_use,
        " not found in SCT data. Skipping."
      )

      summary_list[[length(summary_list) + 1]] <-
        data.frame(
          cluster = cluster_use,
          rank = gene_rank,
          gene = gene_use,
          status = "MISSING_FROM_SCT",
          donor = donor_use,
          assay = "SCT",
          layer = "data",
          n_spots = NA,
          min_expression = NA,
          max_expression = NA,
          q99_expression = NA,
          spot_multiplier = 1.00,
          spot_diameter_hires_px = spot_diameter_hires,
          pdf_file = NA,
          png_file = NA,
          stringsAsFactors = FALSE
        )

      next
    }

    # ========================================================
    # Extract SCT normalized expression
    # ========================================================

    coords <- coords0

    coords$expression <-
      as.numeric(
        expr_mat[
          gene_use,
          coords$cell
        ]
      )

    expr_vals <-
      coords$expression[
        is.finite(coords$expression)
      ]

    if (length(expr_vals) == 0) {

      warning(
        gene_use,
        ": no finite expression values."
      )

      next
    }

    # ========================================================
    # EXPRESSION SCALE
    #
    # Color maximum = 99th percentile
    # Extreme values > q99 are squished to max color.
    # ========================================================

    min_expr <- min(
      expr_vals,
      na.rm = TRUE
    )

    max_expr <- max(
      expr_vals,
      na.rm = TRUE
    )

    q99_expr <- as.numeric(
      quantile(
        expr_vals,
        probs = 0.99,
        na.rm = TRUE
      )
    )

    if (
      !is.finite(q99_expr) ||
      q99_expr <= 0
    ) {
      q99_expr <- max_expr
    }

    if (
      !is.finite(q99_expr) ||
      q99_expr <= 0
    ) {

      warning(
        gene_use,
        ": no positive expression. Skipping."
      )

      next
    }

    mid_expr <-
      q99_expr / 2

    cat(
      "min = ",
      min_expr,
      "\n",
      sep = ""
    )

    cat(
      "max = ",
      max_expr,
      "\n",
      sep = ""
    )

    cat(
      "q99 = ",
      q99_expr,
      "\n",
      sep = ""
    )

    # ========================================================
    # TRUE VISIUM CIRCLES
    # ========================================================

    nspots <-
      nrow(coords)

    circle_df <- data.frame(

      cell = rep(
        coords$cell,
        each = n_circle_points
      ),

      expression = rep(
        coords$expression,
        each = n_circle_points
      ),

      x = rep(
        coords$x_hires,
        each = n_circle_points
      ) +
        rep(
          cos(theta),
          times = nspots
        ) * spot_radius_hires,

      y = rep(
        coords$y_plot,
        each = n_circle_points
      ) +
        rep(
          sin(theta),
          times = nspots
        ) * spot_radius_hires,

      stringsAsFactors = FALSE
    )

    # ========================================================
    # PLOT
    # ========================================================

    p <- ggplot() +

      annotation_raster(
        he,
        xmin = 0,
        xmax = img_w,
        ymin = 0,
        ymax = img_h
      ) +

      geom_polygon(
        data = circle_df,
        aes(
          x = x,
          y = y,
          group = cell,
          fill = expression
        ),
        colour = NA
      ) +

      scale_fill_gradientn(
        colours = hex_cols,

        limits = c(
          0,
          q99_expr
        ),

        oob = scales::squish,

        name = gene_use,

        breaks = c(
          0,
          mid_expr,
          q99_expr
        ),

        labels = c(
          "0",

          format(
            round(
              mid_expr,
              2
            ),
            trim = TRUE
          ),

          format(
            round(
              q99_expr,
              2
            ),
            trim = TRUE
          )
        ),

        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,

          barwidth = unit(
            29,
            "mm"
          ),

          barheight = unit(
            4.2,
            "mm"
          )
        )
      ) +

      coord_fixed(
        xlim = c(cc$xmin, cc$xmax),
        ylim = c(cc$ymin, cc$ymax),
        expand = FALSE
      ) +

      theme_void(
        base_size = 14
      ) +

      theme(
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "inside",
        legend.position.inside = c(0.74, 0.12),
        legend.justification = c(0.5, 0.5),
        legend.direction = "horizontal",

        legend.title = element_text(
          size = 11,
          face = "bold",
          hjust = 0.5
        ),

        legend.text = element_text(
          size = 8.5
        ),

        legend.key = element_rect(
          fill = "white",
          colour = NA
        ),

        legend.background = element_blank(),

        legend.box.background = element_rect(
          fill = "white",
          colour = "grey70",
          linewidth = 0.4
        ),

        legend.margin = margin(2, 0, 2, 0),
        legend.box.margin = margin(4, 16, 4, 16)
      )

    # ========================================================
    # FILE NAMES
    #
    # Prefix rank so markers remain in Top10 order.
    # Example:
    # 01_AQP1_TRUE_Visium_SCT.pdf
    # ========================================================

    rank_string <-
      sprintf(
        "%02d",
        gene_rank
      )

    pdf_file <- file.path(
      cluster_dir,
      paste0(
        rank_string,
        "_",
        gene_use,
        "_TRUE_Visium_SCT.pdf"
      )
    )

    png_file <- file.path(
      cluster_dir,
      paste0(
        rank_string,
        "_",
        gene_use,
        "_TRUE_Visium_SCT.png"
      )
    )

    ggsave(
      filename = pdf_file,
      plot = p,
      width = panel_width_in,
      height = panel_height_in,
      device = cairo_pdf,
      bg = "white"
    )

    ggsave(
      filename = png_file,
      plot = p,
      width = panel_width_in,
      height = panel_height_in,
      dpi = dpi_out,
      device = "png",
      bg = "white"
    )

    plot_counter <-
      plot_counter + 1

    # ========================================================
    # MANIFEST RECORD
    # ========================================================

    summary_list[[length(summary_list) + 1]] <-
      data.frame(
        cluster = cluster_use,
        rank = gene_rank,
        gene = gene_use,
        status = "PLOTTED",
        donor = donor_use,
        assay = "SCT",
        layer = "data",
        n_spots = nrow(coords),
        min_expression = min_expr,
        max_expression = max_expr,
        q99_expression = q99_expr,
        spot_multiplier = 1.00,
        spot_diameter_hires_px = spot_diameter_hires,
        pdf_file = pdf_file,
        png_file = png_file,
        stringsAsFactors = FALSE
      )

    cat("Saved:\n")
    cat(
      "  ",
      pdf_file,
      "\n",
      sep = ""
    )

    cat(
      "  ",
      png_file,
      "\n",
      sep = ""
    )
  }
}

# ============================================================
# COMBINED MANIFEST
# ============================================================

if (length(summary_list) > 0) {

  summary_df <-
    do.call(
      rbind,
      summary_list
    )

} else {

  summary_df <-
    data.frame()
}

write.csv(
  summary_df,
  file.path(
    out_dir,
    "3D_Top10_markers_by_cluster_manifest.csv"
  ),
  row.names = FALSE
)

if (plot_counter != 70) {
  stop(
    "Expected 70 marker plots, but generated ",
    plot_counter
  )
}

object_sha256_after <- digest::digest(
  file = object_path,
  algo = "sha256"
)

if (!identical(object_sha256_before, object_sha256_after)) {
  stop("Source object checksum changed during plotting")
}

# ============================================================
# ONE LIGHTWEIGHT RDS FOR THE WHOLE RUN
# ============================================================

saveRDS(
  list(
    source_object = object_path,
    source_object_sha256 = object_sha256_before,
    donor = donor_use,
    image_key = image_key,
    crop_coordinates = crop_coordinates,
    assay = "SCT",
    layer = "data",
    markers_by_cluster = markers_by_cluster,
    manifest = summary_df,
    spot_geometry = "TRUE Visium / Space Ranger",
    spot_multiplier = 1.00,
    spot_diameter_fullres = spot_diameter_fullres,
    hires_scale_factor = sf_hires,
    spot_diameter_hires_px = spot_diameter_hires,
    expression_color_max = "gene-specific 99th percentile",
    marker_legend = paste(
      "inside bottom-right; opaque white outer legend box;",
      "29-mm bar; 16-point horizontal box margins"
    ),
    cluster_palette_C5red = cluster_colors_red,
    cluster_palette_C5yellow = cluster_colors_yellow,
    H_and_E_scale_bar_um = scale_bar_um,
    H_and_E_source = fullres_he_png,
    H_and_E_source_sha256 = fullres_he_png_sha256,
    H_and_E_registered_source_copy = he_source_copy,
    H_and_E_registered_source_copy_sha256 = he_source_copy_sha256,
    H_and_E_final_manuscript_crop = "manual Adobe Photoshop crop from registered Space Ranger source",
    H_and_E_source_width_px = fullres_img_w,
    H_and_E_source_height_px = fullres_img_h,
    H_and_E_hires_to_fullres_x = hires_to_fullres_x,
    H_and_E_hires_to_fullres_y = hires_to_fullres_y,
    H_and_E_fullres_crop_columns = c(
      fullres_col_start,
      fullres_col_end
    ),
    H_and_E_fullres_crop_rows = c(
      fullres_row_start,
      fullres_row_end
    )
  ),
  file.path(
    out_dir,
    "3D_Top10_markers_by_cluster_plotting_record.rds"
  )
)

# ============================================================
# REPRODUCIBILITY LOG
# ============================================================

writeLines(
  c(
    paste0(
      "source_object: ",
      object_path
    ),

    paste0(
      "source_object_sha256_before: ",
      object_sha256_before
    ),

    paste0(
      "source_object_sha256_after: ",
      object_sha256_after
    ),

    paste0("generator_script: ", script_path),
    paste0("generator_script_sha256: ", script_sha256),

    paste0(
      "donor: ",
      donor_use
    ),

    paste0(
      "image_key: ",
      image_key
    ),

    paste0(
      "hires_image: ",
      hires_png
    ),

    paste0(
      "hires_image_sha256: ",
      hires_png_sha256
    ),

    paste0(
      "H_and_E_registered_source_copy: ",
      he_source_copy
    ),

    paste0(
      "H_and_E_registered_source_copy_sha256: ",
      he_source_copy_sha256
    ),

    "H_and_E_final_manuscript_crop: manual Adobe Photoshop crop from registered Space Ranger source",

    paste0(
      "H_and_E_source_image: ",
      fullres_he_png
    ),

    paste0(
      "H_and_E_source_image_sha256: ",
      fullres_he_png_sha256
    ),

    paste0(
      "H_and_E_source_width_px: ",
      fullres_img_w
    ),

    paste0(
      "H_and_E_source_height_px: ",
      fullres_img_h
    ),

    paste0(
      "H_and_E_coordinate_scale_x: ",
      format(hires_to_fullres_x, digits = 15)
    ),

    paste0(
      "H_and_E_coordinate_scale_y: ",
      format(hires_to_fullres_y, digits = 15)
    ),

    paste0(
      "H_and_E_source_crop_columns: ",
      fullres_col_start,
      ",",
      fullres_col_end
    ),

    paste0(
      "H_and_E_source_crop_rows: ",
      fullres_row_start,
      ",",
      fullres_row_end
    ),

    paste0("output_root: ", output_root),
    paste0("crop_id: ", crop_id),
    paste0("xmin: ", cc$xmin),
    paste0("xmax: ", cc$xmax),
    paste0("ymin: ", cc$ymin),
    paste0("ymax: ", cc$ymax),

    "assay: SCT",
    "layer: data",

    "top_markers_source: Top10_markers_per_final_cluster_dotplot_data.csv",

    paste0(
      "n_cluster_marker_entries: ",
      nrow(marker_table)
    ),

    paste0(
      "n_cluster_marker_entries_found_in_SCT: ",
      sum(marker_table$found_in_SCT)
    ),

    paste0(
      "n_plots_generated: ",
      plot_counter
    ),

    "spot_geometry: TRUE Visium / Space Ranger",

    "spot_multiplier: 1.00",

    paste0(
      "spot_diameter_fullres: ",
      spot_diameter_fullres
    ),

    paste0(
      "hires_scale_factor: ",
      sf_hires
    ),

    paste0(
      "spot_diameter_hires_px: ",
      spot_diameter_hires
    ),

    "expression_source: SCT normalized data layer",

    "expression_color_scale: gene-specific",

    "expression_color_max: 99th percentile",

    "values_above_q99: squished to maximum color",
    "marker_legend_position: inside bottom-right",
    "marker_legend_background: opaque white outer legend box",
    "marker_legend_barwidth_mm: 29",
    "marker_legend_barheight_mm: 4.2",
    "marker_legend_box_horizontal_margin_pt: 16",
    "cluster_legend_position: inside bottom-right",
    paste0("C5_red: ", cluster_colors_red[["C5"]]),
    paste0("C5_yellow: ", cluster_colors_yellow[["C5"]]),
    paste0("H_and_E_scale_bar_um: ", scale_bar_um),
    paste0("panel_width_in: ", panel_width_in),
    paste0("panel_height_in: ", panel_height_in),
    "source_object_modified_or_saved: NO"
  ),
  file.path(
    log_dir,
    "run_parameters.txt"
  )
)

# ============================================================
# SESSION INFO
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    log_dir,
    "sessionInfo.txt"
  )
)

# ============================================================
# OUTPUT VALIDATION AND SHA256 MANIFESTS
# ============================================================

marker_pdfs <- list.files(
  markers_dir,
  pattern = "\\.pdf$",
  recursive = TRUE,
  full.names = TRUE
)

marker_pngs <- list.files(
  markers_dir,
  pattern = "\\.png$",
  recursive = TRUE,
  full.names = TRUE
)

validation <- data.frame(
  check = c(
    "crop_coordinates",
    "marker_pdf_count",
    "marker_png_count",
    "H_and_E_pdf_png",
    "clusters_C5red_pdf_png",
    "clusters_C5yellow_pdf_png",
    "C5_red_hex",
    "C5_yellow_hex",
    "object_checksum_unchanged"
  ),
  observed = c(
    paste(
      cc$xmin,
      cc$xmax,
      cc$ymin,
      cc$ymax,
      sep = ","
    ),
    length(marker_pdfs),
    length(marker_pngs),
    sum(file.exists(he_files)),
    sum(file.exists(cluster_red_files)),
    sum(file.exists(cluster_yellow_files)),
    cluster_colors_red[["C5"]],
    cluster_colors_yellow[["C5"]],
    identical(object_sha256_before, object_sha256_after)
  ),
  expected = c(
    "520,1040,1020,1440",
    "70",
    "70",
    "2",
    "2",
    "2",
    "#E15759",
    "#EDC948",
    "TRUE"
  ),
  stringsAsFactors = FALSE
)

validation$pass <- validation$observed == validation$expected

write.csv(
  validation,
  file.path(records_dir, "generation_validation.csv"),
  row.names = FALSE
)

if (!all(validation$pass)) {
  print(validation)
  stop("Final cropped Figure 1 validation failed")
}

relative_to_output <- function(path) {
  prefix <- paste0(normalizePath(output_root), .Platform$file.sep)
  sub(prefix, "", normalizePath(path), fixed = TRUE)
}

manifest_path <- file.path(output_root, "output_manifest.tsv")
sha_path <- file.path(output_root, "SHA256SUMS.txt")

manifest_files <- list.files(
  output_root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

manifest_files <- manifest_files[
  !file.info(manifest_files)$isdir
]

manifest_files <- setdiff(
  manifest_files,
  c(manifest_path, sha_path)
)

manifest_files <- sort(manifest_files)

manifest_df <- data.frame(
  relative_path = vapply(
    manifest_files,
    relative_to_output,
    character(1)
  ),
  bytes = unname(file.info(manifest_files)$size),
  sha256 = vapply(
    manifest_files,
    function(path) {
      digest::digest(file = path, algo = "sha256")
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.table(
  manifest_df,
  manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

checksum_files <- sort(c(manifest_files, manifest_path))

checksum_lines <- vapply(
  checksum_files,
  function(path) {
    paste0(
      digest::digest(file = path, algo = "sha256"),
      "  ",
      relative_to_output(path)
    )
  },
  character(1)
)

writeLines(checksum_lines, sha_path)

# ============================================================
# DONE
# ============================================================

cat("\n====================================================\n")
cat("FIG1_CROPPED_F3_SAMPLE3D_OK\n")
cat("====================================================\n")

cat(
  "\nCluster-marker entries requested: ",
  nrow(marker_table),
  "\n",
  sep = ""
)

cat(
  "Plots successfully generated: ",
  plot_counter,
  "\n",
  sep = ""
)

if (nrow(missing_table) > 0) {

  cat(
    "Missing cluster-marker entries: ",
    nrow(missing_table),
    "\n",
    sep = ""
  )

} else {

  cat(
    "Missing cluster-marker entries: 0\n"
  )
}

cat("H&E panel pairs generated: 1\n")
cat("C5-red cluster panel pairs generated: 1\n")
cat("C5-yellow cluster panel pairs generated: 1\n")
cat("Source object modified/saved: NO\n")
cat("\nResults:\n")
cat(output_root, "\n")
