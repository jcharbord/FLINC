#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(png)
  library(scales)
  library(grid)
})

repo_dir <- "/Users/jianwu/Documents/FLINC"

# ============================================================
# SOURCE REPRODUCIBLE OBJECT
# ============================================================

object_path <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/results/260701_approach1_cluster_go_marker_deconv",
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)

donor_use <- "3D"
image_key <- "X3D"

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
  "/Users/jianwu/Documents/DATA/spaceranger_output",
  donor_use,
  "spatial",
  "tissue_hires_image.png"
)

cmap_file <- file.path(
  repo_dir,
  "data/st_analysis_pipeline/color_pseudotime.csv"
)

# ============================================================
# OUTPUT FOLDERS
# ============================================================

out_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/results/3D_Top10_markers_by_cluster_TRUE_Visium_SCT"
)

log_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/reproducibility_logs/3D_Top10_markers_by_cluster_TRUE_Visium_SCT"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Make C0-C6 subfolders
for (cluster_use in names(markers_by_cluster)) {

  dir.create(
    file.path(out_dir, cluster_use),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# LOAD OBJECT
# ============================================================

cat("Loading final annotated Seurat object...\n")
cat(object_path, "\n")

object <- readRDS(object_path)

stopifnot("SCT" %in% Assays(object))
stopifnot(image_key %in% names(object@images))
stopifnot(file.exists(hires_png))
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
# LOAD REGISTERED HIRES H&E
# ============================================================

he <- readPNG(hires_png)

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
# EXPORT SETTINGS
# ============================================================

dpi_out <- 300

width_in <-
  img_w / dpi_out

height_in <-
  img_h / dpi_out

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
            35,
            "mm"
          ),

          barheight = unit(
            5,
            "mm"
          )
        )
      ) +

      coord_fixed(
        xlim = c(
          0,
          img_w
        ),
        ylim = c(
          0,
          img_h
        ),
        expand = FALSE
      ) +

      theme_void(
        base_size = 14
      ) +

      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",

        legend.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5
        ),

        legend.text = element_text(
          size = 11
        )
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

    # extra height for horizontal colorbar
    save_height <-
      height_in + 0.7

    ggsave(
      filename = pdf_file,
      plot = p,
      width = width_in,
      height = save_height,
      device = cairo_pdf,
      bg = "white"
    )

    ggsave(
      filename = png_file,
      plot = p,
      width = width_in,
      height = save_height,
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

# ============================================================
# ONE LIGHTWEIGHT RDS FOR THE WHOLE RUN
# ============================================================

saveRDS(
  list(
    source_object = object_path,
    donor = donor_use,
    image_key = image_key,
    assay = "SCT",
    layer = "data",
    markers_by_cluster = markers_by_cluster,
    manifest = summary_df,
    spot_geometry = "TRUE Visium / Space Ranger",
    spot_multiplier = 1.00,
    spot_diameter_fullres = spot_diameter_fullres,
    hires_scale_factor = sf_hires,
    spot_diameter_hires_px = spot_diameter_hires,
    expression_color_max = "gene-specific 99th percentile"
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

    "values_above_q99: squished to maximum color"
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
# DONE
# ============================================================

cat("\n====================================================\n")
cat("DONE: 3D TOP 10 MARKERS PER FINAL CLUSTER\n")
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

cat("\nMain results folder:\n")
cat(out_dir, "\n")

cat("\nSubfolders:\n")

for (cl in names(markers_by_cluster)) {
  cat("  ", file.path(out_dir, cl), "\n", sep = "")
}

cat("\nLogs folder:\n")
cat(log_dir, "\n")
