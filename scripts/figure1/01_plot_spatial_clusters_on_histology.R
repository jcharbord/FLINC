#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(png)
})

repo_dir <- "/Users/jianwu/Documents/FLINC"

# ============================================================
# ORIGINAL reproducible final Seurat object
# We DO NOT modify or overwrite this object
# ============================================================

object_path <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/results/260701_approach1_cluster_go_marker_deconv",
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)

image_root <- "/Users/jianwu/Documents/DATA/spaceranger_output"

out_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/results/all_samples_TRUE_Visium_clusters_C5red_circlelegend"
)

log_dir <- file.path(
  repo_dir,
  "Jeremie/Reproducible_code/reproducibility_logs/all_samples_TRUE_Visium_clusters_C5red_circlelegend"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading final annotated Seurat object...\n")
cat(object_path, "\n")

object <- readRDS(object_path)

# ============================================================
# FINAL FIGURE PALETTE
# C5 changed to red
# ============================================================

cluster_colors <- c(
  C0 = "#56B4E9",  # Portal tracts
  C1 = "#FF9DA7",  # Periportal hepatocytes
  C2 = "#B189A8",  # Centrilobular hepatocytes
  C3 = "#59A14F",  # Immune cells
  C4 = "#F28E2B",  # Myofibroblasts
  C5 = "#E15759",  # Cholangiocytes - RED
  C6 = "#E5C494"   # Matrix-rich fibroblast
)

# Short legend labels only
cluster_labels <- c(
  C0 = "C0",
  C1 = "C1",
  C2 = "C2",
  C3 = "C3",
  C4 = "C4",
  C5 = "C5",
  C6 = "C6"
)

# ============================================================
# Plot settings
# ============================================================

dpi_out <- 300

# polygons approximate each circular Visium spot
n_circle_points <- 48

image_keys <- names(object@images)
donors <- sub("^X", "", image_keys)

cat("\nSamples to process:\n")
print(donors)

manifest_list <- list()
failed_samples <- character(0)

# ============================================================
# LOOP OVER ALL SAMPLES
# ============================================================

for (i in seq_along(image_keys)) {

  image_key <- image_keys[i]
  donor_use <- donors[i]

  cat("\n====================================================\n")
  cat("Processing donor: ", donor_use, "\n", sep = "")
  cat("Image key: ", image_key, "\n", sep = "")
  cat("====================================================\n")

  hires_png <- file.path(
    image_root,
    donor_use,
    "spatial",
    "tissue_hires_image.png"
  )

  if (!file.exists(hires_png)) {
    warning("Missing hires image for donor: ", donor_use)
    failed_samples <- c(failed_samples, donor_use)
    next
  }

  img <- object@images[[image_key]]

  # ==========================================================
  # Spot coordinates + final cluster identity
  # ==========================================================

  coords <- img@coordinates
  coords$cell <- rownames(coords)

  shared_cells <- intersect(
    coords$cell,
    rownames(object@meta.data)
  )

  if (length(shared_cells) == 0) {
    warning("No shared cells for donor: ", donor_use)
    failed_samples <- c(failed_samples, donor_use)
    next
  }

  coords <- coords[shared_cells, , drop = FALSE]

  coords$cluster <- as.character(
    object@meta.data[coords$cell, "cluster"]
  )

  coords$cluster <- factor(
    coords$cluster,
    levels = names(cluster_colors)
  )

  # ==========================================================
  # Load registered Space Ranger HIRES H&E
  # ==========================================================

  he <- readPNG(hires_png)

  img_h <- dim(he)[1]
  img_w <- dim(he)[2]

  sf_hires <- img@scale.factors$hires

  # ==========================================================
  # TRUE VISIUM SPOT GEOMETRY
  #
  # Space Ranger gives spot diameter in full-resolution pixels.
  # Convert to tissue_hires_image.png pixels.
  #
  # multiplier = 1.00 exactly
  # ==========================================================

  spot_diameter_fullres <- img@scale.factors$spot

  spot_diameter_hires <-
    spot_diameter_fullres * sf_hires

  spot_radius_hires <-
    spot_diameter_hires / 2

  cat(
    "True Visium spot diameter on hires image: ",
    spot_diameter_hires,
    " px\n",
    sep = ""
  )

  # ==========================================================
  # Convert coordinates to hires image
  # ==========================================================

  coords$x_hires <- coords$imagecol * sf_hires
  coords$y_hires <- coords$imagerow * sf_hires

  # ggplot origin is bottom-left; image origin is top-left
  coords$y_plot <- img_h - coords$y_hires

  # ==========================================================
  # Build true-size circular polygons
  # ==========================================================

  theta <- seq(
    0,
    2 * pi,
    length.out = n_circle_points + 1
  )[-(n_circle_points + 1)]

  nspots <- nrow(coords)

  circle_df <- data.frame(
    cell = rep(
      coords$cell,
      each = n_circle_points
    ),

    cluster = factor(
      rep(
        as.character(coords$cluster),
        each = n_circle_points
      ),
      levels = names(cluster_colors)
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

  # ==========================================================
  # Dummy points ONLY for circular legend symbols
  # They are placed far outside visible plotting coordinates.
  # ==========================================================

  legend_df <- data.frame(
    x = rep(-1e6, length(cluster_colors)),
    y = rep(-1e6, length(cluster_colors)),
    cluster = factor(
      names(cluster_colors),
      levels = names(cluster_colors)
    )
  )

  # ==========================================================
  # Plot
  # ==========================================================

  p <- ggplot() +

    annotation_raster(
      he,
      xmin = 0,
      xmax = img_w,
      ymin = 0,
      ymax = img_h
    ) +

    # Actual Visium spots
    geom_polygon(
      data = circle_df,
      aes(
        x = x,
        y = y,
        group = cell,
        fill = cluster
      ),
      colour = NA,
      show.legend = FALSE
    ) +

    # Circular legend symbols
    geom_point(
      data = legend_df,
      aes(
        x = x,
        y = y,
        fill = cluster
      ),
      shape = 21,
      size = 5,
      stroke = 0,
      colour = NA,
      show.legend = TRUE
    ) +

    scale_fill_manual(
      values = cluster_colors,
      breaks = names(cluster_colors),
      labels = cluster_labels,
      drop = FALSE,
      name = "Cluster"
    ) +

    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          size = 6,
          stroke = 0,
          colour = NA
        ),
        ncol = 1
      )
    ) +

    coord_fixed(
      xlim = c(0, img_w),
      ylim = c(0, img_h),
      expand = FALSE
    ) +

    theme_void(base_size = 14) +

    theme(
      legend.position = "right",
      legend.title = element_text(
        size = 12,
        face = "bold"
      ),
      legend.text = element_text(
        size = 11
      ),
      legend.key = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.spacing.y = unit(2, "pt")
    )

  # ==========================================================
  # Save PDF + PNG
  # ==========================================================

  width_in <- img_w / dpi_out
  height_in <- img_h / dpi_out

  pdf_file <- file.path(
    out_dir,
    paste0(
      donor_use,
      "_TRUE_Visium_clusters_C5red.pdf"
    )
  )

  png_file <- file.path(
    out_dir,
    paste0(
      donor_use,
      "_TRUE_Visium_clusters_C5red.png"
    )
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    width = width_in,
    height = height_in,
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = width_in,
    height = height_in,
    dpi = dpi_out,
    device = "png",
    bg = "white"
  )

  # ==========================================================
  # Record run information for this donor
  # ==========================================================

  manifest_list[[length(manifest_list) + 1]] <- data.frame(
    donor = donor_use,
    image_key = image_key,
    hires_png = hires_png,
    pdf_file = pdf_file,
    png_file = png_file,
    n_spots = nrow(coords),
    image_width_px = img_w,
    image_height_px = img_h,
    hires_scale_factor = sf_hires,
    spot_diameter_fullres = spot_diameter_fullres,
    spot_diameter_hires = spot_diameter_hires,
    spot_multiplier = 1.00,
    stringsAsFactors = FALSE
  )

  cat("Saved:\n")
  cat("  ", pdf_file, "\n", sep = "")
  cat("  ", png_file, "\n", sep = "")
}

# ============================================================
# Combined manifest
# ============================================================

if (length(manifest_list) > 0) {
  manifest_df <- do.call(rbind, manifest_list)
} else {
  manifest_df <- data.frame()
}

write.csv(
  manifest_df,
  file.path(
    out_dir,
    "all_samples_TRUE_Visium_manifest.csv"
  ),
  row.names = FALSE
)

# ONE lightweight RDS for the entire plotting run
saveRDS(
  list(
    source_object = object_path,
    manifest = manifest_df,
    palette = cluster_colors,
    cluster_labels = cluster_labels,
    spot_multiplier = 1.00,
    geometry = "Space Ranger spot_diameter_fullres x hires scale factor",
    dpi = dpi_out
  ),
  file.path(
    out_dir,
    "all_samples_TRUE_Visium_plotting_record.rds"
  )
)

# ============================================================
# Save palette
# ============================================================

write.csv(
  data.frame(
    cluster = names(cluster_colors),
    color = unname(cluster_colors),
    stringsAsFactors = FALSE
  ),
  file.path(
    log_dir,
    "cluster_palette_C5red.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Run parameters
# ============================================================

writeLines(
  c(
    paste0("source_object: ", object_path),
    paste0("image_root: ", image_root),
    paste0("output_dir: ", out_dir),
    paste0("log_dir: ", log_dir),
    paste0("dpi: ", dpi_out),
    paste0("spot_multiplier: 1.00"),
    paste0("geometry: TRUE Space Ranger / Visium spot diameter"),
    paste0("n_circle_points: ", n_circle_points),
    paste0("n_samples_processed: ", nrow(manifest_df)),
    paste0(
      "failed_samples: ",
      if (length(failed_samples) == 0) {
        "none"
      } else {
        paste(failed_samples, collapse = ", ")
      }
    ),
    "",
    "Final plotting palette:",
    paste(
      names(cluster_colors),
      cluster_colors,
      sep = " = "
    )
  ),
  file.path(
    log_dir,
    "run_parameters.txt"
  )
)

# ============================================================
# Session info
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(
    log_dir,
    "sessionInfo.txt"
  )
)

cat("\n====================================================\n")
cat("DONE: ALL TRUE VISIUM HIGH-RES CLUSTER OVERLAYS\n")
cat("====================================================\n")

cat("\nSamples successfully processed:\n")
cat(nrow(manifest_df), "\n")

cat("\nResults folder:\n")
cat(out_dir, "\n")

cat("\nLogs folder:\n")
cat(log_dir, "\n")

cat("\nOriginal Seurat object was NOT modified.\n")
