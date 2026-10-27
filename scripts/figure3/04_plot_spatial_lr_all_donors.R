#!/usr/bin/env Rscript

# Figure 3C: generate cell-type-weighted local ligand-receptor maps for all
# 24 FLINC donors. Spots are drawn as polygons at the true Space Ranger Visium
# diameter relative to each custom H&E image. Only retained LRlocal PDFs are
# produced; source-detection, target-detection, co-detection, and large derived
# RDS objects are excluded.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(RANN)
  library(ggplot2)
  library(scales)
  library(magick)
})

options(width = 180)

started_at <- Sys.time()

repo_root <- "/Users/jianwu/Documents/FLINC"
repro_root <- file.path(repo_root, "Jeremie", "Reproducible_code")

input_rds <- file.path(
  repro_root,
  "results",
  "260701_approach1_cluster_go_marker_deconv",
  "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
)
selected_rds <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "02_main_and_supplementary_figures",
  "02_selected_interactions_top60.rds"
)
step02_checksums_tsv <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "02_main_and_supplementary_figures",
  "02_OUTPUT_SHA256.tsv"
)
step02_completion_tsv <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "02_main_and_supplementary_figures",
  "02_COMPLETION.tsv"
)
he_image_dir <- "/Users/jianwu/Documents/DATA/H&E"
output_dir <- file.path(
  repro_root,
  "results",
  "Fig3_CCC",
  "03_spatial_LR_all_samples"
)

expected_input_sha256 <-
  "478cbd318aa832381fd9a2cb1f0ecc7e21f7ce55c62b4a5d0de2d522df517fdd"

assay_use <- "Spatial"
input_layer <- "counts"
deconvolution_assay <- "singlecell"
deconvolution_layer <- "counts"
sample_column <- "donor"
condition_column <- "F_score"

cell_types <- c(
  "Bcells", "CentralHep", "Cholangiocyte", "Erythroid",
  "InfMac", "InterHep", "NKTcell", "NonInfMac",
  "PortalEndo", "PortalHep", "Stellate", "cvLSECs"
)

story_edges <- tibble::tribble(
  ~edge_order, ~source, ~target, ~ligand, ~receptor, ~note,
  1L, "Cholangiocyte", "Stellate", "NPNT", "ITGA8", "portal_fibrotic_niche",
  2L, "Cholangiocyte", "InfMac", "MMP7", "CD151", "remodeling_broadcast",
  3L, "Stellate", "PortalEndo", "CCL2", "ACKR1", "trafficking"
) %>%
  dplyr::mutate(
    pair = paste(ligand, receptor, sep = "_"),
    direction = paste(source, "to", target, sep = "_"),
    score_column = paste0(
      "LR__", source, "_to_", target, "__", ligand, "_", receptor, "_log"
    )
  )

donor_stage_map <- tibble::tribble(
  ~stage_order, ~stage, ~stage_folder, ~donor_order, ~donor,
  0L, "F0", "Healthy_F0", 1L, "6B",
  0L, "F0", "Healthy_F0", 2L, "5D",
  0L, "F0", "Healthy_F0", 3L, "4C",
  0L, "F0", "Healthy_F0", 4L, "3B",
  0L, "F0", "Healthy_F0", 5L, "2A",
  0L, "F0", "Healthy_F0", 6L, "1A",
  1L, "F1", "F1", 1L, "6A",
  1L, "F1", "F1", 2L, "5A",
  1L, "F1", "F1", 3L, "4D",
  1L, "F1", "F1", 4L, "3C",
  1L, "F1", "F1", 5L, "2B",
  1L, "F1", "F1", 6L, "1D",
  2L, "F2", "F2", 1L, "6C",
  2L, "F2", "F2", 2L, "5B",
  2L, "F2", "F2", 3L, "1C",
  3L, "F3", "F3", 1L, "4A",
  3L, "F3", "F3", 2L, "3D",
  3L, "F3", "F3", 3L, "2C",
  4L, "F4", "F4", 1L, "6D",
  4L, "F4", "F4", 2L, "5C",
  4L, "F4", "F4", 3L, "4B",
  4L, "F4", "F4", 4L, "3A",
  4L, "F4", "F4", 5L, "2D",
  4L, "F4", "F4", 6L, "1B"
)

pseudocount_multiplier <- 1e4
spatial_normalization_scale_factor <- 1e4
nearest_neighbors <- 30L
include_self <- TRUE
effective_neighbors <- nearest_neighbors + as.integer(include_self)
gaussian_sigma <- 2
spot_multiplier <- 1.00
n_circle_points <- 48L
space_ranger_hires_max_dimension_pixels <- 2000
image_height_pixels <- 1900
pdf_width_inches <- 8
color_quantiles <- c(0.01, 0.99)

expected_spots <- 26056L
expected_input_genes <- 17922L
expected_retained_genes <- 17775L
expected_zero_genes_dropped <- 147L
expected_donors <- 24L
expected_edges <- 3L
expected_pdfs <- expected_donors * expected_edges

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tibble", "RANN",
  "ggplot2", "scales", "magick", "digest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

if (nrow(donor_stage_map) != expected_donors ||
    anyDuplicated(donor_stage_map$donor)) {
  stop("Declared donor-stage map is not a unique 24-donor map.")
}

if (nrow(story_edges) != expected_edges) {
  stop("Expected exactly three retained spatial LR edges.")
}

he_image_paths <- setNames(
  file.path(he_image_dir, paste0(donor_stage_map$donor, ".png")),
  donor_stage_map$donor
)

required_inputs <- c(
  input_rds,
  selected_rds,
  step02_checksums_tsv,
  step02_completion_tsv,
  unname(he_image_paths)
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required spatial-map inputs: ", paste(missing_inputs, collapse = ", "))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", as.character(x))
}

plot_plan <- do.call(
  rbind,
  lapply(seq_len(nrow(donor_stage_map)), function(donor_index) {
    donor_row <- donor_stage_map[donor_index, , drop = FALSE]
    do.call(
      rbind,
      lapply(seq_len(nrow(story_edges)), function(edge_index) {
        edge_row <- story_edges[edge_index, , drop = FALSE]
        filename <- paste0(
          donor_row$donor, "_",
          safe_name(edge_row$source), "_to_", safe_name(edge_row$target), "__",
          safe_name(edge_row$ligand), "_", safe_name(edge_row$receptor), "__",
          safe_name(edge_row$note), "__C_LRlocal.pdf"
        )
        relative_path <- file.path(
          donor_row$stage_folder,
          edge_row$pair,
          edge_row$direction,
          filename
        )
        data.frame(
          stage_order = donor_row$stage_order,
          stage = donor_row$stage,
          stage_folder = donor_row$stage_folder,
          donor_order = donor_row$donor_order,
          donor = donor_row$donor,
          edge_order = edge_row$edge_order,
          source = edge_row$source,
          target = edge_row$target,
          ligand = edge_row$ligand,
          receptor = edge_row$receptor,
          note = edge_row$note,
          pair = edge_row$pair,
          direction = edge_row$direction,
          score_column = edge_row$score_column,
          relative_pdf_path = relative_path,
          stringsAsFactors = FALSE
        )
      })
    )
  })
) %>%
  dplyr::arrange(stage_order, donor_order, edge_order)

if (nrow(plot_plan) != expected_pdfs ||
    anyDuplicated(plot_plan$relative_pdf_path)) {
  stop("Spatial LR plot plan is not a unique 72-PDF plan.")
}

plot_plan$absolute_pdf_path <- file.path(
  output_dir,
  plot_plan$relative_pdf_path
)

metadata_paths <- c(
  story_edges_csv = file.path(output_dir, "04_story_interactions.csv"),
  donor_stage_csv = file.path(output_dir, "04_donor_stage_map.csv"),
  color_limits_csv = file.path(output_dir, "04_spatial_lr_color_limits.csv"),
  plot_manifest_csv = file.path(output_dir, "04_spatial_lr_manifest.csv"),
  parameters_tsv = file.path(output_dir, "04_analysis_parameters.tsv"),
  input_manifest_tsv = file.path(output_dir, "04_input_manifest.tsv"),
  session_info_txt = file.path(output_dir, "04_sessionInfo.txt")
)

planned_outputs <- c(plot_plan$absolute_pdf_path, unname(metadata_paths))
blocked_outputs <- planned_outputs[
  file.exists(planned_outputs) | file.exists(paste0(planned_outputs, ".partial"))
]
if (length(blocked_outputs)) {
  stop(
    "Refusing to replace existing Step 04 outputs: ",
    paste(blocked_outputs, collapse = ", ")
  )
}

atomic_write_csv <- function(object, path, row.names = FALSE) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  utils::write.csv(object, partial, row.names = row.names)
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
}

atomic_write_tsv <- function(object, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  utils::write.table(
    object,
    partial,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
}

atomic_write_lines <- function(lines, path) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing output: ", path)
  }
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) {
    stop("Could not finalize output: ", path)
  }
}

atomic_ggsave_pdf <- function(path, plot, width, height) {
  partial <- paste0(path, ".partial")
  if (file.exists(path) || file.exists(partial)) {
    stop("Refusing to replace existing plot: ", path)
  }
  if (!dir.exists(dirname(path))) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  ggplot2::ggsave(
    filename = partial,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  if (!file.exists(partial) || file.info(partial)$size <= 0) {
    stop("ggsave did not create a non-empty spatial PDF: ", path)
  }
  if (!file.rename(partial, path)) {
    stop("Could not finalize spatial PDF: ", path)
  }
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

relative_to_output <- function(path) {
  substring(path, nchar(output_dir) + 2L)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
} else {
  NA_character_
}

message("Figure 3 spatial LR mapping started at ", format(started_at))
message("Output: ", output_dir)
message("Verifying Step 02 completion and selected-interaction checksum...")

step02_completion <- utils::read.delim(
  step02_completion_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
step02_status <- step02_completion$value[step02_completion$field == "status"]
if (length(step02_status) != 1L || step02_status != "complete") {
  stop("Step 02 completion status is not 'complete'.")
}

step02_checksums <- utils::read.delim(
  step02_checksums_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
selected_expected_sha256 <- step02_checksums$sha256[
  step02_checksums$file == basename(selected_rds)
]
if (length(selected_expected_sha256) != 1L ||
    !identical(sha256_file(selected_rds), selected_expected_sha256)) {
  stop("Selected-interaction RDS does not match the Step 02 checksum manifest.")
}

selected_interactions <- readRDS(selected_rds)
matched_edges <- story_edges %>%
  dplyr::inner_join(
    selected_interactions %>%
      dplyr::transmute(
        source = as.character(source),
        target = as.character(target),
        ligand = as.character(lig_gene),
        receptor = as.character(rec_gene)
      ),
    by = c("source", "target", "ligand", "receptor")
  )
if (nrow(matched_edges) != expected_edges) {
  stop("One or more fixed Figure 3C interactions are absent from the top 60.")
}
rm(selected_interactions, matched_edges)
invisible(gc())

message("Calculating and verifying the 5.1 GB source-object SHA256...")
input_info <- file.info(input_rds)
input_sha256 <- sha256_file(input_rds)
if (!identical(input_sha256, expected_input_sha256)) {
  stop("Validated Figure 1 input-object SHA256 has changed.")
}

message("Calculating H&E image checksums...")
he_image_sha256 <- vapply(unname(he_image_paths), sha256_file, character(1))

message("Loading the validated annotated Seurat object...")
spatial_object <- readRDS(input_rds)

# Compatibility repair for legacy VisiumV1 images serialized before
# SeuratObject required the misc slot. This changes only the in-memory copy.
legacy_visium_misc_slots_added <- character(0)

for (image_name in names(spatial_object@images)) {
  image_object <- spatial_object@images[[image_name]]

  if (inherits(image_object, "VisiumV1")) {
    misc_slot_present <- tryCatch(
      {
        methods::slot(image_object, "misc")
        TRUE
      },
      error = function(e) FALSE
    )

    if (!misc_slot_present) {
      methods::slot(image_object, "misc", check = FALSE) <- list()
      spatial_object@images[[image_name]] <- image_object
      legacy_visium_misc_slots_added <- c(
        legacy_visium_misc_slots_added,
        image_name
      )
    }
  }
}

for (image_name in names(spatial_object@images)) {
  methods::validObject(spatial_object@images[[image_name]])
}

message(
  "Legacy VisiumV1 misc slots added in memory: ",
  length(legacy_visium_misc_slots_added)
)

if (!inherits(spatial_object, "Seurat")) {
  stop("Spatial-map input is not a Seurat object.")
}
if (!all(c(assay_use, deconvolution_assay) %in% names(spatial_object@assays))) {
  stop("Required Spatial and singlecell assays are not both present.")
}
if (!all(c(sample_column, condition_column) %in% names(spatial_object@meta.data))) {
  stop("Required donor and F_score metadata columns are missing.")
}
if (ncol(spatial_object) != expected_spots ||
    nrow(spatial_object[[assay_use]]) != expected_input_genes) {
  stop("Unexpected source-object Spatial dimensions.")
}

spot_names <- colnames(spatial_object)
if (!identical(spot_names, rownames(spatial_object@meta.data))) {
  stop("Source-object spots and metadata are not identically ordered.")
}

deconvolution_values <- SeuratObject::LayerData(
  spatial_object[[deconvolution_assay]],
  layer = deconvolution_layer
)
if (!all(cell_types %in% rownames(deconvolution_values)) ||
    !identical(colnames(deconvolution_values), spot_names)) {
  stop("Deconvolution values do not match the retained cell types and spots.")
}

proportion_matrix <- t(as.matrix(
  deconvolution_values[cell_types, spot_names, drop = FALSE]
))
if (any(!is.finite(proportion_matrix)) ||
    min(proportion_matrix) < -1e-12 ||
    max(proportion_matrix) > 1 + 1e-12) {
  stop("Deconvolution proportions are non-finite or outside [0, 1].")
}
maximum_proportion_sum_error <- max(abs(rowSums(proportion_matrix) - 1))
if (maximum_proportion_sum_error > 1e-8) {
  stop("Deconvolution proportions do not sum to one per spot.")
}

for (cell_type in cell_types) {
  spatial_object@meta.data[[cell_type]] <- proportion_matrix[, cell_type]
}

observed_donor_stage <- spatial_object@meta.data %>%
  tibble::rownames_to_column("spot") %>%
  dplyr::transmute(
    donor = as.character(.data[[sample_column]]),
    stage = as.character(.data[[condition_column]])
  ) %>%
  dplyr::distinct()

if (anyDuplicated(observed_donor_stage$donor) ||
    !setequal(
      paste(observed_donor_stage$donor, observed_donor_stage$stage, sep = "||"),
      paste(donor_stage_map$donor, donor_stage_map$stage, sep = "||")
    )) {
  stop("Observed donor-stage assignments differ from the declared map.")
}

image_keys_before_diet <- names(spatial_object@images)
if (length(image_keys_before_diet) != expected_donors) {
  stop("Expected 24 source-object spatial images.")
}

spatial_counts <- SeuratObject::LayerData(
  spatial_object[[assay_use]],
  layer = input_layer
)
keep_genes <- rownames(spatial_counts)[Matrix::rowSums(spatial_counts) > 0]
zero_genes_dropped <- nrow(spatial_counts) - length(keep_genes)
if (zero_genes_dropped != expected_zero_genes_dropped ||
    length(keep_genes) != expected_retained_genes) {
  stop("Unexpected zero-gene filtering result.")
}

rm(deconvolution_values, proportion_matrix, spatial_counts)
invisible(gc())

message("Thinning the in-memory object to retained Spatial counts and images...")
DefaultAssay(spatial_object) <- assay_use
spatial_object <- Seurat::DietSeurat(
  spatial_object,
  layers = "counts",
  features = keep_genes,
  assays = assay_use,
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)

if (!identical(names(spatial_object@images), image_keys_before_diet) ||
    nrow(spatial_object) != expected_retained_genes ||
    ncol(spatial_object) != expected_spots) {
  stop("In-memory thinning did not preserve the expected Spatial object and images.")
}

message("Running in-memory Spatial LogNormalize with scale.factor = 10000...")
spatial_object <- Seurat::NormalizeData(
  spatial_object,
  assay = assay_use,
  normalization.method = "LogNormalize",
  scale.factor = spatial_normalization_scale_factor,
  verbose = FALSE
)
DefaultAssay(spatial_object) <- assay_use

normalized_expression <- SeuratObject::LayerData(
  spatial_object[[assay_use]],
  layer = "data"
)
required_genes <- unique(c(story_edges$ligand, story_edges$receptor))
if (!all(required_genes %in% rownames(normalized_expression))) {
  stop("A retained spatial LR ligand or receptor is absent after normalization.")
}
rm(normalized_expression)
invisible(gc())

message(
  "Calculating cell-type-weighted local LR scores donor by donor (",
  effective_neighbors, " neighbors including self)..."
)

score_matrix <- matrix(
  NA_real_,
  nrow = expected_spots,
  ncol = expected_edges,
  dimnames = list(spot_names, story_edges$score_column)
)
donor_spot_counts <- setNames(integer(expected_donors), donor_stage_map$donor)

for (donor_id in donor_stage_map$donor) {
  cells <- spot_names[
    as.character(spatial_object@meta.data[[sample_column]]) == donor_id
  ]
  if (!length(cells)) {
    stop("No spots found for donor: ", donor_id)
  }
  donor_spot_counts[[donor_id]] <- length(cells)

  donor_object <- suppressWarnings(spatial_object[, cells])
  donor_expression <- SeuratObject::LayerData(
    donor_object[[assay_use]],
    layer = "data"
  )
  donor_metadata <- donor_object@meta.data

  if (!identical(colnames(donor_expression), rownames(donor_metadata))) {
    stop("Expression and metadata order differ for donor: ", donor_id)
  }

  coordinates <- SeuratObject::GetTissueCoordinates(donor_object)
  required_coordinate_columns <- c("imagerow", "imagecol")
  if (!all(required_coordinate_columns %in% names(coordinates)) ||
      !setequal(rownames(coordinates), colnames(donor_expression))) {
    stop("Unexpected tissue-coordinate structure for donor: ", donor_id)
  }
  coordinates <- coordinates[
    colnames(donor_expression),
    required_coordinate_columns,
    drop = FALSE
  ]

  if (nrow(coordinates) <= effective_neighbors) {
    stop("Too few spots for the requested neighborhood in donor: ", donor_id)
  }

  nearest <- RANN::nn2(coordinates, k = effective_neighbors)
  neighbor_index <- nearest$nn.idx
  neighbor_distance <- nearest$nn.dists
  neighbor_weights <- exp(
    -(neighbor_distance^2) / (2 * gaussian_sigma^2)
  )

  for (edge_index in seq_len(nrow(story_edges))) {
    edge <- story_edges[edge_index, , drop = FALSE]
    ligand_raw <- as.numeric(donor_expression[edge$ligand, ])
    receptor_raw <- as.numeric(donor_expression[edge$receptor, ])
    source_weight <- as.numeric(donor_metadata[[edge$source]])
    target_weight <- as.numeric(donor_metadata[[edge$target]])
    source_weight[is.na(source_weight)] <- 0
    target_weight[is.na(target_weight)] <- 0

    weighted_ligand <- ligand_raw * source_weight
    weighted_receptor <- receptor_raw * target_weight
    neighborhood_receptor <- rowSums(
      neighbor_weights * matrix(
        weighted_receptor[neighbor_index],
        nrow = nrow(neighbor_index)
      ),
      na.rm = TRUE
    )
    score_matrix[colnames(donor_expression), edge$score_column] <-
      weighted_ligand * neighborhood_receptor
  }

  message("  completed donor ", donor_id, " (", length(cells), " spots)")
  rm(
    donor_object,
    donor_expression,
    donor_metadata,
    coordinates,
    nearest,
    neighbor_index,
    neighbor_distance,
    neighbor_weights
  )
  invisible(gc())
}

if (any(!is.finite(score_matrix)) || any(score_matrix < 0)) {
  stop("Computed local LR scores contain non-finite or negative values.")
}

log_score_matrix <- log1p(score_matrix * pseudocount_multiplier)
color_limits <- do.call(
  rbind,
  lapply(seq_len(nrow(story_edges)), function(edge_index) {
    edge <- story_edges[edge_index, , drop = FALSE]
    values <- log_score_matrix[, edge$score_column]
    finite_values <- values[is.finite(values)]
    if (!length(finite_values)) {
      stop("No finite local LR scores for edge: ", edge$pair)
    }
    limits <- as.numeric(stats::quantile(
      finite_values,
      probs = color_quantiles,
      na.rm = TRUE,
      names = FALSE
    ))
    if (!all(is.finite(limits)) || limits[1] >= limits[2]) {
      stop("Invalid 1st-99th percentile color limits for edge: ", edge$pair)
    }
    data.frame(
      edge_order = edge$edge_order,
      source = edge$source,
      target = edge$target,
      ligand = edge$ligand,
      receptor = edge$receptor,
      pair = edge$pair,
      score_column = edge$score_column,
      q01 = limits[1],
      q99 = limits[2],
      finite_spots = length(finite_values),
      stringsAsFactors = FALSE
    )
  })
)

for (edge_index in seq_len(nrow(story_edges))) {
  score_column <- story_edges$score_column[[edge_index]]
  spatial_object@meta.data[[score_column]] <-
    log_score_matrix[rownames(spatial_object@meta.data), score_column]
}

rm(score_matrix, log_score_matrix)
invisible(gc())

image_keys <- names(spatial_object@images)
donors_from_image_keys <- gsub("^X", "", image_keys)
if (!setequal(donors_from_image_keys, donor_stage_map$donor) ||
    anyDuplicated(donors_from_image_keys)) {
  stop("Spatial image keys do not map uniquely to the declared 24 donors.")
}
custom_image_paths <- setNames(
  unname(he_image_paths[donors_from_image_keys]),
  image_keys
)
if (anyNA(custom_image_paths) || any(!file.exists(custom_image_paths))) {
  stop("At least one custom H&E path could not be mapped to a spatial image.")
}

circle_theta <- seq(
  0,
  2 * pi,
  length.out = n_circle_points + 1L
)[seq_len(n_circle_points)]
circle_cos <- cos(circle_theta)
circle_sin <- sin(circle_theta)
image_geometry_rows <- vector("list", expected_donors)

message(
  "Generating 72 LRlocal PDFs with true Space Ranger Visium spot diameters " ,
  "in stage/pair/direction folders..."
)

for (donor_index in seq_len(nrow(donor_stage_map))) {
  donor_id <- donor_stage_map$donor[[donor_index]]
  cells <- spot_names[
    as.character(spatial_object@meta.data[[sample_column]]) == donor_id
  ]
  if (length(cells) != donor_spot_counts[[donor_id]]) {
    stop("Plotting spot count differs from the LR-score spot count for donor: ", donor_id)
  }

  image_key <- image_keys[donors_from_image_keys == donor_id]
  if (length(image_key) != 1L) {
    stop("Could not resolve exactly one source image for donor: ", donor_id)
  }
  image_object <- spatial_object@images[[image_key]]

  image_dimensions <- dim(image_object@image)
  if (length(image_dimensions) < 2L || any(image_dimensions[1:2] <= 0)) {
    stop("Source spatial image has invalid dimensions for donor: ", donor_id)
  }
  registered_image_height_pixels <- as.numeric(image_dimensions[[1]])
  registered_image_width_pixels <- as.numeric(image_dimensions[[2]])

  hires_scale_factor <- as.numeric(image_object@scale.factors$hires)
  spot_diameter_fullres <- as.numeric(image_object@scale.factors$spot)
  if (length(hires_scale_factor) != 1L ||
      length(spot_diameter_fullres) != 1L ||
      !is.finite(hires_scale_factor) ||
      !is.finite(spot_diameter_fullres) ||
      hires_scale_factor <= 0 ||
      spot_diameter_fullres <= 0) {
    stop("Invalid Space Ranger scale factors for donor: ", donor_id)
  }

  raw_image <- magick::image_read(custom_image_paths[[image_key]])
  raw_image_info <- magick::image_info(raw_image)
  if (nrow(raw_image_info) != 1L) {
    stop("Expected one-frame custom H&E image for donor: ", donor_id)
  }

  displayed_image <- magick::image_scale(
    raw_image,
    geometry = paste0("x", image_height_pixels)
  )
  displayed_image_info <- magick::image_info(displayed_image)
  displayed_image_width_pixels <- as.numeric(displayed_image_info$width[[1]])
  displayed_image_height_pixels <- as.numeric(displayed_image_info$height[[1]])
  if (displayed_image_height_pixels != image_height_pixels) {
    stop("Custom H&E image did not scale to the declared height for donor: ", donor_id)
  }

  registered_aspect_ratio <-
    registered_image_width_pixels / registered_image_height_pixels
  custom_aspect_ratio <-
    as.numeric(raw_image_info$width[[1]]) / as.numeric(raw_image_info$height[[1]])
  relative_aspect_ratio_difference <- abs(
    custom_aspect_ratio / registered_aspect_ratio - 1
  )
  if (!is.finite(relative_aspect_ratio_difference) ||
      relative_aspect_ratio_difference > 0.002) {
    stop(
      "Custom and registered H&E aspect ratios differ by more than 0.2% for donor: ",
      donor_id
    )
  }

  hires_to_display_scale <-
    max(displayed_image_width_pixels, displayed_image_height_pixels) /
    space_ranger_hires_max_dimension_pixels
  fullres_to_display_scale <- hires_scale_factor * hires_to_display_scale
  true_spot_diameter_pixels <-
    spot_diameter_fullres * fullres_to_display_scale * spot_multiplier
  true_spot_radius_pixels <- true_spot_diameter_pixels / 2
  donor_pdf_height_inches <- pdf_width_inches *
    displayed_image_height_pixels / displayed_image_width_pixels

  # Plotting requires full-resolution Space Ranger coordinates. The LR-score
  # calculation above intentionally uses GetTissueCoordinates() for the declared
  # low-resolution neighborhood scale; reusing those values here would apply
  # the image scale factor twice.
  tissue_coordinates <- image_object@coordinates
  required_coordinate_columns <- c("imagerow", "imagecol")
  if (!all(required_coordinate_columns %in% names(tissue_coordinates)) ||
      !setequal(rownames(tissue_coordinates), cells)) {
    stop("Unexpected plotting-coordinate structure for donor: ", donor_id)
  }
  tissue_coordinates <- tissue_coordinates[
    cells,
    required_coordinate_columns,
    drop = FALSE
  ]

  spot_centers <- data.frame(
    spot = cells,
    x = as.numeric(tissue_coordinates$imagecol) * fullres_to_display_scale,
    y = displayed_image_height_pixels -
      as.numeric(tissue_coordinates$imagerow) * fullres_to_display_scale,
    stringsAsFactors = FALSE
  )
  maximum_center_overrun_pixels <- max(
    0,
    -spot_centers$x,
    spot_centers$x - displayed_image_width_pixels,
    -spot_centers$y,
    spot_centers$y - displayed_image_height_pixels
  )
  if (!is.finite(maximum_center_overrun_pixels) ||
      maximum_center_overrun_pixels > 2) {
    stop("Transformed spot centers fall outside the custom H&E image for donor: ", donor_id)
  }

  displayed_image_raster <- as.raster(displayed_image)
  image_geometry_rows[[donor_index]] <- data.frame(
    donor = donor_id,
    image_key = image_key,
    registered_image_width_pixels = registered_image_width_pixels,
    registered_image_height_pixels = registered_image_height_pixels,
    custom_he_raw_width_pixels = as.numeric(raw_image_info$width[[1]]),
    custom_he_raw_height_pixels = as.numeric(raw_image_info$height[[1]]),
    displayed_image_width_pixels = displayed_image_width_pixels,
    displayed_image_height_pixels = displayed_image_height_pixels,
    registered_aspect_ratio = registered_aspect_ratio,
    custom_aspect_ratio = custom_aspect_ratio,
    relative_aspect_ratio_difference = relative_aspect_ratio_difference,
    space_ranger_hires_scale_factor = hires_scale_factor,
    space_ranger_spot_diameter_fullres = spot_diameter_fullres,
    hires_to_display_scale = hires_to_display_scale,
    fullres_to_display_scale = fullres_to_display_scale,
    spot_multiplier = spot_multiplier,
    circle_vertices = n_circle_points,
    true_visium_spot_diameter_pixels = true_spot_diameter_pixels,
    true_visium_spot_radius_pixels = true_spot_radius_pixels,
    maximum_center_overrun_pixels = maximum_center_overrun_pixels,
    pdf_width_inches = pdf_width_inches,
    pdf_height_inches = donor_pdf_height_inches,
    stringsAsFactors = FALSE
  )

  for (edge_index in seq_len(nrow(story_edges))) {
    edge <- story_edges[edge_index, , drop = FALSE]
    plan_row <- plot_plan[
      plot_plan$donor == donor_id & plot_plan$edge_order == edge$edge_order,
      ,
      drop = FALSE
    ]
    if (nrow(plan_row) != 1L) {
      stop("Could not resolve one planned output for donor/edge: ", donor_id, "/", edge$pair)
    }
    limits_row <- color_limits[
      color_limits$edge_order == edge$edge_order,
      ,
      drop = FALSE
    ]

    circle_data <- data.frame(
      spot = rep(spot_centers$spot, each = n_circle_points),
      x = rep(spot_centers$x, each = n_circle_points) +
        true_spot_radius_pixels * rep(circle_cos, times = nrow(spot_centers)),
      y = rep(spot_centers$y, each = n_circle_points) +
        true_spot_radius_pixels * rep(circle_sin, times = nrow(spot_centers)),
      score = rep(
        spatial_object@meta.data[cells, edge$score_column],
        each = n_circle_points
      ),
      stringsAsFactors = FALSE
    )

    local_lr_plot <- ggplot2::ggplot() +
      ggplot2::annotation_raster(
        displayed_image_raster,
        xmin = 0,
        xmax = displayed_image_width_pixels,
        ymin = 0,
        ymax = displayed_image_height_pixels
      ) +
      ggplot2::geom_polygon(
        data = circle_data,
        mapping = ggplot2::aes(x = x, y = y, group = spot, fill = score),
        colour = NA
      ) +
      ggplot2::scale_fill_viridis_c(
        limits = c(limits_row$q01, limits_row$q99),
        oob = scales::squish,
        name = "local LR score (log)"
      ) +
      ggplot2::coord_fixed(
        xlim = c(0, displayed_image_width_pixels),
        ylim = c(0, displayed_image_height_pixels),
        expand = FALSE
      ) +
      ggplot2::labs(
        title = paste0(
          donor_id, " local LR score: ", edge$source, " -> ", edge$target,
          " | ", edge$ligand, " -> ", edge$receptor,
          " (log) - ", edge$note
        )
      ) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
        plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
      )

    atomic_ggsave_pdf(
      plan_row$absolute_pdf_path,
      plot = local_lr_plot,
      width = pdf_width_inches,
      height = donor_pdf_height_inches
    )
    rm(circle_data, local_lr_plot)
  }

  message(
    "  plotted donor ", donor_id, " (3 LRlocal PDFs; true diameter ",
    sprintf("%.3f", true_spot_diameter_pixels), " displayed pixels)"
  )
  rm(
    tissue_coordinates,
    spot_centers,
    raw_image,
    raw_image_info,
    displayed_image,
    displayed_image_info,
    displayed_image_raster,
    image_object
  )
  invisible(gc())
}

image_geometry <- do.call(rbind, image_geometry_rows)
if (nrow(image_geometry) != expected_donors ||
    anyDuplicated(image_geometry$donor) ||
    any(!is.finite(image_geometry$true_visium_spot_diameter_pixels)) ||
    any(image_geometry$true_visium_spot_diameter_pixels <= 0)) {
  stop("The per-donor true-Visium geometry record is incomplete or invalid.")
}

rm(spatial_object, image_geometry_rows)
invisible(gc())

if (any(!file.exists(plot_plan$absolute_pdf_path)) ||
    any(file.info(plot_plan$absolute_pdf_path)$size <= 0)) {
  stop("The expected 72 non-empty LRlocal PDFs were not all generated.")
}

plot_plan$spots <- unname(donor_spot_counts[plot_plan$donor])
geometry_match <- match(plot_plan$donor, image_geometry$donor)
if (anyNA(geometry_match)) {
  stop("Could not map all planned PDFs to their donor-level image geometry.")
}
geometry_manifest_columns <- setdiff(names(image_geometry), "donor")
for (geometry_column in geometry_manifest_columns) {
  plot_plan[[geometry_column]] <- image_geometry[[geometry_column]][geometry_match]
}
plot_plan$color_q01 <- color_limits$q01[
  match(plot_plan$edge_order, color_limits$edge_order)
]
plot_plan$color_q99 <- color_limits$q99[
  match(plot_plan$edge_order, color_limits$edge_order)
]
plot_plan$pdf_bytes <- unname(file.info(plot_plan$absolute_pdf_path)$size)
plot_plan$pdf_sha256 <- vapply(
  plot_plan$absolute_pdf_path,
  sha256_file,
  character(1)
)

public_plot_manifest <- plot_plan %>%
  dplyr::select(-absolute_pdf_path)

parameter_table <- data.frame(
  parameter = c(
    "input_assay",
    "input_layer",
    "zero_count_spatial_genes_dropped_in_memory",
    "retained_spatial_genes",
    "spatial_normalization_method",
    "spatial_normalization_scale_factor",
    "local_ligand_signal",
    "local_receptor_signal",
    "nearest_neighbors_excluding_self",
    "include_self",
    "effective_neighbors_including_self",
    "gaussian_sigma",
    "local_lr_score",
    "score_transform",
    "fixed_color_limits",
    "spot_geometry",
    "spot_multiplier",
    "circle_vertices",
    "space_ranger_hires_max_dimension_pixels",
    "custom_he_image_height_pixels",
    "custom_he_aspect_ratio_tolerance",
    "plot_coordinate_transform",
    "pdf_width_inches",
    "pdf_height_inches",
    "output_nesting",
    "source_target_detection_outputs",
    "large_spot_level_score_object_saved",
    "legacy_visium_misc_slots_added_in_memory"
  ),
  value = c(
    assay_use,
    input_layer,
    zero_genes_dropped,
    expected_retained_genes,
    "LogNormalize",
    spatial_normalization_scale_factor,
    "normalized ligand expression * source cell-type proportion",
    "normalized receptor expression * target cell-type proportion",
    nearest_neighbors,
    include_self,
    effective_neighbors,
    gaussian_sigma,
    "weighted ligand at focal spot * Gaussian-weighted receptor neighborhood sum",
    "log1p(local_lr_score * 10000)",
    "per LR edge across all donors: 1st and 99th percentiles",
    paste0(
      "Space Ranger spot diameter fullres * hires scale factor * ",
      "(maximum displayed image dimension / 2000)"
    ),
    spot_multiplier,
    n_circle_points,
    space_ranger_hires_max_dimension_pixels,
    image_height_pixels,
    "relative difference <= 0.002",
    paste0(
      "VisiumV1@coordinates full-resolution; ",
      "x = imagecol * fullres_to_display_scale; ",
      "y = displayed_height - imagerow * fullres_to_display_scale"
    ),
    pdf_width_inches,
    "donor-specific: width * displayed image height / displayed image width",
    "stage/pair/source_to_target/donor_LRlocal.pdf",
    "not generated",
    "FALSE",
    length(legacy_visium_misc_slots_added)
  ),
  stringsAsFactors = FALSE
)

input_roles <- c(
  "validated_annotated_figure1_object",
  "selected_top60_interactions",
  "step02_output_checksum_manifest",
  "step02_completion_record",
  paste0("custom_he_image_", names(he_image_paths))
)
input_sha256_values <- c(
  input_sha256,
  sha256_file(selected_rds),
  sha256_file(step02_checksums_tsv),
  sha256_file(step02_completion_tsv),
  unname(he_image_sha256)
)
input_manifest <- data.frame(
  role = input_roles,
  path = normalizePath(required_inputs, mustWork = TRUE),
  bytes = unname(file.info(required_inputs)$size),
  modified = format(file.info(required_inputs)$mtime, tz = "UTC", usetz = TRUE),
  sha256 = input_sha256_values,
  stringsAsFactors = FALSE
)

session_lines <- capture.output(sessionInfo())

message("Saving Step 04 manifests, parameters, checksums, and software information...")
atomic_write_csv(story_edges, metadata_paths[["story_edges_csv"]])
atomic_write_csv(donor_stage_map, metadata_paths[["donor_stage_csv"]])
atomic_write_csv(color_limits, metadata_paths[["color_limits_csv"]])
atomic_write_csv(public_plot_manifest, metadata_paths[["plot_manifest_csv"]])
atomic_write_tsv(parameter_table, metadata_paths[["parameters_tsv"]])
atomic_write_tsv(input_manifest, metadata_paths[["input_manifest_tsv"]])
atomic_write_lines(session_lines, metadata_paths[["session_info_txt"]])

checksum_targets <- c(
  plot_plan$absolute_pdf_path,
  unname(metadata_paths)
)
output_checksums <- data.frame(
  file = vapply(checksum_targets, relative_to_output, character(1)),
  bytes = unname(file.info(checksum_targets)$size),
  sha256 = vapply(checksum_targets, sha256_file, character(1)),
  stringsAsFactors = FALSE
)

checksum_path <- file.path(output_dir, "04_OUTPUT_SHA256.tsv")
atomic_write_tsv(output_checksums, checksum_path)

finished_at <- Sys.time()

completion_table <- data.frame(
  field = c(
    "status",
    "started_at",
    "finished_at",
    "elapsed_minutes",
    "script_path",
    "script_sha256",
    "source_object_sha256",
    "donors",
    "spatial_lr_edges",
    "lrlocal_pdfs",
    "source_detect_pdfs",
    "target_detect_pdfs",
    "spot_level_score_rds_files",
    "spots",
    "retained_spatial_genes",
    "partial_files"
  ),
  value = c(
    "complete",
    format(started_at, tz = "UTC", usetz = TRUE),
    format(finished_at, tz = "UTC", usetz = TRUE),
    sprintf("%.3f", as.numeric(difftime(finished_at, started_at, units = "mins"))),
    script_path,
    if (is.na(script_path)) NA_character_ else sha256_file(script_path),
    input_sha256,
    nrow(donor_stage_map),
    nrow(story_edges),
    nrow(plot_plan),
    0L,
    0L,
    0L,
    expected_spots,
    expected_retained_genes,
    length(list.files(output_dir, pattern = "\\.partial$", recursive = TRUE))
  ),
  stringsAsFactors = FALSE
)

completion_path <- file.path(output_dir, "04_COMPLETION.tsv")
atomic_write_tsv(completion_table, completion_path)

message("Figure 3 spatial LR mapping completed at ", format(finished_at))
message("LRlocal PDFs generated: ", nrow(plot_plan))
message("Source-detection PDFs generated: 0")
message("Target-detection PDFs generated: 0")
message("Large derived RDS objects saved: 0")
message("Outputs: ", output_dir)
