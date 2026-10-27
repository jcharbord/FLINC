#!/usr/bin/env Rscript

# FLINC Figure 4: reproducible PROGENy spatial maps on H&E.
# This script intentionally uses circles whose radii are expressed in the
# image-coordinate system. It does not use point-size units.

options(stringsAsFactors = FALSE, digits = 8, scipen = 999)
set.seed(1L)

fail <- function(...) stop(paste0(...), call. = FALSE)

required_packages <- c("Seurat", "SeuratObject", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  fail("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L) {
  fail("Arguments must be supplied as --key value pairs")
}
arg_list <- as.list(args[seq(2L, length(args), by = 2L)])
names(arg_list) <- sub("^--", "", args[seq(1L, length(args), by = 2L)])

required_args <- c("input", "activity", "output-dir", "log-dir", "expected-input-sha256")
missing_args <- setdiff(required_args, names(arg_list))
if (length(missing_args) > 0L) {
  fail("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

input_path <- normalizePath(arg_list[["input"]], mustWork = TRUE)
activity_path <- normalizePath(arg_list[["activity"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
expected_input_sha <- tolower(arg_list[["expected-input-sha256"]])

if (dir.exists(output_dir) || file.exists(output_dir)) {
  fail("Output path already exists: ", output_dir)
}
if (dir.exists(log_dir) || file.exists(log_dir)) {
  fail("Log path already exists: ", log_dir)
}

parent_result <- dirname(output_dir)
parent_log <- dirname(log_dir)
dir.create(parent_result, recursive = TRUE, showWarnings = FALSE)
dir.create(parent_log, recursive = TRUE, showWarnings = FALSE)

run_stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
result_tmp <- paste0(output_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
log_tmp <- paste0(log_dir, ".incomplete-", run_stamp, "-pid", Sys.getpid())
dir.create(result_tmp, recursive = FALSE, showWarnings = FALSE)
dir.create(log_tmp, recursive = FALSE, showWarnings = FALSE)

completed <- FALSE
on.exit({
  if (!completed) {
    message("Run did not complete; incomplete directories were retained for inspection:\n",
            result_tmp, "\n", log_tmp)
  }
}, add = TRUE)

atomic_write_table <- function(x, path, row.names = FALSE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(
    x, temporary, sep = "\t", quote = FALSE, row.names = row.names,
    col.names = TRUE, na = ""
  )
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

script_path <- {
  full_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", full_arguments, value = TRUE)
  if (length(file_argument) == 1L) sub("^--file=", "", file_argument) else NA_character_
}
script_sha <- if (!is.na(script_path) && file.exists(script_path)) {
  digest::digest(script_path, algo = "sha256", file = TRUE)
} else {
  NA_character_
}

input_sha <- digest::digest(input_path, algo = "sha256", file = TRUE)
if (!identical(tolower(input_sha), expected_input_sha)) {
  fail("Input SHA256 mismatch. Expected ", expected_input_sha, "; observed ", input_sha)
}
activity_sha <- digest::digest(activity_path, algo = "sha256", file = TRUE)

slot_or_null <- function(x, slot_name) {
  if (slot_name %in% slotNames(x)) slot(x, slot_name) else NULL
}

as_single_numeric <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value)) {
    fail(label, " must be one finite numeric value")
  }
  value
}

make_colour_function <- function(limit_abs) {
  palette <- grDevices::colorRampPalette(c("#FFF5F0", "#FCAE91", "#67000D"))(257L)
  function(values) {
    clipped <- pmax(-limit_abs, pmin(limit_abs, values))
    index <- floor((clipped + limit_abs) / (2 * limit_abs) * 256) + 1L
    palette[pmax(1L, pmin(257L, index))]
  }
}

draw_spatial_map <- function(image_object, coordinates, values, radius, limit_abs,
                             title, draw_y_axis = FALSE) {
  image_array <- slot_or_null(image_object, "image")
  if (is.null(image_array) || length(dim(image_array)) != 3L) {
    fail("Image slot is not a three-channel raster array for ", title)
  }

  image_height <- dim(image_array)[1L]
  image_width <- dim(image_array)[2L]
  if (nrow(coordinates) != length(values)) {
    fail("Coordinate/value length mismatch for ", title)
  }

  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, image_width), ylim = c(image_height, 0), asp = 1,
    xaxs = "i", yaxs = "i"
  )
  graphics::rasterImage(
    as.raster(image_array),
    0, image_height, image_width, 0,
    interpolate = TRUE
  )

  colour_for_value <- make_colour_function(limit_abs)
  colours <- grDevices::adjustcolor(colour_for_value(values), alpha.f = 0.62)
  graphics::symbols(
    x = coordinates$imagecol,
    y = coordinates$imagerow,
    circles = rep(radius, nrow(coordinates)),
    inches = FALSE,
    bg = colours,
    fg = NA,
    add = TRUE
  )

  graphics::box(col = "#666666", lwd = 1.1)
  graphics::title(main = title, line = 0.5, cex.main = 1.1, font.main = 1)
  if (draw_y_axis) {
    graphics::mtext("H&E coordinate frame", side = 2, line = 2.1, cex = 0.65)
  }
}

draw_colour_legend <- function(limit_abs, pathway) {
  palette <- grDevices::colorRampPalette(c("#FFF5F0", "#FCAE91", "#67000D"))(257L)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  positions <- seq(0.15, 0.85, length.out = length(palette) + 1L)
  graphics::rect(0.22, positions[-length(positions)], 0.55, positions[-1L],
                 col = palette, border = NA)
  graphics::rect(0.22, 0.15, 0.55, 0.85, border = "#555555", lwd = 0.7)
  tick_positions <- c(0.15, 0.50, 0.85)
  tick_values <- c(-limit_abs, 0, limit_abs)
  graphics::segments(0.55, tick_positions, 0.59, tick_positions, col = "#333333")
  graphics::text(0.63, tick_positions, labels = format(signif(tick_values, 3), trim = TRUE),
                 adj = c(0, 0.5), cex = 0.72)
  graphics::text(0.5, 0.94, labels = pathway, cex = 0.85, font = 2)
  graphics::text(0.39, 0.06, labels = "PROGENy\nnorm_wmean", cex = 0.65)
}

write_individual_map <- function(path, image_object, coordinates, values, radius,
                                 limit_abs, title, width = 1800L, height = 1500L) {
  grDevices::png(path, width = width, height = height, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(1, 0.17))
  graphics::par(mar = c(0.5, 0.5, 2.6, 0.5))
  draw_spatial_map(image_object, coordinates, values, radius, limit_abs, title)
  graphics::par(mar = c(0.5, 0.1, 2.6, 0.2))
  draw_colour_legend(limit_abs, sub(" .*", "", title))
}

write_individual_map_pdf <- function(path, image_object, coordinates, values, radius,
                                     limit_abs, title) {
  grDevices::pdf(path, width = 9.0, height = 7.0, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(1, 0.17))
  graphics::par(mar = c(0.5, 0.5, 2.6, 0.5))
  draw_spatial_map(image_object, coordinates, values, radius, limit_abs, title)
  graphics::par(mar = c(0.5, 0.1, 2.6, 0.2))
  draw_colour_legend(limit_abs, sub(" .*", "", title))
}

write_panel_pdf <- function(path, pathways_to_draw, height_inches) {
  grDevices::pdf(path, width = 12.0, height = height_inches, useDingbats = FALSE)
  graphics::layout(
    matrix(seq_len(length(pathways_to_draw) * 3L), nrow = length(pathways_to_draw),
           ncol = 3L, byrow = TRUE),
    widths = c(1, 1, 0.17), heights = rep(1, length(pathways_to_draw))
  )

  for (pathway in pathways_to_draw) {
    limit_abs <- colour_limits$symmetric_limit_abs[colour_limits$pathway == pathway]
    for (image_name in sections$image) {
      section <- section_data[[image_name]]
      stage <- sections$stage[match(image_name, sections$image)]
      donor <- sections$donor[match(image_name, sections$image)]
      graphics::par(mar = c(0.5, 0.5, 2.5, 0.5))
      draw_spatial_map(
        image_object = section$image,
        coordinates = section$coordinates,
        values = as.numeric(section$values[pathway, ]),
        radius = section$radius,
        limit_abs = limit_abs,
        title = paste0(pathway, " - ", stage, " donor ", donor)
      )
    }
    graphics::par(mar = c(0.5, 0.1, 2.5, 0.2))
    draw_colour_legend(limit_abs, pathway)
  }
  grDevices::dev.off()
}

started_at <- Sys.time()
message("Loading regenerated Seurat object read-only...")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) fail("Input is not a Seurat object")
if (!"Spatial" %in% names(object@assays)) fail("Input does not contain Spatial assay")
if (!all(c("donor", "F_score") %in% colnames(object[[]]))) {
  fail("Input metadata must contain donor and F_score")
}

message("Loading definitive PROGENy activity object read-only...")
activity_object <- readRDS(activity_path)
if (!is.list(activity_object) || !all(c("activity", "spot_metadata") %in% names(activity_object))) {
  fail("Activity object must contain activity and spot_metadata elements")
}
activity_matrix <- activity_object$activity
spot_metadata <- activity_object$spot_metadata
if (!is.matrix(activity_matrix)) fail("Activity element is not a matrix")
if (!all(c("barcode", "donor", "F_score") %in% colnames(spot_metadata))) {
  fail("spot_metadata must contain barcode, donor, and F_score")
}
if (!identical(colnames(activity_matrix), spot_metadata$barcode)) {
  fail("Activity matrix columns are not identical to activity spot_metadata barcodes")
}
if (!identical(colnames(activity_matrix), colnames(object))) {
  fail("Activity matrix barcode order is not identical to regenerated Seurat object")
}
if (any(!is.finite(activity_matrix))) fail("Activity matrix contains non-finite values")

sections <- data.frame(
  image = c("X1A", "X6D"),
  donor = c("1A", "6D"),
  stage = c("F0", "F4"),
  stringsAsFactors = FALSE
)
pathways <- c("NFkB", "p53", "TNFa", "WNT", "Androgen")
main_pathways <- c("NFkB", "p53", "TNFa")
validation_pathways <- c("WNT", "Androgen")
missing_pathways <- setdiff(pathways, rownames(activity_matrix))
if (length(missing_pathways) > 0L) {
  fail("Missing requested pathways: ", paste(missing_pathways, collapse = ", "))
}

section_data <- list()
geometry_rows <- list()

for (index in seq_len(nrow(sections))) {
  image_name <- sections$image[[index]]
  if (!image_name %in% names(object@images)) fail("Missing image: ", image_name)

  image_object <- object[[image_name]]
  image_cells <- SeuratObject::Cells(image_object)
  metadata_rows <- match(image_cells, rownames(object[[]]))
  observed_donors <- unique(as.character(object[[]]$donor[metadata_rows]))
  observed_stages <- unique(as.character(object[[]]$F_score[metadata_rows]))
  if (!identical(observed_donors, sections$donor[[index]])) {
    fail("Image ", image_name, " does not map uniquely to donor ", sections$donor[[index]])
  }
  if (!identical(observed_stages, sections$stage[[index]])) {
    fail("Image ", image_name, " does not map uniquely to stage ", sections$stage[[index]])
  }

  coordinates <- SeuratObject::GetTissueCoordinates(image_object)
  if (!identical(rownames(coordinates), image_cells)) {
    fail("Coordinate order does not match image-cell order for ", image_name)
  }
  if (!all(c("imagerow", "imagecol") %in% colnames(coordinates))) {
    fail("Expected imagerow/imagecol coordinates for ", image_name)
  }
  if (!all(image_cells %in% colnames(activity_matrix))) {
    fail("Image cells are not all represented in activity matrix for ", image_name)
  }

  scale_factors <- slot_or_null(image_object, "scale.factors")
  if (is.null(scale_factors) || !all(c("spot", "lowres") %in% names(scale_factors))) {
    fail("Image ", image_name, " lacks spot and lowres scale factors")
  }
  lowres_spot_diameter <- as_single_numeric(scale_factors$spot, "spot scale factor") *
    as_single_numeric(scale_factors$lowres, "lowres scale factor")
  lowres_spot_radius <- lowres_spot_diameter / 2
  # Seurat's VisiumV1@spot.radius uses its own normalized coordinate system;
  # it is intentionally not used for a plot drawn in low-resolution H&E pixels.
  slot_radius_native <- as_single_numeric(
    slot_or_null(image_object, "spot.radius"),
    "spot.radius"
  )

  image_array <- slot_or_null(image_object, "image")
  if (is.null(image_array) || length(dim(image_array)) != 3L) {
    fail("Image raster unavailable for ", image_name)
  }
  image_height <- dim(image_array)[1L]
  image_width <- dim(image_array)[2L]
  if (any(coordinates$imagecol < 0 | coordinates$imagecol > image_width) ||
      any(coordinates$imagerow < 0 | coordinates$imagerow > image_height)) {
    fail("Coordinates fall outside the low-resolution image extent for ", image_name)
  }

  selected_values <- activity_matrix[pathways, image_cells, drop = FALSE]
  section_data[[image_name]] <- list(
    image = image_object,
    coordinates = coordinates,
    values = selected_values,
    radius = lowres_spot_radius,
    diameter = lowres_spot_diameter
  )

  geometry_rows[[image_name]] <- data.frame(
    image = image_name,
    donor = sections$donor[[index]],
    stage = sections$stage[[index]],
    n_spots = length(image_cells),
    image_height_lowres_px = image_height,
    image_width_lowres_px = image_width,
    spot_diameter_fullres_px = as_single_numeric(scale_factors$spot, "spot scale factor"),
    lowres_scale_factor = as_single_numeric(scale_factors$lowres, "lowres scale factor"),
    spot_diameter_lowres_px = lowres_spot_diameter,
    radius_used_lowres_px = lowres_spot_radius,
    stored_spot_radius_native = slot_radius_native,
    stored_spot_radius_units = "Seurat VisiumV1 normalized/native units; not low-resolution pixels",
    stringsAsFactors = FALSE
  )
}

geometry_table <- do.call(rbind, geometry_rows)
colour_limits <- do.call(rbind, lapply(pathways, function(pathway) {
  displayed_values <- unlist(lapply(section_data, function(x) x$values[pathway, ]), use.names = FALSE)
  lower_quantile <- as.numeric(stats::quantile(
    displayed_values, probs = 0.01, names = FALSE, type = 8
  ))
  upper_quantile <- as.numeric(stats::quantile(
    displayed_values, probs = 0.99, names = FALSE, type = 8
  ))
  limit_abs <- max(abs(lower_quantile), abs(upper_quantile))
  if (!is.finite(limit_abs) || limit_abs <= 0) fail("Invalid colour limit for ", pathway)
  data.frame(
    pathway = pathway,
    n_displayed_spots = length(displayed_values),
    observed_min = min(displayed_values),
    observed_max = max(displayed_values),
    lower_quantile_0.01 = lower_quantile,
    upper_quantile_0.99 = upper_quantile,
    symmetric_limit_abs = limit_abs,
    colour_scale_min = -limit_abs,
    colour_scale_max = limit_abs,
    n_clipped_low = sum(displayed_values < -limit_abs),
    n_clipped_high = sum(displayed_values > limit_abs),
    stringsAsFactors = FALSE
  )
}))

displayed_spots <- do.call(rbind, lapply(names(section_data), function(image_name) {
  section <- section_data[[image_name]]
  do.call(rbind, lapply(pathways, function(pathway) {
    data.frame(
      image = image_name,
      donor = sections$donor[match(image_name, sections$image)],
      stage = sections$stage[match(image_name, sections$image)],
      pathway = pathway,
      barcode = rownames(section$coordinates),
      imagerow = section$coordinates$imagerow,
      imagecol = section$coordinates$imagecol,
      activity_norm_wmean = as.numeric(section$values[pathway, ]),
      stringsAsFactors = FALSE
    )
  }))
}))

message("Rendering Figure 4B maps with scale-factor-derived Visium radii...")

write_panel_pdf(
  file.path(result_tmp, "Fig4B_PROGENy_spatial_maps_main_F0_vs_F4.pdf"),
  main_pathways,
  height_inches = 8.0
)
write_panel_pdf(
  file.path(result_tmp, "Fig4B_PROGENy_spatial_maps_validation_F0_vs_F4.pdf"),
  validation_pathways,
  height_inches = 5.5
)

for (pathway in pathways) {
  limit_abs <- colour_limits$symmetric_limit_abs[colour_limits$pathway == pathway]
  for (image_name in sections$image) {
    section <- section_data[[image_name]]
    stage <- sections$stage[match(image_name, sections$image)]
    donor <- sections$donor[match(image_name, sections$image)]
    individual_png <- file.path(
      result_tmp,
      paste0("Fig4B_", pathway, "_", stage, "_donor_", donor, ".png")
    )
    write_individual_map(
      path = individual_png,
      image_object = section$image,
      coordinates = section$coordinates,
      values = as.numeric(section$values[pathway, ]),
      radius = section$radius,
      limit_abs = limit_abs,
      title = paste0(pathway, " - ", stage, " donor ", donor)
    )
    individual_pdf <- file.path(
      result_tmp,
      paste0("Fig4B_", pathway, "_", stage, "_donor_", donor, ".pdf")
    )
    write_individual_map_pdf(
      path = individual_pdf,
      image_object = section$image,
      coordinates = section$coordinates,
      values = as.numeric(section$values[pathway, ]),
      radius = section$radius,
      limit_abs = limit_abs,
      title = paste0(pathway, " - ", stage, " donor ", donor)
    )
  }
}

parameters <- data.frame(
  parameter = c(
    "input_object", "input_sha256", "activity_input", "activity_sha256", "script_sha256",
    "displayed_stages", "representative_sections", "displayed_pathways",
    "coordinate_frame", "spot_size_rule", "spot_radius_rule", "colour_scale_rule",
    "raster_interpolation", "H_and_E_alpha", "circle_alpha", "circle_border",
    "individual_section_PDFs", "seed", "figure_role"
  ),
  value = c(
    input_path, input_sha, activity_path, activity_sha, script_sha,
    "F0;F4", "F0=donor_1A:image_X1A;F4=donor_6D:image_X6D",
    paste(pathways, collapse = ";"),
    "low-resolution H&E image pixels; imagecol=x, imagerow=y",
    "Visium low-resolution spot diameter = scale.factors$spot * scale.factors$lowres",
    "radius = (scale.factors$spot * scale.factors$lowres) / 2; image@spot.radius is normalized/native and is not used",
    "sequential red anchors: #FFF5F0 at minimum, #FCAE91 at zero, #67000D at maximum; per-pathway 1st/99th percentile display range across displayed F0/F4 spots; clipped counts recorded",
    "TRUE", "1.00", "0.62", "none", "one PDF and one PNG per pathway-section map",
    "1 (no stochastic plotting operations)",
    "NFkB;p53;TNFa main-panel representatives;WNT;Androgen validation maps"
  ),
  stringsAsFactors = FALSE
)

atomic_write_table(geometry_table, file.path(result_tmp, "04_spot_geometry_audit.tsv"))
atomic_write_table(colour_limits, file.path(result_tmp, "04_map_colour_limits.tsv"))
atomic_write_table(parameters, file.path(result_tmp, "04_run_parameters.tsv"))
utils::write.table(
  displayed_spots,
  gzfile(file.path(result_tmp, "04_displayed_spot_values.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = ""
)

package_versions <- data.frame(
  package = c("R", required_packages),
  version = c(
    R.version.string,
    vapply(required_packages, function(package_name) {
      as.character(utils::packageVersion(package_name))
    }, character(1))
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(package_versions, file.path(log_tmp, "04_package_versions.tsv"))

finished_at <- Sys.time()
execution_summary <- c(
  "FLINC Figure 4 PROGENy spatial-map execution summary",
  paste("Started UTC:", format(started_at, tz = "UTC", usetz = TRUE)),
  paste("Completed UTC:", format(finished_at, tz = "UTC", usetz = TRUE)),
  paste("Elapsed seconds:", as.numeric(difftime(finished_at, started_at, units = "secs"))),
  paste("Input SHA256:", input_sha),
  paste("Activity SHA256:", activity_sha),
  paste("Script SHA256:", script_sha),
  paste("Representative sections:", paste(paste(sections$stage, sections$donor, sections$image, sep = ":"), collapse = ";")),
  paste("True lowres spot diameters:", paste(sprintf("%s=%.8f", geometry_table$image, geometry_table$spot_diameter_lowres_px), collapse = ";")),
  "Spot circles are plotted in image-coordinate units with inches=FALSE.",
  "Main-figure representative maps: NFkB, p53, TNFa.",
  "Validation maps: WNT, Androgen."
)
atomic_write_lines(execution_summary, file.path(log_tmp, "04_execution_summary.txt"))
atomic_write_lines(capture.output(sessionInfo()), file.path(log_tmp, "04_sessionInfo.txt"))

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
manifest <- data.frame(
  file = basename(result_files),
  bytes = file.info(result_files)$size,
  sha256 = vapply(result_files, digest::digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
atomic_write_table(manifest, file.path(result_tmp, "04_output_SHA256.tsv"))
atomic_write_lines("SUCCESS", file.path(result_tmp, "04_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "04_COMPLETED.txt"))

if (!file.rename(result_tmp, output_dir)) fail("Could not atomically publish results to ", output_dir)
if (!file.rename(log_tmp, log_dir)) fail("Could not atomically publish logs to ", log_dir)

completed <- TRUE
message("PROGENy spatial maps completed successfully.")
message("Results: ", output_dir)
message("Compact logs: ", log_dir)
