#!/usr/bin/env Rscript

# FLINC Figure 2 — construct the hdWGCNA network and modules.
#
# This script starts from 01_wgcna_input_prepared.rds and performs:
#   SetupForWGCNA -> MetaspotsByGroups -> NormalizeMetacells ->
#   SetDatExpr -> TestSoftPowers -> ConstructNetwork ->
#   ModuleEigengenes -> ModuleConnectivity -> module names/colors.
#
# It stops before differential module analysis, functional enrichment and
# figure generation.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(WGCNA)
  library(hdWGCNA)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

say <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

require_columns <- function(metadata, columns) {
  missing_columns <- setdiff(columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stopf(
      "Required metadata columns are missing: %s",
      paste(missing_columns, collapse = ", ")
    )
  }
}

move_into_place <- function(partial, path, replace = FALSE) {
  if (file.exists(path)) {
    if (!replace) {
      stopf("Refusing to overwrite existing output: %s", path)
    }
    if (!file.remove(path)) {
      stopf("Could not replace existing output: %s", path)
    }
  }
  if (!file.rename(partial, path)) {
    stopf("Could not move partial file into place: %s", path)
  }
}

write_lines_atomic <- function(lines, path, replace = FALSE) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) {
    stopf("Partial file already exists; inspect it before retrying: %s", partial)
  }
  writeLines(lines, partial, useBytes = TRUE)
  move_into_place(partial, path, replace = replace)
}

write_table_atomic <- function(x, path, replace = FALSE) {
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) {
    stopf("Partial file already exists; inspect it before retrying: %s", partial)
  }
  write.table(
    x,
    file = partial,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  move_into_place(partial, path, replace = replace)
}

save_rds_atomic <- function(object, path, label) {
  partial <- paste0(path, ".partial")
  if (file.exists(path)) {
    stopf("Refusing to overwrite existing RDS: %s", path)
  }
  if (file.exists(partial)) {
    stopf("Partial RDS already exists; inspect it before retrying: %s", partial)
  }
  say("Saving", label, "to", partial)
  saveRDS(object, partial, compress = "gzip", version = 3)
  move_into_place(partial, path, replace = FALSE)
}

save_plot_atomic <- function(plot, path, width, height) {
  partial <- paste0(path, ".partial.pdf")
  if (file.exists(partial)) {
    stopf("Partial plot already exists; inspect it before retrying: %s", partial)
  }
  ggplot2::ggsave(
    filename = partial,
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
  move_into_place(partial, path, replace = TRUE)
}

save_dendrogram_atomic <- function(object, path) {
  partial <- paste0(path, ".partial.pdf")
  if (file.exists(partial)) {
    stopf("Partial plot already exists; inspect it before retrying: %s", partial)
  }
  grDevices::pdf(partial, width = 20, height = 9, useDingbats = FALSE)
  plot_ok <- FALSE
  tryCatch(
    {
      PlotDendrogram(
        object,
        main = "Spatial hdWGCNA dendrogram - Whole slide"
      )
      plot_ok <- TRUE
    },
    finally = grDevices::dev.off()
  )
  if (!plot_ok) {
    stopf("Dendrogram generation failed: %s", partial)
  }
  move_into_place(partial, path, replace = TRUE)
}

check_exact_version <- function(package, expected) {
  observed <- as.character(packageVersion(package))
  if (!identical(observed, expected)) {
    stopf(
      "Package version mismatch for %s: expected %s, observed %s",
      package,
      expected,
      observed
    )
  }
}

validate_base_object <- function(object) {
  if (!inherits(object, "Seurat")) {
    stopf("Input is not a Seurat object; class is: %s", paste(class(object), collapse = ", "))
  }

  if (!identical(Assays(object), "Spatial")) {
    stopf("Expected only the Spatial assay; observed: %s", paste(Assays(object), collapse = ", "))
  }
  if (!identical(DefaultAssay(object), "Spatial")) {
    stopf("Expected Spatial as the default assay; observed: %s", DefaultAssay(object))
  }

  expected_layers <- c("counts", "data", "scale.data")
  observed_layers <- Layers(object[["Spatial"]])
  if (!all(expected_layers %in% observed_layers)) {
    stopf(
      "Spatial assay is missing required layers. Expected [%s], observed [%s]",
      paste(expected_layers, collapse = ", "),
      paste(observed_layers, collapse = ", ")
    )
  }

  if (nrow(object[["Spatial"]]) != 17922L || ncol(object) != 26056L) {
    stopf(
      "Unexpected dimensions: Spatial=%d features x %d spots",
      nrow(object[["Spatial"]]),
      ncol(object)
    )
  }

  required_metadata <- c(
    "F_score", "donor", "F_score_donor", "cluster", "annotation",
    "cluster_named", "row", "col", "imagerow", "imagecol", "barcode"
  )
  require_columns(object[[]], required_metadata)
  for (column_name in required_metadata) {
    if (anyNA(object[[column_name, drop = TRUE]])) {
      stopf("Required metadata column contains missing values: %s", column_name)
    }
  }

  if (length(unique(object$F_score_donor)) != 24L) {
    stopf(
      "Expected 24 F_score_donor samples; observed %d",
      length(unique(object$F_score_donor))
    )
  }
  if (length(Images(object)) != 24L) {
    stopf("Expected 24 spatial images; observed %d", length(Images(object)))
  }
  if (length(levels(Idents(object))) != 7L) {
    stopf("Expected 7 identity classes; observed %d", length(levels(Idents(object))))
  }
  if (!identical(as.character(Idents(object)), as.character(object$cluster_named))) {
    stopf("Seurat identities do not match cluster_named")
  }

  invisible(TRUE)
}

validate_metaspot_object <- function(metaspot_object, source_object) {
  if (!inherits(metaspot_object, "Seurat")) {
    stopf(
      "MetaspotsByGroups did not return a Seurat object; observed class: %s",
      paste(class(metaspot_object), collapse = ", ")
    )
  }

  if (!("Spatial" %in% Assays(metaspot_object))) {
    stopf(
      "Metaspot object is missing the Spatial assay; observed assays: %s",
      paste(Assays(metaspot_object), collapse = ", ")
    )
  }
  if (nrow(metaspot_object[["Spatial"]]) != nrow(source_object[["Spatial"]])) {
    stopf(
      "Metaspot feature count differs from the source object: %d versus %d",
      nrow(metaspot_object[["Spatial"]]),
      nrow(source_object[["Spatial"]])
    )
  }

  # MetaspotsByGroups constructs many spatial aggregates within each sample;
  # group.by defines the 24 slide boundaries, not the number of metaspots.
  # A valid result must therefore contain more than 24 metaspots, fewer than
  # the 26,056 source spots, and retain every F_score_donor group.
  n_metaspots <- ncol(metaspot_object)
  n_source_spots <- ncol(source_object)
  if (n_metaspots <= 24L || n_metaspots >= n_source_spots) {
    stopf(
      paste0(
        "Unexpected metaspot count: %d. Expected more than 24 sample groups ",
        "and fewer than %d source spots"
      ),
      n_metaspots,
      n_source_spots
    )
  }

  require_columns(
    metaspot_object[[]],
    c("F_score_donor", "row", "col", "imagerow", "imagecol")
  )
  if (anyNA(metaspot_object$F_score_donor)) {
    stopf("Metaspot F_score_donor metadata contains missing values")
  }

  expected_groups <- sort(unique(as.character(source_object$F_score_donor)))
  observed_groups <- sort(unique(as.character(metaspot_object$F_score_donor)))
  if (!identical(observed_groups, expected_groups)) {
    stopf(
      "Metaspot sample groups differ from the source object. Expected [%s], observed [%s]",
      paste(expected_groups, collapse = ", "),
      paste(observed_groups, collapse = ", ")
    )
  }

  invisible(TRUE)
}

validate_power_checkpoint <- function(object) {
  validate_base_object(object)
  power_table <- GetPowerTable(object, wgcna_name = "vis")
  if (!is.data.frame(power_table) || nrow(power_table) == 0L) {
    stopf("Soft-power checkpoint does not contain a valid power table")
  }
  required_power_columns <- c("Power", "SFT.R.sq")
  if (!all(required_power_columns %in% colnames(power_table))) {
    stopf(
      "Power table is missing required columns: %s",
      paste(setdiff(required_power_columns, colnames(power_table)), collapse = ", ")
    )
  }
  dat_expr <- GetDatExpr(object, wgcna_name = "vis")
  if (is.null(dim(dat_expr)) || any(dim(dat_expr) == 0L)) {
    stopf("Soft-power checkpoint does not contain a valid datExpr matrix")
  }
  invisible(TRUE)
}

validate_raw_network_checkpoint <- function(object) {
  validate_power_checkpoint(object)
  modules <- GetModules(object, wgcna_name = "vis")
  if (!is.data.frame(modules) || nrow(modules) == 0L) {
    stopf("Raw-network checkpoint does not contain module assignments")
  }
  params <- GetWGCNAParams(object, wgcna_name = "vis")
  if (!is.list(params) || is.null(params$power)) {
    stopf("Raw-network checkpoint does not contain the selected soft power")
  }
  invisible(TRUE)
}

repro_root <- Sys.getenv(
  "FLINC_REPRO_ROOT",
  unset = "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
)

input_file <- Sys.getenv(
  "FIG2_PREPARED_INPUT",
  unset = file.path(repro_root, "objects", "Fig2_WGCNA", "01_wgcna_input_prepared.rds")
)

object_dir <- file.path(repro_root, "objects", "Fig2_WGCNA")
checkpoint_dir <- file.path(object_dir, "checkpoints")
workspace_dir <- file.path(object_dir, "02_network_workspace")
results_dir <- file.path(repro_root, "results", "Fig2_WGCNA", "02_network")
repro_log_dir <- file.path(repro_root, "reproducibility_logs", "Fig2_WGCNA")

softpower_checkpoint <- file.path(checkpoint_dir, "02a_wgcna_softpowers.rds")
raw_network_checkpoint <- file.path(checkpoint_dir, "02b_wgcna_network_raw.rds")
output_file <- Sys.getenv(
  "FIG2_NETWORK_OUTPUT",
  unset = file.path(object_dir, "02_wgcna_network_modules.rds")
)

power_table_file <- file.path(results_dir, "02_soft_power_table.tsv")
softpower_plot_file <- file.path(results_dir, "02_soft_powers.pdf")
dendrogram_file <- file.path(results_dir, "02_module_dendrogram.pdf")
module_table_file <- file.path(results_dir, "02_module_assignments_and_kME.tsv")
module_size_file <- file.path(results_dir, "02_module_sizes.tsv")
hmes_file <- file.path(results_dir, "02_harmonized_module_eigengenes.tsv")
params_file <- file.path(repro_log_dir, "02_wgcna_parameters.txt")
summary_file <- file.path(repro_log_dir, "02_wgcna_network_summary.tsv")
session_file <- file.path(repro_log_dir, "02_wgcna_network_sessionInfo.txt")

dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(workspace_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_log_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(output_file)) {
  stopf("Final network output already exists; do not rerun step 02: %s", output_file)
}
if (file.exists(paste0(output_file, ".partial"))) {
  stopf(
    "A partial final output exists; inspect it before retrying: %s",
    paste0(output_file, ".partial")
  )
}

if (!identical(paste(R.version$major, R.version$minor, sep = "."), "4.4.1")) {
  stopf("R version mismatch: expected 4.4.1, observed %s", R.version.string)
}
check_exact_version("Seurat", "5.3.0")
check_exact_version("SeuratObject", "5.2.0")
check_exact_version("hdWGCNA", "0.4.6")
check_exact_version("WGCNA", "1.73")
check_exact_version("harmony", "1.2.4")

say("R version:", R.version.string)
say("Seurat version:", as.character(packageVersion("Seurat")))
say("SeuratObject version:", as.character(packageVersion("SeuratObject")))
say("hdWGCNA version:", as.character(packageVersion("hdWGCNA")))
say("WGCNA version:", as.character(packageVersion("WGCNA")))
say("harmony version:", as.character(packageVersion("harmony")))
say("Random seed: 12345")

# Do not alter WGCNA's threading state in this script.

if (file.exists(raw_network_checkpoint)) {
  say("Resuming from raw-network checkpoint:", raw_network_checkpoint)
  object <- readRDS(raw_network_checkpoint)
  validate_raw_network_checkpoint(object)
} else {
  if (file.exists(softpower_checkpoint)) {
    say("Resuming from soft-power checkpoint:", softpower_checkpoint)
    object <- readRDS(softpower_checkpoint)
    validate_power_checkpoint(object)
  } else {
    if (!file.exists(input_file)) {
      stopf("Prepared step-01 input does not exist: %s", input_file)
    }

    say("Loading prepared step-01 input:", input_file)
    load_start <- Sys.time()
    object <- readRDS(input_file)
    say(
      "Input loaded in minutes:",
      sprintf("%.2f", as.numeric(difftime(Sys.time(), load_start, units = "mins")))
    )
    validate_base_object(object)

    if ("vis" %in% names(object@misc)) {
      stopf("Input already contains an hdWGCNA experiment named 'vis'")
    }

    say("Running SetupForWGCNA: gene_select=fraction, fraction=0.05, wgcna_name=vis")
    object <- SetupForWGCNA(
      object,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = "vis"
    )

    wgcna_genes <- GetWGCNAGenes(object, wgcna_name = "vis")
    if (length(wgcna_genes) == 0L) {
      stopf("SetupForWGCNA selected no genes")
    }
    say("WGCNA genes selected:", length(wgcna_genes))

    say("Running MetaspotsByGroups by F_score_donor with Spatial counts")
    object <- MetaspotsByGroups(
      object,
      group.by = c("F_score_donor"),
      ident.group = "F_score_donor",
      assay = "Spatial"
    )

    metaspot_object <- GetMetacellObject(object, wgcna_name = "vis")
    validate_metaspot_object(metaspot_object, object)
    say(
      "Metaspots constructed:", ncol(metaspot_object),
      "across", length(unique(metaspot_object$F_score_donor)), "sample groups"
    )

    say("Normalizing metaspots with NormalizeMetacells defaults")
    object <- NormalizeMetacells(object, wgcna_name = "vis")

    say("Setting the whole-slide Spatial expression matrix")
    object <- SetDatExpr(
      object,
      group.by = NULL,
      group_name = NULL,
      assay = "Spatial",
      wgcna_name = "vis"
    )

    dat_expr <- GetDatExpr(object, wgcna_name = "vis")
    if (is.null(dim(dat_expr)) || any(dim(dat_expr) == 0L)) {
      stopf("SetDatExpr produced an invalid expression matrix")
    }
    say("datExpr dimensions:", nrow(dat_expr), "observations x", ncol(dat_expr), "genes")

    say("Testing soft powers with hdWGCNA defaults")
    object <- TestSoftPowers(object, wgcna_name = "vis")
    validate_power_checkpoint(object)

    power_table <- GetPowerTable(object, wgcna_name = "vis")
    write_table_atomic(power_table, power_table_file, replace = TRUE)
    save_rds_atomic(object, softpower_checkpoint, "soft-power checkpoint")
    say("Soft-power checkpoint completed")
  }

  old_working_directory <- getwd()
  on.exit(setwd(old_working_directory), add = TRUE)
  setwd(workspace_dir)

  # Keep tom_outdir relative because hdWGCNA 0.4.6 records the TOM path by
  # prepending getwd(). The controlled working directory makes that path valid.
  dir.create("TOM", recursive = TRUE, showWarnings = FALSE)
  tom_path <- normalizePath(file.path("TOM", "test_TOM.rda"), mustWork = FALSE)
  if (file.exists(tom_path)) {
    say("A TOM from an interrupted attempt exists and will be replaced:", tom_path)
  }

  say("Constructing the network with hdWGCNA 0.4.6 defaults")
  network_start <- Sys.time()
  object <- ConstructNetwork(
    object,
    tom_name = "test",
    overwrite_tom = TRUE,
    wgcna_name = "vis"
  )
  say(
    "Network construction time (minutes):",
    sprintf("%.2f", as.numeric(difftime(Sys.time(), network_start, units = "mins")))
  )

  validate_raw_network_checkpoint(object)
  params <- GetWGCNAParams(object, wgcna_name = "vis")
  selected_power <- as.numeric(params$power)
  if (length(selected_power) != 1L || is.na(selected_power)) {
    stopf("ConstructNetwork did not record exactly one selected soft power")
  }
  say("ConstructNetwork selected soft power:", selected_power)

  network_data <- GetNetworkData(object, wgcna_name = "vis")
  tom_files <- unique(as.character(unlist(network_data$TOMFiles)))
  if (length(tom_files) == 0L || any(!file.exists(tom_files))) {
    stopf(
      "The network object does not point to an existing TOM. Recorded paths: %s",
      paste(tom_files, collapse = ", ")
    )
  }
  say("TOM:", paste(tom_files, collapse = ", "))

  setwd(old_working_directory)
  save_rds_atomic(object, raw_network_checkpoint, "raw-network checkpoint")
  say("Raw-network checkpoint completed")
}

params <- GetWGCNAParams(object, wgcna_name = "vis")
selected_power <- as.numeric(params$power)
raw_modules <- GetModules(object, wgcna_name = "vis")
raw_non_grey_modules <- setdiff(unique(as.character(raw_modules$module)), "grey")
say("Raw non-grey modules:", length(raw_non_grey_modules))

# The Figure 2 analysis contains six biological modules. Preserve the raw
# network checkpoint if this assertion fails so the network can be inspected
# without repeating its construction.
if (length(raw_non_grey_modules) != 6L) {
  raw_size_table <- as.data.frame(table(raw_modules$module), stringsAsFactors = FALSE)
  colnames(raw_size_table) <- c("module", "n_genes")
  write_table_atomic(raw_size_table, module_size_file, replace = TRUE)
  stopf(
    paste0(
      "Expected 6 non-grey modules, but network construction produced %d. ",
      "The raw checkpoint was preserved at %s"
    ),
    length(raw_non_grey_modules),
    raw_network_checkpoint
  )
}

say("Computing harmonized module eigengenes by F_score_donor")
object <- ModuleEigengenes(
  object,
  group.by.vars = "F_score_donor",
  wgcna_name = "vis"
)

say("Computing harmonized module connectivity (kME)")
object <- ModuleConnectivity(
  object,
  harmonized = TRUE,
  wgcna_name = "vis"
)

say("Renaming modules with the SM_ prefix")
object <- ResetModuleNames(
  object,
  new_name = "SM_",
  wgcna_name = "vis"
)

modules_after_names <- GetModules(object, wgcna_name = "vis")
non_grey_after_names <- setdiff(levels(modules_after_names$module), "grey")
expected_module_names <- paste0("SM_", 1:6)
if (!identical(non_grey_after_names, expected_module_names)) {
  stopf(
    "Unexpected renamed modules. Expected [%s], observed [%s]",
    paste(expected_module_names, collapse = ", "),
    paste(non_grey_after_names, collapse = ", ")
  )
}

module_palette <- c(
  "#CAB2D6", "#6A3D9A", "#E69F00", "#7EC0EE", "#993F00", "#023FA5",
  "#117733", "#A52A2A", "#FDB76F", "#E78AC3", "#B3B3B3", "#8B4500",
  "#D3AE7B", "#C39A5B", "#B28634", "#9C7217", "#6B4A00", "#E8A3B8",
  "#BE7E8A", "#B46878", "#A94F64", "#9C3450", "#8E063B", "#8E063B"
)

say("Applying the Figure 2 module-color palette")
object <- ResetModuleColors(
  object,
  module_palette[seq_along(expected_module_names)],
  wgcna_name = "vis"
)

modules <- GetModules(object, wgcna_name = "vis")
observed_module_names <- setdiff(levels(modules$module), "grey")
if (!identical(observed_module_names, expected_module_names)) {
  stopf("Module names changed unexpectedly after ResetModuleColors")
}

hmes <- GetMEs(object, harmonized = TRUE, wgcna_name = "vis")
if (!is.matrix(hmes) && !is.data.frame(hmes)) {
  stopf("GetMEs did not return a matrix/data.frame")
}
hmes <- as.data.frame(hmes, check.names = FALSE)
if (!identical(rownames(hmes), colnames(object))) {
  if (!setequal(rownames(hmes), colnames(object))) {
    stopf("Harmonized module eigengene rows do not match object spots")
  }
  hmes <- hmes[colnames(object), , drop = FALSE]
}
if (!all(expected_module_names %in% colnames(hmes))) {
  stopf(
    "Harmonized module eigengenes are missing expected columns: %s",
    paste(setdiff(expected_module_names, colnames(hmes)), collapse = ", ")
  )
}

duplicate_hme_columns <- intersect(colnames(hmes), colnames(object[[]]))
if (length(duplicate_hme_columns) > 0L) {
  stopf(
    "Refusing to overwrite existing metadata columns with hMEs: %s",
    paste(duplicate_hme_columns, collapse = ", ")
  )
}
object <- AddMetaData(object, metadata = hmes)

module_sizes <- as.data.frame(table(modules$module), stringsAsFactors = FALSE)
colnames(module_sizes) <- c("module", "n_genes")
non_grey_sizes <- module_sizes[module_sizes$module != "grey", , drop = FALSE]
grey_gene_count <- sum(module_sizes$n_genes[module_sizes$module == "grey"])

power_table <- GetPowerTable(object, wgcna_name = "vis")
plot_list <- PlotSoftPowers(object, wgcna_name = "vis")
softpower_plot <- patchwork::wrap_plots(plot_list, ncol = 4)

save_plot_atomic(
  softpower_plot,
  softpower_plot_file,
  width = 14,
  height = 4
)
save_dendrogram_atomic(object, dendrogram_file)

write_table_atomic(power_table, power_table_file, replace = TRUE)
write_table_atomic(modules, module_table_file, replace = TRUE)
write_table_atomic(module_sizes, module_size_file, replace = TRUE)

hmes_output <- data.frame(
  barcode = rownames(hmes),
  hmes,
  row.names = NULL,
  check.names = FALSE
)
write_table_atomic(hmes_output, hmes_file, replace = TRUE)

params_lines <- c(
  "# GetWGCNAParams(object, wgcna_name = 'vis')",
  capture.output(dput(params)),
  "",
  "# formals(hdWGCNA::ConstructNetwork)",
  capture.output(dput(formals(hdWGCNA::ConstructNetwork))),
  "",
  "# formals(hdWGCNA::TestSoftPowers)",
  capture.output(dput(formals(hdWGCNA::TestSoftPowers))),
  "",
  "# formals(hdWGCNA::MetaspotsByGroups)",
  capture.output(dput(formals(hdWGCNA::MetaspotsByGroups)))
)
write_lines_atomic(params_lines, params_file, replace = TRUE)

metaspot_object <- GetMetacellObject(object, wgcna_name = "vis")
dat_expr <- GetDatExpr(object, wgcna_name = "vis")
network_data <- GetNetworkData(object, wgcna_name = "vis")
tom_files <- unique(as.character(unlist(network_data$TOMFiles)))

summary_table <- rbind(
  data.frame(metric = "input_file", value = input_file),
  data.frame(metric = "R_version", value = paste(R.version$major, R.version$minor, sep = ".")),
  data.frame(metric = "Seurat_version", value = as.character(packageVersion("Seurat"))),
  data.frame(metric = "SeuratObject_version", value = as.character(packageVersion("SeuratObject"))),
  data.frame(metric = "hdWGCNA_version", value = as.character(packageVersion("hdWGCNA"))),
  data.frame(metric = "WGCNA_version", value = as.character(packageVersion("WGCNA"))),
  data.frame(metric = "harmony_version", value = as.character(packageVersion("harmony"))),
  data.frame(metric = "random_seed", value = "12345"),
  data.frame(metric = "wgcna_name", value = "vis"),
  data.frame(metric = "gene_select", value = "fraction"),
  data.frame(metric = "gene_fraction", value = "0.05"),
  data.frame(metric = "metaspot_group_by", value = "F_score_donor"),
  data.frame(metric = "n_metaspots", value = ncol(metaspot_object)),
  data.frame(metric = "n_wgcna_genes", value = nrow(modules)),
  data.frame(metric = "datExpr_n_observations", value = nrow(dat_expr)),
  data.frame(metric = "datExpr_n_genes", value = ncol(dat_expr)),
  data.frame(metric = "selected_soft_power", value = selected_power),
  data.frame(metric = "n_non_grey_modules", value = length(expected_module_names)),
  data.frame(metric = "n_grey_genes", value = grey_gene_count),
  data.frame(metric = "tom_files", value = paste(tom_files, collapse = ",")),
  data.frame(metric = "final_output", value = output_file),
  data.frame(
    metric = paste0("genes_", non_grey_sizes$module),
    value = non_grey_sizes$n_genes
  )
)

write_table_atomic(summary_table, summary_file, replace = TRUE)
write_lines_atomic(capture.output(sessionInfo()), session_file, replace = TRUE)
save_rds_atomic(object, output_file, "final step-02 network object")

say("Final object size (GiB):", sprintf("%.2f", file.info(output_file)$size / 1024^3))
say("Selected soft power:", selected_power)
say("Non-grey modules:", paste(expected_module_names, collapse = ", "))
say("Module sizes:")
print(non_grey_sizes, row.names = FALSE)
say("Power table:", power_table_file)
say("Module table:", module_table_file)
say("Summary:", summary_file)
say("Session information:", session_file)
say("Fig. 2 hdWGCNA network construction completed successfully")
