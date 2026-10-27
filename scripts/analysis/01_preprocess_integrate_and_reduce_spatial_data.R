# ============================================================
# 01_preprocess_integrate_and_reduce_spatial_data.R
#
# Purpose:
#   Prepare the 24 FLINC FFPE spatial transcriptomics samples for
#   downstream biological analysis.
#
# Main operations:
#   - load Space Ranger outputs
#   - SCTransform normalization
#   - SCT integration across samples
#   - import previously generated cell-type deconvolution results
#   - PCA
#   - Harmony correction
#   - UMAP
#   - save intermediate and integrated analysis objects
#   - record parameters, package versions and checksums
#
# Reference execution on the validated Mac environment:
#
#   cd "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
#
#   FLINC_THREADS=4 caffeinate -dimsu Rscript \
#     scripts/analysis/01_preprocess_integrate_and_reduce_spatial_data.R \
#     > logs/preprocess_integration_$(date +%Y%m%d_%H%M%S).log 2>&1
# ============================================================

rm(list = ls())
gc()

# ============================================================
# 00. Reproducibility / resources
# ============================================================

N_THREADS <- as.integer(Sys.getenv("FLINC_THREADS", unset = "4"))
if (is.na(N_THREADS) || N_THREADS < 1) {
  N_THREADS <- 4
}

Sys.setenv(
  OMP_NUM_THREADS = as.character(N_THREADS),
  OPENBLAS_NUM_THREADS = as.character(N_THREADS),
  MKL_NUM_THREADS = as.character(N_THREADS),
  VECLIB_MAXIMUM_THREADS = as.character(N_THREADS),
  NUMEXPR_NUM_THREADS = as.character(N_THREADS)
)

set.seed(1234)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 30 * 1024^3)

REPRO_ROOT <- "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"

SCRIPT_DIR <- file.path(REPRO_ROOT, "scripts")
OBJECT_DIR <- file.path(REPRO_ROOT, "objects")
RESULTS_DIR <- file.path(REPRO_ROOT, "results")
REPORTS_DIR <- file.path(REPRO_ROOT, "reports")
LOG_DIR <- file.path(REPRO_ROOT, "logs")
REPRO_LOG_DIR <- file.path(REPRO_ROOT, "reproducibility_logs")

dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPRO_LOG_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(REPRO_ROOT)

# Original inputs
METADATA_FILE <- "/Users/jianwu/Documents/DATA/metadata/metadata_spatial_ffpe0.2.csv"
SPACERANGER_DIR <- "/Users/jianwu/Documents/DATA/spaceranger_output"
DECONV_FILE <- "/Users/jianwu/Documents/FLINC/data/deconvolution/single_nucleus_spot_celltype_frac.txt"

# Analysis outputs
object_list_raw_file <- file.path(OBJECT_DIR, "object.list.ffpe.rds")
object_list_sct_file <- file.path(OBJECT_DIR, "object.list.sct.ffpe.rds")
integration_features_file <- file.path(OBJECT_DIR, "integration.features.rdata")
integration_anchors_file <- file.path(OBJECT_DIR, "integration.anchors.ffpe.rds")
object_integrated_file <- file.path(OBJECT_DIR, "object.ffpe.rds")

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

write_checksum_csv <- function(files, out_file) {
  files <- files[!is.na(files) & file.exists(files)]

  if (length(files) == 0) {
    write.csv(
      data.frame(file = character(), md5 = character()),
      out_file,
      row.names = FALSE
    )
    return(invisible(NULL))
  }

  files_norm <- normalizePath(files, mustWork = TRUE)
  md5 <- tools::md5sum(files_norm)

  write.csv(
    data.frame(
      file = names(md5),
      md5 = unname(md5),
      stringsAsFactors = FALSE
    ),
    out_file,
    row.names = FALSE
  )

  invisible(md5)
}

save_session_info <- function(label) {
  sink(file.path(REPRO_LOG_DIR, paste0("sessionInfo_", label, ".txt")))
  cat("Timestamp:\n")
  print(Sys.time())
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nR version full:\n")
  print(R.version)
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
}

save_installed_packages <- function(label) {
  ip <- installed.packages()
  ip <- ip[, c("Package", "Version", "LibPath", "Built"), drop = FALSE]

  write.csv(
    ip,
    file.path(REPRO_LOG_DIR, paste0("installed_packages_", label, ".csv")),
    row.names = FALSE
  )
}

save_run_parameters <- function() {
  params <- data.frame(
    parameter = c(
      "REPRO_ROOT",
      "OBJECT_DIR",
      "RESULTS_DIR",
      "METADATA_FILE",
      "SPACERANGER_DIR",
      "DECONV_FILE",
      "N_THREADS",
      "sct_seed",
      "sct_variable_features_n",
      "integration_nfeatures",
      "integration_reference_index",
      "integration_anchor_seed",
      "pca_seed",
      "harmony_seed",
      "umap_seed",
      "pca_npcs",
      "harmony_dims",
      "umap_dims"
    ),
    value = c(
      REPRO_ROOT,
      OBJECT_DIR,
      RESULTS_DIR,
      METADATA_FILE,
      SPACERANGER_DIR,
      DECONV_FILE,
      as.character(N_THREADS),
      "5139",
      "8000",
      "8000",
      "1",
      "4751",
      "58291",
      "58291",
      "5147",
      "50",
      "1:40",
      "1:30"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    params,
    file.path(REPRO_LOG_DIR, "run_parameters_01_to_05.csv"),
    row.names = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(Matrix)
})

save_session_info("start_01_to_05")
save_installed_packages("start_01_to_05")
save_run_parameters()

write_checksum_csv(
  c(METADATA_FILE, DECONV_FILE),
  file.path(REPRO_LOG_DIR, "input_text_file_md5_checksums.csv")
)

log_msg("Starting spatial preprocessing, integration and dimensionality reduction")
log_msg("REPRO_ROOT:", REPRO_ROOT)
log_msg("Thread setting:", N_THREADS)

# ============================================================
# 01. Load data
# ============================================================

log_msg("STEP 01: reading metadata")

metadata <- read.csv(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot("linkID_spatial_ffpe" %in% colnames(metadata))

sample.list <- metadata$linkID_spatial_ffpe

if (anyDuplicated(sample.list)) {
  stop("Duplicated sample IDs found in metadata$linkID_spatial_ffpe")
}

sample_dirs <- file.path(SPACERANGER_DIR, sample.list)
missing_dirs <- sample_dirs[!dir.exists(sample_dirs)]

if (length(missing_dirs) > 0) {
  stop(
    "Missing Space Ranger output directories:\n",
    paste(missing_dirs, collapse = "\n")
  )
}

sample_order_table <- data.frame(
  index = seq_along(sample.list),
  sample_id = sample.list,
  spaceranger_dir = sample_dirs,
  is_reference_for_integration = seq_along(sample.list) == 1,
  stringsAsFactors = FALSE
)

write.csv(
  sample_order_table,
  file.path(REPRO_LOG_DIR, "sample_order_and_reference.csv"),
  row.names = FALSE
)

log_msg("Loading one Seurat spatial object per sample")

object.list <- lapply(seq_along(sample.list), function(index) {
  sample_id <- sample.list[index]

  log_msg("Loading sample", index, "of", length(sample.list), ":", sample_id)

  obj <- Load10X_Spatial(
    data.dir = file.path(SPACERANGER_DIR, sample_id)
  )

  # Keep Jian's original indexed cell renaming to preserve compatibility
  # with any downstream files generated from those spot/cell names.
  obj <- RenameCells(
    object = obj,
    add.cell.id = index
  )

  obj <- AddMetaData(
    object = obj,
    metadata = metadata[index, "linkID_spatial_ffpe"],
    col.name = "donor"
  )

  obj
})

names(object.list) <- sample.list

log_msg("Saving raw object.list:", object_list_raw_file)

saveRDS(object.list, file = object_list_raw_file)

write_checksum_csv(
  object_list_raw_file,
  file.path(REPRO_LOG_DIR, "step01_object_list_raw_md5.csv")
)

log_msg("STEP 01 complete")

# ============================================================
# 02. SCTransform normalization
# ============================================================

log_msg("STEP 02: SCTransform normalization")

object.list <- readRDS(object_list_raw_file)

object.list <- lapply(seq_along(object.list), function(i) {
  log_msg("SCTransform sample", i, "of", length(object.list), ":", names(object.list)[i])

  object <- object.list[[i]]

  object <- SCTransform(
    object,
    assay = "Spatial",
    verbose = FALSE,
    return.only.var.genes = FALSE,
    vars.to.regress = c("nFeature_Spatial"),
    variable.features.n = 8000,
    seed.use = 5139
  )

  object
})

log_msg("Saving SCT object.list:", object_list_sct_file)

saveRDS(object.list, file = object_list_sct_file)

write_checksum_csv(
  c(object_list_raw_file, object_list_sct_file),
  file.path(REPRO_LOG_DIR, "step02_sct_objects_md5.csv")
)

save_session_info("after_step02_SCT")

log_msg("STEP 02 complete")

# ============================================================
# 03. Integration
# ============================================================

log_msg("STEP 03: SCT integration")

object.list <- readRDS(object_list_sct_file)

log_msg("Selecting integration features")

integration.features <- SelectIntegrationFeatures(
  object.list = object.list,
  nfeatures = 8000
)

save(
  integration.features,
  file = integration_features_file
)

write.csv(
  data.frame(
    rank = seq_along(integration.features),
    gene = integration.features,
    stringsAsFactors = FALSE
  ),
  file.path(REPRO_LOG_DIR, "integration_features_gene_list.csv"),
  row.names = FALSE
)

log_msg("Running PrepSCTIntegration")

object.list <- PrepSCTIntegration(
  object.list = object.list,
  anchor.features = integration.features,
  verbose = FALSE
)

log_msg("Finding integration anchors; reference index = 1")
set.seed(4751)

integration.anchors <- FindIntegrationAnchors(
  object.list = object.list,
  normalization.method = "SCT",
  anchor.features = integration.features,
  verbose = FALSE,
  reference = 1
)

log_msg("Saving integration anchors:", integration_anchors_file)

saveRDS(
  integration.anchors,
  file = integration_anchors_file
)

write_checksum_csv(
  c(object_list_sct_file, integration_features_file, integration_anchors_file),
  file.path(REPRO_LOG_DIR, "step03_integration_inputs_anchors_md5.csv")
)

log_msg("Integrating data")

plan("sequential")
options(future.globals.maxSize = 30 * 1024^3)

object <- IntegrateData(
  anchorset = integration.anchors,
  normalization.method = "SCT",
  verbose = FALSE
)

VariableFeatures(object) <- integration.features

new_image_names <- paste0("X", unique(object$donor))

if (length(object@images) == length(new_image_names)) {
  names(object@images) <- new_image_names
} else {
  warning(
    "Number of images does not match number of unique donors. ",
    "Keeping original image names."
  )
}

log_msg("Saving integrated object before deconvolution/reduction:", object_integrated_file)

saveRDS(
  object,
  file = object_integrated_file
)

write_checksum_csv(
  c(integration_features_file, integration_anchors_file, object_integrated_file),
  file.path(REPRO_LOG_DIR, "step03_integrated_object_md5.csv")
)

save_session_info("after_step03_integration")

rm(integration.anchors)
gc()

log_msg("STEP 03 complete")

# ============================================================
# 04. Add deconvolution assay
# ============================================================

log_msg("STEP 04: add deconvolution assay")

object <- readRDS(object_integrated_file)

deconv_data <- read.table(
  DECONV_FILE,
  check.names = FALSE
)

log_msg("Deconvolution table dimensions:")
print(dim(deconv_data))

deconv_counts <- as.sparse(t(deconv_data))

if (ncol(deconv_counts) != ncol(object)) {
  stop(
    "Deconvolution matrix has ",
    ncol(deconv_counts),
    " columns after transpose, but object has ",
    ncol(object),
    " spots/cells. Check deconvolution file orientation/order."
  )
}

object@assays[["singlecell"]] <- CreateAssayObject(
  counts = deconv_counts
)

object@assays[["singlecell"]]@key <- "singlecell_"

log_msg("Assays after adding deconvolution:")
print(Assays(object))

log_msg("Saving object with deconvolution assay:", object_integrated_file)

saveRDS(
  object,
  file = object_integrated_file
)

write_checksum_csv(
  c(DECONV_FILE, object_integrated_file),
  file.path(REPRO_LOG_DIR, "step04_deconvolution_object_md5.csv")
)

save_session_info("after_step04_deconvolution")

gc()

log_msg("STEP 04 complete")

# ============================================================
# 05. PCA / Harmony / UMAP reduction
# ============================================================

log_msg("STEP 05: PCA / Harmony / UMAP")

suppressPackageStartupMessages({
  library(Rcpp)
  library(harmony)
})

object <- readRDS(object_integrated_file)
load(integration_features_file)

stopifnot(exists("integration.features"))

log_msg("Running PCA")

object <- RunPCA(
  object,
  assay = "SCT",
  verbose = FALSE,
  features = integration.features,
  npcs = 50,
  reduction.name = "pca_SCT",
  seed.use = 58291
)

log_msg("Running Harmony")

object <- RunHarmony(
  object,
  group.by.vars = c("donor"),
  theta = c(1),
  assay.use = "SCT",
  reduction = "pca_SCT",
  dims.use = 1:40,
  plot_convergence = TRUE,
  verbose = FALSE,
  reduction.save = "harmony_SCT",
  max.iter.harmony = 25,
  seed.use = 58291
)

log_msg("Running UMAP")

object <- RunUMAP(
  object,
  assay.use = "SCT",
  reduction = "harmony_SCT",
  dims = 1:30,
  seed.use = 5147,
  reduction.name = "umap_SCT",
  n.components = 2,
  densmap = FALSE
)

object@reductions[["umap_SCT"]]@cell.embeddings <-
  -object@reductions[["umap_SCT"]]@cell.embeddings

log_msg("Final assays:")
print(Assays(object))

log_msg("Final reductions:")
print(Reductions(object))

stopifnot("SCT" %in% Assays(object))
stopifnot("singlecell" %in% Assays(object))
stopifnot("pca_SCT" %in% Reductions(object))
stopifnot("harmony_SCT" %in% Reductions(object))
stopifnot("umap_SCT" %in% Reductions(object))

log_msg("Saving final post-reduction object:", object_integrated_file)

saveRDS(
  object,
  file = object_integrated_file
)

write_checksum_csv(
  c(
    object_list_raw_file,
    object_list_sct_file,
    integration_features_file,
    integration_anchors_file,
    object_integrated_file
  ),
  file.path(REPRO_LOG_DIR, "ALL_main_object_md5_checksums_after_step05.csv")
)

save_session_info("after_step05_reduction")
save_installed_packages("after_step05_reduction")

cat("\n============================================================\n")
cat("DONE: spatial preprocessing, integration and dimensionality reduction completed\n")
cat("Final object saved here:\n")
cat(object_integrated_file, "\n")
cat("\nAnalysis run records saved here:\n")
cat(REPRO_LOG_DIR, "\n")
cat("============================================================\n")
