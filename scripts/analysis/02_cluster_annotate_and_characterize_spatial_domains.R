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
# Analytical source:
#   The preserved clustering and characterization source is stored
#   under scripts/source_history/ as an R Markdown document.
#
#   This script extracts its R code with knitr::purl() and executes
#   the analysis directly. Document rendering is not required.
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
SOURCE_HISTORY_DIR <- file.path(SCRIPT_DIR, "source_history")
OBJECT_DIR <- file.path(REPRO_ROOT, "objects")
RESULTS_DIR <- file.path(REPRO_ROOT, "results")
LOG_DIR <- file.path(REPRO_ROOT, "logs")
REPRO_LOG_DIR <- file.path(REPRO_ROOT, "reproducibility_logs")

dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_HISTORY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPRO_LOG_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(REPRO_ROOT)

analysis_source_rmd <- file.path(SOURCE_HISTORY_DIR, "clustering_annotation_characterization_source.Rmd")
analysis_extracted_r <- file.path(SOURCE_HISTORY_DIR, "clustering_annotation_characterization_extracted.R")

if (!file.exists(analysis_source_rmd)) {
  stop(
    "Clustering and characterization source not found:\n",
    analysis_source_rmd,
    "\nSee scripts/source_history/README.md for details about the preserved analytical source."
  )
}

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
log_msg("Analytical source:", analysis_source_rmd)
log_msg("Executable R representation:", analysis_extracted_r)
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
  c(analysis_source_rmd, input_rds),
  file.path(REPRO_LOG_DIR, "input_md5_before_06_to_10_as_Rscript.csv")
)

# Extract R chunks from Rmd. This avoids Pandoc entirely.
suppressPackageStartupMessages({
  library(knitr)
})

log_msg("Extracting R chunks with knitr::purl()")

knitr::purl(
  input = analysis_source_rmd,
  output = analysis_extracted_r,
  documentation = 0,
  quiet = FALSE
)

write_checksum_csv(
  analysis_extracted_r,
  file.path(REPRO_LOG_DIR, "extracted_R_script_md5.csv")
)

log_msg("Sourcing extracted R script")
log_msg("This will run clustering, annotation, markers, GO, deconvolution, chord diagrams, proportions, and final checksums.")

# Run in a clean but normal environment.
run_env <- new.env(parent = globalenv())

source(
  file = analysis_extracted_r,
  local = run_env,
  echo = TRUE,
  max.deparse.length = Inf
)

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
