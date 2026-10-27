# ============================================================
# rebuild_clustering_source_from_historical_analysis.R
#
# Prepare the preserved clustering and characterization source from the earlier analysis workflow.
#   /Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code/objects/object.ffpe.rds
#
# and saves outputs under:
#   /Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code/results
#
# Run:
#   cd "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
#   FLINC_THREADS=4 caffeinate -dimsu Rscript scripts/source_history/rebuild_clustering_source_from_historical_analysis.R \
#     > logs/06_to_10_$(date +%Y%m%d_%H%M%S).log 2>&1
#
# If it cannot find the original Rmd automatically, run:
#   Rscript scripts/source_history/rebuild_clustering_source_from_historical_analysis.R "/full/path/to/original.Rmd"
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
REPORTS_DIR <- file.path(REPRO_ROOT, "reports")
LOG_DIR <- file.path(REPRO_ROOT, "logs")
REPRO_LOG_DIR <- file.path(REPRO_ROOT, "reproducibility_logs")

dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPRO_LOG_DIR, recursive = TRUE, showWarnings = FALSE)

input_rds_new <- file.path(OBJECT_DIR, "object.ffpe.rds")
out_dir_new <- file.path(RESULTS_DIR, "260701_approach1_cluster_go_marker_deconv")

if (!file.exists(input_rds_new)) {
  stop("Missing input object from steps 01-05:\n", input_rds_new)
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

# ------------------------------------------------------------
# Locate original Rmd
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1 && nzchar(args[1])) {
  original_rmd <- args[1]
} else {
  code_final <- "/Users/jianwu/Documents/FLINC/Jeremie/Code_final"
  candidates <- list.files(code_final, pattern = "\\.Rmd$", full.names = TRUE, recursive = FALSE)

  if (length(candidates) == 0) {
    stop("No Rmd files found in:\n", code_final, "\nPass the original Rmd path as an argument.")
  }

  nm <- basename(candidates)
  scores <- rep(0, length(candidates))
  scores <- scores + ifelse(grepl("Approach", nm, ignore.case = TRUE), 3, 0)
  scores <- scores + ifelse(grepl("Clustering|cluster", nm, ignore.case = TRUE), 3, 0)
  scores <- scores + ifelse(grepl("Deconvolution|deconv", nm, ignore.case = TRUE), 2, 0)
  scores <- scores + ifelse(grepl("260630|260628|Approach 1", nm, ignore.case = TRUE), 2, 0)

  candidates <- candidates[order(scores, file.info(candidates)$mtime, decreasing = TRUE)]
  original_rmd <- candidates[1]
}

if (!file.exists(original_rmd)) {
  stop("Could not find original Rmd:\n", original_rmd)
}

log_msg("Original Rmd:", original_rmd)

# ------------------------------------------------------------
# Patch Rmd
# ------------------------------------------------------------

x <- paste(readLines(original_rmd, warn = FALSE), collapse = "\n")

# Fix unclosed text block after output path, if present.
x <- sub(
  pattern = "/Users/jianwu/Documents/FLINC/Jeremie/Code_final/260628_approach1_cluster_go_marker_deconv\\n\\s*\\n\\s*\\n\\s*# all together",
  replacement = "/Users/jianwu/Documents/FLINC/Jeremie/Code_final/260628_approach1_cluster_go_marker_deconv\n```\n\n# all together",
  x = x,
  perl = TRUE
)

# Fix duplicated chunk opening before proportion section, if present.
x <- gsub(
  pattern = "```\\{r\\}\\n```\\{r\\}\\n# ============================================================\\n# 10\\.9b",
  replacement = "```{r}\n# ============================================================\n# 10.9b",
  x = x,
  perl = TRUE
)

# Patch paths.
x <- gsub(
  pattern = 'setwd\\("/Users/jianwu/Documents/FLINC/Jeremie/Code_final"\\)',
  replacement = paste0('setwd("', REPRO_ROOT, '")'),
  x = x,
  perl = TRUE
)

x <- gsub(
  pattern = 'input_rds <- "/Users/jianwu/Documents/FLINC/Jeremie/Object/object.ffpe.rds"',
  replacement = paste0('input_rds <- "', input_rds_new, '"'),
  x = x,
  fixed = TRUE
)

x <- gsub(
  pattern = 'out_dir <- "/Users/jianwu/Documents/FLINC/Jeremie/Code_final/260628_approach1_cluster_go_marker_deconv"',
  replacement = paste0('out_dir <- "', out_dir_new, '"'),
  x = x,
  fixed = TRUE
)

# Add reproducibility patch after out_dir creation.
target_line <- 'dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)'

repro_patch <- paste0(
  target_line, "\n\n",
  'repro_dir <- file.path(out_dir, "reproducibility_logs")\n',
  'dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)\n\n',
  'write_checksum_csv <- function(files, out_file) {\n',
  '  files <- files[!is.na(files) & file.exists(files)]\n',
  '  if (length(files) == 0) {\n',
  '    write.csv(data.frame(file = character(), md5 = character()), out_file, row.names = FALSE)\n',
  '    return(invisible(NULL))\n',
  '  }\n',
  '  files_norm <- normalizePath(files, mustWork = TRUE)\n',
  '  md5 <- tools::md5sum(files_norm)\n',
  '  write.csv(data.frame(file = names(md5), md5 = unname(md5), stringsAsFactors = FALSE), out_file, row.names = FALSE)\n',
  '  invisible(md5)\n',
  '}\n\n',
  'sink(file.path(repro_dir, "sessionInfo_before_06_to_10.txt"))\n',
  'cat("Timestamp:\\n"); print(Sys.time())\n',
  'cat("\\nR version:\\n"); print(R.version.string)\n',
  'cat("\\n.libPaths():\\n"); print(.libPaths())\n',
  'cat("\\nSys.info():\\n"); print(Sys.info())\n',
  'cat("\\nsessionInfo():\\n"); print(sessionInfo())\n',
  'sink()\n\n',
  'write.csv(installed.packages()[, c("Package", "Version", "LibPath", "Built")], file.path(repro_dir, "installed_packages_before_06_to_10.csv"), row.names = FALSE)\n\n',
  'write_checksum_csv(input_rds, file.path(repro_dir, "input_object_ffpe_md5_before_06_to_10.csv"))\n'
)

if (grepl(target_line, x, fixed = TRUE)) {
  x <- sub(target_line, repro_patch, x, fixed = TRUE)
} else {
  warning("Could not find target line for reproducibility patch: ", target_line)
}

final_chunk <- '
```{r}
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
cat("Timestamp:\\n"); print(Sys.time())
cat("\\nR version:\\n"); print(R.version.string)
cat("\\n.libPaths():\\n"); print(.libPaths())
cat("\\nSys.info():\\n"); print(Sys.info())
cat("\\nsessionInfo():\\n"); print(sessionInfo())
sink()

write.csv(
  installed.packages()[, c("Package", "Version", "LibPath", "Built")],
  file.path(repro_dir, "installed_packages_after_06_to_10.csv"),
  row.names = FALSE
)

cat("\\nReproducibility checksums and sessionInfo saved in:\\n")
cat(repro_dir, "\\n")
```
'

x <- paste0(x, "\n\n", final_chunk)

patched_rmd <- file.path(SCRIPT_DIR, "source_history", "clustering_annotation_characterization_source.Rmd")
writeLines(x, patched_rmd)

log_msg("Patched Rmd written to:", patched_rmd)

patch_info <- data.frame(
  field = c("timestamp", "original_rmd", "patched_rmd", "input_rds_new", "out_dir_new", "N_THREADS"),
  value = c(as.character(Sys.time()), original_rmd, patched_rmd, input_rds_new, out_dir_new, as.character(N_THREADS)),
  stringsAsFactors = FALSE
)

write.csv(
  patch_info,
  file.path(REPRO_LOG_DIR, "patch_and_render_06_to_10_info.csv"),
  row.names = FALSE
)

write_checksum_csv(
  c(original_rmd, patched_rmd, input_rds_new),
  file.path(REPRO_LOG_DIR, "patch_and_render_06_to_10_input_md5.csv")
)

# ------------------------------------------------------------
# Render
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(rmarkdown)
})

log_msg("Rendering patched Rmd")

rendered <- rmarkdown::render(
  input = patched_rmd,
  output_dir = REPORTS_DIR,
  envir = new.env(parent = globalenv()),
  clean = TRUE
)

log_msg("Rendered report:", rendered)
log_msg("Main outputs:", out_dir_new)
log_msg("Done.")
