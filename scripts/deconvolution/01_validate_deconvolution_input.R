args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L || any(!startsWith(args[seq(1L, length(args), 2L)], "--"))) {
  stop("Arguments must be supplied as --key value pairs", call. = FALSE)
}
arg_list <- setNames(args[seq(2L, length(args), 2L)], sub("^--", "", args[seq(1L, length(args), 2L)]))
required <- c("fractions", "output-dir")
missing <- setdiff(required, names(arg_list))
if (length(missing) > 0L) stop("Missing argument(s): ", paste(missing, collapse = ", "), call. = FALSE)

fractions_path <- normalizePath(arg_list[["fractions"]], mustWork = TRUE)
output_dir <- normalizePath(arg_list[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_sha256 <- "aa2e9aef96fb51afeb5488a61ff374a03eb6ea9a96b724e96f04ce3cdd24815e"
expected_spots <- 26056L
expected_celltypes <- 12L
tolerance <- 1e-6

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L || !grepl("^[0-9a-f]{64}[[:space:]]", output)) {
    stop("Could not calculate SHA-256 with shasum", call. = FALSE)
  }
  sub("[[:space:]].*$", "", output)
}

lines <- readLines(fractions_path, warn = FALSE)
if (length(lines) != expected_spots + 1L) stop("Unexpected number of lines", call. = FALSE)
fields <- strsplit(lines, "	", fixed = TRUE)
field_counts <- lengths(fields)
if (!all(field_counts == expected_celltypes + 1L)) {
  stop("Rows do not contain one spot ID plus 12 fractions", call. = FALSE)
}

spot_ids <- vapply(fields[-1L], `[[`, character(1), 1L)
if (any(!nzchar(spot_ids)) || anyDuplicated(spot_ids)) {
  stop("Spot IDs are missing or duplicated", call. = FALSE)
}
fraction_values <- unlist(lapply(fields[-1L], function(x) as.numeric(x[-1L])), use.names = FALSE)
if (any(!is.finite(fraction_values))) stop("Fractions contain non-finite values", call. = FALSE)
if (any(fraction_values < -tolerance | fraction_values > 1 + tolerance)) {
  stop("Fractions fall outside [0, 1]", call. = FALSE)
}
fraction_matrix <- matrix(fraction_values, ncol = expected_celltypes, byrow = TRUE)
row_sums <- rowSums(fraction_matrix)
maximum_sum_error <- max(abs(row_sums - 1))
if (maximum_sum_error > tolerance) {
  stop("Fraction rows do not sum to one within tolerance", call. = FALSE)
}

actual_sha256 <- sha256(fractions_path)
if (!identical(actual_sha256, expected_sha256)) {
  stop("Fraction matrix SHA-256 differs from canonical input", call. = FALSE)
}

summary <- c(
  "FLINC deconvolution input validation",
  "Status: PASS",
  paste("Fraction matrix:", fractions_path),
  paste("SHA-256:", actual_sha256),
  paste("Spots:", length(spot_ids)),
  paste("Fraction columns:", expected_celltypes),
  paste("Maximum absolute row-sum error:", format(maximum_sum_error, scientific = TRUE)),
  paste("Tolerance:", tolerance),
  "Interpretation: STRIDE relative cell-type fractions; not raw counts."
)
writeLines(summary, file.path(output_dir, "01_deconvolution_input_validation.txt"))
cat("DECONVOLUTION_INPUT_VALIDATION_OK\n")
