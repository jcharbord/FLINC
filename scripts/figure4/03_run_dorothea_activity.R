#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 10, scipen = 999)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
set.seed(42L)

fail <- function(...) stop(paste0(...), call. = FALSE)

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "decoupleR", "dorothea",
  "digest", "dplyr", "ggplot2", "ggrepel"
)
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

required_args <- c("input", "output-dir", "log-dir", "expected-input-sha256")
missing_args <- setdiff(required_args, names(arg_list))
if (length(missing_args) > 0L) {
  fail("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

input_path <- normalizePath(arg_list[["input"]], mustWork = TRUE)
output_dir <- arg_list[["output-dir"]]
log_dir <- arg_list[["log-dir"]]
expected_input_sha <- tolower(arg_list[["expected-input-sha256"]])

if (!grepl("^[0-9a-f]{64}$", expected_input_sha)) {
  fail("--expected-input-sha256 must be a lowercase 64-character SHA256")
}
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
    message(
      "Run did not complete; incomplete directories were retained for inspection:\n",
      result_tmp, "\n", log_tmp
    )
  }
}, add = TRUE)

atomic_write_table <- function(x, path, row.names = FALSE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(
    x,
    temporary,
    sep = "\t",
    quote = FALSE,
    row.names = row.names,
    col.names = TRUE,
    na = ""
  )
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

atomic_write_lines <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  writeLines(x, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

atomic_save_rds <- function(x, path, compress = TRUE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, temporary, compress = compress)
  if (!file.rename(temporary, path)) fail("Could not atomically write ", path)
}

script_path <- {
  full_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", full_arguments, value = TRUE)
  if (length(file_argument) == 1L) {
    normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
  } else {
    NA_character_
  }
}
script_sha <- if (!is.na(script_path)) {
  digest::digest(script_path, algo = "sha256", file = TRUE)
} else {
  NA_character_
}

started_at <- Sys.time()

message("Computing SHA256 for the regenerated input object...")
input_sha <- digest::digest(input_path, algo = "sha256", file = TRUE)
if (!identical(tolower(input_sha), expected_input_sha)) {
  fail(
    "Input SHA256 mismatch. Expected ", expected_input_sha,
    "; observed ", input_sha
  )
}

message("Loading regenerated Seurat object read-only...")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) fail("Input is not a Seurat object")
if (!"Spatial" %in% names(object@assays)) {
  fail("Input does not contain the Spatial assay")
}

metadata <- object[[]]
required_metadata <- c("donor", "F_score")
if (!all(required_metadata %in% colnames(metadata))) {
  fail(
    "Input metadata is missing: ",
    paste(setdiff(required_metadata, colnames(metadata)), collapse = ", ")
  )
}

counts <- SeuratObject::LayerData(object[["Spatial"]], layer = "counts")
if (!inherits(counts, "sparseMatrix")) {
  counts <- methods::as(counts, "dgCMatrix")
}
if (!identical(colnames(counts), rownames(metadata))) {
  fail("Spatial count columns are not identical to metadata row names")
}
if (any(counts@x < 0) || any(!is.finite(counts@x))) {
  fail("Spatial counts contain negative or non-finite values")
}

input_feature_count <- nrow(counts)
input_spot_count <- ncol(counts)
gene_totals <- Matrix::rowSums(counts)
zero_count_genes <- names(gene_totals)[gene_totals == 0]
observed_genes <- names(gene_totals)[gene_totals > 0]
library_sizes <- Matrix::colSums(counts)
if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
  fail("At least one spot has a non-positive or non-finite library size")
}

spot_metadata <- data.frame(
  barcode = colnames(counts),
  donor = as.character(metadata[colnames(counts), "donor"]),
  F_score = as.character(metadata[colnames(counts), "F_score"]),
  stringsAsFactors = FALSE
)
if (anyNA(spot_metadata) || any(spot_metadata$donor == "") ||
    any(spot_metadata$F_score == "")) {
  fail("Donor or fibrosis-stage metadata contains missing/empty values")
}
if (!all(spot_metadata$F_score %in% paste0("F", 0:4))) {
  fail(
    "Unexpected fibrosis stages: ",
    paste(sort(unique(spot_metadata$F_score)), collapse = ", ")
  )
}

donor_stage <- unique(spot_metadata[c("donor", "F_score")])
if (anyDuplicated(donor_stage$donor)) {
  fail("At least one donor maps to more than one fibrosis stage")
}
donor_stage$stage_number <- match(donor_stage$F_score, paste0("F", 0:4)) - 1L
donor_stage <- donor_stage[order(donor_stage$stage_number, donor_stage$donor), ]
rownames(donor_stage) <- NULL

observed_donor_counts <- table(
  factor(donor_stage$F_score, levels = paste0("F", 0:4))
)
expected_donor_counts <- c(F0 = 6L, F1 = 6L, F2 = 3L, F3 = 3L, F4 = 6L)
if (!identical(as.integer(observed_donor_counts), as.integer(expected_donor_counts))) {
  fail(
    "Unexpected donor counts by stage. Observed: ",
    paste(names(observed_donor_counts), observed_donor_counts, sep = "=", collapse = ";")
  )
}

message("Loading and validating DoRothEA human A/B regulons...")
resource_environment <- new.env(parent = emptyenv())
utils::data("dorothea_hs", package = "dorothea", envir = resource_environment)
if (!exists("dorothea_hs", envir = resource_environment, inherits = FALSE)) {
  fail("Could not load dorothea_hs from the dorothea package")
}
dorothea_hs <- get("dorothea_hs", envir = resource_environment)

required_resource_columns <- c("tf", "confidence", "target", "mor")
if (!all(required_resource_columns %in% colnames(dorothea_hs))) {
  fail(
    "DoRothEA resource is missing columns: ",
    paste(setdiff(required_resource_columns, colnames(dorothea_hs)), collapse = ", ")
  )
}

regulon_ab_resource <- dorothea_hs[
  dorothea_hs$confidence %in% c("A", "B"),
  ,
  drop = FALSE
]
regulon_ab_resource <- dplyr::arrange(
  regulon_ab_resource,
  tf,
  target,
  confidence,
  mor
)
resource_sha <- digest::digest(
  regulon_ab_resource,
  algo = "sha256",
  serialize = TRUE
)
expected_resource_sha <-
  "7bf9a0895d3737354ab1d777c25165d57ba583b43e24daf0d0b0589c6bf73117"

resource_summary <- c(
  rows = nrow(regulon_ab_resource),
  TFs = length(unique(regulon_ab_resource$tf)),
  targets = length(unique(regulon_ab_resource$target)),
  confidence_A = sum(regulon_ab_resource$confidence == "A"),
  confidence_B = sum(regulon_ab_resource$confidence == "B")
)
expected_resource_summary <- c(
  rows = 6408L,
  TFs = 118L,
  targets = 3254L,
  confidence_A = 5434L,
  confidence_B = 974L
)
if (!identical(as.integer(resource_summary), as.integer(expected_resource_summary))) {
  fail(
    "Unexpected DoRothEA A/B resource structure: ",
    paste(names(resource_summary), resource_summary, sep = "=", collapse = ";")
  )
}
if (!identical(resource_sha, expected_resource_sha)) {
  fail(
    "DoRothEA A/B resource SHA256 mismatch. Expected ", expected_resource_sha,
    "; observed ", resource_sha
  )
}

network <- data.frame(
  source = as.character(regulon_ab_resource$tf),
  target = as.character(regulon_ab_resource$target),
  mor = as.numeric(regulon_ab_resource$mor),
  confidence = as.character(regulon_ab_resource$confidence),
  stringsAsFactors = FALSE
)
if (anyNA(network) || any(!is.finite(network$mor))) {
  fail("DoRothEA A/B network contains missing or non-finite values")
}
network <- unique(network)
network <- network[
  order(network$source, network$target, network$confidence, network$mor),
  ,
  drop = FALSE
]

network_observed <- network[network$target %in% observed_genes, , drop = FALSE]
coverage <- aggregate(
  network_observed$target,
  by = list(TF = network_observed$source),
  FUN = function(x) length(unique(x))
)
colnames(coverage)[2L] <- "observed_unique_targets"
all_tfs <- sort(unique(network$source))
coverage <- merge(
  data.frame(TF = all_tfs, stringsAsFactors = FALSE),
  coverage,
  by = "TF",
  all.x = TRUE,
  sort = TRUE
)
coverage$observed_unique_targets[is.na(coverage$observed_unique_targets)] <- 0L
coverage$total_unique_targets <- vapply(
  coverage$TF,
  function(tf_name) length(unique(network$target[network$source == tf_name])),
  integer(1)
)
coverage$coverage_fraction <-
  coverage$observed_unique_targets / coverage$total_unique_targets

minsize <- 5L
eligible_tfs <- coverage$TF[coverage$observed_unique_targets >= minsize]
if (!identical(eligible_tfs, all_tfs)) {
  fail(
    "Not all 118 A/B TFs meet minsize=5 after feature matching. Ineligible: ",
    paste(setdiff(all_tfs, eligible_tfs), collapse = ", ")
  )
}

network_scored <- network_observed[
  network_observed$source %in% eligible_tfs,
  c("source", "target", "mor"),
  drop = FALSE
]
scored_targets <- sort(unique(network_scored$target))
message("Normalizing ", length(scored_targets), " observed DoRothEA target genes...")

target_counts <- counts[scored_targets, , drop = FALSE]
rm(object, counts, metadata, dorothea_hs)
invisible(gc())

scale_factor <- 10000
normalization_multiplier <- scale_factor / library_sizes
target_log <- target_counts %*% Matrix::Diagonal(x = normalization_multiplier)
target_log@x <- log1p(target_log@x)
rownames(target_log) <- rownames(target_counts)
colnames(target_log) <- colnames(target_counts)

gene_means <- Matrix::rowMeans(target_log)
gene_mean_squares <- Matrix::rowMeans(target_log ^ 2)
gene_variances <- (input_spot_count / (input_spot_count - 1)) *
  pmax(0, gene_mean_squares - gene_means ^ 2)
gene_sds <- sqrt(gene_variances)
zero_variance_targets <- names(gene_sds)[!is.finite(gene_sds) | gene_sds <= 0]

if (length(zero_variance_targets) > 0L) {
  network_scored <- network_scored[
    !network_scored$target %in% zero_variance_targets,
    ,
    drop = FALSE
  ]
  scored_targets <- setdiff(scored_targets, zero_variance_targets)
  target_log <- target_log[scored_targets, , drop = FALSE]
  gene_means <- gene_means[scored_targets]
  gene_sds <- gene_sds[scored_targets]
}

post_variance_coverage <- aggregate(
  network_scored$target,
  by = list(TF = network_scored$source),
  FUN = function(x) length(unique(x))
)
colnames(post_variance_coverage)[2L] <- "scored_unique_targets"
coverage <- merge(
  coverage,
  post_variance_coverage,
  by = "TF",
  all.x = TRUE,
  sort = TRUE
)
coverage$scored_unique_targets[is.na(coverage$scored_unique_targets)] <- 0L
if (any(coverage$scored_unique_targets < minsize)) {
  fail("At least one TF has fewer than five targets after variance filtering")
}

target_z <- as.matrix(target_log)
target_z <- sweep(target_z, 1L, gene_means, FUN = "-")
target_z <- sweep(target_z, 1L, gene_sds, FUN = "/")
if (any(!is.finite(target_z))) fail("Target z-score matrix contains non-finite values")

rm(target_counts, target_log, gene_mean_squares, gene_variances)
invisible(gc())

times <- 1000L
seed <- 42L
spot_batch_size <- 512L
batch_equivalence_spots <- 12L
message(
  "Validating batched decoupleR::run_wmean equivalence at production settings..."
)
expected_statistics <- c("corr_wmean", "norm_wmean", "wmean")
run_wmean_batch <- function(batch_matrix) {
  batch_result <- decoupleR::run_wmean(
    mat = batch_matrix,
    network = network_scored,
    .source = "source",
    .target = "target",
    .mor = "mor",
    times = times,
    seed = seed,
    sparse = FALSE,
    randomize_type = "rows",
    minsize = minsize
  )
  observed_statistics <- sort(unique(as.character(batch_result$statistic)))
  if (!identical(observed_statistics, sort(expected_statistics))) {
    fail(
      "Unexpected run_wmean statistic set: ",
      paste(observed_statistics, collapse = ", ")
    )
  }
  batch_result[
    batch_result$statistic == "norm_wmean",
    c("source", "condition", "score", "p_value"),
    drop = FALSE
  ]
}

spots <- colnames(target_z)
tfs <- sort(unique(network_scored$source))
audit_spots <- spots[seq_len(min(batch_equivalence_spots, length(spots)))]
audit_split <- split(
  audit_spots,
  rep(seq_len(2L), each = ceiling(length(audit_spots) / 2L))[
    seq_along(audit_spots)
  ]
)
audit_full <- run_wmean_batch(target_z[, audit_spots, drop = FALSE])
audit_batched <- do.call(
  rbind,
  lapply(
    audit_split,
    function(batch_spots) {
      run_wmean_batch(target_z[, batch_spots, drop = FALSE])
    }
  )
)
audit_order_full <- order(audit_full$source, audit_full$condition)
audit_order_batched <- order(audit_batched$source, audit_batched$condition)
audit_full <- audit_full[audit_order_full, , drop = FALSE]
audit_batched <- audit_batched[audit_order_batched, , drop = FALSE]
rownames(audit_full) <- NULL
rownames(audit_batched) <- NULL
if (
  any(!is.finite(audit_full$score)) ||
  any(!is.finite(audit_batched$score)) ||
  any(!is.finite(audit_full$p_value)) ||
  any(!is.finite(audit_batched$p_value))
) {
  fail("Batched run_wmean equivalence audit returned non-finite values")
}
batch_equivalence_max_score_difference <- max(
  abs(audit_full$score - audit_batched$score)
)
batch_equivalence_max_p_difference <- max(
  abs(audit_full$p_value - audit_batched$p_value)
)
if (
  !identical(audit_full[c("source", "condition")],
             audit_batched[c("source", "condition")]) ||
  batch_equivalence_max_score_difference > 1e-12 ||
  batch_equivalence_max_p_difference > 1e-12
) {
  fail(
    "Batched run_wmean equivalence audit failed. Maximum score difference=",
    batch_equivalence_max_score_difference,
    "; maximum p-value difference=", batch_equivalence_max_p_difference
  )
}
message(
  "Batched equivalence audit passed: maximum norm_wmean score difference=",
  signif(batch_equivalence_max_score_difference, 6),
  "; maximum p-value difference=",
  signif(batch_equivalence_max_p_difference, 6)
)
rm(audit_full, audit_batched)
invisible(gc())

activity_matrix <- matrix(
  NA_real_,
  nrow = length(tfs),
  ncol = length(spots),
  dimnames = list(tfs, spots)
)

batch_starts <- seq.int(1L, length(spots), by = spot_batch_size)
message(
  "Running definitive decoupleR::run_wmean in ", length(batch_starts),
  " deterministic spot batches; times=", times, ", seed=", seed, "..."
)
for (batch_number in seq_along(batch_starts)) {
  batch_start <- batch_starts[[batch_number]]
  batch_end <- min(batch_start + spot_batch_size - 1L, length(spots))
  batch_spots <- spots[batch_start:batch_end]
  message(
    "DoRothEA batch ", batch_number, "/", length(batch_starts),
    " (spots ", batch_start, "-", batch_end, ")"
  )
  activity_long <- run_wmean_batch(
    target_z[, batch_spots, drop = FALSE]
  )
  expected_batch_rows <- length(tfs) * length(batch_spots)
  if (nrow(activity_long) != expected_batch_rows) {
    fail(
      "Unexpected norm_wmean row count in batch ", batch_number,
      ". Expected ", expected_batch_rows, "; observed ", nrow(activity_long)
    )
  }
  if (anyDuplicated(activity_long[c("source", "condition")])) {
    fail("Duplicate TF/spot rows in batch ", batch_number)
  }
  tf_index <- match(activity_long$source, tfs)
  spot_index <- match(activity_long$condition, spots)
  if (anyNA(tf_index) || anyNA(spot_index)) {
    fail("Could not map all norm_wmean rows in batch ", batch_number)
  }
  activity_matrix[cbind(tf_index, spot_index)] <- activity_long$score
  rm(activity_long)
  invisible(gc())
}
if (any(!is.finite(activity_matrix))) {
  fail("Definitive TF activity matrix contains missing or non-finite values")
}

rm(target_z)
invisible(gc())

message("Summarizing spot-level TF activity by donor median...")
donor_levels <- donor_stage$donor
donor_activity <- matrix(
  NA_real_,
  nrow = length(tfs),
  ncol = length(donor_levels),
  dimnames = list(tfs, donor_levels)
)
for (donor_name in donor_levels) {
  donor_spots <- spot_metadata$barcode[spot_metadata$donor == donor_name]
  donor_activity[, donor_name] <- apply(
    activity_matrix[, donor_spots, drop = FALSE],
    1L,
    stats::median,
    na.rm = TRUE
  )
}
if (any(!is.finite(donor_activity))) {
  fail("Donor-median TF activity matrix contains non-finite values")
}

donor_activity_long <- data.frame(
  TF = rep(rownames(donor_activity), times = ncol(donor_activity)),
  donor = rep(colnames(donor_activity), each = nrow(donor_activity)),
  activity_median = as.vector(donor_activity),
  stringsAsFactors = FALSE
)
donor_stage_match <- match(donor_activity_long$donor, donor_stage$donor)
donor_activity_long$F_score <- donor_stage$F_score[donor_stage_match]
donor_activity_long <- donor_activity_long[
  order(
    match(donor_activity_long$F_score, paste0("F", 0:4)),
    donor_activity_long$donor,
    donor_activity_long$TF
  ),
]
rownames(donor_activity_long) <- NULL

test_one_tf <- function(tf_name, stage_k) {
  f0_donors <- donor_stage$donor[donor_stage$F_score == "F0"]
  fk_donors <- donor_stage$donor[donor_stage$F_score == stage_k]
  y <- as.numeric(donor_activity[tf_name, f0_donors])
  x <- as.numeric(donor_activity[tf_name, fk_donors])
  ties_present <- anyDuplicated(c(x, y)) > 0L

  p_asymptotic <- suppressWarnings(
    stats::wilcox.test(
      x,
      y,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE
    )$p.value
  )

  if (!ties_present) {
    p_primary <- suppressWarnings(
      stats::wilcox.test(
        x,
        y,
        alternative = "two.sided",
        exact = TRUE,
        correct = FALSE
      )$p.value
    )
    test_method <- "two-sided exact Wilcoxon rank-sum"
  } else {
    p_primary <- p_asymptotic
    test_method <- "two-sided asymptotic Wilcoxon rank-sum; tie fallback"
  }

  data.frame(
    TF = tf_name,
    contrast = paste0(stage_k, "_vs_F0"),
    n_F0 = length(y),
    n_FK = length(x),
    mean_donor_median_F0 = mean(y),
    mean_donor_median_FK = mean(x),
    median_donor_median_F0 = stats::median(y),
    median_donor_median_FK = stats::median(x),
    delta_mean_donor_medians = mean(x) - mean(y),
    delta_median_donor_medians = stats::median(x) - stats::median(y),
    p_primary = p_primary,
    p_asymptotic = p_asymptotic,
    ties_present = ties_present,
    test_method = test_method,
    stringsAsFactors = FALSE
  )
}

contrast_stages <- paste0("F", 1:4)
statistics <- do.call(
  rbind,
  lapply(contrast_stages, function(stage_k) {
    do.call(rbind, lapply(tfs, test_one_tf, stage_k = stage_k))
  })
)
statistics$FDR_BH <- ave(
  statistics$p_primary,
  statistics$contrast,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)
statistics$figure_role <- ifelse(
  statistics$contrast == "F4_vs_F0",
  "main",
  ifelse(
    statistics$contrast == "F1_vs_F0",
    "supplementary",
    "validation_only"
  )
)
statistics <- statistics[
  order(
    match(statistics$contrast, paste0(contrast_stages, "_vs_F0")),
    statistics$FDR_BH,
    -abs(statistics$delta_mean_donor_medians),
    statistics$TF
  ),
]
rownames(statistics) <- NULL

make_volcano <- function(contrast_name, figure_title, output_stem) {
  plot_data <- statistics[statistics$contrast == contrast_name, , drop = FALSE]
  plot_data$significant <- is.finite(plot_data$FDR_BH) & plot_data$FDR_BH <= 0.05
  plot_data$label <- ifelse(plot_data$significant, plot_data$TF, "")
  plot_data$minus_log10_FDR <- -log10(pmax(plot_data$FDR_BH, .Machine$double.xmin))

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = delta_mean_donor_medians,
      y = minus_log10_FDR
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.45,
      colour = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "black"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(alpha = significant),
      colour = "black",
      size = 1.8
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      seed = seed,
      size = 2.7,
      box.padding = 0.30,
      point.padding = 0.18,
      min.segment.length = 0,
      max.overlaps = Inf,
      segment.alpha = 0.65,
      segment.size = 0.25
    ) +
    ggplot2::scale_alpha_manual(
      values = c(`TRUE` = 1.0, `FALSE` = 0.45),
      guide = "none"
    ) +
    ggplot2::labs(
      title = figure_title,
      x = "Delta TF activity (stage - F0) [mean donor medians]",
      y = expression(-log[10]("FDR"))
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(result_tmp, paste0(output_stem, ".pdf")),
    plot_object,
    width = 7.2,
    height = 6.2,
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    file.path(result_tmp, paste0(output_stem, ".png")),
    plot_object,
    width = 7.2,
    height = 6.2,
    dpi = 300,
    bg = "white"
  )
}

message("Rendering definitive DoRothEA volcanoes...")
make_volcano(
  "F4_vs_F0",
  "F4 vs F0",
  "Fig4C_DoRothEA_activity_F4_vs_F0"
)
make_volcano(
  "F1_vs_F0",
  "F1 vs F0",
  "SupFig9B_DoRothEA_activity_F1_vs_F0"
)
make_volcano(
  "F2_vs_F0",
  "F2 vs F0 - validation only",
  "Validation_DoRothEA_activity_F2_vs_F0"
)
make_volcano(
  "F3_vs_F0",
  "F3 vs F0 - validation only",
  "Validation_DoRothEA_activity_F3_vs_F0"
)

activity_output <- list(
  activity = activity_matrix,
  spot_metadata = spot_metadata,
  donor_activity = donor_activity,
  donor_stage = donor_stage,
  regulon_resource_sha256 = resource_sha,
  selected_statistic = "norm_wmean",
  times = times,
  seed = seed
)
atomic_save_rds(
  activity_output,
  file.path(result_tmp, "05_dorothea_norm_wmean_activity.rds"),
  compress = TRUE
)
atomic_write_table(
  donor_activity_long,
  file.path(result_tmp, "05_dorothea_donor_medians.tsv")
)
atomic_write_table(
  statistics,
  file.path(result_tmp, "05_dorothea_all_contrasts.tsv")
)
atomic_write_table(
  donor_stage,
  file.path(result_tmp, "05_donor_stage_map.tsv")
)
atomic_write_table(
  coverage,
  file.path(result_tmp, "05_dorothea_model_coverage.tsv")
)
atomic_write_table(
  data.frame(
    gene = zero_count_genes,
    reason = rep("zero_total_count_in_regenerated_input", length(zero_count_genes)),
    stringsAsFactors = FALSE
  ),
  file.path(result_tmp, "05_removed_zero_count_genes.tsv")
)
atomic_write_table(
  data.frame(
    gene = zero_variance_targets,
    reason = rep("zero_variance_after_log_normalization", length(zero_variance_targets)),
    stringsAsFactors = FALSE
  ),
  file.path(result_tmp, "05_removed_zero_variance_targets.tsv")
)

parameters <- data.frame(
  parameter = c(
    "input_object", "input_sha256", "assay", "input_layer",
    "zero_count_gene_rule", "normalization", "scale_factor",
    "gene_standardization", "regulon_database", "regulon_confidence",
    "regulon_resource_sha256", "decoupler_method", "selected_statistic",
    "minsize", "times", "seed", "randomize_type", "sparse_argument",
    "spot_batch_size", "batch_equivalence_spots",
    "batch_equivalence_tolerance", "batch_equivalence_max_score_difference",
    "batch_equivalence_max_p_value_difference",
    "donor_summary", "primary_test", "tie_fallback", "multiple_testing",
    "effect_size", "main_contrast", "supplementary_contrast",
    "validation_only_contrasts", "blas_threads", "volcano_label_rule"
  ),
  value = c(
    input_path, input_sha, "Spatial", "counts",
    "remove genes with rowSums(counts) == 0",
    "log1p(count / spot_library_size * scale_factor)", as.character(scale_factor),
    "per-gene z-score across all 26056 spots; exclude zero-variance targets",
    "dorothea::dorothea_hs", "A;B", resource_sha,
    "decoupleR::run_wmean", "norm_wmean", as.character(minsize),
    as.character(times), as.character(seed), "rows", "FALSE",
    as.character(spot_batch_size), as.character(length(audit_spots)), "1e-12",
    as.character(batch_equivalence_max_score_difference),
    as.character(batch_equivalence_max_p_difference),
    "median across spots per TF and donor",
    "two-sided exact Wilcoxon rank-sum when no ties",
    "two-sided asymptotic Wilcoxon rank-sum with continuity correction",
    "Benjamini-Hochberg within each 118-TF contrast",
    "mean(stage donor medians) - mean(F0 donor medians)",
    "F4_vs_F0", "F1_vs_F0", "F2_vs_F0;F3_vs_F0", "1",
    "label all TFs with FDR_BH <= 0.05"
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(parameters, file.path(result_tmp, "05_run_parameters.tsv"))

significance_summary <- do.call(
  rbind,
  lapply(paste0(contrast_stages, "_vs_F0"), function(contrast_name) {
    contrast_data <- statistics[
      statistics$contrast == contrast_name,
      ,
      drop = FALSE
    ]
    significant_names <- sort(contrast_data$TF[contrast_data$FDR_BH <= 0.05])
    data.frame(
      contrast = contrast_name,
      figure_role = unique(contrast_data$figure_role),
      TFs_tested = nrow(contrast_data),
      significant_TFs_FDR_le_0.05 = length(significant_names),
      significant_TF_names = paste(significant_names, collapse = ";"),
      RFX5_delta_mean_donor_medians = contrast_data$delta_mean_donor_medians[
        contrast_data$TF == "RFX5"
      ],
      RFX5_FDR_BH = contrast_data$FDR_BH[contrast_data$TF == "RFX5"],
      stringsAsFactors = FALSE
    )
  })
)
atomic_write_table(
  significance_summary,
  file.path(log_tmp, "05_significance_summary.tsv")
)
atomic_write_table(
  coverage,
  file.path(log_tmp, "05_model_coverage.tsv")
)

package_versions <- data.frame(
  package = c("R", required_packages),
  version = c(
    R.version.string,
    vapply(
      required_packages,
      function(package_name) as.character(utils::packageVersion(package_name)),
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(
  package_versions,
  file.path(log_tmp, "05_package_versions.tsv")
)

finished_at <- Sys.time()
execution_summary <- c(
  "FLINC Figure 4 DoRothEA execution summary",
  paste("Started UTC:", format(started_at, tz = "UTC", usetz = TRUE)),
  paste("Completed UTC:", format(finished_at, tz = "UTC", usetz = TRUE)),
  paste(
    "Elapsed seconds:",
    as.numeric(difftime(finished_at, started_at, units = "secs"))
  ),
  paste("Input object:", input_path),
  paste("Input SHA256:", input_sha),
  paste("Script SHA256:", script_sha),
  paste("Input spots:", input_spot_count),
  paste("Input Spatial features before zero-count filtering:", input_feature_count),
  paste("Removed zero-count genes:", length(zero_count_genes)),
  paste("Observed DoRothEA model targets:", length(unique(network_observed$target))),
  paste("Scored DoRothEA targets:", length(scored_targets)),
  paste("Removed zero-variance model targets:", length(zero_variance_targets)),
  paste("Scored TFs:", length(tfs)),
  paste("Selected decoupleR statistic:", "norm_wmean"),
  paste("Deterministic spot batch size:", spot_batch_size),
  paste("Number of definitive spot batches:", length(batch_starts)),
  paste("Batch-equivalence audit spots:", length(audit_spots)),
  paste(
    "Batch-equivalence maximum norm_wmean score difference:",
    batch_equivalence_max_score_difference
  ),
  paste(
    "Batch-equivalence maximum p-value difference:",
    batch_equivalence_max_p_difference
  ),
  paste("DoRothEA A/B resource serialized SHA256:", resource_sha),
  paste("Primary main contrast:", "F4_vs_F0"),
  paste("Supplementary contrast:", "F1_vs_F0"),
  paste("Validation-only contrasts:", "F2_vs_F0;F3_vs_F0"),
  paste("RFX5 main-contrast FDR:", signif(
    statistics$FDR_BH[
      statistics$contrast == "F4_vs_F0" & statistics$TF == "RFX5"
    ],
    8
  )),
  paste("Command:", paste(commandArgs(), collapse = " "))
)
atomic_write_lines(
  execution_summary,
  file.path(log_tmp, "05_execution_summary.txt")
)
atomic_write_lines(
  capture.output(sessionInfo()),
  file.path(log_tmp, "05_sessionInfo.txt")
)

result_files <- sort(list.files(result_tmp, full.names = TRUE, recursive = FALSE))
manifest <- data.frame(
  file = basename(result_files),
  bytes = file.info(result_files)$size,
  sha256 = vapply(
    result_files,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
atomic_write_table(manifest, file.path(result_tmp, "05_output_SHA256.tsv"))
atomic_write_lines("SUCCESS", file.path(result_tmp, "05_COMPLETED.txt"))
atomic_write_lines("SUCCESS", file.path(log_tmp, "05_COMPLETED.txt"))

if (!file.rename(result_tmp, output_dir)) {
  fail("Could not atomically publish results to ", output_dir)
}
if (!file.rename(log_tmp, log_dir)) {
  fail("Could not atomically publish logs to ", log_dir)
}

completed <- TRUE
message("DoRothEA analysis completed successfully.")
message("Results: ", output_dir)
message("Compact logs: ", log_dir)
