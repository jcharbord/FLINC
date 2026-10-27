#!/usr/bin/env Rscript

# FLINC Figure 2 supplementary stagewise GSEA.
#
# Independent biological unit: donor.
# Expression input: raw counts from the Spatial assay of the final annotated
# Seurat object. SCT values and SCT models are intentionally not used.
# Differential-expression model: edgeR filtering/TMM followed by limma-voom
# and robust empirical Bayes moderation.
# GSEA input: the complete expression-filtered gene list ranked by the limma
# moderated t-statistic; no DEG significance threshold is applied.

options(stringsAsFactors = FALSE, warn = 1)

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) {
  stop(
    paste(
      "Usage: 08_run_supplementary_stage_gsea.R",
      "[input_rds] [results_root] [reproducibility_log_root]"
    ),
    call. = FALSE
  )
}

repro_root <- normalizePath(
  path.expand(env_or(
    "FLINC_REPRO_ROOT",
    "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
  )),
  mustWork = TRUE
)

input_rds <- if (length(args) >= 1L) {
  args[[1]]
} else {
  env_or(
    "FIG2_GSEA_INPUT_OBJECT",
    file.path(
      repro_root,
      "results",
      "260701_approach1_cluster_go_marker_deconv",
      "object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds"
    )
  )
}
results_root <- if (length(args) >= 2L) {
  args[[2]]
} else {
  env_or(
    "FIG2_GSEA_RESULTS_DIR",
    file.path(repro_root, "results", "Fig2_GSEA", "08_stagewise_gsea")
  )
}
logs_root <- if (length(args) >= 3L) {
  args[[3]]
} else {
  env_or(
    "FIG2_GSEA_LOG_DIR",
    file.path(
      repro_root,
      "reproducibility_logs",
      "Fig2_GSEA",
      "08_stagewise_gsea"
    )
  )
}

input_rds <- normalizePath(path.expand(input_rds), mustWork = TRUE)
results_root <- path.expand(results_root)
logs_root <- path.expand(logs_root)
if (identical(normalizePath(dirname(results_root), mustWork = FALSE),
              normalizePath(dirname(logs_root), mustWork = FALSE)) &&
    identical(basename(results_root), basename(logs_root))) {
  stop("Results and reproducibility-log roots must be different.", call. = FALSE)
}

seed_base <- 260814L
stage_levels <- c("F0", "F1", "F2", "F3", "F4")
comparisons <- c(F1_vs_F0 = "F1", F2_vs_F0 = "F2", F3_vs_F0 = "F3", F4_vs_F0 = "F4")
collections <- c("H", "C2", "C3", "C5")

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "edgeR", "limma", "statmod",
  "msigdbr", "clusterProfiler", "fgsea", "AnnotationDbi",
  "org.Hs.eg.db", "ggplot2", "stringr", "digest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE)
}

rng_kind <- RNGkind()
set.seed(seed_base)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_root, recursive = TRUE, showWarnings = FALSE)
results_root <- normalizePath(results_root, mustWork = TRUE)
logs_root <- normalizePath(logs_root, mustWork = TRUE)

run_stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
run_name <- paste0("run_", run_stamp, "Z")
temporary_results_dir <- file.path(
  results_root,
  paste0(".", run_name, ".partial.", Sys.getpid())
)
temporary_logs_dir <- file.path(
  logs_root,
  paste0(".", run_name, ".partial.", Sys.getpid())
)
final_results_dir <- file.path(results_root, run_name)
final_logs_dir <- file.path(logs_root, run_name)

output_candidates <- c(
  temporary_results_dir,
  temporary_logs_dir,
  final_results_dir,
  final_logs_dir
)
if (any(file.exists(output_candidates))) {
  stop(
    "Refusing to overwrite an existing final or partial run: ",
    paste(output_candidates[file.exists(output_candidates)], collapse = "; "),
    call. = FALSE
  )
}

dir.create(temporary_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temporary_logs_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(temporary_results_dir) || !dir.exists(temporary_logs_dir)) {
  stop("Could not create temporary output directories.", call. = FALSE)
}

result_subdirs <- c(
  "de", "rankings", "gsea_full", "gsea_significant",
  "figures", "qc", "manifests", "objects"
)
invisible(vapply(
  file.path(temporary_results_dir, result_subdirs),
  dir.create,
  logical(1),
  recursive = TRUE,
  showWarnings = FALSE
))

completed <- FALSE
on.exit({
  if (!completed) {
    message(
      "Operation 08 did not complete. Partial files, if any, remain in: ",
      temporary_results_dir,
      " and ",
      temporary_logs_dir
    )
  }
}, add = TRUE)

message("Loading input object: ", input_rds)
obj <- readRDS(input_rds)

if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)
if (!"Spatial" %in% names(obj@assays)) stop("Spatial assay is missing.", call. = FALSE)

meta <- obj[[]]
required_metadata <- c("donor", "F_score")
missing_metadata <- setdiff(required_metadata, colnames(meta))
if (length(missing_metadata)) {
  stop("Missing metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}

counts <- tryCatch(
  SeuratObject::LayerData(obj[["Spatial"]], layer = "counts"),
  error = function(e) Seurat::GetAssayData(obj, assay = "Spatial", layer = "counts")
)

if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
if (!identical(colnames(counts), rownames(meta))) {
  if (!all(colnames(counts) %in% rownames(meta))) {
    stop("Spatial count columns do not match metadata row names.", call. = FALSE)
  }
  meta <- meta[colnames(counts), , drop = FALSE]
}

if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
  stop("Spatial counts contain non-finite or negative values.", call. = FALSE)
}
if (any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("Spatial counts are not integer-valued raw counts.", call. = FALSE)
}

n_sct_models <- if (
  "SCT" %in% names(obj@assays) && inherits(obj[["SCT"]], "SCTAssay")
) {
  length(obj[["SCT"]]@SCTModel.list)
} else {
  NA_integer_
}
input_object_summary <- data.frame(
  item = c(
    "input_path",
    "seurat_object_version",
    "default_assay",
    "assays_present",
    "spatial_assay_class",
    "spatial_features",
    "spots",
    "required_metadata",
    "sct_models_present_not_used",
    "expression_input"
  ),
  value = c(
    input_rds,
    as.character(obj@version),
    SeuratObject::DefaultAssay(obj),
    paste(names(obj@assays), collapse = ","),
    paste(class(obj[["Spatial"]]), collapse = ","),
    nrow(counts),
    ncol(counts),
    paste(required_metadata, collapse = ","),
    n_sct_models,
    "raw Spatial counts; SCT values and SCT models are not used"
  ),
  stringsAsFactors = FALSE
)
write.table(
  input_object_summary,
  file.path(temporary_logs_dir, "input_object_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

donor <- trimws(as.character(meta$donor))
stage <- trimws(as.character(meta$F_score))
if (anyNA(donor) || any(!nzchar(donor))) stop("Missing donor labels.", call. = FALSE)
if (anyNA(stage) || any(!stage %in% stage_levels)) {
  stop("Unexpected or missing F_score labels.", call. = FALSE)
}

donor_stage <- unique(data.frame(donor = donor, stage = stage))
donor_stage_counts <- table(donor_stage$donor)
if (any(donor_stage_counts != 1L)) {
  stop("At least one donor maps to more than one fibrosis stage.", call. = FALSE)
}

donor_stage$stage <- factor(donor_stage$stage, levels = stage_levels)
donor_stage <- donor_stage[order(donor_stage$stage, donor_stage$donor), , drop = FALSE]
donor_order <- donor_stage$donor
observed_replicates <- table(donor_stage$stage)
expected_replicates <- c(F0 = 6L, F1 = 6L, F2 = 3L, F3 = 3L, F4 = 6L)
if (!identical(as.integer(observed_replicates), as.integer(expected_replicates))) {
  stop(
    "Unexpected donor replication: ",
    paste(names(observed_replicates), observed_replicates, collapse = ", "),
    call. = FALSE
  )
}
if (any(observed_replicates < 2L)) {
  stop("At least one fibrosis stage has fewer than two donors.", call. = FALSE)
}

message("Aggregating raw Spatial counts to donor-level pseudobulk profiles.")
donor_factor <- factor(donor, levels = donor_order)
membership <- Matrix::sparseMatrix(
  i = seq_along(donor_factor),
  j = as.integer(donor_factor),
  x = 1,
  dims = c(length(donor_factor), length(donor_order)),
  dimnames = list(colnames(counts), donor_order)
)

pseudobulk_counts <- counts %*% membership
colnames(pseudobulk_counts) <- donor_order
rownames(pseudobulk_counts) <- rownames(counts)

sample_metadata <- donor_stage
sample_metadata$n_spots <- as.integer(table(donor_factor)[donor_order])
sample_metadata$library_size <- as.numeric(Matrix::colSums(pseudobulk_counts))
sample_metadata$stage <- factor(sample_metadata$stage, levels = stage_levels)
rownames(sample_metadata) <- sample_metadata$donor

write.table(
  transform(sample_metadata, stage = as.character(stage)),
  file.path(temporary_logs_dir, "sample_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
saveRDS(
  pseudobulk_counts,
  file.path(temporary_results_dir, "objects", "donor_pseudobulk_raw_counts.rds"),
  compress = "xz"
)

rm(obj, meta, counts, membership)
invisible(gc())

message("Filtering, TMM-normalizing and fitting limma-voom model.")
dge <- edgeR::DGEList(counts = pseudobulk_counts, samples = sample_metadata)
design <- stats::model.matrix(~ 0 + stage, data = sample_metadata)
colnames(design) <- sub("^stage", "", colnames(design))
rownames(design) <- sample_metadata$donor

keep <- edgeR::filterByExpr(dge, design = design)
if (sum(keep) < 1000L) stop("Fewer than 1,000 genes passed filterByExpr.", call. = FALSE)

filter_manifest <- data.frame(
  gene_id = rownames(dge),
  retained = keep,
  total_count = Matrix::rowSums(dge$counts),
  stringsAsFactors = FALSE
)
write_tsv_gz(
  filter_manifest,
  file.path(
    temporary_results_dir,
    "manifests",
    "gene_filter_manifest.tsv.gz"
  )
)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")

grDevices::pdf(
  file.path(temporary_results_dir, "qc", "voom_mean_variance.pdf"),
  width = 7,
  height = 6
)
voom_object <- limma::voom(dge, design = design, plot = TRUE)
grDevices::dev.off()

stage_colors <- c(
  F0 = "#7F7F7F", F1 = "#1B9E77", F2 = "#7570B3",
  F3 = "#E6AB02", F4 = "#D95F02"
)
grDevices::pdf(
  file.path(temporary_results_dir, "qc", "pseudobulk_MDS.pdf"),
  width = 8,
  height = 7
)
limma::plotMDS(
  voom_object$E,
  labels = sample_metadata$donor,
  col = unname(stage_colors[as.character(sample_metadata$stage)]),
  main = "Donor-level pseudobulk MDS"
)
graphics::legend(
  "topright", legend = stage_levels, col = stage_colors[stage_levels],
  pch = 16, bty = "n", title = "Fibrosis stage"
)
grDevices::dev.off()

fit <- limma::lmFit(voom_object, design)
contrast_matrix <- limma::makeContrasts(
  F1_vs_F0 = F1 - F0,
  F2_vs_F0 = F2 - F0,
  F3_vs_F0 = F3 - F0,
  F4_vs_F0 = F4 - F0,
  levels = design
)
fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
fit_contrasts <- limma::eBayes(fit_contrasts, robust = TRUE)
saveRDS(
  fit_contrasts,
  file.path(temporary_results_dir, "objects", "limma_voom_contrast_fit.rds"),
  compress = "xz"
)

message("Loading and freezing MSigDB collections.")
msigdb <- tryCatch(
  msigdbr::msigdbr(db_species = "HS", species = "Homo sapiens"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens")
)

collection_column <- if ("gs_collection" %in% colnames(msigdb)) {
  "gs_collection"
} else if ("gs_cat" %in% colnames(msigdb)) {
  "gs_cat"
} else {
  stop("Could not identify the MSigDB collection column.", call. = FALSE)
}

if (!all(c("gs_name", "gene_symbol") %in% colnames(msigdb))) {
  stop("MSigDB data lack gs_name or gene_symbol.", call. = FALSE)
}

msigdb <- msigdb[msigdb[[collection_column]] %in% collections, , drop = FALSE]
term2gene_all <- unique(data.frame(
  collection = msigdb[[collection_column]],
  term = msigdb$gs_name,
  gene = msigdb$gene_symbol,
  stringsAsFactors = FALSE
))
term2gene_all <- term2gene_all[order(term2gene_all$collection, term2gene_all$term, term2gene_all$gene), ]
write_tsv_gz(
  term2gene_all,
  file.path(
    temporary_results_dir,
    "manifests",
    "msigdb_term2gene_H_C2_C3_C5.tsv.gz"
  )
)

msigdb_version <- if ("db_version" %in% colnames(msigdb)) {
  paste(sort(unique(msigdb$db_version)), collapse = ",")
} else {
  "not_exposed_by_msigdbr"
}

collection_manifest <- do.call(
  rbind,
  lapply(collections, function(coll) {
    x <- term2gene_all[term2gene_all$collection == coll, ]
    data.frame(
      collection = coll,
      n_gene_sets = length(unique(x$term)),
      n_term_gene_pairs = nrow(x),
      stringsAsFactors = FALSE
    )
  })
)
collection_manifest$msigdbr_package_version <- as.character(packageVersion("msigdbr"))
collection_manifest$msigdb_version <- msigdb_version
write.table(
  collection_manifest,
  file.path(
    temporary_results_dir,
    "manifests",
    "msigdb_collection_manifest.tsv"
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)

clean_gene_ids <- function(x) {
  trimws(sub("\\..*$", "", as.character(x)))
}

build_symbol_ranking <- function(gene_ids, statistics) {
  cleaned <- clean_gene_ids(gene_ids)
  ensembl_fraction <- mean(grepl("^ENSG", cleaned))

  if (ensembl_fraction > 0.5) {
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(cleaned),
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    )
    mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL), , drop = FALSE]
    source <- data.frame(ENSEMBL = cleaned, statistic = statistics, gene_id = gene_ids)
    mapped <- merge(source, mapping, by = "ENSEMBL", all.x = FALSE, sort = FALSE)
  } else {
    mapped <- data.frame(
      gene_id = gene_ids,
      SYMBOL = cleaned,
      statistic = statistics,
      stringsAsFactors = FALSE
    )
  }

  mapped <- mapped[
    is.finite(mapped$statistic) & !is.na(mapped$SYMBOL) & nzchar(mapped$SYMBOL),
    , drop = FALSE
  ]
  mapped <- mapped[order(-abs(mapped$statistic), mapped$SYMBOL, mapped$gene_id), , drop = FALSE]
  mapped <- mapped[!duplicated(mapped$SYMBOL), , drop = FALSE]

  ranking <- mapped$statistic
  names(ranking) <- mapped$SYMBOL
  ranking <- sort(ranking, decreasing = TRUE, method = "radix")
  list(ranking = ranking, mapping = mapped, ensembl_fraction = ensembl_fraction)
}

gsea_full_tables <- list()
gsea_significant_tables <- list()
significance_counts <- list()
ranking_manifests <- list()
hallmark_significant_tables <- list()

for (comparison_index in seq_along(comparisons)) {
  comparison_name <- names(comparisons)[[comparison_index]]
  stage_name <- unname(comparisons[[comparison_index]])
  message("DE and GSEA: ", comparison_name)

  de <- limma::topTable(
    fit_contrasts,
    coef = comparison_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  de$gene_id <- rownames(de)
  de$comparison <- comparison_name
  de$stage_1 <- stage_name
  de$stage_2 <- "F0"
  leading_columns <- c("gene_id", "comparison", "stage_1", "stage_2")
  de <- de[, c(leading_columns, setdiff(colnames(de), leading_columns))]
  write_tsv_gz(
    de,
    file.path(
      temporary_results_dir,
      "de",
      paste0(comparison_name, "_limma_voom_DE.tsv.gz")
    )
  )

  ranking_result <- build_symbol_ranking(de$gene_id, de$t)
  gene_list <- ranking_result$ranking
  if (length(gene_list) < 1000L) stop("Too few mapped genes for ", comparison_name, call. = FALSE)
  if (anyDuplicated(names(gene_list))) stop("Duplicated symbols remain in ranking.", call. = FALSE)

  ranking_table <- data.frame(
    gene_symbol = names(gene_list),
    moderated_t = unname(gene_list),
    rank = seq_along(gene_list),
    stringsAsFactors = FALSE
  )
  write_tsv_gz(
    ranking_table,
    file.path(
      temporary_results_dir,
      "rankings",
      paste0(comparison_name, "_moderated_t_ranking.tsv.gz")
    )
  )
  write_tsv_gz(
    ranking_result$mapping,
    file.path(
      temporary_results_dir,
      "rankings",
      paste0(comparison_name, "_gene_id_to_symbol_mapping.tsv.gz")
    )
  )
  ranking_manifests[[comparison_name]] <- data.frame(
    comparison = comparison_name,
    expression_filtered_genes = nrow(de),
    ranked_unique_symbols = length(gene_list),
    ensembl_identifier_fraction = ranking_result$ensembl_fraction,
    ranking_statistic = "limma moderated t",
    DEG_significance_filter_applied = FALSE,
    stringsAsFactors = FALSE
  )

  for (collection_index in seq_along(collections)) {
    collection_name <- collections[[collection_index]]
    term2gene <- term2gene_all[term2gene_all$collection == collection_name, c("term", "gene")]
    colnames(term2gene) <- c("term", "gene")

    set.seed(seed_base + comparison_index * 100L + collection_index)
    gsea_object <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      by = "fgsea",
      exponent = 1,
      nPermSimple = 20000,
      minGSSize = 15,
      maxGSSize = 5000,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      eps = 0,
      seed = TRUE,
      verbose = FALSE
    )

    gsea_table <- as.data.frame(gsea_object@result)
    gsea_table$comparison <- rep(comparison_name, nrow(gsea_table))
    gsea_table$stage_1 <- rep(stage_name, nrow(gsea_table))
    gsea_table$stage_2 <- rep("F0", nrow(gsea_table))
    gsea_table$collection <- rep(collection_name, nrow(gsea_table))
    gsea_table$direction <- ifelse(
      gsea_table$NES > 0,
      "Activated",
      "Suppressed"
    )
    if (nrow(gsea_table)) {
      gsea_table <- gsea_table[
        order(gsea_table$p.adjust, gsea_table$pvalue, gsea_table$ID),
        ,
        drop = FALSE
      ]
    }

    key <- paste(comparison_name, collection_name, sep = "__")
    gsea_full_tables[[key]] <- gsea_table
    significant <- gsea_table[!is.na(gsea_table$p.adjust) & gsea_table$p.adjust < 0.05, , drop = FALSE]
    gsea_significant_tables[[key]] <- significant

    write_tsv_gz(
      gsea_table,
      file.path(
        temporary_results_dir,
        "gsea_full",
        paste0(comparison_name, "_", collection_name, "_GSEA_full.tsv.gz")
      )
    )
    write_tsv_gz(
      significant,
      file.path(
        temporary_results_dir,
        "gsea_significant",
        paste0(
          comparison_name,
          "_",
          collection_name,
          "_GSEA_FDR_lt_0.05.tsv.gz"
        )
      )
    )

    significance_counts[[key]] <- data.frame(
      comparison = comparison_name,
      collection = collection_name,
      tested_gene_sets = nrow(gsea_table),
      significant_FDR_lt_0.05 = nrow(significant),
      activated = sum(significant$NES > 0),
      suppressed = sum(significant$NES < 0),
      stringsAsFactors = FALSE
    )

    if (collection_name == "H") {
      hallmark_significant_tables[[comparison_name]] <- significant
    }
  }
}

ranking_manifest <- do.call(rbind, ranking_manifests)
rownames(ranking_manifest) <- NULL
write.table(
  ranking_manifest,
  file.path(
    temporary_results_dir,
    "manifests",
    "ranking_manifest.tsv"
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)

all_significant_hallmarks <- do.call(rbind, hallmark_significant_tables)
if (nrow(all_significant_hallmarks)) {
  hallmark_size_limits <- range(all_significant_hallmarks$setSize, na.rm = TRUE)
  hallmark_nes_limit <- ceiling(
    max(abs(all_significant_hallmarks$NES), na.rm = TRUE) * 10
  ) / 10
} else {
  hallmark_size_limits <- c(15, 500)
  hallmark_nes_limit <- 1
}

message("Generating four Hallmark plots from all significant Hallmark sets.")
for (comparison_name in names(comparisons)) {
  stage_name <- unname(comparisons[[comparison_name]])
  plot_data <- hallmark_significant_tables[[comparison_name]]
  plot_path <- file.path(
    temporary_results_dir,
    "figures",
    paste0(comparison_name, "_H_dotplot.pdf")
  )

  if (nrow(plot_data)) {
    plot_data$pathway <- stringr::str_to_title(
      stringr::str_replace_all(
        stringr::str_remove(plot_data$Description, "^HALLMARK_"),
        "_",
        " "
      )
    )
    plot_data$negative_log10_FDR <- pmin(
      -log10(pmax(plot_data$p.adjust, .Machine$double.xmin)),
      20
    )
    plot_data <- plot_data[order(plot_data$NES, plot_data$pathway), , drop = FALSE]
    plot_data$pathway <- factor(
      plot_data$pathway,
      levels = unique(plot_data$pathway)
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = NES,
        y = pathway,
        size = setSize,
        color = negative_log10_FDR
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::scale_x_continuous(
        limits = c(-hallmark_nes_limit, hallmark_nes_limit)
      ) +
      ggplot2::scale_size_continuous(
        name = "Gene set size",
        limits = hallmark_size_limits,
        range = c(2.5, 7)
      ) +
      ggplot2::scale_color_gradient(
        name = "-log10(FDR)\n(capped at 20)",
        low = "#D9EAF7",
        high = "#B2182B",
        limits = c(0, 20)
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title.position = "plot"
      ) +
      ggplot2::labs(
        title = paste0("Donor-pseudobulk Hallmark GSEA: ", stage_name, " vs F0"),
        subtitle = paste(
          "All Hallmark gene sets with BH FDR < 0.05;",
          "positive NES = activated, negative NES = suppressed"
        ),
        x = "Normalized enrichment score (NES)",
        y = NULL
      )
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No Hallmark gene sets at BH FDR < 0.05",
        size = 5
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = paste0("Donor-pseudobulk Hallmark GSEA: ", stage_name, " vs F0")
      )
  }

  plot_height <- max(7, 2.8 + 0.28 * nrow(plot_data))
  ggplot2::ggsave(
    plot_path,
    p,
    width = 9,
    height = plot_height,
    units = "in",
    limitsize = FALSE
  )
}

gsea_full_summary <- do.call(rbind, gsea_full_tables)
gsea_significant_summary <- do.call(rbind, gsea_significant_tables)
significance_summary <- do.call(rbind, significance_counts)
rownames(gsea_full_summary) <- NULL
rownames(gsea_significant_summary) <- NULL
rownames(significance_summary) <- NULL

write_tsv_gz(
  gsea_full_summary,
  file.path(
    temporary_results_dir,
    "gsea_full",
    "all_comparisons_all_collections_GSEA_full.tsv.gz"
  )
)
write_tsv_gz(
  gsea_significant_summary,
  file.path(
    temporary_results_dir,
    "gsea_significant",
    "all_comparisons_all_collections_GSEA_FDR_lt_0.05.tsv.gz"
  )
)
write.table(
  significance_summary,
  file.path(temporary_results_dir, "GSEA_significance_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

package_versions <- rbind(
  data.frame(
    package = c("R", "R_platform"),
    version = c(R.version.string, R.version$platform),
    stringsAsFactors = FALSE
  ),
  data.frame(
    package = required_packages,
    version = vapply(
      required_packages,
      function(x) as.character(packageVersion(x)),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
)
write.table(
  package_versions,
  file.path(temporary_logs_dir, "package_versions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  NA_character_
}

sha256_file <- function(path) {
  digest::digest(object = path, algo = "sha256", serialize = FALSE, file = TRUE)
}

input_manifest <- data.frame(
  item = c("input_rds", "analysis_script"),
  path = c(input_rds, script_path),
  sha256 = c(
    sha256_file(input_rds),
    if (!is.na(script_path)) sha256_file(script_path) else NA_character_
  ),
  stringsAsFactors = FALSE
)
write.table(
  input_manifest,
  file.path(temporary_logs_dir, "input_checksums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

run_parameters <- c(
  paste0("run_name=", run_name),
  paste0("seed_base=", seed_base),
  paste0("RNGkind=", paste(rng_kind, collapse = ",")),
  "computational_threads=1",
  paste0("comparisons=", paste(names(comparisons), collapse = ",")),
  paste0("collections=", paste(collections, collapse = ",")),
  "independent_biological_unit=donor",
  "de_method=donor_pseudobulk_limma_voom",
  "ranking_statistic=limma_moderated_t",
  "ranking_gene_universe=all_expression_filtered_genes",
  "DEG_significance_filter_before_GSEA=none",
  "raw_count_assay=Spatial",
  "SCT_values_or_models_used=FALSE",
  "filter=edgeR_filterByExpr",
  "normalization=edgeR_TMM",
  "empirical_Bayes=robust",
  "gsea_engine=clusterProfiler_GSEA_by_fgsea",
  "minGSSize=15",
  "maxGSSize=5000",
  "nPermSimple=20000",
  "gsea_exponent=1",
  "clusterProfiler_GSEA_nproc_argument=not_overridden",
  "gsea_FDR_threshold=0.05",
  "multiple_testing_adjustment=BH_within_comparison_and_collection",
  "Hallmark_plot_selection=all_significant_Hallmark_gene_sets",
  "Hallmark_plot_common_NES_size_and_FDR_scales=TRUE",
  paste0("msigdb_version=", msigdb_version)
)
writeLines(run_parameters, file.path(temporary_logs_dir, "run_parameters.txt"))
writeLines(
  capture.output(sessionInfo()),
  file.path(temporary_logs_dir, "sessionInfo.txt")
)
writeLines(
  c(
    "FLINC Figure 2 supplementary stagewise GSEA.",
    "The independent biological unit is the donor.",
    "Raw Spatial counts are aggregated by donor before modelling.",
    "GSEA uses every expression-filtered gene ranked by the moderated t-statistic.",
    "No individual-gene significance cutoff is applied before GSEA.",
    "F2 and F3 each contain three independent donors and should be interpreted accordingly."
  ),
  file.path(temporary_logs_dir, "analysis_scope.txt")
)

completion_lines <- c(
  paste0("completed_utc=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("run_name=", run_name),
  "status=complete"
)
writeLines(
  completion_lines,
  file.path(temporary_results_dir, "RESULTS_COMPLETE.txt")
)
writeLines(completion_lines, file.path(temporary_logs_dir, "COMPLETED.txt"))

result_files <- list.files(
  temporary_results_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)
result_files <- result_files[
  file.info(result_files)$isdir %in% FALSE
]
result_relative_paths <- substring(
  result_files,
  nchar(temporary_results_dir) + 2L
)
result_checksums <- data.frame(
  sha256 = vapply(result_files, sha256_file, character(1)),
  path = result_relative_paths,
  stringsAsFactors = FALSE
)
result_checksums <- result_checksums[order(result_checksums$path), ]
write.table(
  result_checksums,
  file.path(temporary_logs_dir, "result_SHA256SUMS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

execution_record_files <- list.files(
  temporary_logs_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)
execution_record_files <- execution_record_files[
  file.info(execution_record_files)$isdir %in% FALSE
]
execution_record_relative_paths <- substring(
  execution_record_files,
  nchar(temporary_logs_dir) + 2L
)
execution_record_checksums <- data.frame(
  sha256 = vapply(execution_record_files, sha256_file, character(1)),
  path = execution_record_relative_paths,
  stringsAsFactors = FALSE
)
execution_record_checksums <- execution_record_checksums[
  order(execution_record_checksums$path),
]
write.table(
  execution_record_checksums,
  file.path(temporary_logs_dir, "execution_record_SHA256SUMS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (!file.rename(temporary_results_dir, final_results_dir)) {
  stop("Could not atomically publish the completed results directory.", call. = FALSE)
}
if (!file.rename(temporary_logs_dir, final_logs_dir)) {
  stop(
    "Results were published, but execution records could not be published: ",
    final_results_dir,
    call. = FALSE
  )
}

completed <- TRUE
message("Completed donor-pseudobulk GSEA results: ", final_results_dir)
message("Execution records: ", final_logs_dir)
message("Significance summary:")
print(significance_summary, row.names = FALSE)
