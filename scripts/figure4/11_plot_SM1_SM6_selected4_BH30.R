#!/usr/bin/env Rscript
# Script 11; renumbered from local script 24 (2026-10-24), derived from script 23 v2.
# Display-only change: q < 0.001 / q < 0.01; exact BH30 values remain in tables.
# FLINC Fig4 H/I: donor-level Spearman profiles, BH across 30 genes PER module.
# No WGCNA rerun. Reuses the validated top-30 lists and Spatial counts.
# Both plots: original stage colours, 2x2 panels, names-only and rho/q variants.
# Finds the original v6 input/topgenes in its run_manifest.tsv. Optional overrides:
# FLINC_SECRETED_MANIFEST, FLINC_SM6_SOURCE_DIR (exact original source directories).
# Run with Rscript --vanilla in hdWGCNA_pipeline-env; FLINC_REPRO_ROOT as usual.

main <- function() {
  options(stringsAsFactors = FALSE, digits = 10)
  needed <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "digest")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  library(ggplot2)
  root <- Sys.getenv("FLINC_REPRO_ROOT",
    "/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code")
  parent <- file.path(root, "results/Fig4_regulatory_and_module_networks")
  read_tsv <- function(p) read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  fail <- function(...) stop(..., call. = FALSE)
  sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  read_manifest <- function(p) {
    z <- tryCatch(read.delim(p, header = FALSE, sep = "\t", quote = "",
      fill = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(z) || ncol(z) < 2L) return(NULL)
    setNames(as.character(z[[2L]]), z[[1L]])
  }
  message("1/5 Finding the validated secreted-factor inputs...")
  override <- Sys.getenv("FLINC_SECRETED_MANIFEST", "")
  if (nzchar(override)) {
    manifests <- override
  } else {
    dirs <- c(file.path(root, "reproducibility_logs"), parent)
    manifests <- unique(unlist(lapply(dirs[dir.exists(dirs)], list.files,
      pattern = "^run_manifest[.]tsv$", recursive = TRUE, full.names = TRUE)))
    manifests <- manifests[!grepl("incomplete|partial_failed", manifests)]
  }
  candidates <- lapply(manifests, function(p) {
    x <- read_manifest(p)
    req <- c("input_path", "input_sha256", "topgenes_path", "topgenes_sha256",
      "stage_test", "main_heatmap_summary")
    if (is.null(x) || !all(req %in% names(x))) return(NULL)
    if (!grepl("four mapped genes", x[["stage_test"]], fixed = TRUE)) return(NULL)
    if (!all(file.exists(x[c("input_path", "topgenes_path")]))) return(NULL)
    list(path = p, parameters = x)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) fail("Cannot locate the original secreted-factor run_manifest.tsv. ",
    "Set FLINC_SECRETED_MANIFEST to its exact path; do not rerun the analysis pipeline.")
  keys <- vapply(candidates, function(x) paste(x$parameters[c("input_sha256",
    "topgenes_sha256")], collapse = ":"), character(1))
  if (length(unique(keys)) != 1L) fail("Multiple different input versions found. Set ",
    "FLINC_SECRETED_MANIFEST to the manifest for the final v6 figure. Candidates:\n",
    paste(vapply(candidates, `[[`, character(1), "path"), collapse = "\n"))
  manifest <- candidates[[1L]]$parameters
  input <- manifest[["input_path"]]
  topfile <- manifest[["topgenes_path"]]
  if (!identical(tolower(sha(input)), tolower(manifest[["input_sha256"]]))) fail("Input object checksum mismatch")
  if (!identical(tolower(sha(topfile)), tolower(manifest[["topgenes_sha256"]]))) fail("Top-30 table checksum mismatch")
  # v6 may point to a 16-gene selection audit, not the original full top-30
  # export. Validate that audit, then recover a complete list consistent with it.
  auditfile <- topfile
  normalize_top <- function(z) {
    if (!all(c("gene", "module", "rank") %in% names(z))) return(NULL)
    z$gene <- toupper(trimws(as.character(z$gene)))
    z$module <- gsub("-", "_", trimws(as.character(z$module)))
    z$rank <- suppressWarnings(as.numeric(z$rank))
    z <- z[z$module %in% c("SM_1", "SM_6") & !is.na(z$rank) & z$rank >= 1 & z$rank <= 30, ]
    z <- z[order(z$module, z$rank), , drop = FALSE]
    rownames(z) <- NULL
    z
  }
  audit <- normalize_top(read_tsv(auditfile))
  if (is.null(audit) || !nrow(audit)) fail("Unexpected original top-gene audit format")
  complete_top <- function(z) {
    !is.null(z) && nrow(z) == 60L && !anyDuplicated(z$gene) &&
      setequal(z$module, c("SM_1", "SM_6")) &&
      all(vapply(c("SM_1", "SM_6"), function(m)
        identical(sort(z$rank[z$module == m]), as.numeric(1:30)), logical(1)))
  }
  agrees_with_audit <- function(z) {
    key <- function(x) paste(x$module, x$gene, sep = "::")
    idx <- match(key(audit), key(z))
    if (anyNA(idx) || !identical(audit$rank, z$rank[idx])) return(FALSE)
    if ("mean_normalized_expression" %in% names(audit)) {
      if (!("mean_normalized_expression" %in% names(z))) return(FALSE)
      if (!isTRUE(all.equal(as.numeric(audit$mean_normalized_expression),
        as.numeric(z$mean_normalized_expression[idx]), tolerance = 1e-8))) return(FALSE)
    }
    TRUE
  }
  top <- audit
  if (!complete_top(top)) {
    message("Original file is a selected-gene audit; locating the complete top-30 export...")
    requested <- Sys.getenv("FLINC_SECRETED_TOP30", "")
    paths <- if (nzchar(requested)) requested else list.files(parent,
      pattern = "[.](tsv|csv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    paths <- paths[!grepl("incomplete|partial_failed", paths)]
    # Only small summary tables are candidates, never expression matrices.
    paths <- paths[file.exists(paths) & !is.na(file.info(paths)$size) & file.info(paths)$size < 5e6]
    found <- list()
    for (p in paths) {
      sep <- if (grepl("[.]csv$", p, ignore.case = TRUE)) "," else "\t"
      header <- tryCatch(read.table(p, header = TRUE, sep = sep, nrows = 0L,
        check.names = FALSE, comment.char = ""), error = function(e) NULL)
      if (is.null(header) || !all(c("gene", "module", "rank") %in% names(header))) next
      z <- tryCatch(normalize_top(read.table(p, header = TRUE, sep = sep,
        check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)), error = function(e) NULL)
      if (complete_top(z) && agrees_with_audit(z)) found[[p]] <- z
    }
    if (!length(found)) fail("No complete top-30 table matches the recorded gene ranks and expression values. ",
      "Please provide the original full top-30 export (gene, module, rank, mean_normalized_expression). ",
      "Original selected-gene audit: ", auditfile)
    signatures <- vapply(found, function(z) digest::digest(z[c("gene", "module", "rank")]), character(1))
    if (length(unique(signatures)) != 1L) fail("Different complete lists match the audit. ",
      "Set FLINC_SECRETED_TOP30 to the correct full table:\n", paste(names(found), collapse = "\n"))
    topfile <- names(found)[1L]
    top <- found[[1L]]
  }
  message("Complete top-30 source: ", topfile)
  selected <- list(SM_1 = c("TIMP1", "IGFBP7", "ADAMTSL2", "DCN"),
    SM_6 = c("AHSG", "FN1", "ITIH2", "SERPINC1"))
  for (m in names(selected)) if (!all(selected[[m]] %in% top$gene[top$module == m]))
    fail("Selected genes absent from top-30 table for ", m)

  sm6dir <- Sys.getenv("FLINC_SM6_SOURCE_DIR", "")
  if (!nzchar(sm6dir)) {
    dirs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[grepl("^20_SM6_top30_stage_profiles_", basename(dirs)) &
      !grepl("incomplete", dirs)]
    dirs <- dirs[file.exists(file.path(dirs, "SM6_top30_donor_stage_expression.tsv")) &
      file.exists(file.path(dirs, "SM6_top30_stage_spearman_BH30.tsv"))]
    if (length(dirs) != 1L) fail("Set FLINC_SM6_SOURCE_DIR to the original 20_SM6_top30_stage_profiles_* folder.")
    sm6dir <- dirs[[1L]]
  }
  old6 <- read_tsv(file.path(sm6dir, "SM6_top30_donor_stage_expression.tsv"))
  old6stats <- read_tsv(file.path(sm6dir, "SM6_top30_stage_spearman_BH30.tsv"))

  message("2/5 Reading saved Spatial counts and calculating donor means for 60 genes...")
  object <- readRDS(input)
  if (!inherits(object, "Seurat")) fail("Input is not a Seurat object")
  counts <- SeuratObject::LayerData(object[["Spatial"]], layer = "counts")
  md <- object[[]]
  if (!all(c("donor", "F_score") %in% names(md))) fail("Missing donor or F_score metadata")
  if (anyDuplicated(rownames(md)) || !setequal(rownames(md), colnames(counts))) fail("Spot metadata mismatch")
  md <- md[colnames(counts), , drop = FALSE]
  genes <- toupper(rownames(counts))
  if (anyDuplicated(genes) || !all(top$gene %in% genes)) fail("Missing or duplicated gene identifiers")
  rownames(counts) <- genes
  stage <- as.character(md$F_score)
  stage[stage == "Healthy"] <- "F0"
  if (all(stage %in% as.character(0:4))) stage <- paste0("F", stage)
  donor <- as.character(md$donor)
  if (anyNA(donor) || anyNA(stage) || any(!nzchar(donor))) fail("Missing donor/stage")
  donors <- unique(data.frame(donor = donor, F_score = stage))
  stages <- paste0("F", 0:4)
  if (anyDuplicated(donors$donor) || !setequal(stage, stages) ||
      !identical(as.integer(table(factor(donors$F_score, levels = stages))),
        c(6L, 6L, 3L, 3L, 6L))) fail("Expected 24 donors with stage counts 6/6/3/3/6")
  donors <- donors[order(match(donors$F_score, stages), donors$donor), ]
  totals <- Matrix::colSums(counts)
  if (any(!is.finite(totals)) || any(totals <= 0)) fail("Invalid library sizes")
  norm <- counts[unique(top$gene), , drop = FALSE] %*% Matrix::Diagonal(x = 10000 / totals)
  norm <- log1p(norm)
  dimnames(norm) <- list(unique(top$gene), colnames(counts))
  expr <- do.call(rbind, lapply(seq_len(nrow(donors)), function(i) {
    idx <- which(donor == donors$donor[i])
    data.frame(donor = donors$donor[i], F_score = donors$F_score[i],
      gene = rownames(norm), mean_normalized_expression = as.numeric(Matrix::rowMeans(norm[, idx, drop = FALSE])))
  }))
  expr <- merge(expr, top[c("gene", "module")], by = "gene", sort = FALSE)
  rm(object, counts, norm); invisible(gc())

  check_expression <- function(reference, current, label) {
    cols <- c("gene", "donor", "F_score", "mean_normalized_expression")
    if (!all(cols %in% names(reference)) || anyDuplicated(reference[c("gene", "donor")]))
      fail("Invalid reference table: ", label)
    key <- function(x) paste(x$gene, x$donor, sep = "::")
    idx <- match(key(reference), key(current))
    if (anyNA(idx) || !identical(as.character(reference$F_score), as.character(current$F_score[idx])) ||
      !isTRUE(all.equal(as.numeric(reference$mean_normalized_expression),
        as.numeric(current$mean_normalized_expression[idx]), tolerance = 1e-8)))
      fail("Donor expression disagrees with original export: ", label)
    data.frame(check = label, n_values = nrow(reference), result = "PASS")
  }
  checks <- list(check_expression(old6, expr[expr$module == "SM_6", ], "SM6 donor values"))
  # Verify SM1 against matching v6 exports wherever present; do not assume that
  # every old/backup version belongs to the same input object.
  refs <- list.files(parent, pattern = "^selected_secreted_factor_donor_stage_statistics[.]tsv$",
    recursive = TRUE, full.names = TRUE)
  refs <- refs[!grepl("incomplete|partial_failed", refs)]
  matching <- list()
  for (p in refs) {
    z <- read_tsv(p)
    if (all(c("gene", "donor", "F_score", "mean_normalized_expression") %in% names(z))) {
      z <- z[z$gene %in% selected$SM_1, ]
      if (nrow(z) == 96L) {
        ok <- tryCatch(check_expression(z, expr[expr$module == "SM_1", ], p), error = function(e) NULL)
        if (!is.null(ok)) matching[[length(matching) + 1L]] <- ok
      }
    }
  }
  if (!length(matching)) fail("Could not verify the four SM1 donor profiles against a previous export. ",
    "Keep the existing figures and inspect the input paths; no final outputs were written.")
  checks <- c(checks, matching)
  message("3/5 Spearman tests and BH correction separately across 30 genes per module...")
  stats <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    z <- expr[expr$gene == top$gene[i] & expr$module == top$module[i], ]
    if (length(unique(z$mean_normalized_expression)) < 2L) {
      rho <- p <- NA_real_
    } else {
      tt <- suppressWarnings(stats::cor.test(match(z$F_score, stages) - 1L,
        z$mean_normalized_expression, method = "spearman", exact = FALSE))
      rho <- unname(tt$estimate); p <- tt$p.value
    }
    data.frame(module = top$module[i], gene = top$gene[i], spearman_rho = rho, p_value = p)
  }))
  stats$q_bh_30_genes <- NA_real_
  for (m in names(selected)) {
    idx <- which(stats$module == m)
    stats$q_bh_30_genes[idx] <- p.adjust(stats$p_value[idx], method = "BH", n = 30L)
  }
  s6 <- stats[stats$module == "SM_6", ]
  if (!all(c("gene", "spearman_rho", "p_value", "q_bh_30_genes") %in% names(old6stats)) ||
      nrow(old6stats) != 30L || anyDuplicated(old6stats$gene) || !setequal(old6stats$gene, s6$gene))
    fail("Invalid original SM6 statistics table")
  s6 <- s6[match(old6stats$gene, s6$gene), ]
  for (col in c("spearman_rho", "p_value", "q_bh_30_genes"))
    if (!isTRUE(all.equal(as.numeric(old6stats[[col]]), as.numeric(s6[[col]]), tolerance = 1e-8)))
      fail("SM6 statistics changed unexpectedly: ", col)
  shown <- do.call(rbind, lapply(names(selected), function(m)
    stats[match(selected[[m]], stats$gene), ]))
  shown$q_bh_4_reference <- NA_real_
  for (m in names(selected)) {
    idx <- which(shown$module == m)
    shown$q_bh_4_reference[idx] <- p.adjust(shown$p_value[idx], method = "BH", n = 4L)
  }
  shown$significant_BH30 <- !is.na(shown$q_bh_30_genes) & shown$q_bh_30_genes < 0.05

  message("4/5 Exporting matched 2x2 panels into a new folder...")
  out <- file.path(parent, paste0("11_SM1_SM6_selected4_BH30_q_thresholds_",
    format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_", Sys.getpid()))
  tmp <- paste0(out, ".incomplete")
  if (dir.exists(out) || !dir.create(tmp, recursive = TRUE)) fail("Cannot create output folder")
  palette <- c(F0 = "#BDBDBD", F1 = "#F8C8D8", F2 = "#F06292", F3 = "#FF7043", F4 = "#C62828")
  q_label <- function(q) {
    if (is.na(q)) return("q = NA")
    if (q < 0.001) return("q < 0.001")
    if (q < 0.01) return("q < 0.01")
    paste0("q = ", format(signif(q, 3), trim = TRUE))
  }
  for (m in names(selected)) {
    chosen <- selected[[m]]
    d <- expr[expr$module == m & expr$gene %in% chosen, ]
    d <- d[order(match(d$gene, chosen), match(d$F_score, stages), d$donor), ]
    d$gene <- factor(d$gene, levels = chosen)
    d$F_score <- factor(d$F_score, levels = stages)
    set.seed(42)
    d$plot_x <- as.integer(d$F_score) + runif(nrow(d), -0.10, 0.10)
    s <- shown[match(chosen, shown$gene), ]
    labs_stats <- setNames(paste0(chosen, "\nSpearman ρ = ", sprintf("%.2f", s$spearman_rho),
      "; ", vapply(s$q_bh_30_genes, q_label, character(1))), chosen)
    for (variant in c("names_only", "Spearman_BH30")) {
      labels <- if (variant == "names_only") setNames(chosen, chosen) else labs_stats
      p <- ggplot(d, aes(x = as.integer(F_score), y = mean_normalized_expression)) +
        geom_boxplot(aes(group = F_score, fill = F_score), width = 0.58,
          outlier.shape = NA, colour = "#444444", alpha = 0.88) +
        geom_point(aes(x = plot_x, fill = F_score), shape = 21, size = 2,
          stroke = 0.35, alpha = 0.88, colour = "#444444") +
        scale_fill_manual(values = palette, guide = "none", drop = FALSE) +
        scale_x_continuous(breaks = 1:5, labels = c("Healthy", "F1", "F2", "F3", "F4"), limits = c(0.5, 5.5)) +
        facet_wrap(~gene, scales = "free_y", ncol = 2L, labeller = as_labeller(labels)) +
        labs(x = NULL, y = "Mean log-normalized expression") + theme_classic(base_size = 12) +
        theme(strip.background = element_rect(fill = "#EEEEEE", colour = NA),
          strip.text = element_text(colour = "black", size = 12), axis.text = element_text(colour = "black"),
          panel.spacing = grid::unit(1.2, "lines"), plot.background = element_rect(fill = "white", colour = NA))
      stem <- paste0(gsub("_", "", m), "_selected4_donor_stage_profiles_", variant)
      ggsave(file.path(tmp, paste0(stem, ".pdf")), p, width = 9, height = 8,
        units = "in", device = grDevices::cairo_pdf, bg = "white")
      ggsave(file.path(tmp, paste0(stem, ".png")), p, width = 9, height = 8, units = "in", dpi = 300, bg = "white")
    }
  }
  write_tsv <- function(x, name) write.table(x, file.path(tmp, name), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(expr, "SM1_SM6_top30_donor_expression.tsv")
  write_tsv(stats, "SM1_SM6_top30_spearman_BH30.tsv")
  write_tsv(shown, "selected4_statistics_BH30_and_BH4_comparison.tsv")
  write_tsv(top, "top30_gene_lists_used.tsv")
  write_tsv(do.call(rbind, checks), "expression_concordance_checks.tsv")
  write_tsv(donors, "donor_stage_map.tsv")
  source_files <- unique(c(input, auditfile, topfile))
  write_tsv(data.frame(file = source_files, sha256 = vapply(source_files, sha, character(1))), "input_checksums.tsv")
  writeLines(c(
    "Methods: For each of the 30 highest-expressed secreted-factor genes in SM_1 and SM_6,",
    "log-normalized expression was averaged across all spots within each donor, including zeros.",
    "Two-sided Spearman correlations with ordinal stage (Healthy=0 through F4=4) were computed",
    "using exact=FALSE. Benjamini-Hochberg correction was applied separately across 30 genes",
    "within each module. The resulting q values were retained for the four displayed factors.",
    "No stage-versus-Healthy gene tests or network reconstruction were performed.",
    "BH4 values are exported for comparison only; all new plots display BH30 values.",
    "Script 11 changes label formatting only: q < 0.001 or q < 0.01 when applicable.",
    "Exact numeric q values are retained unchanged in the exported statistics tables.",
    "Check the exported statistics before retaining q-value thresholds in manuscript Results.",
    paste("Source manifest:", candidates[[1L]]$path), paste("SM6 reference:", sm6dir)
  ), file.path(tmp, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script) == 1L && !file.copy(sub("^--file=", "", script), tmp)) fail("Could not archive source script")
  if (!file.rename(tmp, out)) fail("Could not finalize output; inspect ", tmp)
  message("5/5 Done. SM6 statistics verified unchanged. Review SM1 BH30 q values below.")
  print(shown, row.names = FALSE)
  message("Saved figures, tables and source script: ", out)
  if (Sys.info()[["sysname"]] == "Darwin") system2("open", shQuote(out), wait = FALSE)
}
main()
