#!/usr/bin/env Rscript
# FLINC used hdWGCNA 0.4.06 with WGCNA 1.73. In this specific saved
# consTomDS object, the off-diagonal entries are similarities (see README).
options(stringsAsFactors = FALSE, digits = 15)

# R's dist layout stores the lower triangle by column. Use double arithmetic.
tom_entry <- function(x, i, j) {
  n <- as.double(attr(x, "Size"))
  stopifnot(!anyNA(i), !anyNA(j), all(i != j),
            all(i >= 1 & i <= n), all(j >= 1 & j <= n))
  lower <- as.double(pmin(i, j))
  upper <- as.double(pmax(i, j))
  position <- n * (lower - 1) - lower * (lower - 1) / 2 + upper - lower
  as.numeric(x[position])
}

rank_neighbors <- function(index, module, genes, tom) {
  i <- match(index, genes$gene_name)
  stopifnot(!is.na(i), genes$module[i] == module)
  j <- which(genes$module == module & genes$gene_name != index)
  similarity <- tom_entry(tom, i, j)
  stopifnot(length(j) >= 20L, all(is.finite(similarity)),
            all(similarity >= 0 & similarity <= 1))
  # Descending similarity is equivalent to ascending (1 - similarity).
  ord <- order(-similarity, genes$gene_name[j], method = "radix")
  data.frame(index_gene = index, module = module, rank = seq_along(ord),
             neighbor = genes$gene_name[j[ord]], TOM_similarity = similarity[ord],
             TOM_dissimilarity = 1 - similarity[ord])
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: Rscript 12_rfx5_sm1_neighbours.R [test_TOM.rda] [new_output_folder]")
  arg_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  script_dir <- dirname(normalizePath(arg_file, mustWork = TRUE))
  project <- dirname(dirname(script_dir))
  default_tom <- file.path(project, "objects/Fig2_WGCNA/02_network_workspace/TOM/test_TOM.rda")
  default_out <- file.path(project, "results/Fig4_regulatory_and_module_networks/261102_RFX5_SM1")
  tom_path <- normalizePath(if (length(args)) args[1] else default_tom, mustWork = TRUE)
  out <- path.expand(if (length(args) == 2L) args[2] else default_out)
  if (file.exists(out)) stop("Output folder already exists. Supply a new output folder as the second argument: ", out)
  source(file.path(script_dir, "12_plot_rfx5_sm1.R"), local = TRUE)
  fig2 <- file.path(project, "resources/figure2/results")
  fig4 <- file.path(project, "resources/figure4")
  read_input <- function(path) read.delim(path, check.names = FALSE,
                                          stringsAsFactors = FALSE, quote = "", comment.char = "")
  genes <- read_input(file.path(fig2, "02_module_assignments_and_kME.tsv"))
  genes <- genes[, c("gene_name", "module")]
  pairs <- read_input(file.path(fig2, "11_ego_top20_neighbors.tsv"))
  # The historical column labelled TOM_dissimilarity contains the raw saved entries.
  pairs$saved_TOM_value <- pairs$TOM_dissimilarity
  targets <- read_input(file.path(fig4, "results/07_rfx5_all_observed_target_effects.tsv"))
  annotation <- read_input(file.path(fig4, "reference_inputs/SM1_function.tsv.gz"))
  flags <- c("hla_gene", "mhc_ii", "ecm", "actin")
  stopifnot(nrow(genes) == 11515L, !anyDuplicated(genes$gene_name), !anyDuplicated(annotation$gene),
            nrow(targets) == 14L, !anyDuplicated(targets$gene),
            all(c("gene", "gene_name_full", flags) %in% names(annotation)),
            nrow(annotation) == 1652L, setequal(annotation$gene, genes$gene_name[genes$module == "SM_1"]))
  for (field in flags) {
    flag <- annotation[[field]]
    stopifnot(all(is.na(flag) | toupper(as.character(flag)) %in% c("TRUE","FALSE")))
    annotation[[field]] <- as.logical(toupper(as.character(flag)))
  }

  message("Checking the saved FLINC TOM against 40 archived gene pairs...")
  holder <- new.env(parent = emptyenv())
  load(tom_path, envir = holder)
  tom <- holder$consTomDS
  n <- nrow(genes)
  stopifnot(inherits(tom, "dist"), attr(tom,"Size") == n, length(tom) == as.double(n)*(n-1)/2)
  lab <- attr(tom, "Labels")
  if (!is.null(lab) && !identical(as.character(lab), genes$gene_name) &&
      !identical(as.character(lab), as.character(seq_len(n)))) stop("TOM labels do not match the FLINC gene order.")
  measured <- tom_entry(tom, match(pairs$index_gene, genes$gene_name),
                        match(pairs$neighbor, genes$gene_name))
  stopifnot(length(measured) == 40L, all(abs(measured - pairs$saved_TOM_value) < 1e-12))

  # Complete the ranking before joining any functional or regulatory annotation.
  rfx <- rank_neighbors("RFX5", "SM_1", genes, tom)
  stopifnot(nrow(rfx) == 1651L)
  hit <- match(rfx$neighbor, targets$gene)
  rfx$RFX5_DoRothEA_mode <- targets$mor_sign[hit]
  rfx$RFX5_DoRothEA_confidence <- targets$confidence[hit]
  top <- head(rfx, 20)
  top <- cbind(top, annotation[match(top$neighbor, annotation$gene), setdiff(names(annotation),"gene"), drop=FALSE])
  top <- top[, c("rank","neighbor","TOM_similarity","gene_name_full", flags,
                  "RFX5_DoRothEA_mode","RFX5_DoRothEA_confidence",
                  setdiff(names(top),c("rank","neighbor","TOM_similarity","gene_name_full",flags,
                                       "RFX5_DoRothEA_mode","RFX5_DoRothEA_confidence")))]
  lookup <- match(targets$gene, rfx$neighbor)
  target_ranks <- data.frame(target = targets$gene, module = genes$module[match(targets$gene,genes$gene_name)],
                             DoRothEA_mode = targets$mor_sign, confidence = targets$confidence,
                             SM1_neighbor_rank = rfx$rank[lookup], TOM_similarity = rfx$TOM_similarity[lookup],
                             in_top20 = rfx$rank[lookup] <= 20L)
  stopifnot(sum(!is.na(lookup)) == 9L)
  target_ranks <- target_ranks[order(target_ranks$SM1_neighbor_rank, target_ranks$target, na.last=TRUE), ]
  counts <- data.frame(category = c(flags, "FLINC_RFX5_DoRothEA_target"),
                       matches_in_top20 = c(colSums(top[flags], na.rm=TRUE), sum(!is.na(top$RFX5_DoRothEA_mode))),
                       unknown_in_top20 = c(colSums(is.na(top[flags])), 0L), total_neighbors = nrow(top))

  parent <- dirname(out)
  if (!dir.exists(parent) && !dir.create(parent, recursive=TRUE)) stop("Cannot create: ", parent)
  tmp <- tempfile(".RFX5-", tmpdir = parent)
  if (!dir.create(tmp)) stop("Cannot create temporary output directory.")
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  save_tsv <- function(x, name) write.table(x, file.path(tmp, name), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  plotted <- annotate_neighbours(top)
  top$figure_category <- plotted$group
  save_tsv(top, "RFX5_top20_annotated.tsv")
  save_tsv(counts, "RFX5_top20_counts.tsv")
  save_tsv(target_ranks, "RFX5_DoRothEA_target_ranks.tsv")
  save_tsv(rfx, "RFX5_all_SM1_neighbor_ranks.tsv")
  draw_table(plotted, file.path(tmp, "RFX5_top20_functions_hubs.pdf"))
  draw_network(plotted, file.path(tmp, "RFX5_top20_network_by_function.pdf"), tom_scale(plotted$TOM_similarity))
  writeLines(c(paste("TOM filename:", basename(tom_path)), paste("TOM MD5:", tools::md5sum(tom_path)),
               "40/40 archived values matched; 1651 within-module neighbours ranked by decreasing raw TOM.",
               capture.output(sessionInfo())), file.path(tmp,"run_info.txt"))
  if (file.exists(out) || !file.rename(tmp,out)) stop("Could not publish results. Choose a new output directory.")
  print(top[,c("rank","neighbor","TOM_similarity",flags,"RFX5_DoRothEA_mode")], row.names=FALSE)
  print(counts, row.names=FALSE)
  message("Results: ", normalizePath(out))
}
if (sys.nframe() == 0L) main()
