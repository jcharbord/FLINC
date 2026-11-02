# Plotting functions for 12_rfx5_sm1_neighbours.R.
# Frozen rows of resources/figure2/results/06_hubs_with_functional_subcategories.tsv
# in jcharbord/FLINC, revision 184b505abc561bc3ce456a5c62e66cbb4f26b99f.
# Original hub categories are retained unless the MHC class II target subset takes priority.
hubs <- read.delim(text = "gene	hub_rank	class
IGFBP7	1	ECM
TMSB4X	2	Cytoskeleton
VIM	3	Cytoskeleton
CD74	4	Immune process
HLA-DRA	5	Oxidative stress
COL1A2	6	ECM
LYZ	7	Immune process
S100A6	8	Cytoskeleton
COL3A1	9	ECM
AEBP1	10	ECM
FLNA	11	Cytoskeleton
HLA-DPB1	12	Immune process
DCN	13	ECM
COL1A1	14	ECM
LUM	15	ECM
TIMP1	16	ECM
ACTB	18	Cytoskeleton
LSP1	25	Immune process
UBD	27	Immune process
LAPTM5	31	Immune process", stringsAsFactors=FALSE)

annotate_neighbours <- function(d) {
  required <- c("rank", "neighbor", "RFX5_DoRothEA_mode")
  if (!all(required %in% names(d))) stop("The annotated TSV must contain: ", paste(required, collapse = ", "))
  if ("index_gene" %in% names(d) && !all(d$index_gene == "RFX5")) stop("Expected RFX5 as the index gene.")
  if ("module" %in% names(d) && !all(d$module == "SM_1")) stop("Expected SM_1 as the module.")
  d <- d[order(d$rank), , drop = FALSE]
  if (!identical(as.integer(d$rank), 1:20) || !setequal(d$neighbor, hubs$gene))
    stop("Expected the 20 distinct genes in the FLINC RFX5 figure and ranks 1 to 20.")
  rows <- match(d$neighbor, hubs$gene)
  d$hub_rank <- hubs$hub_rank[rows]
  d$group <- unname(c(ECM = "ecm", Cytoskeleton = "cytoskeleton",
                      "Immune process" = "immune")[hubs$class[rows]])
  mhc <- d$neighbor %in% c("CD74", "HLA-DRA", "HLA-DPB1")
  d$group[mhc] <- "mhc"
  mode <- d$RFX5_DoRothEA_mode
  stopifnot(!anyNA(d$group), all(mode[mhc] == "activates"),
            identical(mode[d$neighbor == "COL1A2"], "represses"), sum(!is.na(mode)) == 4L)
  d
}

open_pdf <- function(path, width, height, title) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width, height = height, family = "sans")
  } else {
    grDevices::pdf(path, width = width, height = height, family = "Helvetica",
                   useDingbats = FALSE, title = title)
  }
}

# SM_1 colours from Figure 2A; darker orange distinguishes the MHC class II subset.
rfx5_colours <- c(index = "#F29797", mhc = "#C35A0D", immune = "#F4C243",
                  ecm = "#1432F3", cytoskeleton = "#4FAC5B")

draw_table <- function(d, path) {
  open_pdf(path, 10.7, 8.7, "RFX5 neighbours and module hubs in SM_1")
  on.exit(dev.off())
  par(mar=c(0,0,0,0), xaxs="i", yaxs="i", family="sans")
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0,100))
  label <- function(x,y,txt,align=0,...) text(x,y,txt,adj=c(align,0.5),...)
  ink <- "#172B38"
  muted <- "#475A63"
  label(3,95.4,"The 20 nearest RFX5 neighbours in SM_1",font=2,cex=1.45,col=ink)
  label(3,91.7,"Index gene = RFX5     |     SM_1 kME rank = 712 of 1,652 (kME = 0.239).",cex=1.0,col=muted)
  rect(2.5,82.5,97.5,88.5,col="#E5EBEE",border=NA)
  x <- c(7.6,18.0,24.6,43.4,56.6,67.9,78.0,91.2)
  labs <- c("RFX5\nTOM rank", "SM_1\nkME rank", "Gene", "MHC class II\nRFX5 target", "Other immune\nprocess", "ECM", "Cytoskeleton", "DoRothEA\nregulation")
  for (j in seq_along(x)) label(x[j],85.6,labs[j],align=if (j==3) 0 else .5,font=2,cex=.83,col=ink)
  row_y <- 80.8 - (seq_len(nrow(d))-1)*3.22
  for (i in seq_len(nrow(d))) {
    y <- row_y[i]
    if (i %% 2 == 1) rect(2.5,y-1.58,97.5,y+1.58,col="#F2F5F6",border=NA)
    label(x[1],y,d$rank[i],align=.5,cex=.85,col=muted)
    label(x[2],y,d$hub_rank[i],align=.5,cex=.85,col=muted)
    label(x[3],y,d$neighbor[i],cex=.94,font=2,col=ink)
    if (d$group[i] == "mhc") points(x[4],y,pch=18,col=rfx5_colours["mhc"],cex=1.18)
    if (d$group[i] == "immune") points(x[5],y,pch=16,col=rfx5_colours["immune"],cex=.87)
    if (d$group[i] == "ecm") points(x[6],y,pch=16,col=rfx5_colours["ecm"],cex=.87)
    if (d$group[i] == "cytoskeleton") points(x[7],y,pch=16,col=rfx5_colours["cytoskeleton"],cex=.87)
    if (!is.na(d$RFX5_DoRothEA_mode[i])) label(x[8],y,if (d$RFX5_DoRothEA_mode[i]=="activates") "Activation" else "Repression",
                               align=.5, cex=.86, col=if (d$RFX5_DoRothEA_mode[i]=="activates") "#A23B36" else "#2549BC")
  }
  segments(2.5,17.6,97.5,17.6,col="#D2DCE1",lwd=.9)
  label(3,15.2,"TOM ranks genes by their overlap with RFX5. kME ranks genes by correlation with the SM_1 eigengene.",cex=.79,col=muted)
  label(3,12.4,"Filled symbols show Figure 2 functional groups, with HLA-DRA reclassified under MHC class II antigen presentation.",cex=.79,col=muted)
  label(3,9.6,"DoRothEA supplies the regulatory direction. The three MHC class II-related targets include the invariant chain CD74.",cex=.79,col=muted)
}

tom_scale <- function(scores) {
  # Inputs are TOM similarities in neighbour-rank order, largest first.
  if (length(scores) != 20L || !is.numeric(scores) || any(!is.finite(scores)) ||
      any(diff(scores) > 0)) stop("Expected 20 finite TOM similarities in decreasing rank order.")
  bounds <- range(scores)
  if (diff(bounds) == 0) return(list(width = rep(2, 20), key = 2, labels = sprintf("%.4g", bounds[1])))
  list(width = 0.8 + 3.2 * (scores - bounds[1]) / diff(bounds),
       key = c(0.8, 4), labels = sprintf("%.4g", bounds))
}

draw_network <- function(d, path, scale) {
  stopifnot(length(scale$width) == nrow(d), length(scale$key) == length(scale$labels))
  open_pdf(path, 11.5, 7.0, "RFX5 co-expression neighbourhood in SM_1")
  on.exit(dev.off())
  par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans")
  plot.new()
  plot.window(xlim = c(0, 11.5), ylim = c(0, 7.0), asp = 1)
  ink <- "#18292F"
  muted <- "#4E6066"
  edge <- "#90999B"
  dot <- function(x, y, radius, fill) {
    symbols(x, y, circles = rep_len(radius, length(x)), inches = FALSE,
            bg = fill, fg = NA, add = TRUE)
  }
  left <- function(x, y, label, ...) text(x, y, label, adj = c(0, 0.5), ...)
  left(0.40, 6.65, "RFX5 co-expression neighbourhood in SM_1", cex = 1.52, font = 2, col = ink)
  left(0.40, 6.31, "The 20 nearest genes, ordered clockwise by decreasing TOM similarity", cex = 0.94, col = muted)

  angle <- pi / 2 - (seq_len(nrow(d)) - 0.5) * 2 * pi / nrow(d)
  cx <- 3.63
  cy <- 3.16
  radius <- 2.43
  x <- cx + radius * cos(angle)
  y <- cy + radius * sin(angle)
  segments(cx, cy, x, y, col = edge, lwd = scale$width)
  dot(x, y, 0.109, rfx5_colours[d$group])
  # Centre the top and bottom labels vertically outside their nodes.
  side <- sign(cos(angle))
  for (i in seq_along(x)) {
    if (d$neighbor[i] %in% c("TMSB4X", "UBD")) {
      text(x[i], y[i] + 0.245, d$neighbor[i], adj = c(0.5, 0.5), cex = 0.98, col = ink)
    } else if (d$neighbor[i] %in% c("ACTB", "HLA-DRA")) {
      text(x[i], y[i] - 0.245, d$neighbor[i], adj = c(0.5, 0.5), cex = 0.98, col = ink)
    } else {
      text(x[i] + 0.17 * side[i], y[i] + 0.145 * sin(angle[i]), d$neighbor[i],
           adj = c(if (side[i] > 0) 0 else 1, 0.5), cex = 0.98, col = ink)
    }
  }
  dot(cx, cy, 0.315, rfx5_colours["index"])
  text(cx, cy, "RFX5", cex = 1.0, font = 2, col = ink)

  lx <- 7.62
  left(lx, 5.45, "Functional annotation", cex = 1.10, font = 2, col = ink)
  key_names <- c("RFX5 (index gene)", "MHC class II-related RFX5 target",
                 "Other immune process", "Extracellular matrix (ECM)", "Cytoskeleton")
  key_y <- 5.02 - (0:4) * 0.405
  for (i in seq_along(key_y)) {
    dot(lx + 0.10, key_y[i], 0.092, rfx5_colours[i])
    left(lx + 0.35, key_y[i], key_names[i], cex = 0.89, col = ink)
  }
  left(lx, 2.67, "TOM similarity to RFX5", cex = 1.10, font = 2, col = ink)
  for (i in seq_along(scale$key)) {
    ky <- 2.27 - (i - 1) * 0.36
    segments(lx, ky, lx + 0.45, ky, col = edge, lwd = scale$key[i])
    left(lx + 0.60, ky, scale$labels[i], cex = 0.89, col = ink)
  }
}

# Standalone redraw from existing results; 12_rfx5_sm1_neighbours.R sources the functions above.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript 12_plot_rfx5_sm1.R <RFX5 results folder>")
  out <- normalizePath(args[1], mustWork = TRUE)
  d <- annotate_neighbours(read.delim(file.path(out, "RFX5_top20_annotated.tsv"),
                                     stringsAsFactors = FALSE, check.names = FALSE))
  path <- file.path(out, "RFX5_top20_network_by_function.pdf")
  draw_network(d, path, tom_scale(d$TOM_similarity))
  table_path <- file.path(out, "RFX5_top20_functions_hubs.pdf")
  draw_table(d, table_path)
  message("Saved ", path, "\nSaved ", table_path)
}
