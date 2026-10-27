# Figure 4 source workflow: regulatory, network, and secreted-factor analyses

## Scope

This document describes the final, public Figure 4 workflow for the FLINC spatial transcriptomics study.  It is deliberately organised as one reproducible analysis sequence rather than as a record of exploratory plotting and troubleshooting.

The workflow uses the regenerated final Seurat object and the regenerated Figure 2 WGCNA outputs.  Large input objects and generated figures are intentionally excluded from Git; their identities are recorded by SHA-256 checksums in the local execution records and in `FIGURE4_SHA256.txt`.

## Canonical public workflow

| Step | Script | Figure content |
|---:|---|---|
| 01 | `scripts/figure4/01_run_progeny_activity.R` | PROGENy pathway activities and donor-level F4 versus F0 inference (Figure 4A). |
| 02 | `scripts/figure4/02_plot_progeny_spatial_maps.R` | Representative-section PROGENy maps (Figure 4B). |
| 03 | `scripts/figure4/03_run_dorothea_activity.R` | DoRothEA A/B transcription-factor activity and donor-level inference (Figure 4C). |
| 04 | `scripts/figure4/04_build_rfx5_regulon.R` | RFX5 A/B regulon, restricted to targets observed in the regenerated object (Figure 4D). |
| 05 | `scripts/figure4/05_rfx5_celltype_associations.R` | Within-stage spot-level RFX5/cell-type fraction associations (Figure 4E). |
| 06 | `scripts/figure4/06_build_ego_networks.R` | PCSK7 and PNPLA3 strict top-20 same-module TOM-neighbour ego networks (Figure 4F/G). |
| 07 | `scripts/figure4/07_run_secreted_factor_analysis.R` | Strict secreted-factor selection and mean-expression ranking (Figure 4H/I and supplementary counts). |
| 08 | `scripts/figure4/08_plot_secreted_enrichment.R` | Readable supplementary enrichment rendering from the versioned enrichment reference table. |
| 09 | `scripts/figure4/09_validate_figure4_outputs.R` | Public-package validation and SHA-256 manifest generation. |

## Main analytical definitions

### PROGENy

`01_run_progeny_activity.R` uses human PROGENy top-100 target models and `decoupleR::run_wmean`, retaining `norm_wmean`. Counts are library-size normalised with `log1p(count / spot-library-size * 10000)` and target genes are z-scored across spots. The definitive run uses 1,000 permutations and seed 42. Spot activities are summarised as donor medians; inference is performed on donors, not on individual spots. The primary contrast is F4 versus F0; F1 versus F0 is supplementary, and F2/F3 contrasts are validation-only.

### Spatial pathway maps

`02_plot_progeny_spatial_maps.R` draws Visium spots in low-resolution H&E image coordinates. Spot radius is calculated from the Visium scale factors rather than from Seurat's native normalised radius. The main representative maps are NFkB, p53, and TNFa; WNT and Androgen are validation maps. Colour scales and section choices are recorded in the local output manifest.

### DoRothEA and RFX5

`03_run_dorothea_activity.R` uses the human DoRothEA A/B regulon with `norm_wmean`, 1,000 permutations, and seed 42. Deterministic spot batching is used solely to stay within memory limits; an equivalence audit verifies agreement with an unbatched subset. `04_build_rfx5_regulon.R` uses only RFX5 A/B targets that are actually present and non-zero in the regenerated expression object. `05_rfx5_celltype_associations.R` reports descriptive within-stage, spot-level Spearman correlations with deconvolution fractions; these are not interpreted as cell-type-specific causal effects.

### WGCNA ego networks

`06_build_ego_networks.R` uses the regenerated Figure 2 network and TOM workspace. `GetTOM()` returns TOM **dissimilarity**, so all edge weights and rankings use:

```r
TOM_similarity <- 1 - GetTOM_output
```

For each index gene, the plot contains the genuine top 20 genes ranked by descending TOM similarity within the same regenerated module. No functional gene is force-included. Node roles are controlled by the versioned `ego_network_functional_context_annotations.tsv` table; unannotated neighbours are labelled `Other`.

### Strict secreted factors and enrichment rendering

`07_run_secreted_factor_analysis.R` applies the predefined strict secretion rule to the versioned UniProt snapshot: signal peptide present, no transmembrane feature, and no membrane localisation. Genes are then ranked by mean log-normalised spatial expression within their regenerated WGCNA module.

The public resource directory also contains compact reference tables used to verify the historical-equivalent count/ranking results and to render the supplementary enrichment panel. The enrichment panel is rendered from a versioned reference table rather than from a live GO/KEGG web request; this avoids silently changing database mappings between reruns. The input secreted-gene lists and WGCNA universes were checked to be identical to the reference inputs before rendering.

## Versioned small resources

The following are retained in Git because they are small, explicit inputs to the public workflow:

```text
resources/Fig4_secreted_factor/
  uniprot_human_secreted_only_KW-0964_reviewed.tsv
  reference_secreted_counts_by_module.csv
  reference_top30_secreted_factors_by_module.csv
  reference_secreted_enrichment_top2.csv

scripts/figure4/
  ego_network_functional_context_annotations.tsv
```

## Inputs intentionally excluded from Git

- Regenerated Seurat object (multi-gigabyte RDS).
- Regenerated Figure 2 WGCNA network RDS and TOM workspace.
- Broad `results/` and `reproducibility_logs/` directories.
- Rendered figure PDFs and PNGs.

These files remain available in the protected analysis environment. The final scripts accept their paths explicitly and record their SHA-256 identities during execution.

## Reproduction sequence

1. Recreate or obtain the final Seurat object and Figure 2 WGCNA outputs.
2. Run scripts 01 through 08 in numerical order, supplying the documented input/output paths.
3. Run script 09 to validate the public package and generate `FIGURE4_SHA256.txt`.
4. Compare regenerated local figures/tables with the compact execution records before manuscript use.

The exact command-line arguments and required files for each script are listed in [`scripts/figure4/README.md`](scripts/figure4/README.md).
