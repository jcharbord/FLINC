# Figure 2 spatial coexpression and supplementary GSEA workflow

This folder contains the analysis used to construct and characterize the six
spatial coexpression modules shown in Figure 2, together with the independent
donor-level stagewise GSEA used as a supplementary analysis supporting the
Figure 2 interpretation.

The numbered scripts should be run sequentially. Operations 01 through 07
perform the WGCNA analysis and plotting (including the same-module ego networks
in step `05_bis`), operation 08 performs the
supplementary GSEA, and operation 09 performs the final reproducibility audit.

Reader-facing companion files at the reproducibility root are:

- [`FIGURE2_SOURCE_WORKFLOW.md`](../../FIGURE2_SOURCE_WORKFLOW.md), which maps
  the WGCNA, supplementary GSEA and final-audit sources
- [`FIGURE2_SHA256.txt`](../../FIGURE2_SHA256.txt), which verifies the public
  Figure 2 documentation, scripts and static reference table

## Workflow

### 1. `01_prepare_wgcna_input.R`

Prepares the final annotated spatial-expression object for hdWGCNA.

Main operations:

- validates the 24-sample Seurat object and its seven annotated spatial domains
- imports and checks Visium coordinates
- retains the Spatial assay required for network construction
- writes `objects/Fig2_WGCNA/01_wgcna_input_prepared.rds`
- records a run summary and `sessionInfo()`

### 2. `02_build_wgcna_network.R`

Constructs the spatial coexpression network and calculates module properties.

Main operations:

- selects genes detected in at least 5% of spots
- constructs metaspots by `F_score_donor`
- normalizes metaspot expression
- tests soft-thresholding powers
- constructs the signed coexpression network and topological overlap matrix
- identifies six non-grey modules, named `SM_1` through `SM_6`
- calculates harmonized module eigengenes by `F_score_donor`
- calculates harmonized module connectivity (`kME`)
- writes the network object, module assignments, hMEs and network diagnostics

The validated environment for network construction uses R 4.4.1, Seurat
5.3.0, SeuratObject 5.2.0, hdWGCNA 0.4.6, WGCNA 1.73 and harmony 1.2.4.

### 3. `03_analyze_differential_module_eigengenes.R` (historical spot-level analysis)

Retained for traceability; not used for the final manuscript significance claims.
The donor-level analysis is independent and does not require running this script.

### 3 bis. `03_bis_analyze_differential_module_eigengenes_donor_level.R`

Tests donor-mean harmonized module eigengenes across fibrosis stages (updated 2026-10-27).

Main operations:

- retrieves hMEs with `GetMEs(..., harmonized = TRUE)`
- averages spot-level hMEs within each donor, using all available spots
- compares independent donors in F1, F2, F3 and F4 with Healthy (F0) using
  two-sided unpaired Wilcoxon rank-sum tests (`exact = NULL`, `correct = TRUE`)
- applies Benjamini-Hochberg correction across the six modules within each contrast
- reports the difference between mean donor scores, giving each donor equal weight
- exports complete results, donor means, donor/spot QC and module assignments
- generates lollipop and volcano summaries with Healthy display labels

This supplies the final donor-level tests; historical spot-level testing remains in `03_analyze_differential_module_eigengenes.R`. The
network is not reconstructed. The default input is the completed step-02 object,
with `GetMEs(..., harmonized = TRUE, wgcna_name = "vis")`; no automatic fallback
to another score matrix is allowed. The summary and result tables record the
selected ME source. Alternative sources require an explicit `FIG2_ME_SOURCE`
override and do not pass the harmonized-input check in the final Figure 2 audit.

Outputs are written to `results/Fig2_WGCNA/03_DME_donor_level_2026-10-13/`.
Existing outputs are protected against overwrite. The expected cohort is
Healthy/F1/F2/F3/F4 = 6/6/3/3/6 donors; confirm the exported QC for the actual run.

The final manuscript lollipop formatting is built into step 03_bis: 7 mm
asterisks, right-aligned module IDs and descriptions, uniform 12 pt labels with
wrapping, and **Interferon** for SM_2. The lollipop is exported at 11.5 x 7 inches.
No separate replot script is needed; the statistical calculations are unchanged.

### 4. `04_export_module_tables_and_run_enrichr.R`

Exports module membership and hub-gene tables and runs functional enrichment.

Main operations:

- validates the six module assignments, colors, sizes and harmonized kME columns
- ranks genes by signed module-specific kME
- exports all genes, top-70% gene sets and displayed hub-gene sets
- queries EnrichR using every gene assigned to each non-grey module
- saves raw results, significant top terms, figures, query metadata and session
  information

EnrichR databases:

- GO Biological Process 2023
- GO Molecular Function 2023
- GO Cellular Component 2023
- KEGG 2019 Human
- WikiPathways 2024 Human

### 5. `05_generate_functionally_colored_hub_networks.R`

Generates the six functionally colored hub-gene networks.

Main operations:

- selects 50 hubs for `SM_1` and 30 hubs for each of `SM_2` through `SM_6`
- retrieves hub-to-hub topological overlap values
- assigns the accepted module-specific functional categories from
  `resources/figure2/functional_subcategories.tsv`
- exports individual and combined PDF/PNG network figures
- exports the annotated hub table, category legend and category counts

The reference table is part of the plotting specification and is versioned in
`resources/figure2/`; paths in this guide are relative to the repository root.

### 5 bis. `05_bis_build_ego_networks.R`

Sources:

- `05_bis_build_ego_networks.R`
- `resources/figure2/ego_network_functional_context_annotations.tsv`

This analysis uses the regenerated Figure 2 network and TOM workspace. Since
`GetTOM()` returns TOM dissimilarity, edge weights and rankings use
`TOM_similarity <- 1 - GetTOM_output`. Each index gene retains the genuine top
20 neighbours ranked by descending TOM similarity within its own module
(PCSK7: SM_1; PNPLA3: SM_3). No functional gene is force-included.

The annotation table assigns biological-context categories to selected true
neighbours; unannotated genes are labelled `Other`. It does not alter rankings.
This analysis was moved from Figure 4 step 06 to Figure 2 step `05_bis`.
Calculations and plotting settings are unchanged; new figure filenames use
`Fig2_` and gene names. Legacy `11_` table/log prefixes and F/G export identifiers
are retained for traceability and do not assign current manuscript panels.

Supply the completed network, matching TOM workspace, module table and annotation
table through `--network`, `--tom`, `--module-table` and `--annotation`; use new
directories for `--output-dir` and `--log-dir`. The script validates matching TOM
values and writes input/output checksums and completion records. See the command
under Running the workflow. Required packages are hdWGCNA, digest, ggplot2,
ggrepel and scales.

### 6. `06_plot_harmonized_module_eigengenes_and_celltype_associations.R`

Generates the remaining hME-based Figure 2 outputs.

Main operations:

- retrieves `SM_1` through `SM_6` using
  `GetMEs(..., harmonized = TRUE)`
- generates the hME dot plot across the seven final spatial domains
- generates hME spatial maps across all 24 samples
- calculates 72 spot-level Spearman correlations between the six hMEs and 12
  deconvolved cell-type fractions
- generates the hME-cell-type association heatmap
- exports correlation tables, validation results and `sessionInfo()`

The dot plot uses 45-degree cluster labels for readability. This formatting
changes only label presentation; it does not change hME values, spatial-domain
assignments or statistical calculations.

### 7. `07_generate_supplementary_wgcna_umaps.R`

Generates two supplementary UMAP outputs from the completed WGCNA object.

Main operations:

- runs `RunModuleUMAP` on the non-grey network genes using five hub features
  per module, 15 neighbors, minimum distance 0.3, spread 2 and seed 3210
- sizes each gene by its module-specific kME and colors it by module
- exports `hubgene_umap_ggplot.pdf` without individual gene labels or hub circles
- uses `ModuleFeaturePlot(features = "hMEs")` to plot the six harmonized module
  eigengenes on the retained `umap_SCT` spot embedding
- exports `umap.MEs.pdf`, PNG previews, the gene-UMAP coordinates, parameters,
  completion record and `sessionInfo()`

The operation reads the validated network object and its existing topological
overlap matrix. It does not overwrite the network RDS or save another large
Seurat object.

### 8. `08_run_supplementary_stage_gsea.R`

Runs whole-transcriptome donor-level differential expression and preranked GSEA
for F1, F2, F3 and F4 versus F0.

Main operations:

- reads raw counts from the `Spatial` assay of the final annotated Seurat
  object; SCT values and SCT models are not used
- sums spot counts within each of the 24 donors, making donor the independent
  biological unit
- applies `edgeR::filterByExpr()`, TMM normalization and `limma::voom()`
- fits F1, F2, F3 and F4 versus F0 contrasts and applies robust empirical Bayes
  moderation
- ranks every expression-filtered gene by the limma moderated t-statistic,
  without applying an individual-gene significance cutoff
- runs seeded GSEA through clusterProfiler's fgsea backend for MSigDB H, C2,
  C3 and C5 with gene-set
  sizes 15-5,000, `nPermSimple = 20000` and BH FDR reporting
- exports complete and significant GSEA tables, differential-expression tables,
  rankings, pseudobulk objects and the frozen TERM2GENE database snapshot
- generates MDS and voom mean-variance diagnostics
- generates four horizontal Hallmark plots displaying every Hallmark set with
  BH FDR below 0.05 on common NES, gene-set-size and capped `-log10(FDR)` scales
- records the input and result SHA-256 checksums, package/MSigDB versions,
  parameters, completion status and `sessionInfo()`

The GSEA is analytically independent of hdWGCNA and does not read the WGCNA
network object. F2 and F3 each contain three independent donors; their results
should therefore be interpreted with the corresponding limitation in precision.

### 9. `09_validate_figure2_outputs.R`

Performs a final, lightweight audit after operations 01 through 08.

Main operations:

- confirms that the public analysis and replot scripts, documentation and reference table are
  present
- confirms that the prepared and completed WGCNA objects and execution records
  exist
- checks the six-module assignment table and harmonized kME columns
- checks the 24 donor-level differential hME tests, donor counts, equal-donor
  effect sizes and within-contrast BH correction
- checks 30 module-database EnrichR queries, 200 displayed hubs and 72
  hME-cell-type correlations
- checks the expected hub-network, dot-plot, spatial-map and heatmap files
- checks the 2,206-gene module UMAP, the absence of gene-label requests and the
  harmonized-ME feature-map specification
- checks donor replication, the whole-list moderated-t rankings, all four
  comparisons and all four MSigDB collections in the supplementary GSEA
- verifies the GSEA result and execution-record SHA-256 manifests
- detects incomplete `.partial` outputs
- records code, reference and final-output checksums, component status and
  `sessionInfo()` under
  `reproducibility_logs/Fig2_WGCNA/Fig2_final_audit/`

This audit checks output completeness and internal consistency and performs no
additional biological analysis. For step `05_bis`, it checks the source script
and annotation table; the ego-network script records and validates its own
outputs during execution.

## Harmonized module eigengenes

Harmonized module eigengenes are used consistently for:

- differential module eigengene analysis
- the spatial-domain dot plot
- the 24-sample spatial maps
- module-cell-type correlations

Hub networks and enrichment analyses use module assignments, harmonized kME
values and topological overlap, as appropriate for those analyses.

## Module names

| Module | Biological interpretation |
| --- | --- |
| SM_1 | Fibro-immune ECM remodeling |
| SM_2 | Interferon |
| SM_3 | Cholesterol metabolism and lipogenesis |
| SM_4 | Hepatocyte homeostatic function |
| SM_5 | Stress and protein homeostasis |
| SM_6 | Hepatic plasma |

## Required inputs

The default paths are relative to the workflow root unless an environment
variable overrides them.

| Input | Default location | Used by |
| --- | --- | --- |
| Final annotated Seurat object | `results/260701_approach1_cluster_go_marker_deconv/object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds` | WGCNA input preparation and supplementary donor-level GSEA |
| Prepared WGCNA object | `objects/Fig2_WGCNA/01_wgcna_input_prepared.rds` | network construction |
| Completed network object | `objects/Fig2_WGCNA/02_wgcna_network_modules.rds` | DME, hub networks, hME plots and module-gene UMAP |
| Module assignment and kME table | `results/Fig2_WGCNA/02_network/02_module_assignments_and_kME.tsv` | enrichment |
| Functional hub categories | `resources/figure2/functional_subcategories.tsv` | colored hub networks |
| Deconvolution fractions | network-object metadata or the object specified by `FIG2_DECONVOLUTION_OBJECT_INPUT` | hME-cell-type associations |

The deconvolution calculation is upstream of this workflow. Figure 2 imports
the resulting spot-level fractions and does not recalculate them.

## Software requirements

The validated network-construction environment uses R 4.4.1 with Seurat 5.3.0,
SeuratObject 5.2.0, hdWGCNA 0.4.6, WGCNA 1.73, harmony 1.2.4 and patchwork.
The downstream WGCNA operations additionally use ggplot2, ggrepel, writexl,
enrichR, igraph, pheatmap, viridis and ggrastr.

The validated supplementary GSEA run uses R 4.4.1, Seurat 5.3.0,
SeuratObject 5.2.0, Matrix 1.7.4, edgeR 4.4.2, limma 3.62.2, statmod 1.5.1,
msigdbr 25.1.1, clusterProfiler 4.14.6, fgsea 1.32.4, AnnotationDbi 1.68.0,
org.Hs.eg.db 3.20.0, ggplot2 4.0.0, stringr 1.5.2 and digest 0.6.37. The
execution record freezes the exact MSigDB version and H/C2/C3/C5 TERM2GENE
snapshot used by the run. Operation 09 additionally requires digest for
SHA-256 verification.

## Running the workflow

From the workflow root, activate the validated R environment
(`conda activate hdWGCNA_pipeline-env` on Jian's Mac) and set:

```bash
export FLINC_REPRO_ROOT="$PWD"
```

Run the scripts sequentially:

```bash
Rscript --vanilla scripts/figure2/01_prepare_wgcna_input.R
Rscript --vanilla scripts/figure2/02_build_wgcna_network.R
Rscript --vanilla scripts/figure2/03_bis_analyze_differential_module_eigengenes_donor_level.R
Rscript --vanilla scripts/figure2/04_export_module_tables_and_run_enrichr.R
Rscript --vanilla scripts/figure2/05_generate_functionally_colored_hub_networks.R
Rscript --vanilla scripts/figure2/05_bis_build_ego_networks.R \
  --network objects/Fig2_WGCNA/02_wgcna_network_modules.rds \
  --tom objects/Fig2_WGCNA/02_network_workspace/TOM/test_TOM.rda \
  --module-table results/Fig2_WGCNA/02_network/02_module_assignments_and_kME.tsv \
  --annotation resources/figure2/ego_network_functional_context_annotations.tsv \
  --output-dir results/Fig2_WGCNA/05_bis_ego_networks \
  --log-dir reproducibility_logs/Fig2_WGCNA/05_bis_ego_networks
Rscript --vanilla scripts/figure2/06_plot_harmonized_module_eigengenes_and_celltype_associations.R
Rscript --vanilla scripts/figure2/07_generate_supplementary_wgcna_umaps.R
Rscript --vanilla scripts/figure2/08_run_supplementary_stage_gsea.R
Rscript --vanilla scripts/figure2/09_validate_figure2_outputs.R
```

The analysis scripts refuse to overwrite existing final or partial results.
Existing outputs should be reviewed or moved before intentionally repeating a
stage. The GSEA script accepts optional explicit positional paths for the input
RDS, results root and reproducibility-log root; the defaults shown above are
used when these are omitted. Set `FIG2_GSEA_RUN_NAME` to a specific completed
run before operation 09 when more than one GSEA run is present. Otherwise the
latest paired completed result/log run is audited. The final audit can be rerun
to refresh its compact audit record.

## Outputs and execution records

Generated objects and results are written under:

- `objects/Fig2_WGCNA/`
- `results/Fig2_WGCNA/`
- `reproducibility_logs/Fig2_WGCNA/`
- `results/Fig2_GSEA/08_stagewise_gsea/`
- `reproducibility_logs/Fig2_GSEA/08_stagewise_gsea/`

The `results/` directory contains generated figures, tables and objects and is
ignored by Git. These generated GSEA results must not be committed. Compact
parameters, completion reports, input manifests, checksums and session-information
files under `reproducibility_logs/` provide the execution record and may be
staged selectively using explicit paths.


## Repository snapshot note

Only software-version and session-information records are bundled under
`reproducibility_logs/`. Other run records and generated outputs described above
are local analysis artifacts and are not included in this checkout. See the
[root README](../../README.md) for migration scope and required external inputs.

## Published numerical source tables

Selected completed-run exports are available in [source_data/figure2/](../../source_data/figure2/README.md).
The [source-data manifest](../../source_data/manifest.tsv) identifies the original
run directories and checks performed. This is a partial source-data release;
see the linked guide for included and missing outputs. It does not replace the
external inputs required to run the analysis.
