# Figure 3 cell-cell communication workflow

This folder contains the retained LIANA cell-cell communication analysis for
Figure 3 and Supplementary Figure 8. The numbered scripts should be run in
order. Operations 01 and 02 perform the retained analysis, operations 03 and 04
generate the main and supplementary plots, operation 05 records the fixed
Figure 3C panel selection, and operation 06 performs the final reproducibility
audit.

Reader-facing companion files at the reproducibility root are:

- [`FIGURE3_SOURCE_WORKFLOW.md`](../../FIGURE3_SOURCE_WORKFLOW.md), which maps
  the main and supplementary panels to their retained sources
- [`FIGURE3_SHA256.txt`](../../FIGURE3_SHA256.txt), which verifies the public
  Figure 3 documentation and scripts

Historical notebooks, previous outputs and read-only equivalence checks were
used internally to identify and validate the retained workflow. They are not
part of this public workflow and are not required to run it.

## Workflow

### 1. `01_run_liana_inference.R`

Reconstructs the donor-by-cell-type expression input and runs LIANA.

Main operations:

- reads the validated final annotated Figure 1 Seurat object
- uses raw `Spatial` counts and the existing `singlecell` deconvolution assay
- removes zero-count Spatial genes from the in-memory copy
- constructs 288 weighted pseudobulks for 12 cell types across 24 donors
- applies `LogNormalize` with scale factor 1,000,000
- converts the LIANA input to a `SingleCellExperiment` for compatibility with
  the validated Seurat 5 environment
- runs LIANA with the OmniPath resource and the CellChat, NATMI, SCA, logFC and
  Connectome methods
- aggregates method-specific ranks and retains interactions with
  `aggregate_rank <= 0.01`
- records parameters, input and output checksums, completion status and
  `sessionInfo()`

The validated run contains 609,984 aggregate interactions, of which 29,205
meet the aggregate-rank threshold. It retains 17,775 expressed Spatial genes.

### 2. `02_analyze_stage_activity.R`

Calculates the retained stage-activity summaries used by Figure 3B and
Supplementary Figure 8.

Main operations:

- selects the top 3,000 aggregate-ranked single-gene ligand-receptor candidates
- calculates donor-level LR activity and summarizes it by fibrosis stage
- tests monotonic F0-to-F4 activity trends with Spearman correlation
- applies Benjamini-Hochberg correction
- retains 60 interactions for the stage heatmap and dot plots
- exports the row-z-scored heatmap table
- counts significant source-to-target interactions for the chord diagram
- counts incoming and outgoing cholangiocyte interactions
- creates the four F1, F2, F3 and F4 versus F0 dot-plot input tables
- records the complete statistical-test counts and checksums

The validated run performs 2,903 finite Spearman tests among 3,000 candidates,
identifies 349 BH-significant trends and 315 interactions that additionally
meet the recorded F4-versus-F0 effect-size threshold. The public workflow
retains 60 interactions for visualization.

### 3. `03_plot_main_and_supplementary_ccc.R`

Generates seven deterministic plotting outputs without performing additional
biological analysis:

- Figure 3A: significant-interaction chord diagram
- Figure 3B: stage-activity row-z heatmap and LR annotation table
- Supplementary Figure 8A: cholangiocyte incoming/outgoing partner counts
- Supplementary Figure 8B-E: F1, F2, F3 and F4 versus F0 dot plots

The plot manifest records the panel identity, content and PDF dimensions. The
plotting operation also records its input checksums, parameters and software
information.

### 4. `04_plot_spatial_lr_all_donors.R`

Generates the spatial LRlocal maps for all 24 donors and three retained
interactions:

| Source | Target | Ligand | Receptor |
| --- | --- | --- | --- |
| Cholangiocyte | Stellate | NPNT | ITGA8 |
| Cholangiocyte | InfMac | MMP7 | CD151 |
| Stellate | PortalEndo | CCL2 | ACKR1 |

For spot `i`, the untransformed local score is

```text
S_i = (E_L,i * p_source,i) *
      sum[j in {i plus 30 nearest neighbours}]
      exp(-d_ij^2 / (2 * sigma^2)) * (E_R,j * p_target,j)
```

where `sigma = 2`, `E_L` and `E_R` are LogNormalized Spatial expression and
`p_source` and `p_target` are deconvolution-derived cell-type fractions. The
displayed score is `log1p(S_i * 10000)`.

For each LR pair, identical color limits are used for every donor. The limits
are the first and 99th percentiles across all 26,056 spots, with values outside
the limits clipped. Colors are therefore comparable between stages within one
LR pair, but absolute color values must not be compared between different LR
pairs. These spatial maps are descriptive and are not used for statistical
testing.

The operation uses the custom H&E images while deriving the displayed spot
diameter from the Space Ranger scale factors. Spots are drawn as 48-vertex
polygons at the true Visium diameter with multiplier 1.00. The manifest records
the registered and custom-image geometry, coordinate transformation, diameter,
color limits, file sizes and SHA-256 checksums for every map.

Outputs are organized as:

```text
results/Fig3_CCC/03_spatial_LR_all_samples/
  Healthy_F0|F1|F2|F3|F4/
    LR_PAIR/
      SOURCE_to_TARGET/
        donor_...__C_LRlocal.pdf
```

Exactly 72 LRlocal PDFs are generated. Source-detection, target-detection and
co-detection PDFs are not generated, and no derived spot-level RDS object is
saved.

### 5. `05_select_main_figure3C_panels.R`

Copies six byte-identical PDFs into the ordered Figure 3C panel folder. It
performs no scientific analysis and does not recalculate spatial scores.

The panel identities were fixed from the accepted manuscript composition
before visual inspection of the regenerated maps:

| Position | Stage | Donor | Interaction |
| --- | --- | --- | --- |
| top left | F0 | 6B | NPNT-ITGA8, Cholangiocyte to Stellate |
| top middle | F0 | 2A | MMP7-CD151, Cholangiocyte to InfMac |
| top right | F0 | 4C | CCL2-ACKR1, Stellate to PortalEndo |
| bottom left | F4 | 4B | NPNT-ITGA8, Cholangiocyte to Stellate |
| bottom middle | F4 | 3A | MMP7-CD151, Cholangiocyte to InfMac |
| bottom right | F4 | 2D | CCL2-ACKR1, Stellate to PortalEndo |

The operation verifies the exact LIANA aggregate ranks, shared within-pair
color limits and source-PDF checksums before copying.

### 6. `06_validate_figure3_outputs.R`

Performs the final audit after operations 01 through 05.

Main checks:

- confirms the six public scripts and Figure 3 documentation
- verifies `FIGURE3_SHA256.txt`
- verifies every Step 01-05 input and output checksum manifest
- checks the 609,984 aggregate and 29,205 significant LIANA interactions
- checks the 288 pseudobulks and the 24-donor stage structure
- checks all stage-activity dimensions and recorded statistical-test counts
- checks the 60 selected interactions and row-z-score consistency
- checks all seven Figure 3A-B and Supplementary Figure 8A-E PDFs
- checks the 24-donor by three-interaction spatial manifest, true-Visium
  geometry, shared color limits and all 72 PDF checksums
- checks the six fixed Figure 3C panels, exact aggregate ranks and byte-identical
  copies
- confirms that detection PDFs, derived spatial-score RDS files and incomplete
  `.partial` outputs are absent
- copies only compact verified execution records into
  `reproducibility_logs/Fig3_CCC/Fig3_final_audit/`

The audit checks completeness and internal consistency and performs no
additional biological analysis. It does not run historical comparisons.

## Required inputs

The validated Mac locations are recorded in each operation's input manifest.
The principal inputs are:

| Input | Validated location or role | Used by |
| --- | --- | --- |
| Final annotated Seurat object | `results/260701_approach1_cluster_go_marker_deconv/object_full_with_final_annotation_markers_deconv_proportions_ready_C3green_C5yellow.rds` | LIANA inference and spatial maps |
| Spatial expression | raw counts in the object's `Spatial` assay | LIANA pseudobulks and spatial maps |
| Cell-type fractions | object's `singlecell` assay | weighted pseudobulks and LRlocal maps |
| Sample/stage information | `donor` and `F_score` metadata, validated against the recorded 24-donor map | all donor- and stage-level operations |
| Custom H&E images | `/Users/jianwu/Documents/DATA/H&E/` | spatial maps |
| Step 01 normalized pseudobulks and LIANA outputs | `results/Fig3_CCC/01_liana_inference/` | stage activity |
| Step 02 retained analysis tables | `results/Fig3_CCC/02_main_and_supplementary_figures/` | plotting and spatial interaction specification |

No source RDS is overwritten. The Figure 2 WGCNA object is not modified or
converted into a combined object by this workflow.

## Validated software environment

The reference run used R 4.4.3 with Seurat 5.4.0, SeuratObject 5.3.0, Matrix
1.7.3, dplyr 1.1.4, tidyr 1.3.2, purrr 1.2.0, tibble 3.3.0, LIANA 0.1.14,
OmnipathR 3.17.4, circlize 0.4.17, ComplexHeatmap 2.22.0, RColorBrewer 1.1.3,
scales 1.4.0, semla 1.4.0, ggplot2 4.0.1, patchwork 1.3.2, magick 2.9.0 and
RANN 2.6.2. The operation-specific `sessionInfo()` files are authoritative.

LIANA 0.1.14 calls Seurat's retired `slot` interface when given a Seurat 5.4
object directly. Operation 01 therefore constructs an equivalent
`SingleCellExperiment` with `counts`, natural-log `logcounts` and the same
cell-type identities before calling LIANA. This is a software-compatibility
adaptation and does not alter the retained expression values or LIANA methods.

## Running the workflow

From the reproducibility root, activate the validated `r44-liana` environment
and set:

```bash
export FLINC_REPRO_ROOT="$PWD"
RSCRIPT_BIN="/Users/jianwu/miniforge3/envs/r44-liana/bin/Rscript"
```

Run operations 01 through 05 sequentially:

```bash
"$RSCRIPT_BIN" --vanilla scripts/figure3/01_run_liana_inference.R
"$RSCRIPT_BIN" --vanilla scripts/figure3/02_analyze_stage_activity.R
"$RSCRIPT_BIN" --vanilla scripts/figure3/03_plot_main_and_supplementary_ccc.R
"$RSCRIPT_BIN" --vanilla scripts/figure3/04_plot_spatial_lr_all_donors.R
"$RSCRIPT_BIN" --vanilla scripts/figure3/05_select_main_figure3C_panels.R
```

After the documentation and `FIGURE3_SHA256.txt` are present, run the final
audit last:

```bash
"$RSCRIPT_BIN" --vanilla scripts/figure3/06_validate_figure3_outputs.R
```

The analysis scripts protect existing outputs. The final audit also refuses to
replace a non-empty audit directory; preserve any previous audit in a
timestamped backup before rerunning it.

## Generated outputs and public records

Generated results are written under `results/Fig3_CCC/` and run logs under
`logs/Fig3_CCC/`. These directories are ignored and are not committed to
GitHub. This includes LIANA RDS files, complete result tables and the 78
generated spatial PDFs.

The final audit copies verified compact records—including parameters, input and
output manifests, selected interaction tables, panel manifests, completion
records, checksums and `sessionInfo()`—to
`reproducibility_logs/Fig3_CCC/Fig3_final_audit/`. These compact records may be
staged explicitly after the audit passes.

## Integrity manifest

`FIGURE3_SHA256.txt` records the SHA-256 checksums of the public Figure 3
documentation and scripts. From the reproducibility root, run:

```bash
shasum -a 256 -c FIGURE3_SHA256.txt
```

## Repository snapshot note

Only software-version and session-information records are bundled under
`reproducibility_logs/`. Other run records and generated outputs described above
are local analysis artifacts and are not included in this checkout. See the
[root README](../../README.md) for migration scope and required external inputs.

## Published numerical source tables

[resources/figure3/](../../resources/figure3/README.md) contains the
selected completed-run tables and available reference inputs for this figure.
The [resource manifest](../../resources/manifest.tsv) records exact original paths,
checksums and verification details. The figure resource guide distinguishes
final tests, descriptive summaries and external inputs.
