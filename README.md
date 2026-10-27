# FLINC spatial transcriptomics analysis

This directory contains the spatial transcriptomics workflow used to analyze
24 human liver FFPE samples across healthy liver and MASLD fibrosis stages.

The current repository structure describes the validated reference analysis
performed on the Mac environment on which the final results were reproduced.

The workflow is organized by scientific operation rather than by manuscript
panel letters or historical script numbers.

## Workflow overview

The validated workflow links the 24-sample cohort and Space Ranger-derived
count matrices to shared preprocessing in R, STRIDE deconvolution, and the
downstream analyses used for Figures 1–4.

**Workflow overview figure:** [View the original FLINC workflow diagram](https://github.com/jcharbord/FLINC_archive/blob/jeremie-jian-reproducible-rerun/Jeremie/Reproducible_code/workflow/FLINC_workflow_overview.png).

The original PNG is awaiting transfer into this repository; the diagram remains
available in the private archive at the link above.

## Analysis overview

The principal workflow consists of two sequential analysis stages:

1. `scripts/analysis/01_preprocess_integrate_and_reduce_spatial_data.R`
2. `scripts/analysis/02_cluster_annotate_and_characterize_spatial_domains.R`

The first stage performs preprocessing, normalization, integration and
dimensionality reduction. The second stage performs spatial-domain clustering,
biological annotation and downstream characterization.

Figure 1 visualization scripts are stored under `scripts/figure1/`. The
Figure 2 spatial coexpression and supplementary donor-level GSEA workflow is
stored under `scripts/figure2/`. The Figure 3 LIANA cell-cell communication,
stage-activity and spatial LR-mapping workflow is stored under
`scripts/figure3/`. The Figure 4 regulatory, network and strict secreted-factor
workflow is stored under `scripts/figure4/`. Historical analytical source files for the main spatial
clustering workflow are retained under `scripts/source_history/`.

## Final manuscript updates

The final donor-level DME script is
`scripts/figure2/03_bis_analyze_differential_module_eigengenes_donor_level.R`,
including the final plot formatting. Spot-level module eigengenes are averaged
within each donor; stage-versus-Healthy contrasts use independent donors and
BH correction across the six modules within each contrast. The companion `03`
script is the historical spot-level analysis, not the final inferential result.

The final secreted-factor script is
`scripts/figure4/11_plot_SM1_SM6_selected4_BH30.R`. It calculates donor-level
Spearman correlations for all 30 genes per module and applies BH correction
separately within SM_1 and SM_6. Four selected factors per module are plotted.
This guarded rerun also requires saved manifests and validation exports listed
in [the Figure 4 guide](scripts/figure4/README.md#validation-inputs-required-for-reproducibility);
a fresh checkout alone is not sufficient to run it.

## Required input data

### 1. Sample metadata

A metadata table describing the 24 spatial transcriptomics samples.

Expected sample identifier column:

`linkID_spatial_ffpe`

Validated Mac location:

`/Users/jianwu/Documents/DATA/metadata/metadata_spatial_ffpe0.2.csv`

### 2. Space Ranger outputs

One Space Ranger output directory is required for each spatial sample.

Validated Mac location:

`/Users/jianwu/Documents/DATA/spaceranger_output/`

Each sample directory is loaded with `Load10X_Spatial()`.

Registered high-resolution histology images are read from:

`<sample>/spatial/tissue_hires_image.png`

### 3. Cell-type deconvolution results

Spot-level cell-type compositions are imported from the validated STRIDE output
`data/deconvolution/single_nucleus_spot_celltype_frac.txt` (26,056 spots × 12
cell-type fractions). The matrix is stored in the final Seurat object as the
`singlecell` assay and used downstream as relative composition data.

The principal workflow imports this validated input rather than silently
rerunning STRIDE. The upstream reference export, spatial export, STRIDE command,
and Seurat import scripts are retained in `scripts/deconvolution/jian_original/`. Full provenance,
validation requirements, and the reproducibility boundary are in
[`DECONVOLUTION_SOURCE_WORKFLOW.md`](DECONVOLUTION_SOURCE_WORKFLOW.md).

## Main analysis workflow

### Stage 1 — preprocessing, integration and dimensionality reduction

Run:

`scripts/analysis/01_preprocess_integrate_and_reduce_spatial_data.R`

Main operations:

- load the 24 Space Ranger datasets
- create Seurat spatial objects
- rename spots and attach sample metadata
- perform SCTransform normalization
- select SCT integration features
- identify integration anchors
- integrate the 24 samples
- import the precomputed cell-type deconvolution matrix
- perform PCA
- perform Harmony correction using sample identity
- calculate the UMAP representation
- save intermediate and integrated analysis objects
- record parameters, checksums and software information

### Stage 2 — clustering, annotation and biological characterization

Run:

`scripts/analysis/02_cluster_annotate_and_characterize_spatial_domains.R`

Main operations:

- graph-based clustering
- hierarchical subclustering where required
- final spatial-domain annotation
- cluster marker identification
- Gene Ontology enrichment
- deconvolution-based cluster characterization
- cluster-cell-type chord relationships
- stage-associated cluster proportion analyses
- generation of the final annotated Seurat object
- generation of analytical results used in Figure 1

Final annotated spatial domains:

- C0 — Portal tracts
- C1 — Periportal hepatocytes
- C2 — Centrilobular hepatocytes
- C3 — Immune cells
- C4 — Myofibroblasts
- C5 — Cholangiocytes
- C6 — Matrix-rich fibroblasts

## Figure 1 source material

Figure-specific source scripts are located in `scripts/figure1/`.

They generate or validate:

- spatial cluster overlays on registered H&E images
- spatial maps of top cluster markers
- the final annotated UMAP
- the selected histology crop
- Figure 1 source-output validation

Several Figure 1 results are generated directly during the main biological
analysis rather than by separate plotting scripts:

- top cluster marker summaries
- Gene Ontology enrichment
- cell-type deconvolution summaries
- cluster-cell-type chord diagrams
- stage-associated cluster proportions

See `FIGURE1_SOURCE_WORKFLOW.md` for detailed source information.

## Figure 2 spatial coexpression workflow

The Figure 2 workflow is located in `scripts/figure2/` and contains eight
numbered analytical and plotting operations followed by a final reproducibility
audit:

1. preparation of the final annotated spatial-expression object for hdWGCNA
2. spatial coexpression network construction and hME calculation
3. donor-level differential harmonized module eigengene analysis (`03_bis`)
4. module-table export and EnrichR functional enrichment
5. generation of functionally colored hub-gene networks
6. hME dot plots, spatial maps and cell-type associations
7. supplementary module-gene UMAP and harmonized-ME feature maps
8. independent donor-pseudobulk stagewise GSEA for F1, F2, F3 and F4 versus F0
9. final output, execution-record and checksum audit

Harmonized module eigengenes are used consistently for differential module
analysis, spatial-domain dot plots, 24-sample spatial maps and spot-level
module-cell-type correlations.

The supplementary GSEA uses raw `Spatial` counts aggregated to the 24 donors,
limma-voom moderated-t whole-gene-list rankings and MSigDB H, C2, C3 and C5.
It is analytically independent of hdWGCNA and does not read the WGCNA network
object.

The final audit checks output completeness and internal consistency and performs
no additional biological analysis.

The six modules are:

- SM_1 — Fibro-immune ECM remodeling
- SM_2 — Interferon
- SM_3 — Cholesterol metabolism and lipogenesis
- SM_4 — Hepatocyte homeostatic function
- SM_5 — Stress and protein homeostasis
- SM_6 — Hepatic plasma

See `scripts/figure2/README.md` for inputs, software requirements, execution
order and output details.

## Figure 3 cell-cell communication workflow

The Figure 3 and Supplementary Figure 8 workflow is located in
`scripts/figure3/` and contains five numbered analytical or plotting operations
followed by a final reproducibility audit:

1. LIANA inference from 12 cell-type by 24-donor weighted pseudobulks
2. stage-progressive LR-activity analysis and retained plotting tables
3. Figure 3A-B and Supplementary Figure 8A-E plotting
4. all-donor spatial mapping of three retained LR interactions
5. checksum-verified selection of the six accepted Figure 3C panels
6. final output, execution-record and checksum audit

LIANA uses the OmniPath resource and aggregates the CellChat, NATMI, SCA,
logFC and Connectome methods. Significant interactions are retained at
`aggregate_rank <= 0.01`.

The spatial maps use normalized ligand and receptor expression weighted by
deconvolved source- and target-cell fractions. Receptor signal is aggregated
over the focal spot and its 30 nearest neighbours with a Gaussian kernel. Each
LR pair uses one common first-to-99th-percentile color scale across all 24
donors, allowing within-pair stage comparisons. The maps are descriptive and
are not used for statistical testing.

True Visium spot diameters are derived from the Space Ranger scale factors and
recorded for every donor. The workflow generates only LRlocal maps; separate
source-detection, target-detection and co-detection PDFs are excluded.

The final audit checks completeness and internal consistency and performs no
additional biological analysis. Historical notebooks and internal equivalence
checks are not part of the public workflow.

See `scripts/figure3/README.md` for inputs, software requirements, execution
order, formulas and output details.

## Figure 4 regulatory, network and secreted-factor workflow

The Figure 4 workflow is located in `scripts/figure4/` and contains eight
core analytical or plotting operations, a public-package validation step,
and the spatial-context and final donor-profile scripts:

1. PROGENy pathway activity and donor-level F4 versus F0 inference
2. true-Visium-spot PROGENy spatial maps
3. DoRothEA A/B transcription-factor activity
4. RFX5 regulon construction from observed DoRothEA targets
5. RFX5/cell-type deconvolution-fraction associations
6. strict same-module PCSK7 and PNPLA3 TOM-similarity ego networks
7. strict secreted-factor selection and expression ranking
8. versioned enrichment-panel rendering
9. public-file checksum validation
10. secreted-factor spatial-cluster summaries and maps
11. final SM_1/SM_6 donor-level stage profiles with BH correction across 30 genes per module

The pathway and transcription-factor analyses use donor-level inference; the
RFX5/fraction correlations and spatial maps are descriptive. Ego-network edges
use TOM similarity defined explicitly as `1 - GetTOM()` because `GetTOM()`
returns TOM dissimilarity. Generated objects, broad result folders and rendered
figures are retained locally rather than committed to Git.

See `FIGURE4_SOURCE_WORKFLOW.md` and `scripts/figure4/README.md` for the
complete source map, methods, required versioned resources and execution order.

## Figure-specific documentation and integrity manifests

Each figure has the same reader-facing documentation pattern:

| Figure | Source-workflow map | Script-folder guide | SHA-256 manifest |
| --- | --- | --- | --- |
| Figure 1 | [`FIGURE1_SOURCE_WORKFLOW.md`](FIGURE1_SOURCE_WORKFLOW.md) | [`scripts/figure1/README.md`](scripts/figure1/README.md) | [`FIGURE1_SHA256.txt`](FIGURE1_SHA256.txt) |
| Figure 2 | [`FIGURE2_SOURCE_WORKFLOW.md`](FIGURE2_SOURCE_WORKFLOW.md) | [`scripts/figure2/README.md`](scripts/figure2/README.md) | [`FIGURE2_SHA256.txt`](FIGURE2_SHA256.txt) |
| Figure 3 | [`FIGURE3_SOURCE_WORKFLOW.md`](FIGURE3_SOURCE_WORKFLOW.md) | [`scripts/figure3/README.md`](scripts/figure3/README.md) | [`FIGURE3_SHA256.txt`](FIGURE3_SHA256.txt) |
| Figure 4 | [`FIGURE4_SOURCE_WORKFLOW.md`](FIGURE4_SOURCE_WORKFLOW.md) | [`scripts/figure4/README.md`](scripts/figure4/README.md) | [`FIGURE4_SHA256.txt`](FIGURE4_SHA256.txt) |

From this directory, verify the tracked reader-facing source files with:

```bash
shasum -a 256 -c FIGURE1_SHA256.txt
shasum -a 256 -c FIGURE2_SHA256.txt
shasum -a 256 -c FIGURE3_SHA256.txt
shasum -a 256 -c FIGURE4_SHA256.txt
```

These public manifests cover figure-specific documentation, scripts and static
reference files. Generated outputs and execution records have their own
checksums in the local analysis folders. Only software-environment records are
bundled in this repository.

## Analysis objects and outputs

Large generated objects and result directories are not committed to GitHub.

Generated locations include:

- `objects/`
- `results/`
- `reports/`
- `logs/`

Important intermediate objects include:

- `object.list.ffpe.rds`
- `object.list.sct.ffpe.rds`
- `integration.features.rdata`
- `integration.anchors.ffpe.rds`
- `object.ffpe.rds`

## Reproducibility records

Software-version and session-information records are included under
`reproducibility_logs/`. The complete historical run records remain in the
private archive. Donor metadata, sample-level run tables and generated result
folders were not copied into this code repository.

The scripts generate local execution records, including input/output checksums,
parameters and session information. References to those output locations in the
figure guides do not imply that all generated records are bundled here.

## Software environment

The exact software environment used for the validated execution is recorded in
the session-information and installed-package files under
`reproducibility_logs/`.

These records are the authoritative description of the reference R environment.

The current workflow should first be considered a validated reference
implementation on the original Mac environment.

## Running the reference analysis

From the workflow root:

```bash
FLINC_THREADS=4 caffeinate -dimsu Rscript \
  scripts/analysis/01_preprocess_integrate_and_reduce_spatial_data.R
```

After completion of the first stage:

```bash
FLINC_THREADS=4 caffeinate -dimsu Rscript \
  scripts/analysis/02_cluster_annotate_and_characterize_spatial_domains.R
```

Figure-specific source scripts can subsequently be run from `scripts/figure1/`.

The Figure 2 workflow can subsequently be run in the order documented in
`scripts/figure2/README.md`.

The Figure 3 workflow can subsequently be run in the order documented in
`scripts/figure3/README.md`.

## Source history

`scripts/source_history/` contains historical analytical source files retained
for traceability. The reader-facing Stage 2 script uses the preserved R Markdown
analytical source internally.

## Portability

This version documents the validated reference implementation on the original
Mac analysis environment.

Cross-platform environment reconstruction and Linux portability are separate
steps. Future portability work may add environment-management or
containerization files and replace machine-specific paths without changing the
underlying scientific analysis.

## Repository layout and license

The workflow now sits directly at this repository's root. Old commands referring
to `/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code` describe the original
Mac folder. They do not require that nested structure on GitHub. Scripts retain
their original path assumptions; use the documented environment/path settings
when working on a different machine.

The migration preserved the analysis scripts and scientific reference files.
It did not rerun the biological analyses or import the old repository's Git
history. Human-participant input data and large objects are external inputs.

Jian Wu's original STRIDE scripts and copyright attribution are preserved.
See the [MIT license](LICENSE).
