# Figure 1 source workflow

This document identifies the analytical and visualization sources used to
generate the current Figure 1 material.

The workflow is described by scientific content rather than by fixed panel
letters because the final figure layout may change.

## Spatial clustering and biological characterization

Primary analytical source:

`scripts/analysis/02_cluster_annotate_and_characterize_spatial_domains.R`

This analysis provides:

- spatial clustering and hierarchical subclustering
- final biological cluster annotation
- top cluster marker identification
- Gene Ontology enrichment
- cell-type deconvolution summaries
- cluster-cell-type chord relationships
- stage-associated cluster proportion analyses
- final annotated analysis objects

The final annotated spatial domains are:

- C0 Portal tracts
- C1 Periportal hepatocytes
- C2 Centrilobular hepatocytes
- C3 Immune cells
- C4 Myofibroblasts
- C5 Cholangiocytes
- C6 Matrix-rich fibroblasts

## Final annotated UMAP

Source:

`scripts/figure1/03_plot_final_annotated_umap.R`

The final UMAP uses the existing `umap_SCT` coordinates and final cluster
annotations.

It does not recalculate dimensionality reduction or clustering.

Final cluster palette:

- C0 `#56B4E9`
- C1 `#FF9DA7`
- C2 `#B189A8`
- C3 `#59A14F`
- C4 `#F28E2B`
- C5 `#E15759`
- C6 `#E5C494`

## Spatial cluster overlays on histology

Source:

`scripts/figure1/01_plot_spatial_clusters_on_histology.R`

This generates spatial cluster overlays using registered Space Ranger
high-resolution H&E images and true Visium capture-spot geometry.

## Spatial maps of top cluster markers

Source:

`scripts/figure1/02_plot_top_cluster_markers_spatially.R`

This generates spatial expression maps of the top markers of each final
cluster in donor 3D.

Expression source:

- assay: `SCT`
- layer: `data`

Visualization uses:

- registered high-resolution H&E
- true Visium spot geometry
- gene-specific display scaling using the 99th percentile

## Representative donor-3D Figure 1 material

Supporting visualization source:

`scripts/figure1/04_prepare_selected_histology_crop.R`

This workflow generates the representative donor-3D presentation material
using the final annotated spatial object and the registered Space Ranger
high-resolution H&E image.

It provides:

- the fixed presentation crop
- the registered H&E source used for spatial alignment
- C5-red and historical C5-yellow cluster-overlay versions
- the top 10 spatial marker maps for each of C0-C6
- true Visium capture-spot geometry
- gene-specific expression scaling at the 99th percentile
- source-object and image checksums
- reproducibility records and output manifests

The registered H&E source is also copied into the output structure so that
the final manuscript H&E panel may be cropped manually in Adobe Photoshop
from the same registered image.

The presentation workflow does not recalculate clustering, expression,
deconvolution or biological annotation and does not modify the source
Seurat object.

## Supplementary spatial QC metrics

Source: `scripts/figure1/07_summarize_spatial_QC.R`.

This function uses `nFeature_Spatial`, `nCount_Spatial` and `donor` metadata in
the final annotated Seurat object to summarize 26,056 spots from 24 biopsies.
The manuscript values, 4,845 detected genes and 20,064 UMIs per spot, are pooled
medians. The function also calculates means and alternative sample-weighted
summaries and generates the Supplementary Figure 1D-E violin plots.

The [published QC tables](resources/figure1/results/qc/README.md) are cohort-level
summaries. Sample-specific and spot-level tables and rendered plots remain
local. The script neither modifies the source object nor applies additional
filtering, normalization or clustering. Run versions are recorded in
`reproducibility_logs/Fig1_QC_20261029/sessionInfo.txt`.

## Figure 1 output validation

Validation source:

`scripts/figure1/05_validate_figure1_outputs.R`

This checks:

- expected analysis and figure source files
- critical analytical inputs
- run records and session information
- checksums
- expected output counts
- cluster-palette consistency
- supporting selected-crop material

The validation step does not rerun the biological analysis.

## Integrity manifest

`FIGURE1_SHA256.txt` records the SHA-256 checksums of the public Figure 1
documentation, main Stage 2 analysis and figure-specific scripts. From the reproducibility root, run:

```bash
shasum -a 256 -c FIGURE1_SHA256.txt
```

## Published numerical source tables

[resources/figure1/](resources/figure1/README.md) contains the
selected completed-run tables and available reference inputs for this figure.
The [resource manifest](resources/manifest.tsv) records exact original paths,
checksums and verification details. The figure resource guide distinguishes
final tests, descriptive summaries and external inputs.

The final Figure 1G C4 plot uses Healthy versus F1–F4, with Holm correction across
four tests and only significant comparisons annotated. Other cluster plots
retain selected nominal comparisons. The separate all-pairwise table uses
Holm correction across ten pairs within each cluster. The final plotter is
`scripts/figure1/06_plot_final_cluster_proportions.R`; it reads the published
per-cluster donor proportions and does not rerun clustering.

