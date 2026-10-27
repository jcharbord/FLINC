# Figure 1 source scripts

This folder contains scripts used to generate the current Figure 1 source
visualizations and validate their provenance.

The scripts are organized by scientific content rather than by fixed panel
letters because the final figure layout may change.

Reader-facing companion files at the reproducibility root are:

- [`FIGURE1_SOURCE_WORKFLOW.md`](../../FIGURE1_SOURCE_WORKFLOW.md), which maps
  the analytical and visualization sources
- [`FIGURE1_SHA256.txt`](../../FIGURE1_SHA256.txt), which verifies the public
  Figure 1 documentation and scripts

## 01_plot_spatial_clusters_on_histology.R

Generates spatial cluster overlays using registered Space Ranger
high-resolution H&E images and true Visium spot geometry.

## 02_plot_top_cluster_markers_spatially.R

Generates spatial expression maps of top markers for the final spatial
clusters in donor 3D.

Expression source:

- assay: SCT
- layer: data

Visualization uses registered high-resolution histology, true Visium spot
geometry and gene-specific display scaling.

## 03_plot_final_annotated_umap.R

Generates the final annotated UMAP using the existing UMAP coordinates and
final cluster assignments.

It does not recalculate dimensionality reduction or clustering.

## 04_prepare_selected_histology_crop.R

Generates the representative donor-3D material used for Figure 1 presentation.

The workflow includes:

- the fixed donor-3D presentation region
- registered Space Ranger high-resolution H&E
- cluster overlays with both the final C5-red palette and historical C5-yellow palette
- spatial expression maps for the top 10 markers of each final cluster
  (70 marker maps in total)
- true Visium capture-spot geometry
- gene-specific expression display scaling
- run records, source-object checksums and output manifests

The registered Space Ranger H&E source is also copied to the output directory
so that the final manuscript H&E crop can be prepared manually in Adobe
Photoshop from the same registered image used for the spatial plots.

Cropping and plotting do not modify expression values, cluster identities,
spatial coordinates, spot geometry or the source Seurat object.

## 05_validate_figure1_outputs.R

Checks Figure 1 analytical and plotting sources, critical inputs,
software/session records, checksums, output counts and palette consistency.

It does not rerun the biological analysis.

## Analytical results used directly in Figure 1

Several Figure 1 results are produced directly by:

`../analysis/02_cluster_annotate_and_characterize_spatial_domains.R`

These include:

- top cluster marker summaries
- Gene Ontology enrichment
- cell-type deconvolution summaries
- cluster-cell-type chord relationships
- stage-associated cluster proportion analyses

These analyses are not duplicated here because they belong to the main
biological analysis.
