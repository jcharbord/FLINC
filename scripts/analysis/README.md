# Main spatial transcriptomics analysis

This folder contains the two principal analysis scripts.

They should be read and run in numerical order.

## 01_preprocess_integrate_and_reduce_spatial_data.R

Performs the initial spatial transcriptomics processing:

- loading of Space Ranger outputs for 24 FFPE samples
- creation of Seurat spatial objects
- SCTransform normalization
- SCT integration
- import of previously generated cell-type deconvolution results
- PCA
- Harmony correction
- UMAP
- saving of intermediate and integrated analysis objects
- recording of parameters, software versions and checksums

## 02_cluster_annotate_and_characterize_spatial_domains.R

Performs the principal biological analysis of the integrated spatial dataset:

- graph-based clustering
- hierarchical subclustering
- final biological annotation
- cluster marker identification
- Gene Ontology enrichment
- cell-type deconvolution summaries
- cluster-cell-type chord relationships
- stage-associated cluster proportion analyses
- generation of final annotated objects
- generation of analytical results used in Figure 1

The final annotated spatial domains are:

- C0 Portal tracts
- C1 Periportal hepatocytes
- C2 Centrilobular hepatocytes
- C3 Immune cells
- C4 Myofibroblasts
- C5 Cholangiocytes
- C6 Matrix-rich fibroblasts

The second script requires the integrated object produced by the first script.

Historical preparation files are retained separately in:

`../source_history/`
