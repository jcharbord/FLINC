# Figure 1 resources

The [results directory](results/) contains the source tables for cluster
characterization and donor-level abundance comparisons.

| Files | Purpose |
| --- | --- |
| [cluster_marker_final_clusters.csv](results/cluster_marker_final_clusters.csv) | All 13,306 positive marker rows across seven clusters. |
| [cluster_marker_final_clusters_top20_per_cluster.csv](results/cluster_marker_final_clusters_top20_per_cluster.csv) | Top 20 markers per cluster. |
| [Top10_markers_per_final_cluster_dotplot_data.csv](results/Top10_markers_per_final_cluster_dotplot_data.csv) | Numerical input for the marker dotplot. |
| `GO_BP_all_terms_final_clusters.csv`, `GO_BP_top_terms_final_clusters.csv` | GO Biological Process results and displayed terms. |
| `deconvolution_*_FIXED.csv`, `chordDiagram_*_data.csv` | Cluster-level cell-type summaries and exported chord data. |
| `final_cluster*_counts.csv`, `approach1_raw_to_final_cluster_mapping.csv` | Final cluster sizes and annotation mapping. |
| [cluster_proportions/](results/cluster_proportions/) | Donor proportions, stage mean/SEM and statistical comparisons. |

The marker/GO/deconvolution exports come from the Stage 2 analysis script.
The deconvolution tables retain the original cell-type categories; downstream
display grouping (for example, pooling immune populations) is a plotting step.
Not every later display-specific regrouped table is included.

## Cluster-proportion comparisons

Donor-level cluster proportions were compared using two-sided Wilcoxon rank-sum
tests. For C4 (Figure 1G), Healthy was compared with F1–F4, with Holm correction
across the four tests. Other displayed cluster comparisons use nominal P-values.
The separate all-pairwise table includes Holm-adjusted P-values within each cluster.

## Replot from the published donor tables

Generate the cluster-proportion plots with
[`06_plot_final_cluster_proportions.R`](../../scripts/figure1/06_plot_final_cluster_proportions.R).
It uses ggplot2, dplyr and scales and takes input and output directories:

```bash
Rscript scripts/figure1/06_plot_final_cluster_proportions.R \
  resources/figure1/results/cluster_proportions \
  results/Fig1_final_cluster_proportions
```

Use a new output directory. The script reads the supplied donor proportions.
