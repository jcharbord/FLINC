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

## Final C4 comparison in Figure 1G

The author confirmed the manuscript panel basename
`C4_Myofibroblasts_Healthy_vs_F1_F2_F3_and_F4_wilcox`, from local directory
`260824_cluster_proportion_by_F_score_healthy_sig_C4_F1_F2_F3_F4/`.
Its matching four-test CSV and selected-comparison table are included under
`results/cluster_proportions/`. The older two-comparison C4 variant is excluded.

C4 donor proportions are compared between Healthy and each of F1–F4 with
two-sided Wilcoxon rank-sum tests (`exact = FALSE`, `correct = FALSE`). Holm
correction is applied across these four tests. Adjusted P values are 0.2971,
0.0604, 0.0604 and 0.0416, respectively. Only significant comparisons are
annotated, so the main panel shows Healthy versus F4.

Other cluster plots retain their selected nominal comparisons: C0 uses F1 versus
F4; C3 uses Healthy versus F3; C1, C2, C5, C6 and combined C4+C6 use Healthy versus
F4. Separately, `ALL_pairwise_Wilcoxon_tests_Holm_corrected.csv` contains all ten
stage pairs per cluster, with Holm correction within each ten-test family.
These are distinct adjustment families and must not be interchanged.

The eight-cluster-group summary includes the seven individual clusters plus
combined C4+C6. There are still seven spatial clusters and 24 independent donors.

## Replot from the published donor tables

The confirmed local plotting script is preserved without calculation changes as
[`06_plot_final_cluster_proportions.R`](../../scripts/figure1/06_plot_final_cluster_proportions.R).
It uses ggplot2, dplyr and scales and takes input and output directories:

```bash
Rscript scripts/figure1/06_plot_final_cluster_proportions.R \
  resources/figure1/results/cluster_proportions \
  results/Fig1_final_cluster_proportions
```

Use a new output directory. This generates plots from the supplied donor
proportions; it does not rerun clustering or marker detection. The source script
was preserved exactly; the full R plotting run was not repeated during packaging.
