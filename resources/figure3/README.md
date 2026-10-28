# Figure 3 resources

The [results directory](results/) contains retained LIANA interactions and the
numerical tables used for stage-activity, chord, partner-count and dotplot panels.

| Files | Purpose |
| --- | --- |
| `01_liana_significant_rank_le_0.01.csv.gz` | All 29,205 interactions retained at LIANA aggregate rank <= 0.01. |
| `01_pseudobulk_metadata.csv`, `01_donor_stage_map.csv` | Metadata for the donor/cell-type pseudobulks and stage mapping. |
| `02_stage_activity.csv` | All 3,000 tested candidates, their stage summaries, Spearman statistics and BH values. |
| `02_selected_interactions_top60.csv`, `02_selected_interactions_row_z_scores.csv` | Selected interaction list and heatmap values. |
| `02_chord_*`, `02_cholangiocyte_*`, `02_dotplot_input_*` | Chord, cholangiocyte-partner and stage-comparison plot inputs. |
| [spatial_panels/](results/spatial_panels/) | Spatial interaction definitions, colour limits and final panel-selection manifests. |

The retained LIANA table exactly matches filtering the supplied 609,984-row full
aggregate table. That approximately 99 MB uncompressed full table remains
external; the retained table is supplied once as lossless gzip.

Stage trends are Spearman correlations across the five stage medians, not tests
across individual spots or independent donor observations. Of 3,000 candidates,
2,903 yield finite tests and 349 have BH q <= 0.05. Selection additionally requires
an absolute F4-minus-F0 difference >= 0.5, then retains the top 60 by that
difference. The exported correlations and BH values were independently checked.

Spatial panel manifests describe the selected plots but do not contain the
original images, RDS objects or complete per-spot spatial scores. Figure images
and the frozen OmniPath input network are not included here. Refer to the
[Figure 3 guide](../../scripts/figure3/README.md) for the workflow and external inputs.
