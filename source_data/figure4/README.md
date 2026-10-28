# Figure 4 source tables

These exports come from local run
`24_SM1_SM6_selected4_BH30_q_thresholds_20261025T095730Z_26706`.
The generating local script 24 corresponds to repository
[`11_plot_SM1_SM6_selected4_BH30.R`](../../scripts/figure4/11_plot_SM1_SM6_selected4_BH30.R).
Its statistical calculations are unchanged by that renumbering.

| File | Contents |
| --- | --- |
| `SM1_SM6_top30_spearman_BH30.tsv` | All 60 gene–stage Spearman correlations; BH correction separately across 30 genes in SM_1 and SM_6. |
| `SM1_SM6_top30_donor_expression.tsv` | 1,440 donor/gene values: 60 genes × 24 donors. |
| `selected4_statistics_BH30_and_BH4_comparison.tsv` | The eight displayed factors with final BH30 values and BH4 values for historical comparison only. |
| `top30_gene_lists_used.tsv` | Exact 60-gene lists, module membership, expression ranking and ranks used by the analysis. |
| `donor_stage_map.tsv` | Donor codes and ordinal fibrosis stages used in the analysis. |

Log-normalized expression is averaged over all spots within each donor,
including zeros. Two-sided Spearman correlations use ordinal stage F0/Healthy=0
through F4=4 and `exact=FALSE`. Use `q_bh_30_genes` for manuscript inference;
`q_bh_4_reference` is not the final correction.

The displayed factors are TIMP1, IGFBP7, ADAMTSL2 and DCN in SM_1, and AHSG,
FN1, ITIH2 and SERPINC1 in SM_6. All four SM_1 factors have BH30 q < 0.001;
all four SM_6 factors have BH30 q < 0.01. Exact values are retained in the tables.

These five tables are byte-identical between local runs 23 and 24. The later
run changes q-value label formatting only. The exact 60-gene list here differs
in scope and schema from the older cross-module reference-ranking CSV under
[`resources/figure4/`](../../resources/figure4/); the two files serve different roles.

Publishing these exports does not supply all original object and historical
validation inputs required by the guarded rerun. See the
[Figure 4 guide](../../scripts/figure4/README.md#validation-inputs-required-for-reproducibility).
PROGENy, DoRothEA, RFX5 and additional spatial/cluster numerical exports are not
included in this source-data release.

See the [source-data manifest](../manifest.tsv) for original paths and checksums.
