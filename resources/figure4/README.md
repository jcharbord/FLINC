# Figure 4 resources

| Location | Contents |
| --- | --- |
| [reference_inputs/](reference_inputs/) | UniProt secreted-protein snapshot, frozen secreted-factor references, the exported PROGENy top-100 model/network and RFX5 resource targets. |
| [results/](results/) | PROGENy, DoRothEA, RFX5 and secreted-factor numerical results. |
| [results/spatial_context/](results/spatial_context/) | v6 secreted-factor donor/cluster summaries, donor-balanced heatmap values and exported map data. |

## Regulatory analyses

`01_progeny_contrast_statistics.tsv` contains 56 tests (14 pathways × four stages
versus Healthy), with BH correction across 14 pathways within each contrast.
`05_dorothea_all_contrasts.tsv` contains 472 tests (118 TFs × four contrasts),
with BH correction across 118 TFs within each contrast. Their exported donor
medians, mean differences, rank-sum tests and BH values were checked independently.
The full DoRothEA AB resource network and complete per-spot activity matrices
remain external; coverage summaries and the RFX5 target subset are included.

The `07_rfx5_*` result tables contain target donor medians, observed effects and
displayed regulatory effects. Regulatory effects were checked against the
donor summaries and mode-of-regulation weights. `08_rfx5_celltype_spearman_by_stage.tsv`
contains descriptive within-stage spot-level correlations, not donor-level tests.

## Final secreted-factor stage profiles

`SM1_SM6_top30_spearman_BH30.tsv` contains all 60 tests, and
`SM1_SM6_top30_donor_expression.tsv` contains 1,440 donor/gene values. Expression
is averaged over all spots within each donor, including zeros. Two-sided
Spearman tests use ordinal stage Healthy/F0=0 through F4=4 and `exact=FALSE`.
BH correction is applied separately across 30 genes in each module.

Use `q_bh_30_genes` for inference. `q_bh_4_reference` in the selected-factor
comparison table is historical only. The four displayed SM_1 factors are TIMP1,
IGFBP7, ADAMTSL2 and DCN (all BH30 q < 0.001); the four SM_6 factors are AHSG,
FN1, ITIH2 and SERPINC1 (all BH30 q < 0.01). Exact values are unchanged.

The final run24 tables match run23 byte-for-byte; run24 changed label formatting
only. Local script 24 is repository script 11. `top30_gene_lists_used.tsv` contains
the exact 60-gene statistical input; the older reference-ranking CSV under
`reference_inputs/` covers multiple modules and has a different purpose/schema.
`12_top30_strict_secreted_factors_by_module.tsv` matches the input SHA recorded
by the final BH30 run.

## Spatial-cluster summaries

The v6 spatial-context export contains 16 selected genes. Donor-balanced cluster
means were independently checked, and 384 shared donor/gene expression values
agree with the final BH30 bundle. Its historical BH4 stage-correlation file is
excluded in favour of the final BH30 results. The spatial tables are descriptive;
their presence does not change the final stage-profile correction.

Saved historical validation records and original expression objects are still
needed for the guarded full rerun. See the
[Figure 4 guide](../../scripts/figure4/README.md#validation-inputs-required-for-reproducibility).
