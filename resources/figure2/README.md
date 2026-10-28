# Figure 2 resources

| Location | Contents |
| --- | --- |
| [reference_inputs/](reference_inputs/) | Functional categories for kME-selected hubs and PCSK7/PNPLA3 neighbourhood annotations. |
| [results/](results/) | Module assignments, kME, final donor-level DME, module–cell-type correlations, hubs, neighbours and module-gene UMAP coordinates. |
| [results/module_enrichment/](results/module_enrichment/) | Ranked module gene lists, module metadata and archived Enrichr results. |
| [results/GSEA/](results/GSEA/) | Independent donor-pseudobulk GSEA, limma-voom differential-expression tables, rankings and database-version records. |

## Differential module eigengenes and cell-type correlations

`03_DiffModuleEigengenes_donor_level.csv` contains the final 24 tests: six modules
for each of F1–F4 versus Healthy/F0. Spot-level harmonized eigengenes were averaged
within donors; independent-donor Wilcoxon tests use BH adjustment across six
modules separately within each contrast. `p_adj` is the final BH6 value and
`delta_mean` is stage minus Healthy. Donor means and QC counts are supplied.
These are the `03_bis` results, not the historical spot-level `03` analysis.
The final bigger-stars replot retained exactly the same statistical table.

`02_module_assignments_and_kME.tsv` includes all 11,515 network genes, including
9,309 grey/unassigned genes. `03_module_assignments_used.csv` contains only the
2,206 genes assigned to SM_1–SM_6; their assignments agree exactly.

`07_hME_celltype_spearman_long.csv` contains 72 descriptive spot-level
correlations across 26,056 spots, with raw P values and both global and
within-module corrections. These are not independent-donor tests. The matrix
matches the final heatmap without annotation bars. The historical word `driver`
in the strongest-association filename does not establish causality.

## Hub genes and neighbourhoods

`06_hubs_with_functional_subcategories.tsv` contains 200 kME-selected hubs:
50 in SM_1 and 30 each in SM_2–SM_6. Their categories match the frozen reference.
The `11_ego_*` tables instead contain top-20 same-module neighbours for PCSK7
(SM_1) and PNPLA3 (SM_3), ranked by TOM similarity. Their original checksums pass;
run 11 preserves run 09 ranks/weights and adds extended functional categories.
Legacy F/G identifiers reflect the previous Figure 4 placement; the current
script is Figure 2 `05_bis`. The two selection procedures are distinct.
Circular-network edge weights are not part of this release. No network
calculation or plot was changed during packaging.

## Enrichment and supplementary GSEA

Archived Enrichr outputs are supplied without querying live databases again.
The GSEA tables come from `run_20260814T185515Z`, using MSigDB 2025.1.Hs and
collections H, C2, C3 and C5. All 16 full contrast/collection tables are supplied
once under `GSEA/gsea_full/`; they jointly contain 67,500 rows. The combined
and significant-only duplicate exports are omitted. Select `p.adjust < 0.05`
to recover the 11,537 significant rows. BH correction is within each
contrast/collection. GSEA is based on donor-pseudobulk limma-voom rankings,
independently of WGCNA.

The full MSigDB term-to-gene snapshot and spatial plotting matrices are not
included. Model objects and original spatial inputs remain external.
See the [manifest](../manifest.tsv) for exact provenance and compression details.
