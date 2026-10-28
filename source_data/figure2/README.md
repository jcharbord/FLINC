# Figure 2 source tables

## Differential module eigengenes

| File | Contents |
| --- | --- |
| `03_DiffModuleEigengenes_donor_level.csv` | Final 24 comparisons: six modules for each of F1–F4 versus F0/Healthy. |
| `03_donor_mean_module_eigengenes.csv` | Six module eigengene means for each of 24 independent donors. |
| `03_QC_spots_per_donor_and_module.csv` | Spot counts contributing to donor/module summaries. |
| `03_QC_donors_per_stage.csv` | Donor and spot counts by stage. |
| `03_module_assignments_used.csv` | Gene-to-module assignments used by the DME run. |

Spot-level harmonized module eigengenes were averaged within each donor.
Stage-versus-Healthy tests use independent donors and BH adjustment across six
modules separately within each contrast. Donor counts are F0=6, F1=6, F2=3,
F3=3 and F4=6. `delta_mean` is stage minus Healthy. `p_adj` is the final BH6
value. These are the donor-level `03_bis` results, not the historical spot-level
`03` analysis. The October 15 bigger-stars replot retained identical statistics.

## Module–cell-type associations

| File | Contents |
| --- | --- |
| `07_hME_celltype_spearman_long.csv` | All 72 spot-level correlations (six modules × 12 cell types), with raw P values, global adjustment and within-module adjustment. |
| `07_hME_celltype_spearman_matrix.csv` | Correlation matrix used by the heatmap. |
| `07_strongest_positive_driver_by_module.csv` | Exported strongest positive association per module. |

These are descriptive spot-level associations across 26,056 spots, not
independent-donor tests. The October 15 heatmap without annotation bars used
the identical correlation matrix. The historical word `driver` in the filename
does not establish causality.

## Hub lists and PCSK7/PNPLA3 neighbourhoods

| File | Contents |
| --- | --- |
| `06_hubs_with_functional_subcategories.tsv` | 200 kME-selected hubs: 50 for SM_1 and 30 each for SM_2–SM_6, with ranks and curated categories. |
| `06_functional_subcategory_legend.tsv` | Category labels, order and colours. |
| `06_functional_subcategory_counts.tsv` | Category counts in the selected hub lists. |
| `11_ego_top20_neighbors.tsv` | Ranked top 20 same-module neighbours for each of PCSK7 and PNPLA3, with dissimilarity and similarity. |
| `11_ego_network_nodes.tsv` | The two index genes and their neighbours with functional-context annotations. |
| `11_ego_network_edges.tsv` | Index-to-neighbour edges weighted by TOM similarity. |

Hub selection is based on kME. PCSK7/PNPLA3 neighbour selection is instead based
on descending TOM similarity within the index gene's module (PCSK7: SM_1;
PNPLA3: SM_3). These are distinct selections. The `11_` tables are copied from
the extended-category run; their ranks and numerical weights match run 09.
Legacy F/G panel identifiers in these exports reflect their original Figure 4
location; the current analysis is Figure 2 script `05_bis`.

The six-module circular-network edge weights are not included in this bundle.
No network calculation or plot was changed in preparing these source tables.
Frozen annotation inputs are under [`resources/figure2/`](../../resources/figure2/).

See the [source-data manifest](../manifest.tsv) for exact original paths and checksums.
