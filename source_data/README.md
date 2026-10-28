# Numerical source data

Selected numerical exports from completed FLINC analyses are provided here as
CSV/TSV files. They were copied without changing their contents from the local
analysis bundle `FLINC_review_20261028_185935.zip` on 28 October 2026.

| Figure | Included source data | Coverage still incomplete |
| --- | --- | --- |
| [Figure 2](figure2/README.md) | Donor-level DME and QC, module assignments, module–cell-type correlations, annotated hub lists, PCSK7/PNPLA3 neighbours and edges. | Module enrichment, supplementary GSEA and other spatial/plotting exports. |
| [Figure 4](figure4/README.md) | SM_1/SM_6 donor expression, all 60 Spearman tests with BH30, selected-factor statistics and gene lists. | PROGENy, DoRothEA, RFX5 and additional spatial/cluster summaries. |

Figure 1 and Figure 3 numerical exports and manuscript figure PDFs are not yet
included. This is a partial source-data release, not a complete set of outputs
for every manuscript panel.

## Provenance and verification

[`manifest.tsv`](manifest.tsv) records each table's original relative path, row
count, SHA-256, current repository script and checks performed. Original paths
are relative to the local `Reproducible_code/` directory. Local script numbers
and folder names can differ from their public repository counterparts.

From the repository root, verify this source-data bundle with:

```bash
shasum -a 256 -c source_data/SHA256SUMS
```

The donor-level DME table was identical in the original run and final bigger-stars
replot. Its available original MD5 records passed. Rank-sum tests, mean differences
and BH correction were independently checked from the exported donor means.
The final correlation heatmap reused the same matrix as the original export.
PCSK7/PNPLA3 tables passed their original SHA-256 checks; neighbourhoods and
weights agree between runs 09 and 11. Secreted-factor tables in runs 23 and 24
were identical, and Spearman correlations and BH30 were independently checked
from exported donor expression. Numerical differences in the independent checks
were below 1e-12.

These checks did not rerun preprocessing, WGCNA, LIANA or regulatory-activity
analysis. Original RDS/TOM inputs and manuscript PDFs were not supplied in this
bundle. Donor codes are the original analytical sample identifiers; no additional
clinical metadata has been added. Required external inputs and saved validation
records remain described in the figure guides.
