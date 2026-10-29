# Figure resources and source tables

Start with your figure below. Each figure has one directory containing its
numerical results and, where available, the frozen reference inputs used by its
scripts. There is no separate top-level source-data directory.

| Figure | Contents |
| --- | --- |
| [Figure 1](figure1/README.md) | Marker genes, GO enrichment, deconvolution, cluster counts, final donor-level cluster-proportion comparisons and pooled/per-sample spatial QC summaries. |
| [Figure 2](figure2/README.md) | Module assignments and kME, DME, hub and neighbour lists, cell-type correlations, module enrichment and donor-pseudobulk GSEA. |
| [Figure 3](figure3/README.md) | Retained LIANA interactions, stage-activity tests, selected interactions, chord/dotplot data and spatial-panel selections. |
| [Figure 4](figure4/README.md) | PROGENy, DoRothEA, RFX5, secreted-factor rankings, donor expression, final BH30 tests and spatial-cluster summaries. |

Within each figure, `results/` contains exports from completed analyses.
`reference_inputs/`, when present, contains curated annotations or frozen
database/model references used by the scripts. Similar names do not imply
identical roles: for example, the UniProt input list and the final gene–stage
correlation results are different tables.

## Provenance and integrity

[manifest.tsv](manifest.tsv) records the 135 selected numerical/model exports,
their original run paths, row counts, generating scripts and checks performed.
The six pre-existing curated/reference inputs are retained unchanged and covered
by the figure manifests and the checksum file below. The original 131 exports
were supplied in `FLINC_review_20261028_185935.zip` and
`FLINC_remaining_tables_20261028_192056.zip`. Three pooled spatial QC tables
and one per-sample technical QC table were exported separately on 29 October 2026 and checked against the matching
spot-level QC metadata. Their source and summary definitions are described in
[the QC guide](figure1/results/qc/README.md).
Original paths are relative to the local `Reproducible_code/` directory.

No table values were edited. Selected larger text tables were gzip-compressed;
their decompressed bytes match the originals. Files already supplied as gzip
were copied byte-for-byte. The manifest records both original and stored hashes.
From the repository root:

```bash
shasum -a 256 -c resources/SHA256SUMS
```

All 141 available original checksum comparisons in the second bundle passed.
Independent checks also reproduced the donor-level DME and regulatory tests,
cluster-proportion tests, final secreted-factor Spearman/BH30 statistics and
several derived summaries. GSEA tables match the combined export, including
BH correction within each contrast/collection. Files lacking supplied historical
checksums are explicitly identified in the manifest.

These checks do not replace a full rerun from original expression objects.
Large Seurat/TOM objects, original histology, complete per-spot activity matrices,
the full unfiltered LIANA aggregate table, the MSigDB term-to-gene snapshot and
assembled manuscript figures are not included. The figure guides describe
external dependencies and remaining plot-specific inputs. Donor codes are the
original analytical identifiers; no additional clinical metadata was added.

