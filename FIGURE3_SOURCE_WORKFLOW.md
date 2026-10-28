# Figure 3 source workflow

This document identifies the retained analysis and visualization sources used
to generate Figure 3 and Supplementary Figure 8. The workflow is described by
scientific content rather than by historical notebook names.

## LIANA cell-cell communication inference

Source:

`scripts/figure3/01_run_liana_inference.R`

This operation reconstructs the 12-cell-type by 24-donor weighted pseudobulk
input from the validated final annotated Figure 1 object, applies
LogNormalization and runs LIANA with OmniPath. Results from CellChat, NATMI,
SCA, logFC and Connectome are aggregated, and interactions with
`aggregate_rank <= 0.01` are retained.

## Stage-progressive interaction activity

Source:

`scripts/figure3/02_analyze_stage_activity.R`

This operation calculates donor- and stage-level activity for the top 3,000
single-gene LR candidates, performs the 2,903 finite Spearman trend tests,
applies Benjamini-Hochberg correction and retains the 60 interactions used in
the stage-activity heatmap and supplementary dot plots. It also prepares the
significant-interaction chord counts and cholangiocyte incoming/outgoing
partner counts.

## Figure 3A-B and Supplementary Figure 8A-E

Source:

`scripts/figure3/03_plot_main_and_supplementary_ccc.R`

This deterministic plotting operation generates:

- the Figure 3A significant-interaction chord diagram
- the Figure 3B stage-activity row-z heatmap and annotation table
- the Supplementary Figure 8A cholangiocyte communication-hub counts
- the Supplementary Figure 8B-E F1, F2, F3 and F4 versus F0 dot plots

It reads only the retained Step 02 tables and performs no new biological
analysis.

## Figure 3C all-donor spatial LR maps

Source:

`scripts/figure3/04_plot_spatial_lr_all_donors.R`

This operation calculates and plots cell-type-abundance-weighted local LR
scores for NPNT-ITGA8 (Cholangiocyte to Stellate), MMP7-CD151
(Cholangiocyte to InfMac) and CCL2-ACKR1 (Stellate to PortalEndo) in all 24
donors. It generates 72 LRlocal PDFs, arranged by stage, LR pair and direction.

The score uses normalized ligand signal at the focal spot and a
Gaussian-weighted receptor signal in that spot plus its 30 nearest neighbours
(`sigma = 2`), followed by `log1p(score * 10000)`. Each pair uses a common
first-to-99th-percentile scale across all donors. The maps are descriptive and
are not used for statistical testing.

The custom H&E images are combined with full-resolution Space Ranger spot
coordinates. The true Visium spot diameter is derived from the recorded Space
Ranger scale factors and drawn with multiplier 1.00. Geometry and PDF
checksums are recorded for every donor and interaction.

## Fixed Figure 3C manuscript selection

Source:

`scripts/figure3/05_select_main_figure3C_panels.R`

This operation verifies and copies the six accepted Figure 3C PDFs:

- F0: donor 6B NPNT-ITGA8, donor 2A MMP7-CD151 and donor 4C CCL2-ACKR1
- F4: donor 4B NPNT-ITGA8, donor 3A MMP7-CD151 and donor 2D CCL2-ACKR1

The accepted donor/pair composition was fixed before visual inspection of the
regenerated maps. The operation performs no biological analysis, verifies the
three exact LIANA aggregate ranks and shared color limits, and creates
byte-identical copies in manuscript order.

## Final Figure 3 audit

Validation source:

`scripts/figure3/06_validate_figure3_outputs.R`

This checks the public source files, all Step 01-05 completion records, expected
dimensions and statistical-test counts, donor/stage structure, plot manifests,
true-Visium geometry, fixed panel identities, input and output SHA-256
manifests, software information and absence of incomplete or excluded outputs.
It performs no additional biological analysis.

## Generated outputs and records

Generated results are written under `results/Fig3_CCC/` and are not committed
to GitHub. Compact verified execution records are copied by the final audit to
`reproducibility_logs/Fig3_CCC/Fig3_final_audit/`.

Historical notebooks, historical result folders and internal equivalence
checks are provenance material and are not part of the public reproducible
workflow.

## Integrity manifest

`FIGURE3_SHA256.txt` records the SHA-256 checksums of the public Figure 3
documentation and scripts. From the reproducibility root, run:

```bash
shasum -a 256 -c FIGURE3_SHA256.txt
```

## Repository snapshot note

Only software-version and session-information records are bundled under
`reproducibility_logs/`. Other run records and generated outputs described above
are local analysis artifacts and are not included in this checkout. See the
[root README](README.md) for migration scope and required external inputs.

## Published numerical source tables

[resources/figure3/](resources/figure3/README.md) contains the
selected completed-run tables and available reference inputs for this figure.
The [resource manifest](resources/manifest.tsv) records exact original paths,
checksums and verification details. The figure resource guide distinguishes
final tests, descriptive summaries and external inputs.
