# Figure 4 scripts

This directory contains the final public scripts for Figure 4. Steps 01–05 and 07–10 use explicit command-line paths. The later manuscript update, script 11, uses the environment variables and saved validation inputs described below. No large input object or generated output is stored in Git.

| Script | Required upstream material | Main output |
|---|---|---|
| `01_run_progeny_activity.R` | Final Seurat object | PROGENy activities, donor summaries, F4/F0 and validation contrasts. |
| `02_plot_progeny_spatial_maps.R` | Final Seurat object; output from step 01 | Spatial pathway maps. |
| `03_run_dorothea_activity.R` | Final Seurat object | DoRothEA TF activities and contrasts. |
| `04_build_rfx5_regulon.R` | Final Seurat object; output from step 03 | RFX5 observed-target regulon. |
| `05_rfx5_celltype_associations.R` | Final Seurat object; output from step 03 | RFX5/deconvolution-fraction associations. |
| `07_run_secreted_factor_analysis.R` | Final Seurat object, module table, versioned UniProt snapshot, reference count/ranking tables | Strict secreted-factor summaries and rankings. |
| `08_plot_secreted_enrichment.R` | Versioned enrichment reference table | Supplementary enrichment panel. |
| `09_validate_figure4_outputs.R` | This script directory and the Figure 4 resource directory | Public-file checksum manifest. |
| `10_plot_secreted_factor_spatial_context_v6.R` | Final Seurat object, strict top-30 table and expected input SHA-256 | Spatial-cluster summaries, maps and historical SM_1 BH4 profiles. |
| `11_plot_SM1_SM6_selected4_BH30.R` | Validated v6 manifest, complete top-30 table, earlier SM_1 donor export and SM_6 top-30 validation exports | Final Fig. 4H/I donor-level profiles with BH30 correction. |

The former step 06 and its functional-context annotation table now belong to
[the Figure 2 ego-network workflow](../figure2/README.md#5-bis-05_bis_build_ego_networksr).
The remaining Figure 4 script numbers are preserved.

## Important interpretation safeguards

- DoRothEA and PROGENy infer pathway/TF activities from expression targets; they do not measure protein activity directly.
- RFX5/fraction correlations are descriptive, within-stage spot-level associations.
- The secreted-factor enrichment panel uses a versioned reference table, not a live GO/KEGG query. This keeps the rendered database mappings stable across reruns.

## Validation

After installing the public workflow, run:

```bash
Rscript scripts/figure4/09_validate_figure4_outputs.R \
  --script-dir scripts/figure4 \
  --resource-dir resources/figure4 \
  --manifest FIGURE4_SHA256.txt \
  --record-dir reproducibility_records/Fig4
```

The validator writes a checksum manifest and confirms that the expected final scripts and small resources are present. It does not access or alter the large local result directories.


## Secreted-factor donor-level stage profiles

Script 11 is the successful final SM_1/SM_6 panel script, renumbered from the
local script 24. Only its filename, script-number wording and new-output folder
prefix changed; the calculations and plotting settings are unchanged.
It supersedes the four-gene BH correction and stage-profile figure produced by
script 10; script 10 remains the source for spatial-cluster summaries and maps.

- SM_1: TIMP1, IGFBP7, ADAMTSL2, DCN.
- SM_6: AHSG, FN1, ITIH2, SERPINC1.
- Normalize raw Spatial counts as `log1p(count / spot library size * 10000)`.
- Average across all spots within each donor, including zero expression.
- Test two-sided Spearman associations across the 24 independent donors using
  ordinal stage 0–4 and `exact = FALSE`.
- Apply BH correction separately across all 30 tested genes in each module.
  The four displayed genes retain these BH30 q values; BH4 is exported only
  as a historical comparison.
- Produce matched 2×2 panels with Healthy labels, the established stage colors,
  and both names-only and rho/q variants. Labels show `q < 0.001` or `q < 0.01`
  only when the numeric q value meets that threshold. Exact values stay in tables.

### Run against the existing validated outputs on Jian's Mac

```bash
conda activate hdWGCNA_pipeline-env
export FLINC_REPRO_ROOT="/Users/jianwu/Documents/FLINC/Jeremie/Reproducible_code"
cd "$FLINC_REPRO_ROOT"
Rscript --vanilla scripts/figure4/11_plot_SM1_SM6_selected4_BH30.R
```

The script locates the matching v6 `run_manifest.tsv`, verifies its recorded
object and top-gene checksums, and recovers the complete top-30 export when the
manifest points to the smaller selected-gene audit. The validated full table was
`results/Fig4_regulatory_and_module_networks/12_secreted_factor_analysis/12_top30_strict_secreted_factors_by_module.tsv`.

### Validation inputs required for reproducibility

This is a guarded rerun of an existing analysis, not a standalone fresh-run
pipeline. In addition to the Seurat object and top-30 lists it requires:

- a matching v6 `run_manifest.tsv` with input/top-gene paths and SHA-256 values;
- `SM6_top30_donor_stage_expression.tsv` and
  `SM6_top30_stage_spearman_BH30.tsv` from the earlier
  `20_SM6_top30_stage_profiles_*` result folder;
- a matching `selected_secreted_factor_donor_stage_statistics.tsv` with the
  four SM_1 donor profiles (96 rows), under the Figure 4 result tree.

The historical SM_6 validation tables are not generated by this repository's
numbered steps and are not included in this code update. A fresh checkout needs
these saved artifacts as well as the large input object; their availability in
any external reproducibility archive must be checked separately. Do not disable
the concordance checks simply to make a fresh checkout run.

If auto-discovery is ambiguous, set `FLINC_SECRETED_MANIFEST`,
`FLINC_SECRETED_TOP30`, and `FLINC_SM6_SOURCE_DIR` to the matching files/directory.
On another computer, manifest paths must resolve to the same checksummed files.

Each successful run creates a new
`results/Fig4_regulatory_and_module_networks/11_SM1_SM6_selected4_BH30_q_thresholds_*`
folder containing PDF/PNG panels, donor expression, full 30-gene statistics,
selected-gene statistics, the gene lists, concordance results, input checksums,
`sessionInfo()` and a copy of the source script. It does not overwrite old figures
or reconstruct WGCNA modules.

## Published numerical source tables

Selected completed-run exports are available in [source_data/figure4/](../../source_data/figure4/README.md).
The [source-data manifest](../../source_data/manifest.tsv) identifies the original
run directories and checks performed. This is a partial source-data release;
see the linked guide for included and missing outputs. It does not replace the
external inputs required to run the analysis.
