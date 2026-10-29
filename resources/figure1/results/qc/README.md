# Spatial QC summaries

These tables summarize the final annotated dataset: 26,056 spatial spots from
24 liver biopsies. They support the QC sentence and Supplementary Figure 1D-E.

| Metric | Pooled spot median | Pooled spot mean |
| --- | ---: | ---: |
| Detected genes (`nFeature_Spatial`) | 4,845 | 4,666.37 |
| Spatial UMIs (`nCount_Spatial`) | 20,064 | 20,916.38 |

The manuscript values 4,845 and 20,064 are medians across all retained spots.
They are not means or the median of the 24 sample medians. Samples with more
spots contribute more observations to these pooled per-spot summaries.

## Files

- [QC_overall_spots.csv](QC_overall_spots.csv): spot counts, means, medians,
  quartiles and ranges for both QC metrics.
- [QC_manuscript_value_check.csv](QC_manuscript_value_check.csv): the two
  manuscript values alongside the calculated pooled means and medians.
- [QC_summary_definitions.csv](QC_summary_definitions.csv): pooled summaries
  and alternative equal-sample-weight summaries, with their definitions in
  the column names. These alternatives are distinct from the manuscript values.

## Source and verification

The tables were generated on 29 October 2026 by
[`07_summarize_spatial_QC.R`](../../../../scripts/figure1/07_summarize_spatial_QC.R),
using `nFeature_Spatial`, `nCount_Spatial` and `donor` in the final annotated
Seurat object. Recalculation from the corresponding spot-level QC export
reproduced the published means, medians, quartiles, ranges and alternative
summary definitions. CSV bytes are unchanged from the completed run.

Original local output directory, relative to the workflow root:

`results/261029_QC_final_object_20261029T174325_30456/`

Software versions for that run are in
[`sessionInfo.txt`](../../../../reproducibility_logs/Fig1_QC_20261029/sessionInfo.txt).
The public script uses the same QC calculations as the completed run; only its
commented usage example was adjusted to the repository path.

## Regenerating the plots

With approved access to the final Seurat object loaded as `obj`, run from the
repository root:

```r
source("scripts/figure1/07_summarize_spatial_QC.R")
qc_result <- run_flinc_qc(obj, output_parent = "results")
```

This calculates QC summaries and writes per-sample violin plots to a new local
output folder. Black dots show the median within each sample; these differ
from the pooled medians reported in the manuscript. The function also exports
local sample-level and spot-level QC tables. It does not modify the Seurat
object, apply additional spot filtering, or rerun normalization or clustering.

Only cohort-level summary tables are distributed in this folder. Sample-level
QC tables, spot metadata and the Seurat object are not included.
