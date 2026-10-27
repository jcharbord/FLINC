# STRIDE deconvolution source and validation

The four original analysis scripts and readme by Jian Wu are preserved unchanged
in `jian_original/`, copied from `jcharbord/FLINC`, `main`, `Deconvolution/`.
They are also retained in the repository-root `code/deconvolution/` on the
reproducibility branch. These are analysis sources, not just validation code.

| Order | Original script | Purpose |
| --- | --- | --- |
| 1 | `jian_original/1_write_single_nucleus_count_matrix.R` | Export reference expression matrix and cell-type labels. |
| 2 | `jian_original/2_write_ST_count_matrix.R` | Export spatial expression matrix. |
| 3 | `jian_original/3_deconvolution_stride.sh` | Run STRIDE deconvolution. |
| 4 | `jian_original/4_add_assay_deconvolution.R` | Import fractions as the Seurat `singlecell` assay. |

`01_validate_deconvolution_input.R` independently checks the saved fraction table
used downstream; it does not calculate deconvolution.

## Before a fresh upstream rerun

The originals retain Jian's working directory, relative paths (including
`.../../` output paths), direct Seurat slot access and the `stride` environment.
Scripts 1 and 2 export the `data` slot, despite referring to counts in their names.
Confirm the contents of those slots in the original objects before treating the
exports as raw counts; do not silently substitute a different expression matrix.
Script 4 saves over its input object. Reconstruct the intended paths, inputs and
software environment before execution. These original scripts were not rerun as
part of this reorganization, and copying them does not establish a portable rerun.

See ../../DECONVOLUTION_SOURCE_WORKFLOW.md for the canonical fraction-table hash,
input validation and the documented reproducibility boundary.
