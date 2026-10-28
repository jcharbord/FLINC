# Reference resources

This directory holds frozen reference inputs and curated annotations used by the
analysis scripts. Numerical exports from completed analyses are under
[`source_data/`](../source_data/README.md).

| Directory | Contents |
| --- | --- |
| [`figure2/`](figure2/) | Functional categories for the 200 displayed kME-selected hubs and functional-context annotations for PCSK7/PNPLA3 neighbours. |
| [`figure4/`](figure4/) | Reviewed human UniProt secreted-protein snapshot and frozen secreted-factor count, ranking and enrichment references. |

The Figure 2 annotations were previously stored under `scripts/figure2/reference/`;
the Figure 4 resources were previously under `resources/Fig4_secreted_factor/`.
Their contents are unchanged. Script defaults, documented commands and integrity
manifests use the paths above.

There are no Figure 1 or Figure 3 resource bundles in this release. Additional
resources should be added when their exact analysis versions have been identified.
Do not substitute a current database download for a frozen reference silently.
