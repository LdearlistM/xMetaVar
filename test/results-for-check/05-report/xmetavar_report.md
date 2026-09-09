# xMetaVar Workflow Report

## Run overview

- Samples: 2
- Sequencing mode: PE
- Ready output files: 11
- Existing output files: 11
- Standardized features: 7758
- Native matrix dimensions: 2 samples x 7758 features
- Native non-zero fraction: 88.45%
- Scaled positive fraction: 76.06%
- vSV binary threshold: abs(value) >= 2.0
- CNV binary threshold: copy number >= 0.5
- Partial deliverables mode: enabled. Modules without available source files are reported as Not run.


## Main output directories

- Webserver-ready files: `03-webserver-inputs/`
- Summary tables: `04-summary/`
- Type-specific native matrices: `04-summary/type_matrices/`
- HTML report: `05-report/xmetavar_report.html`
- Original tool outputs: `02-variant-calling/`

## Core summary files

- `04-summary/output_manifest.tsv`
- `04-summary/variant_feature_manifest.tsv`
- `04-summary/xmetavar_feature_matrix.tsv`
- `04-summary/xmetavar_scaled_matrix.tsv`
- `04-summary/sample_variant_burden.tsv`
- `04-summary/matrix_overview.tsv`
