# InSituViewer

Local desktop app for Xenium-style spatial transcriptomics analysis.

`InSituViewer` lets you run the pipeline and inspect results in one place:
- QC plots
- Spatial maps (KaroSpace)
- UMAP
- MANA-based compartment maps

## Quick start

Install core dependencies:

```bash
pip install scanpy pandas numpy scipy matplotlib seaborn PySide6
```

Optional dependencies:

```bash
pip install squidpy
pip install PySide6-QtWebEngine
```

- `squidpy` is used when MANA needs to build spatial neighbors.
- `PySide6-QtWebEngine` enables embedded KaroSpace in-app (otherwise it opens in browser).

Run the app:

```bash
python3 -m app.main
```

or:

```bash
python3 run_app.py
```

## App workflow

1. Set `Data dir` and `Output dir`.
2. Choose run discovery:
- `Run depth = Direct folders only` for `/dataset/output-*`.
- `Run depth = One level below samples` for `/dataset/sample_x/output-*`.
3. Choose sample ID source (`Auto`, `From run label`, `From parent folder`).
4. Optional: enable `MANA weighted aggregation`.
5. Optional: enable `Export KaroSpace HTML`.
6. Click `Run`.
7. Use top-bar actions: `Load Outputs`, `Generate UMAP`, `Generate Compartments`.

## Supported input layouts

Direct runs:

```text
/dataset/
  output-.../
    cell_feature_matrix.h5
    cells.csv.gz
```

Nested sample folders:

```text
/dataset/
  sample_A/
    output-.../
      cell_feature_matrix.h5
      cells.csv.gz
  sample_B/
    output-.../
      cell_feature_matrix.h5
      cells.csv.gz
```

## Outputs

Written under `--out-dir`:

- `data/raw.h5ad`
- `data/clustered.h5ad`
- `data/cluster_info.json`
- `data/markers_by_cluster.csv`
- `xenium_qc/summary_by_run.csv`
- `xenium_qc/gene_detection_overall.csv`
- `xenium_qc/*.png`
- `plots/umap.png` (generated from app action)
- `plots/compartments.png` (generated from app action)
- `karospace.html` (if export enabled)

## CLI fallback

You can still run everything without the app:

```bash
python3 run_xenium_analysis.py \
  --data-dir /absolute/path/to/dataset \
  --out-dir /absolute/path/to/output
```

Useful options:

- `--run-prefix`
- `--run-search-depth` (`1` or `2`)
- `--sample-id-source` (`auto`, `run`, `parent`)
- `--mana-aggregate`
- `--karospace-html /absolute/path/to/karospace.html`

## Visual themes

- App name: `InSituViewer`
- Themes: `app/theme_light.qss` and `app/theme_dark.qss`
- Toggle in top bar: `Dark` / `Light`

## Project notes

- Recent runs are stored at `~/.spatial-analysis-for-dummies/recent.json`.
- Utility modules are vendored under `utils/` (`MANA` + `KaroSpace` helpers).
