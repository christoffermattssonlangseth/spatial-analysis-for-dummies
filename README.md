![InSituCore](assets/logo.png)

Local desktop app for Xenium-style spatial transcriptomics analysis.

`InSituCore` lets you run the pipeline and inspect results in one place:
- QC plots
- Spatial maps (static plot + optional KaroSpace)
- UMAP
- MANA-based compartment maps
- Gene-expression dotplots

## Quick start

Create and activate an environment (recommended):

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

Install core dependencies:

```bash
pip install -r requirements.txt
```

Optional dependencies:

```bash
pip install -r requirements-optional.txt
```

- `squidpy` is used when MANA needs to build spatial neighbors from coordinates.
- `louvain` is required only when using Louvain clustering.
- `scvi-tools` is required when using scVI latent representation for MANA.
- `cellcharter` is optional and used to remove long spatial links after graph construction.
- `pyarrow` (or another parquet backend) is recommended when using transcript-level count matrix mode (`nucleus_or_distance`).
- MANA spatial graph is built per sample (`library_key=sample_id` when available) to avoid cross-sample edges.
- Long spatial edges are pruned with `cellcharter.gr.remove_long_links` by default after graph construction.
- Embedded KaroSpace in-app requires `PySide6.QtWebEngineWidgets`.
  Depending on Python and platform, this may come from `PySide6` directly or require:

```bash
pip install PySide6-QtWebEngine
```

If `PySide6-QtWebEngine` is unavailable for your Python version, the app still works and opens KaroSpace in your browser.

Check your environment:

```bash
python3 check_env.py
```

Strict mode (fails if optional packages are missing):

```bash
python3 check_env.py --require-optional
```

Run the app:

```bash
python3 -m app.main
```

or:

```bash
python3 run_app.py
```

Generate/refresh macOS icon assets:

```bash
python3 scripts/generate_macos_icon.py
```

## Open as app (macOS)

Build a clickable app launcher (no Terminal window):

```bash
bash scripts/build_macos_app.sh
```

This creates:

- `dist/InSituCore.app`

Double-click `dist/InSituCore.app` in Finder to launch the app.
Logs are written to:

- `~/Library/Logs/InSituCore.log`

## Portable bundle (macOS)

Create a portable folder you can zip/share and run without opening Terminal:

```bash
bash scripts/build_portable_macos_bundle.sh
```

Outputs:

- `dist/InSituCore-portable/InSituCore.app`
- `dist/InSituCore-portable/InSituCore/` (project + local runtime)
- `dist/InSituCore-portable/environment.yml` (conda environment file)
- `dist/InSituCore-portable/InSituCore.log`

Create a zip in one step:

```bash
bash scripts/build_portable_macos_bundle.sh dist --zip
```

This will produce `dist/InSituCore-portable.zip`.

## App workflow

1. Set `Data dir` and `Output dir`.
2. Choose run discovery:
- `Run depth = Direct folders only` for `/dataset/output-*`.
- `Run depth = One level below samples` for `/dataset/sample_x/output-*`.
3. Choose sample ID source (`Auto`, `From run label`, `From parent folder`).
4. Choose `Count matrix mode`:
- `Standard Xenium cell matrix` uses `cell_feature_matrix.h5`.
- `Nucleus OR distance-filtered transcripts` builds counts from `transcripts.parquet` using:
  `overlaps_nucleus OR nucleus_distance <= Tx max distance (um)`.
5. Optional: enable `MANA weighted aggregation`.
6. Optional: enable `Export KaroSpace HTML`.
7. Open `Analysis Controls` in the workspace navigator (expanded by default) to configure:
- neighbors / PCs / UMAP min_dist
- clustering method (`Leiden`, `Louvain`, or `KMeans`)
- method-specific sweep values (resolutions or K list)
8. In `Run` -> `Optional Steps`, choose MANA representation (`scVI` recommended).
9. Click `Run` (top bar).
10. Use top-bar `Load Outputs` to refresh an existing output directory.
11. Use the workspace navigator on the left to switch pages:
- `Spatial Static` for `Generate Spatial Map` (fast native plot)
- `Spatial Interactive` for embedded KaroSpace
- `UMAP` for cluster-key-aware UMAP generation
- `Compartment Map` for compartment key selection and map generation
- `Gene Expression` for top-N marker dotplots
12. For compartments, choose a compartment key in `Compartment Map`:
- `Auto (primary)` uses the pipeline-selected default.
- Or choose a specific key like `compartment_gmm_k6` / `compartment_leiden_1.0`.
13. For dotplots, choose a group key and top-gene count in `Gene Expression`.

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

Transcript-derived count matrix mode (`nucleus_or_distance`):

```text
/dataset/
  output-.../
    transcripts.parquet          # file or parquet shard directory
    cells.parquet                # preferred
    # or cells.csv.gz
```

The transcript filter matches your notebook workflow:

- keep transcript if `overlaps_nucleus == True`
- OR if `nucleus_distance <= --tx-max-distance-um`
- and keep only configured categories (default: `predesigned_gene,custom_gene`)

## Outputs

Written under `--out-dir`:

- `data/raw.h5ad`
- `data/clustered.h5ad`
- `data/cluster_info.json`
- `data/markers_by_cluster.csv`
- `data/compartment_models.csv` (model summary for MANA compartments, including AIC/BIC for GMM)
- `xenium_qc/summary_by_run.csv`
- `xenium_qc/gene_detection_overall.csv`
- `xenium_qc/*.png`
- `plots/spatial.png` (generated from `Spatial Static`)
- `plots/umap.png` (generated from `UMAP`)
- `plots/compartments.png` (generated from `Compartment Map`)
- `plots/gene_expression_dotplot.png` (generated from `Gene Expression`)
- `karospace.html` (if export enabled)

## CLI fallback

You can still run everything without the app:

```bash
python3 run_xenium_analysis.py \
  --data-dir /absolute/path/to/dataset \
  --out-dir /absolute/path/to/output
```

Example using transcript-derived count matrix (nucleus OR <= 5 um):

```bash
python3 run_xenium_analysis.py \
  --data-dir /absolute/path/to/dataset \
  --out-dir /absolute/path/to/output \
  --count-matrix-mode nucleus_or_distance \
  --tx-max-distance-um 5 \
  --tx-nucleus-distance-key nucleus_distance
```

Useful options:

- `--run-prefix`
- `--run-search-depth` (`1` or `2`)
- `--sample-id-source` (`auto`, `run`, `parent`)
- `--count-matrix-mode` (`cell_feature_matrix`, `nucleus_or_distance`)
- `--tx-max-distance-um` (default `5.0`)
- `--tx-nucleus-distance-key` (default `nucleus_distance`; auto-infers common names if missing)
- `--tx-allowed-categories` (default `predesigned_gene,custom_gene`)
- `--cluster-method` (`leiden`, `louvain`, `kmeans`)
- `--n-neighbors`, `--n-pcs`, `--umap-min-dist`
- `--leiden-resolutions`, `--louvain-resolutions`, `--kmeans-clusters`
- `--spatial-long-links-percentile` (default `99.0`)
- `--spatial-no-remove-long-links` (disables long-link pruning)
- `--mana-aggregate`
- `--mana-representation-mode` (`scvi`, `pca`, `custom`, `auto`; default `scvi`)
- `--scvi-latent-key`, `--scvi-n-latent`, `--scvi-max-epochs` (default `30`)
- `--scvi-hvg-top-genes` (default `500`) and `--scvi-hvg-flavor` (`seurat_v3`)
- notebook-like MANA defaults: `--mana-distance-kernel gaussian`, `--mana-hop-decay 0.2`, `--mana-out-key X_mana_gauss`
- `--mana-compartment-method` (`gmm`, `leiden`, `both`; default `gmm`)
- `--mana-gmm-components` (e.g. `6,10,14`)
- `--mana-gmm-max-dims` (PCA cap before GMM; default `30`)
- `--mana-gmm-covariance-type` (default `diag`, safer for large data)
- `--mana-compartment-resolutions` (for Leiden compartments)
- `--karospace-html /absolute/path/to/karospace.html`

## Visual themes

- App name: `InSituCore`
- Themes: `app/theme_light.qss` and `app/theme_dark.qss`
- Startup behavior: detects system appearance (dark/light) automatically
- Toggle in top bar: `Dark` / `Light` (manual override)

## Project notes

- Recent runs are stored at `~/.insitucore/recent.json`.
- Utility modules are vendored under `utils/` (`MANA` + `KaroSpace` helpers).
- Local build/output artifacts are git-ignored (`dist/`, caches, generated icons).
- Run `bash scripts/repo_hygiene.sh` to inspect large files.
- Run `bash scripts/repo_hygiene.sh --clean-local` to remove local build/cache artifacts.
