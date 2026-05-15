# GTEx Protein Ratio

This project prompts you for two proteins (gene symbols), pulls GTEx tissue median expression data for both, computes the `protein1/protein2` ratio per tissue, and creates a plot.

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 gtex_ldh_ratio.py
```

When prompted, enter your two gene symbols for a `protein1/protein2` ratio.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
python3 gtex_ldh_ratio.py
```

You will be prompted:
- `Enter protein 1 (gene symbol, numerator):`
- `Enter protein 2 (gene symbol, denominator):`

Outputs:
- `<protein1>_<protein2>_ratio_by_tissue.csv`
- `<protein1>_<protein2>_ratio_by_tissue.png`

## Optional arguments

```bash
python3 gtex_ldh_ratio.py \
  --protein-1 LDHA \
  --protein-2 LDHB \
  --dataset-id gtex_v8 \
  --gencode-version v26 \
  --genome-build "GRCh38/hg38" \
  --out-csv results.csv \
  --out-plot results.png
```
