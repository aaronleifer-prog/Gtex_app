#!/usr/bin/env python3
"""Compute and plot a user-selected gene expression ratio per tissue from GTEx API v2.

This script:
1. Accepts two user-entered proteins (gene symbols).
2. Resolves versioned GENCODE IDs for both genes.
3. Pulls median GTEx tissue expression (TPM) for both genes.
4. Computes protein1/protein2 ratio per tissue.
5. Writes CSV output and a PNG chart.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import requests

BASE_URL = "https://gtexportal.org/api/v2"


def _get(endpoint: str, params: List[Tuple[str, str]], timeout: int) -> dict:
    url = f"{BASE_URL}/{endpoint}"
    try:
        response = requests.get(url, params=params, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(f"Failed to fetch GTEx endpoint '{endpoint}': {exc}") from exc
    payload = response.json()
    if not isinstance(payload, dict) or "data" not in payload:
        raise RuntimeError(f"Unexpected GTEx response format from {endpoint}")
    return payload


def _extract_gencode_id(row: dict) -> str | None:
    candidates = []
    for key in ("gencodeId", "gencodeIdVersion", "gencode_id", "gencode_id_version"):
        value = row.get(key)
        if isinstance(value, str) and value.startswith("ENSG"):
            candidates.append(value)
    if not candidates:
        return None
    # Prefer a versioned ID if available.
    versioned = [value for value in candidates if "." in value]
    return versioned[0] if versioned else candidates[0]


def fetch_gene_id_map(
    gene_symbols: Iterable[str],
    gencode_version: str,
    genome_build: str,
    timeout: int,
) -> Dict[str, str]:
    symbols = [symbol.upper() for symbol in gene_symbols]
    params: List[Tuple[str, str]] = [("geneId", symbol) for symbol in symbols]
    params.extend(
        [
            ("gencodeVersion", gencode_version),
            ("genomeBuild", genome_build),
            ("itemsPerPage", "1000"),
            ("page", "0"),
        ]
    )

    payload = _get("reference/gene", params=params, timeout=timeout)

    id_map: Dict[str, str] = {}
    for row in payload.get("data", []):
        if not isinstance(row, dict):
            continue

        symbol = row.get("geneSymbol")
        if not isinstance(symbol, str):
            continue
        symbol = symbol.upper()
        if symbol not in symbols:
            continue

        gencode_id = _extract_gencode_id(row)
        if gencode_id:
            id_map[symbol] = gencode_id

    missing = [symbol for symbol in symbols if symbol not in id_map]
    if missing:
        missing_text = ", ".join(missing)
        raise RuntimeError(
            "Could not resolve versioned GENCODE IDs for: "
            f"{missing_text}. Check gencodeVersion/genomeBuild."
        )

    return id_map


def fetch_median_gene_expression(
    gencode_ids: Dict[str, str],
    dataset_id: str,
    timeout: int,
) -> pd.DataFrame:
    id_to_symbol = {gencode_id: symbol for symbol, gencode_id in gencode_ids.items()}
    params: List[Tuple[str, str]] = [("gencodeId", gencode_id) for gencode_id in id_to_symbol]
    params.extend(
        [
            ("datasetId", dataset_id),
            ("itemsPerPage", "100000"),
            ("page", "0"),
        ]
    )

    payload = _get("expression/medianGeneExpression", params=params, timeout=timeout)

    rows = []
    for row in payload.get("data", []):
        if not isinstance(row, dict):
            continue

        gencode_id = row.get("gencodeId")
        symbol = row.get("geneSymbol")
        tissue = row.get("tissueSiteDetail")
        if not isinstance(tissue, str):
            tissue = row.get("tissueSiteDetailId")
        median = row.get("median")

        if not isinstance(gencode_id, str):
            continue
        if isinstance(symbol, str):
            symbol = symbol.upper()
        else:
            symbol = id_to_symbol.get(gencode_id, "")
        if symbol not in gencode_ids:
            symbol = id_to_symbol.get(gencode_id, symbol)

        if not isinstance(tissue, str):
            continue

        if median is None:
            continue
        median_value = float(median)

        rows.append(
            {
                "tissue": tissue,
                "gene_symbol": symbol,
                "gencode_id": gencode_id,
                "median_tpm": median_value,
            }
        )

    if not rows:
        raise RuntimeError("No expression rows returned from GTEx medianGeneExpression endpoint.")

    return pd.DataFrame(rows)


def build_ratio_table(expression_df: pd.DataFrame, protein_1: str, protein_2: str) -> pd.DataFrame:
    pivot = expression_df.pivot_table(
        index="tissue",
        columns="gene_symbol",
        values="median_tpm",
        aggfunc="first",
    )

    required = {protein_1, protein_2}
    missing = required.difference(set(pivot.columns))
    if missing:
        raise RuntimeError(f"Missing required genes in expression table: {', '.join(sorted(missing))}")

    ratio_df = pivot.reset_index().rename(
        columns={protein_1: "protein_1_median_tpm", protein_2: "protein_2_median_tpm"}
    )
    ratio_df.insert(1, "protein_1_symbol", protein_1)
    ratio_df.insert(2, "protein_2_symbol", protein_2)
    ratio_df["protein1_over_protein2_ratio"] = ratio_df.apply(
        lambda row: (row["protein_1_median_tpm"] / row["protein_2_median_tpm"])
        if row["protein_2_median_tpm"] > 0
        else math.nan,
        axis=1,
    )
    ratio_df["log2_ratio"] = ratio_df["protein1_over_protein2_ratio"].apply(
        lambda x: math.log2(x) if pd.notna(x) and x > 0 else math.nan
    )

    return ratio_df.sort_values("protein1_over_protein2_ratio", ascending=False).reset_index(drop=True)


def plot_ratio_table(ratio_df: pd.DataFrame, output_path: Path, protein_1: str, protein_2: str) -> None:
    plot_df = ratio_df.dropna(subset=["protein1_over_protein2_ratio"]).copy()
    plot_df = plot_df[plot_df["protein1_over_protein2_ratio"] > 0]

    if plot_df.empty:
        raise RuntimeError(f"No positive {protein_1}/{protein_2} ratios available to plot.")

    plot_df = plot_df.sort_values("protein1_over_protein2_ratio", ascending=True)

    fig_height = max(7, len(plot_df) * 0.22)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    ax.barh(plot_df["tissue"], plot_df["protein1_over_protein2_ratio"], color="#1f77b4", edgecolor="none")
    ax.set_xscale("log")
    ax.axvline(1.0, linestyle="--", linewidth=1.2, color="#333333", alpha=0.9)
    ax.set_xlabel(f"{protein_1} / {protein_2} median TPM ratio (log scale)")
    ax.set_ylabel("Tissue")
    ax.set_title(f"GTEx {protein_1}/{protein_2} Median Expression Ratio by Tissue")
    ax.grid(axis="x", linestyle=":", linewidth=0.8, alpha=0.5)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prompt for two proteins (gene symbols), pull GTEx tissue median expression, compute ratio per tissue, "
            "and output CSV + plot."
        )
    )
    parser.add_argument(
        "--protein-1",
        "--gene-1",
        dest="protein_1",
        help="Protein 1 gene symbol (numerator). If omitted, you will be prompted.",
    )
    parser.add_argument(
        "--protein-2",
        "--gene-2",
        dest="protein_2",
        help="Protein 2 gene symbol (denominator). If omitted, you will be prompted.",
    )
    parser.add_argument("--dataset-id", default="gtex_v8", help="GTEx dataset ID (default: gtex_v8).")
    parser.add_argument(
        "--gencode-version",
        default="v26",
        help="GENCODE version for gene symbol lookup (default: v26).",
    )
    parser.add_argument(
        "--genome-build",
        default="GRCh38/hg38",
        help="Genome build for gene symbol lookup (default: GRCh38/hg38).",
    )
    parser.add_argument(
        "--out-csv",
        default=None,
        help="Output CSV path (default: <protein1>_<protein2>_ratio_by_tissue.csv).",
    )
    parser.add_argument(
        "--out-plot",
        default=None,
        help="Output plot path (default: <protein1>_<protein2>_ratio_by_tissue.png).",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="HTTP timeout in seconds for GTEx requests (default: 30).",
    )
    return parser.parse_args()


def _prompt_symbol(prompt_text: str) -> str:
    value = input(prompt_text).strip().upper()
    if not value:
        raise RuntimeError("Protein symbol cannot be empty.")
    return value


def resolve_protein_symbols(args: argparse.Namespace) -> Tuple[str, str]:
    protein_1 = args.protein_1.strip().upper() if args.protein_1 else _prompt_symbol(
        "Enter protein 1 (gene symbol, numerator): "
    )
    protein_2 = args.protein_2.strip().upper() if args.protein_2 else _prompt_symbol(
        "Enter protein 2 (gene symbol, denominator): "
    )

    if protein_1 == protein_2:
        raise RuntimeError("Protein 1 and protein 2 must be different symbols.")
    return protein_1, protein_2


def _slugify_symbol(symbol: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", symbol).strip("_").lower()
    return slug or "gene"


def main() -> None:
    args = parse_args()

    try:
        protein_1, protein_2 = resolve_protein_symbols(args)
        gene_map = fetch_gene_id_map(
            gene_symbols=[protein_1, protein_2],
            gencode_version=args.gencode_version,
            genome_build=args.genome_build,
            timeout=args.timeout,
        )

        expression_df = fetch_median_gene_expression(
            gencode_ids=gene_map,
            dataset_id=args.dataset_id,
            timeout=args.timeout,
        )

        ratio_df = build_ratio_table(expression_df, protein_1=protein_1, protein_2=protein_2)

        default_stem = f"{_slugify_symbol(protein_1)}_{_slugify_symbol(protein_2)}_ratio_by_tissue"
        out_csv = Path(args.out_csv) if args.out_csv else Path(f"{default_stem}.csv")
        out_plot = Path(args.out_plot) if args.out_plot else Path(f"{default_stem}.png")

        ratio_df.to_csv(out_csv, index=False)
        plot_ratio_table(ratio_df, out_plot, protein_1=protein_1, protein_2=protein_2)

    except Exception as exc:
        raise SystemExit(f"Error: {exc}")

    print(f"Saved ratio table: {out_csv.resolve()}")
    print(f"Saved plot:        {out_plot.resolve()}")
    print(f"Rows: {len(ratio_df)} tissues")


if __name__ == "__main__":
    main()
