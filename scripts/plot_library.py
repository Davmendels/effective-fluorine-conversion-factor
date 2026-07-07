from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def load_json(path: str | Path) -> pd.DataFrame:
    path = Path(path)

    with path.open("r", encoding="utf-8") as f:
        records = json.load(f)

    df = pd.DataFrame(records)

    if "fluorine_fraction" not in df.columns:
        df["fluorine_fraction"] = (
            df["fluorine_atoms"].astype(float) * 18.998403163 / df["mw"].astype(float)
        )

    df["fluorine_percent"] = 100.0 * df["fluorine_fraction"].astype(float)

    return df


def add_chemistry_flags(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    name = df["name"].fillna("").str.lower()
    smiles = df.get("smiles", pd.Series("", index=df.index)).fillna("").str.lower()
    formula = df.get("formula", pd.Series("", index=df.index)).fillna("")

    df["has_ether"] = name.str.contains("ether|oxa|hfpo|adona", regex=True) | smiles.str.contains("o", regex=False)
    df["has_aromatic_name"] = name.str.contains("benz|phenyl|aromatic|tolyl|xyl", regex=True)
    df["has_sulfon"] = name.str.contains("sulfonic|sulfonate|sulfonamide|sulfone", regex=True)
    df["has_carboxyl"] = name.str.contains("carboxylic|carboxylate", regex=True)
    df["has_phosph"] = name.str.contains("phosphonic|phosphonate|phosphate", regex=True)
    df["has_chlorine"] = formula.str.contains("Cl", regex=False)

    return df


def save_histogram(df: pd.DataFrame, column: str, output_path: Path, title: str, xlabel: str, bins: int = 80) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(df[column].dropna(), bins=bins)
    ax.axvline(df[column].mean(), linestyle="--", linewidth=1.5, label=f"Mean = {df[column].mean():.3g}")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_scatter(df: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(df["mw"], df["fluorine_percent"], s=5, alpha=0.25)
    ax.set_title("Molecular weight versus fluorine mass fraction")
    ax.set_xlabel("Molecular weight (g/mol)")
    ax.set_ylabel("Fluorine mass fraction (%)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def summarize_flags(df: pd.DataFrame) -> pd.DataFrame:
    flag_cols = [
        "has_ether",
        "has_aromatic_name",
        "has_sulfon",
        "has_carboxyl",
        "has_phosph",
        "has_chlorine",
    ]

    rows = []

    for col in flag_cols:
        subset = df[df[col]]
        rows.append(
            {
                "flag": col,
                "count": int(len(subset)),
                "percent_of_library": 100.0 * len(subset) / len(df),
                "mean_mw": float(subset["mw"].mean()) if len(subset) else None,
                "mean_f_percent": float(subset["fluorine_percent"].mean()) if len(subset) else None,
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("library_json", help="Path to PFAS JSON library")
    parser.add_argument("--label", default=None, help="Short label for output folder")
    args = parser.parse_args()

    library_path = Path(args.library_json)
    label = args.label or library_path.stem

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("outputs/runs") / f"library_{label}_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_json(library_path)
    df = add_chemistry_flags(df)

    summary = {
        "library": str(library_path),
        "n_compounds": int(len(df)),
        "mw_mean": float(df["mw"].mean()),
        "mw_median": float(df["mw"].median()),
        "fluorine_percent_mean": float(df["fluorine_percent"].mean()),
        "fluorine_percent_median": float(df["fluorine_percent"].median()),
        "fluorine_percent_min": float(df["fluorine_percent"].min()),
        "fluorine_percent_max": float(df["fluorine_percent"].max()),
    }

    with (output_dir / "library_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    df.to_csv(output_dir / "library_table_with_flags.csv", index=False)

    flag_summary = summarize_flags(df)
    flag_summary.to_csv(output_dir / "chemistry_flag_summary.csv", index=False)

    save_histogram(
        df,
        "mw",
        output_dir / "mw_histogram.png",
        title=f"Molecular weight distribution: {label}",
        xlabel="Molecular weight (g/mol)",
        bins=80,
    )

    save_histogram(
        df,
        "fluorine_percent",
        output_dir / "fluorine_fraction_histogram.png",
        title=f"Fluorine mass fraction distribution: {label}",
        xlabel="Fluorine mass fraction (%)",
        bins=80,
    )

    save_histogram(
        df,
        "fluorine_atoms",
        output_dir / "fluorine_atoms_histogram.png",
        title=f"Fluorine atom count distribution: {label}",
        xlabel="Number of fluorine atoms",
        bins=60,
    )

    save_scatter(df, output_dir / "mw_vs_fluorine_fraction.png")

    print(f"Output written to: {output_dir}")
    print(json.dumps(summary, indent=2))
    print()
    print(flag_summary.to_string(index=False))


if __name__ == "__main__":
    main()
