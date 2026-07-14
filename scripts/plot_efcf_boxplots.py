from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


RUNS_DIR = Path("outputs/runs")
OUTPUT_DIR = Path("figures")


# Folder name -> (display name, compartment key)
from pfb_simulator.config import (
    COMPARTMENT_COLORS,
    COMPARTMENT_LABELS,
    DATASET_ORDER,
    DATASET_REGISTRY,
)

COMPARTMENT_ORDER = [
    "soil_extract",
    "sea_water",
    "surface_water",
    "industrial_discharge",
    "groundwater",
    "mixed_water",
    "surface_water_with_TFA",
    "tap_water",
    "drinking_water",
    "environmental_hotspot_water",
]

def get_run_dir(run_name: str) -> Path:
    """
    Return the deterministic output directory for one dataset.

    Individual analyses must first be generated with
    scripts/run_all_individual_experiments.py.
    """
    run_dir = RUNS_DIR / run_name
    csv_path = run_dir / "sample_burdens.csv"

    if not run_dir.is_dir():
        raise FileNotFoundError(
            f"Run directory not found:\n  {run_dir}\n"
            "Run scripts/run_all_individual_experiments.py first."
        )

    if not csv_path.is_file():
        raise FileNotFoundError(
            f"Individual-analysis output not found:\n  {csv_path}\n"
            "Run scripts/run_all_individual_experiments.py first."
        )

    return run_dir

def load_all_runs() -> pd.DataFrame:
    """
    Load all deterministic dataset-level analysis outputs.

    Each dataset is expected at:

        outputs/runs/<run_name>/sample_burdens.csv
    """
    frames: list[pd.DataFrame] = []

    for run_name, metadata in DATASET_REGISTRY.items():
        dataset_name, _country, compartment = metadata

        csv_path = RUNS_DIR / run_name / "sample_burdens.csv"

        if not csv_path.is_file():
            raise FileNotFoundError(
                f"Missing sample file for {dataset_name}:\n"
                f"  {csv_path}\n"
                "Run scripts/run_all_individual_experiments.py first."
            )

        df = pd.read_csv(csv_path)

        if "conversion_factor" not in df.columns:
            raise ValueError(
                f"'conversion_factor' missing from:\n  {csv_path}"
            )

        df["conversion_factor"] = pd.to_numeric(
            df["conversion_factor"],
            errors="coerce",
        )

        df = df[
            df["conversion_factor"].notna()
            & np.isfinite(df["conversion_factor"])
            & (df["conversion_factor"] > 0)
        ].copy()

        df["dataset"] = dataset_name
        df["compartment"] = compartment

        frames.append(
            df[
                [
                    "dataset",
                    "compartment",
                    "conversion_factor",
                ]
            ]
        )

        print(
            f"Loaded {dataset_name:24s}: "
            f"{len(df):>7,} valid samples from {csv_path}"
        )

    if not frames:
        raise RuntimeError("No datasets were loaded.")

    return pd.concat(
        frames,
        ignore_index=True,
        sort=False,
    )


def box_stats(values: pd.Series, label: str) -> dict:
    """
    Explicit boxplot statistics:
      box     = Q1 to Q3
      line    = median
      whisker = P05 to P95
    """
    values = pd.to_numeric(values, errors="coerce").dropna()

    if values.empty:
        raise ValueError(f"No valid EFCF values for {label}")

    return {
        "label": label,
        "med": float(values.median()),
        "q1": float(values.quantile(0.25)),
        "q3": float(values.quantile(0.75)),
        "whislo": float(values.quantile(0.05)),
        "whishi": float(values.quantile(0.95)),
        "fliers": [],
    }


def draw_box_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    group_column: str,
    order: list[str],
    labels: dict[str, str],
    colors: dict[str, str],
    title: str,
    global_median: float,
) -> None:
    stats: list[dict] = []
    used_colors: list[str] = []

    for key in order:
        group = df.loc[df[group_column] == key, "conversion_factor"]

        if group.empty:
            print(f"Warning: no observations for {group_column}={key}")
            continue

        stats.append(box_stats(group, labels.get(key, key)))
        used_colors.append(colors[key])

    artists = ax.bxp(
        stats,
        showfliers=False,
        patch_artist=True,
        widths=0.55,
        medianprops={
            "color": "#b66a3c",
            "linewidth": 1.4,
        },
        boxprops={
            "edgecolor": "#222222",
            "linewidth": 1.0,
        },
        whiskerprops={
            "color": "#222222",
            "linewidth": 1.0,
        },
        capprops={
            "color": "#222222",
            "linewidth": 1.0,
        },
    )

    for box, color in zip(artists["boxes"], used_colors, strict=True):
        box.set_facecolor(color)
        box.set_alpha(0.82)

    ax.axhline(
        global_median,
        linestyle="--",
        linewidth=1.2,
        color="#3f79a8",
        zorder=0,
    )

    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.set_ylabel("EFCF")
    ax.grid(axis="y", linewidth=0.5, alpha=0.18)
    ax.set_axisbelow(True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_all_runs()

    global_median = float(df["conversion_factor"].median())
    global_p05 = float(df["conversion_factor"].quantile(0.05))
    global_p95 = float(df["conversion_factor"].quantile(0.95))

    print(f"Loaded valid samples: {len(df):,}")
    print(f"Global median EFCF: {global_median:.6f}")
    print(f"Global P05–P95: {global_p05:.6f}–{global_p95:.6f}")

    # Dataset colours inherited from their compartment.
    dataset_to_compartment = {
        dataset_name: compartment
        for _, (dataset_name, _country, compartment) in DATASET_REGISTRY.items()
    }

    dataset_colors = {
        dataset_name: COMPARTMENT_COLORS[compartment]
        for dataset_name, compartment in dataset_to_compartment.items()
    }

    dataset_labels = {name: name for name in DATASET_ORDER}

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(12.5, 8.2),
        sharey=True,
        gridspec_kw={
            "height_ratios": [1.0, 1.05],
            "hspace": 0.37,
        },
    )

    draw_box_panel(
        ax=axes[0],
        df=df,
        group_column="compartment",
        order=COMPARTMENT_ORDER,
        labels=COMPARTMENT_LABELS,
        colors=COMPARTMENT_COLORS,
        title="(a) EFCF by environmental compartment",
        global_median=global_median,
    )

    draw_box_panel(
        ax=axes[1],
        df=df,
        group_column="dataset",
        order=DATASET_ORDER,
        labels=dataset_labels,
        colors=dataset_colors,
        title="(b) EFCF across independent monitoring programs",
        global_median=global_median,
    )

    # Shared y-range chosen to display the P05–P95 whiskers clearly.
    # ADES has P95 ≈ 2.00, so the upper limit must include it.
    axes[0].set_ylim(1.38, 1.80)

    axes[0].tick_params(axis="x", labelrotation=0, labelsize=8.5)
    axes[1].tick_params(axis="x", labelrotation=43, labelsize=8.3)

    for tick in axes[1].get_xticklabels():
        tick.set_horizontalalignment("right")

    legend_items = [
        mpatches.Patch(
            facecolor=COMPARTMENT_COLORS["soil_extract"],
            edgecolor="#222222",
            label="Soil extracts",
        ),
        mpatches.Patch(
            facecolor=COMPARTMENT_COLORS["mixed_water"],
            edgecolor="#222222",
            label="Natural aquatic environments",
        ),
        mpatches.Patch(
            facecolor=COMPARTMENT_COLORS["industrial_discharge"],
            edgecolor="#222222",
            label="Industrial discharges",
        ),
        mpatches.Patch(
            facecolor=COMPARTMENT_COLORS["drinking_water"],
            edgecolor="#222222",
            label="Drinking and tap water",
        ),
        mpatches.Patch(
            facecolor=COMPARTMENT_COLORS["surface_water_with_TFA"],
            edgecolor="#222222",
            label="TFA-focused dataset",
        ),
        mpatches.Patch(
            facecolor=COMPARTMENT_COLORS[
                "environmental_hotspot_water"
            ],
            edgecolor="#222222",
            label="Highly contaminated sites",
        ),
    ]

    fig.legend(
        handles=legend_items,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=3,
        frameon=True,
        fontsize=8.5,
    )

    fig.text(
        0.5,
        0.073,
        f"Dashed line: global median EFCF = {global_median:.3f}",
        ha="center",
        va="center",
        fontsize=9,
    )

    fig.subplots_adjust(
        left=0.075,
        right=0.99,
        top=0.97,
        bottom=0.19,
    )

    png_path = OUTPUT_DIR / "figure_1_efcf_boxplots.png"
    pdf_path = OUTPUT_DIR / "figure_1_efcf_boxplots.pdf"

    fig.savefig(png_path, dpi=400, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Written: {png_path}")
    print(f"Written: {pdf_path}")


if __name__ == "__main__":
    main()
