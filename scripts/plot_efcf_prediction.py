from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


RUNS_DIR = Path("outputs/runs")
OUTPUT_DIR = Path("figures")


# Prefix of output folder -> (display name, compartment)
#
# Prefix matching avoids hard-coding timestamps. If several matching runs exist,
# the most recently modified valid run is used.
from pfb_simulator.config import (
    COMPARTMENT_COLORS,
    DATASET_ORDER,
    DATASET_REGISTRY,
)

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

def load_all_samples() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []

    required_columns = {
        "total_pfas_ng_l",
        "fluorine_burden_ng_l",
        "conversion_factor",
    }

    for run_name, metadata in DATASET_REGISTRY.items():
        dataset, _country, compartment = metadata
        run_dir = get_run_dir(run_name)
        csv_path = run_dir / "sample_burdens.csv"

        df = pd.read_csv(csv_path)

        missing = required_columns.difference(df.columns)
        if missing:
            raise ValueError(
                f"{csv_path} is missing required columns: {sorted(missing)}"
            )

        for column in required_columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

        df = df.dropna(subset=list(required_columns)).copy()

        df = df[
            np.isfinite(df["total_pfas_ng_l"])
            & np.isfinite(df["fluorine_burden_ng_l"])
            & np.isfinite(df["conversion_factor"])
            & (df["total_pfas_ng_l"] > 0)
            & (df["fluorine_burden_ng_l"] > 0)
            & (df["conversion_factor"] > 0)
        ].copy()

        df["dataset"] = dataset
        df["compartment"] = compartment
        df["source_run"] = run_dir.name

        frames.append(
            df[
                [
                    "dataset",
                    "compartment",
                    "source_run",
                    "total_pfas_ng_l",
                    "fluorine_burden_ng_l",
                    "conversion_factor",
                ]
            ]
        )

        print(
            f"Loaded {dataset:24s}: "
            f"{len(df):>7,} valid samples from {run_dir.name}"
        )

    if not frames:
        raise RuntimeError("No valid datasets were loaded.")

    return pd.concat(frames, ignore_index=True)


def summarize_datasets(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []

    for dataset in DATASET_ORDER:
        group = df[df["dataset"] == dataset]

        if group.empty:
            print(f"Warning: no observations found for {dataset}")
            continue

        ratio = group["estimation_ratio"]

        rows.append(
            {
                "dataset": dataset,
                "compartment": group["compartment"].iloc[0],
                "n": len(group),
                "median_reported_pfas_ng_l":
                    group["total_pfas_ng_l"].median(),
                "median_estimated_pfas_ng_l":
                    group["estimated_pfas_ng_l"].median(),
                "median_ratio": ratio.median(),
                "p05_ratio": ratio.quantile(0.05),
                "p95_ratio": ratio.quantile(0.95),
                "median_relative_error_percent":
                    100.0 * (ratio.median() - 1.0),
                "p05_relative_error_percent":
                    100.0 * (ratio.quantile(0.05) - 1.0),
                "p95_relative_error_percent":
                    100.0 * (ratio.quantile(0.95) - 1.0),
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_all_samples()

    # Global conversion statistics.
    global_efcf = float(df["conversion_factor"].median())
    global_efcf_p05 = float(df["conversion_factor"].quantile(0.05))
    global_efcf_p95 = float(df["conversion_factor"].quantile(0.95))

    # Estimate total PFAS solely from fluorine burden and the global EFCF.
    df["estimated_pfas_ng_l"] = (
        global_efcf * df["fluorine_burden_ng_l"]
    )

    # Ratio = estimated PFAS / reported PFAS.
    df["estimation_ratio"] = (
        df["estimated_pfas_ng_l"] / df["total_pfas_ng_l"]
    )

    ratio_median = float(df["estimation_ratio"].median())
    ratio_p05 = float(df["estimation_ratio"].quantile(0.05))
    ratio_p95 = float(df["estimation_ratio"].quantile(0.95))

    dataset_summary = summarize_datasets(df)

    print()
    print(f"Total valid samples: {len(df):,}")
    print(f"Global median EFCF: {global_efcf:.6f}")
    print(
        "Global EFCF P05–P95: "
        f"{global_efcf_p05:.6f}–{global_efcf_p95:.6f}"
    )
    print(f"Median estimate / reported ratio: {ratio_median:.6f}")
    print(
        "Estimate / reported P05–P95: "
        f"{ratio_p05:.6f}–{ratio_p95:.6f}"
    )
    print(
        "Equivalent relative-error interval: "
        f"{100 * (ratio_p05 - 1):+.1f}% to "
        f"{100 * (ratio_p95 - 1):+.1f}%"
    )

    # Save numerical results used in the figure.
    sample_output = OUTPUT_DIR / "figure_3_sample_predictions.csv.gz"
    summary_output = OUTPUT_DIR / "figure_3_dataset_summary.csv"

    df.to_csv(sample_output, index=False, compression="gzip")
    dataset_summary.to_csv(summary_output, index=False)

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(13.5, 7.2))

    grid = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[1.45, 1.0],
        wspace=0.28,
    )

    ax_main = fig.add_subplot(grid[0, 0])
    ax_ratio = fig.add_subplot(grid[0, 1])

    # ------------------------------------------------------------------
    # Panel (a): sample-level estimated versus reported PFAS
    # ------------------------------------------------------------------
    x = df["total_pfas_ng_l"].to_numpy()
    y = df["estimated_pfas_ng_l"].to_numpy()

    # Display range only. All valid observations remain included in the
    # global statistics and in panel (b).
    plot_min = 1.0
    plot_max = 1.0e5

    mask = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x >= plot_min)
        & (x <= plot_max)
        & (y >= plot_min)
        & (y <= plot_max)
    )

    x_plot = x[mask]
    y_plot = y[mask]

    print(
        f"Panel (a): displaying {len(x_plot):,} of {len(x):,} valid samples "
        f"within {plot_min:g}–{plot_max:g} ng/L."
    )

    identity_x = np.logspace(
        np.log10(plot_min),
        np.log10(plot_max),
        400,
    )

    # Individual environmental samples.
    #
    # rasterized=True keeps the PDF manageable while retaining vector
    # axes, labels, identity line, and prediction boundaries.
    ax_main.scatter(
        x_plot,
        y_plot,
        s=1.8,
        alpha=0.012,
        color="#4f88b8",
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )

    # Dataset-level median positions.
    #
    # These are the same monitoring-program summaries represented in
    # panel (b), projected into reported-versus-estimated concentration
    # space. Their alignment with the identity line provides the clearest
    # visual summary of the global EFCF performance.
    for _, row in dataset_summary.iterrows():
        x_median = float(row["median_reported_pfas_ng_l"])
        y_median = float(row["median_estimated_pfas_ng_l"])

        # Show only medians within the displayed axis range.
        if not (
            plot_min <= x_median <= plot_max
            and plot_min <= y_median <= plot_max
        ):
            continue

        color = COMPARTMENT_COLORS[row["compartment"]]

        ax_main.scatter(
            x_median,
            y_median,
            s=72,
            facecolor=color,
            edgecolor="#222222",
            linewidth=0.8,
            zorder=5,
        )

    # Empirical 5th and 95th percentile prediction boundaries.
    #
    # These delimit the range containing 90% of sample-level estimates.
    ax_main.plot(
        identity_x,
        ratio_p05 * identity_x,
        linestyle="--",
        linewidth=1.15,
        color="#58758d",
        label="90% empirical prediction interval",
        zorder=3,
    )

    ax_main.plot(
        identity_x,
        ratio_p95 * identity_x,
        linestyle="--",
        linewidth=1.15,
        color="#58758d",
        zorder=3,
    )

    # Identity line: perfect agreement.
    ax_main.plot(
        identity_x,
        identity_x,
        color="#303030",
        linewidth=1.55,
        label="Identity",
        zorder=4,
    )

    ax_main.set_xscale("log")
    ax_main.set_yscale("log")

    ax_main.set_xlim(plot_min, plot_max)
    ax_main.set_ylim(plot_min, plot_max)
    ax_main.set_aspect("equal", adjustable="box")

    ax_main.set_xlabel(
        r"Reported total PFAS concentration (ng L$^{-1}$)"
    )
    ax_main.set_ylabel(
        r"EFCF-estimated total PFAS concentration (ng L$^{-1}$)"
    )
    ax_main.set_title(
        "(a) Estimated versus reported total PFAS",
        loc="left",
        fontweight="bold",
    )

    ax_main.grid(
        which="major",
        linewidth=0.55,
        alpha=0.20,
    )
    ax_main.set_axisbelow(True)

    ax_main.legend(
        loc="upper left",
        frameon=True,
        fontsize=8.5,
    )

    ax_main.text(
        0.97,
        0.04,
        (
            f"Global median EFCF = {global_efcf:.3f}\n"
            f"90% prediction interval:\n"
            f"{100 * (ratio_p05 - 1):+.1f}% to "
            f"{100 * (ratio_p95 - 1):+.1f}%"
        ),
        transform=ax_main.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#aaaaaa",
            "alpha": 0.92,
        },
        zorder=6,
    )

    # ------------------------------------------------------------------
    # Panel (b): dataset-level ratios and error bars
    # ------------------------------------------------------------------
    display_summary = dataset_summary.copy()

    # Reverse order so first item appears at the top.
    display_summary = display_summary.iloc[::-1].reset_index(drop=True)

    positions = np.arange(len(display_summary))

    for position, (_, row) in zip(
        positions,
        display_summary.iterrows(),
        strict=True,
    ):
        median = row["median_ratio"]
        lower = median - row["p05_ratio"]
        upper = row["p95_ratio"] - median

        color = COMPARTMENT_COLORS[row["compartment"]]

        ax_ratio.errorbar(
            median,
            position,
            xerr=np.array([[lower], [upper]]),
            fmt="o",
            markersize=5.8,
            markerfacecolor=color,
            markeredgecolor="#222222",
            markeredgewidth=0.65,
            ecolor="#555555",
            elinewidth=1.0,
            capsize=3.0,
            zorder=4,
        )

    # Global empirical prediction interval.
    ax_ratio.axvspan(
        ratio_p05,
        ratio_p95,
        facecolor="#9fb6ca",
        alpha=0.22,
        zorder=0,
    )

    ax_ratio.axvline(
        1.0,
        color="#303030",
        linewidth=1.35,
        zorder=2,
    )

    ax_ratio.axvline(
        ratio_p05,
        color="#58758d",
        linewidth=0.85,
        linestyle="--",
        zorder=1,
    )

    ax_ratio.axvline(
        ratio_p95,
        color="#58758d",
        linewidth=0.85,
        linestyle="--",
        zorder=1,
    )

    ax_ratio.set_yticks(positions)
    ax_ratio.set_yticklabels(display_summary["dataset"], fontsize=8.3)

    # Choose a range based on dataset intervals, with modest padding.
    ratio_axis_min = min(
        0.85,
        float(display_summary["p05_ratio"].min()) - 0.025,
    )
    ratio_axis_max = max(
        1.15,
        float(display_summary["p95_ratio"].max()) + 0.025,
    )

    ax_ratio.set_xlim(ratio_axis_min, ratio_axis_max)

    ax_ratio.set_xlabel(
        r"Estimated / reported total PFAS"
    )
    ax_ratio.set_title(
        "(b) Accuracy across monitoring programs",
        loc="left",
        fontweight="bold",
    )

    ax_ratio.grid(
        axis="x",
        linewidth=0.55,
        alpha=0.20,
    )
    ax_ratio.set_axisbelow(True)

    # Global interval annotation.
    ax_ratio.text(
        0.02,
        0.02,
        (
            f"Global median EFCF = {global_efcf:.3f}\n"
            f"90% estimate ratio: {ratio_p05:.3f}–{ratio_p95:.3f}\n"
            f"Relative error: "
            f"{100 * (ratio_p05 - 1):+.1f}% to "
            f"{100 * (ratio_p95 - 1):+.1f}%"
        ),
        transform=ax_ratio.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.2,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#aaaaaa",
            "alpha": 0.92,
        },
    )

    # Legend for compartment colors.
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
        bbox_to_anchor=(0.5, -0.015),
        ncol=3,
        fontsize=8.2,
        frameon=True,
    )

    fig.subplots_adjust(
        left=0.075,
        right=0.985,
        top=0.94,
        bottom=0.17,
    )

    png_path = OUTPUT_DIR / "figure_3_efcf_prediction.png"
    pdf_path = OUTPUT_DIR / "figure_3_efcf_prediction.pdf"

    fig.savefig(
        png_path,
        dpi=400,
        bbox_inches="tight",
    )
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )
    plt.close(fig)

    print()
    print(f"Written: {png_path}")
    print(f"Written: {pdf_path}")
    print(f"Written: {sample_output}")
    print(f"Written: {summary_output}")


if __name__ == "__main__":
    main()
