from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RUNS_DIR = Path("outputs/runs")
GLOBAL_RESULTS_DIR = Path("global_results")


# Deterministic run-directory name:
#     (display name, country, environmental compartment)
from pfb_simulator.config import DATASET_REGISTRY

REQUIRED_COLUMNS = {
    "total_pfas_ng_l",
    "fluorine_burden_ng_l",
    "effective_f_fraction",
    "conversion_factor",
    "effective_k",
    "dominant_pfas",
}


NUMERIC_COLUMNS = [
    "total_pfas_ng_l",
    "fluorine_burden_ng_l",
    "effective_f_fraction",
    "conversion_factor",
    "effective_k",
    "dominant_fraction",
    "effective_molecular_weight",
]


def load_dataset(
    run_name: str,
    dataset: str,
    country: str,
    compartment: str,
) -> pd.DataFrame:
    """
    Load one deterministic individual-analysis result.
    """
    run_dir = RUNS_DIR / run_name
    csv_path = run_dir / "sample_burdens.csv"

    if not run_dir.exists():
        raise FileNotFoundError(
            f"Run directory does not exist:\n  {run_dir}\n"
            "Run scripts/run_all_individual_experiments.py first."
        )

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Expected result file does not exist:\n  {csv_path}"
        )

    df = pd.read_csv(csv_path)

    missing = REQUIRED_COLUMNS.difference(df.columns)

    if missing:
        raise ValueError(
            f"{csv_path} is missing required columns:\n"
            f"  {sorted(missing)}"
        )

    for column in NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(
                df[column],
                errors="coerce",
            )

    df["dataset"] = dataset
    df["country"] = country
    df["compartment"] = compartment
    df["source_run"] = run_name

    print(
        f"Loaded {dataset:24s}: "
        f"{len(df):>8,} samples from {csv_path}"
    )

    return df


def load_all_datasets() -> pd.DataFrame:
    """
    Load and concatenate all registered individual analyses.
    """
    frames: list[pd.DataFrame] = []

    for run_name, metadata in DATASET_REGISTRY.items():
        dataset, country, compartment = metadata

        frame = load_dataset(
            run_name=run_name,
            dataset=dataset,
            country=country,
            compartment=compartment,
        )

        frames.append(frame)

    if not frames:
        raise RuntimeError("No environmental datasets were loaded.")

    return pd.concat(
        frames,
        ignore_index=True,
        sort=False,
    )


def select_valid_samples(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select samples having a finite, positive EFCF.
    """
    valid = df[
        df["conversion_factor"].notna()
        & np.isfinite(df["conversion_factor"])
        & (df["conversion_factor"] > 0)
        & df["total_pfas_ng_l"].notna()
        & np.isfinite(df["total_pfas_ng_l"])
        & (df["total_pfas_ng_l"] > 0)
        & df["fluorine_burden_ng_l"].notna()
        & np.isfinite(df["fluorine_burden_ng_l"])
        & (df["fluorine_burden_ng_l"] > 0)
    ].copy()

    return valid


def most_common_non_null(values: pd.Series) -> str | None:
    """
    Return the most frequent non-null string value.
    """
    clean = values.dropna().astype(str)

    if clean.empty:
        return None

    counts = clean.value_counts()

    if counts.empty:
        return None

    return str(counts.index[0])


def summarize_datasets(
    all_samples: pd.DataFrame,
    valid_samples: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build one summary row per environmental dataset.
    """
    rows: list[dict] = []

    for run_name, metadata in DATASET_REGISTRY.items():
        dataset, country, compartment = metadata

        all_group = all_samples[
            all_samples["dataset"] == dataset
        ]

        valid_group = valid_samples[
            valid_samples["dataset"] == dataset
        ]

        rows.append(
            {
                "dataset": dataset,
                "country": country,
                "compartment": compartment,
                "n_samples": len(all_group),
                "n_valid_cf": len(valid_group),
                "median_total_pfas_ng_l":
                    valid_group["total_pfas_ng_l"].median(),
                "median_efcf":
                    valid_group["conversion_factor"].median(),
                "mean_efcf":
                    valid_group["conversion_factor"].mean(),
                "p05_efcf":
                    valid_group["conversion_factor"].quantile(0.05),
                "p95_efcf":
                    valid_group["conversion_factor"].quantile(0.95),
                "median_effective_k":
                    valid_group["effective_k"].median(),
                "dominant_pfas":
                    most_common_non_null(
                        valid_group["dominant_pfas"]
                    ),
            }
        )

    summary = pd.DataFrame(rows)

    return summary.sort_values(
        by=["country", "dataset"],
        kind="stable",
    ).reset_index(drop=True)


def summarize_global(valid_samples: pd.DataFrame) -> dict[str, float | int]:
    """
    Calculate statistics for the combined environmental population.
    """
    conversion_factor = valid_samples["conversion_factor"]

    return {
        "n_samples": int(len(valid_samples)),
        "median_efcf": float(conversion_factor.median()),
        "mean_efcf": float(conversion_factor.mean()),
        "std_efcf": float(conversion_factor.std()),
        "p05_efcf": float(conversion_factor.quantile(0.05)),
        "p25_efcf": float(conversion_factor.quantile(0.25)),
        "p75_efcf": float(conversion_factor.quantile(0.75)),
        "p95_efcf": float(conversion_factor.quantile(0.95)),
        "median_effective_f_fraction": float(
            valid_samples["effective_f_fraction"].median()
        ),
        "mean_effective_f_fraction": float(
            valid_samples["effective_f_fraction"].mean()
        ),
        "median_effective_k": float(
            valid_samples["effective_k"].median()
        ),
    }


def plot_global_efcf_histogram(
    valid_samples: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Plot the global EFCF distribution.
    """
    values = valid_samples["conversion_factor"].dropna()

    median = values.median()
    p05 = values.quantile(0.05)
    p95 = values.quantile(0.95)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))

    ax.hist(
        values,
        bins=120,
        density=True,
        alpha=0.75,
    )

    ax.axvline(
        median,
        linewidth=1.5,
        label=f"Median = {median:.3f}",
    )

    ax.axvline(
        p05,
        linewidth=1.0,
        linestyle="--",
        label=f"P05 = {p05:.3f}",
    )

    ax.axvline(
        p95,
        linewidth=1.0,
        linestyle="--",
        label=f"P95 = {p95:.3f}",
    )

    ax.set_xlabel("Effective Fluorine Conversion Factor")
    ax.set_ylabel("Probability density")
    ax.set_title("Global EFCF distribution")
    ax.grid(alpha=0.2)
    ax.legend()

    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_efcf_boxplot_by_dataset(
    valid_samples: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Exploratory dataset-level EFCF boxplot.

    Publication Figure 1 is generated by plot_efcf_boxplots.py.
    """
    datasets = list(DATASET_REGISTRY.keys())

    groups: list[np.ndarray] = []
    labels: list[str] = []

    for run_name in datasets:
        dataset = DATASET_REGISTRY[run_name][0]

        values = valid_samples.loc[
            valid_samples["dataset"] == dataset,
            "conversion_factor",
        ].dropna()

        if values.empty:
            continue

        groups.append(values.to_numpy())
        labels.append(dataset)

    fig, ax = plt.subplots(figsize=(13.0, 6.5))

    ax.boxplot(
        groups,
        tick_labels=labels,
        showfliers=False,
        whis=(5, 95),
    )

    global_median = valid_samples[
        "conversion_factor"
    ].median()

    ax.axhline(
        global_median,
        linestyle="--",
        linewidth=1.2,
        label=f"Global median = {global_median:.3f}",
    )

    ax.set_ylabel("Effective Fluorine Conversion Factor")
    ax.set_title("EFCF by monitoring dataset")
    ax.tick_params(axis="x", rotation=60)
    ax.grid(axis="y", alpha=0.2)
    ax.legend()

    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_efcf_boxplot_by_compartment(
    valid_samples: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Exploratory compartment-level EFCF boxplot.

    Publication Figure 1 is generated by plot_efcf_boxplots.py.
    """
    compartments = sorted(
        valid_samples["compartment"].dropna().unique()
    )

    groups: list[np.ndarray] = []
    labels: list[str] = []

    for compartment in compartments:
        values = valid_samples.loc[
            valid_samples["compartment"] == compartment,
            "conversion_factor",
        ].dropna()

        if values.empty:
            continue

        groups.append(values.to_numpy())
        labels.append(compartment)

    fig, ax = plt.subplots(figsize=(11.0, 6.5))

    ax.boxplot(
        groups,
        tick_labels=labels,
        showfliers=False,
        whis=(5, 95),
    )

    global_median = valid_samples[
        "conversion_factor"
    ].median()

    ax.axhline(
        global_median,
        linestyle="--",
        linewidth=1.2,
        label=f"Global median = {global_median:.3f}",
    )

    ax.set_ylabel("Effective Fluorine Conversion Factor")
    ax.set_title("EFCF by environmental compartment")
    ax.tick_params(axis="x", rotation=55)
    ax.grid(axis="y", alpha=0.2)
    ax.legend()

    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_dataset_centroids(
    dataset_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Plot median PFAS concentration against median EFCF.
    """
    plot_df = dataset_summary.dropna(
        subset=[
            "median_total_pfas_ng_l",
            "median_efcf",
        ]
    ).copy()

    fig, ax = plt.subplots(figsize=(9.0, 6.0))

    for _, row in plot_df.iterrows():
        ax.scatter(
            row["median_total_pfas_ng_l"],
            row["median_efcf"],
            s=45,
        )

        ax.annotate(
            row["dataset"],
            (
                row["median_total_pfas_ng_l"],
                row["median_efcf"],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7.5,
        )

    ax.set_xscale("log")
    ax.set_xlabel(
        r"Median total PFAS concentration (ng L$^{-1}$)"
    )
    ax.set_ylabel("Median EFCF")
    ax.set_title("Dataset centroids")
    ax.grid(alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_sample_ashby_space(
    valid_samples: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Plot reconstructed sample fluorine fraction against effective
    molecular weight when that quantity is present in the individual
    analysis outputs.
    """
    required = {
        "effective_molecular_weight",
        "effective_f_fraction",
    }

    if not required.issubset(valid_samples.columns):
        print(
            "Skipping sample Ashby plot: "
            "'effective_molecular_weight' is not present."
        )
        return

    plot_df = valid_samples.dropna(
        subset=[
            "effective_molecular_weight",
            "effective_f_fraction",
        ]
    ).copy()

    plot_df = plot_df[
        np.isfinite(plot_df["effective_molecular_weight"])
        & np.isfinite(plot_df["effective_f_fraction"])
        & (plot_df["effective_molecular_weight"] > 0)
        & (plot_df["effective_f_fraction"] > 0)
    ]

    if plot_df.empty:
        print(
            "Skipping sample Ashby plot: "
            "no valid effective molecular-weight values."
        )
        return

    fig, ax = plt.subplots(figsize=(9.0, 6.0))

    ax.scatter(
        plot_df["effective_molecular_weight"],
        plot_df["effective_f_fraction"],
        s=3,
        alpha=0.08,
        edgecolors="none",
        rasterized=True,
    )

    ax.set_xscale("log")
    ax.set_xlabel(
        r"Effective molecular weight (g mol$^{-1}$)"
    )
    ax.set_ylabel("Effective fluorine mass fraction")
    ax.set_title(
        "Environmental PFAS mixtures in molecular-property space"
    )
    ax.grid(alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    output_dir = (
        GLOBAL_RESULTS_DIR
        / f"global_analysis_{timestamp}"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    all_samples = load_all_datasets()
    valid_samples = select_valid_samples(all_samples)

    dataset_summary = summarize_datasets(
        all_samples=all_samples,
        valid_samples=valid_samples,
    )

    global_summary = summarize_global(valid_samples)

    all_samples.to_csv(
        output_dir / "all_samples.csv",
        index=False,
    )

    valid_samples.to_csv(
        output_dir / "valid_samples.csv",
        index=False,
    )

    dataset_summary.to_csv(
        output_dir / "dataset_summary.csv",
        index=False,
    )

    pd.DataFrame(
        [global_summary]
    ).to_csv(
        output_dir / "global_summary.csv",
        index=False,
    )

    pd.Series(
        global_summary
    ).to_json(
        output_dir / "global_summary.json",
        indent=2,
    )

    plot_global_efcf_histogram(
        valid_samples,
        output_dir / "efcf_histogram_global.png",
    )

    plot_efcf_boxplot_by_dataset(
        valid_samples,
        output_dir / "efcf_boxplot_by_dataset.png",
    )

    plot_efcf_boxplot_by_compartment(
        valid_samples,
        output_dir / "efcf_boxplot_by_compartment.png",
    )

    plot_dataset_centroids(
        dataset_summary,
        output_dir / "dataset_centroids.png",
    )

    plot_sample_ashby_space(
        valid_samples,
        output_dir
        / "ashby_samples_effective_mw_vs_f_fraction.png",
    )

    print()
    print(f"Output written to: {output_dir}")
    print()

    print(
        dataset_summary.to_string(
            index=False,
            float_format=lambda value: f"{value:.6f}",
        )
    )

    print()
    print("Global EFCF:")
    print(
        f"  median: {global_summary['median_efcf']:.6f}"
    )
    print(
        f"  mean  : {global_summary['mean_efcf']:.6f}"
    )
    print(
        f"  p05   : {global_summary['p05_efcf']:.6f}"
    )
    print(
        f"  p95   : {global_summary['p95_efcf']:.6f}"
    )


if __name__ == "__main__":
    main()