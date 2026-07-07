from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from pfb_simulator.compounds import load_compounds
from pfb_simulator.monte_carlo import run_subset_size_experiment
from pfb_simulator.validation import validate_compounds


def main() -> None:
    compounds = load_compounds("data/pfas_extended_environmental.json")
    validate_compounds(compounds)

    k_values = [1, 2, 3, 5, 10, 20, 50, 100, 250, 500]
    n_iterations_per_k = 50_000
    alpha = 1.0

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("outputs/runs") / f"subset_k_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = run_subset_size_experiment(
        compounds,
        k_values=k_values,
        n_iterations_per_k=n_iterations_per_k,
        alpha=alpha,
        seed=42,
    )

    df = pd.DataFrame(rows)
    df.to_csv(output_dir / "subset_size_summary.csv", index=False)

    print(df.to_string(index=False))

    # Plot fluorine fraction versus k
    fig, ax = plt.subplots(figsize=(9, 5))

    x = df["k"]
    y = 100.0 * df["p50_f"]
    yerr_lower = 100.0 * (df["p50_f"] - df["p05_f"])
    yerr_upper = 100.0 * (df["p95_f"] - df["p50_f"])

    ax.errorbar(
        x,
        y,
        yerr=[yerr_lower, yerr_upper],
        fmt="o-",
        capsize=4,
        label="5–95% interval",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Effective number of PFAS contributors, k")
    ax.set_ylabel("Weighted fluorine fraction (%)")
    ax.set_title("Uncertainty collapse with increasing PFAS mixture complexity")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "fluorine_fraction_vs_k.png", dpi=200)
    plt.close(fig)

    # Plot conversion factor versus k
    fig, ax = plt.subplots(figsize=(9, 5))

    y = df["p50_cf"]
    yerr_lower = df["p50_cf"] - df["p05_cf"]
    yerr_upper = df["p95_cf"] - df["p50_cf"]

    ax.errorbar(
        x,
        y,
        yerr=[yerr_lower, yerr_upper],
        fmt="o-",
        capsize=4,
        label="5–95% interval",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Effective number of PFAS contributors, k")
    ax.set_ylabel("TFB-to-PFAS conversion factor")
    ax.set_title("Conversion-factor uncertainty versus mixture complexity")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "conversion_factor_vs_k.png", dpi=200)
    plt.close(fig)

    print()
    print(f"Output written to: {output_dir}")


if __name__ == "__main__":
    main()
