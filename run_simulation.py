from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

from pfb_simulator.compounds import load_compounds, filter_compounds
from pfb_simulator.monte_carlo import run_dirichlet_monte_carlo
from pfb_simulator.plotting import plot_histogram
from pfb_simulator.statistics import (
    compute_mixture_bounds,
    summarize_distribution,
)
from pfb_simulator.validation import validate_compounds


def main() -> None:
    n_iterations = 1_000_000
    alpha = 1.0
    seed = 42

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("outputs/runs") / f"run_{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)

    compounds = load_compounds("data/pfas_reference.json")
    validate_compounds(compounds)

    compounds = filter_compounds(compounds, exclude_abbrs={"TFA"})

    result = run_dirichlet_monte_carlo(
        compounds,
        n_iterations=n_iterations,
        alpha=alpha,
        seed=seed,
        batch_size=100_000,
    )

    f_summary = summarize_distribution(result.fluorine_fractions)
    cf_summary = summarize_distribution(result.conversion_factors)
    bounds = compute_mixture_bounds(compounds)

    config = {
        "n_iterations": n_iterations,
        "alpha": alpha,
        "seed": seed,
        "excluded_abbrs": ["TFA"],
        "n_compounds": len(compounds),
    }

    summary = {
        "config": config,
        "fluorine_fraction": f_summary.to_dict(),
        "conversion_factor": cf_summary.to_dict(),
        "hard_bounds": bounds.to_dict(),
    }

    with (output_dir / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    pd.DataFrame(
        {
            "fluorine_fraction": result.fluorine_fractions,
            "conversion_factor": result.conversion_factors,
        }
    ).to_csv(output_dir / "samples.csv", index=False)

    plot_histogram(
        100.0 * result.fluorine_fractions,
        output_dir / "fluorine_fraction_histogram.png",
        title="Monte Carlo distribution of weighted fluorine fraction",
        xlabel="Weighted fluorine fraction (%)",
        bins=100,
    )

    plot_histogram(
        result.conversion_factors,
        output_dir / "conversion_factor_histogram.png",
        title="Monte Carlo distribution of TFB-to-PFAS conversion factor",
        xlabel="Conversion factor: PFAS mass / fluorine burden",
        bins=100,
    )

    print(f"Output written to: {output_dir}")
    print()
    print("Fluorine fraction")
    print(json.dumps(f_summary.to_dict(), indent=2))
    print()
    print("Conversion factor")
    print(json.dumps(cf_summary.to_dict(), indent=2))
    print()
    print("Hard bounds")
    print(json.dumps(bounds.to_dict(), indent=2))


if __name__ == "__main__":
    main()
