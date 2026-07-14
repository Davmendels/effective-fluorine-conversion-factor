from __future__ import annotations

import argparse
import shutil
import tempfile
from datetime import datetime
from pathlib import Path

import pandas as pd

from pfb_simulator.compounds import load_compounds
from pfb_simulator.validation.compute import (
    compute_sample_burden,
    make_compound_lookup,
    make_fraction_lookup,
)
from pfb_simulator.validation.pdh import load_pdh_samples


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("input_path")

    parser.add_argument(
        "--format",
        choices=["pdh", "ucmr5", "wqp"],
        default="pdh",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Deterministic output directory. If omitted, a timestamped "
            "development directory is created."
        ),
    )

    args = parser.parse_args()

    input_path = Path(args.input_path)

    if args.format == "pdh":
        samples = load_pdh_samples(input_path)

    elif args.format == "ucmr5":
        from pfb_simulator.validation.ucmr5 import load_ucmr5_samples

        samples = load_ucmr5_samples(input_path)

    elif args.format == "wqp":
        from pfb_simulator.validation.wqp import load_wqp_samples

        samples = load_wqp_samples(input_path)

    else:
        raise ValueError(args.format)

    compounds = load_compounds("data/pfas_reference.json")

    fluorine_fraction_lookup = make_fraction_lookup(compounds)
    compound_lookup = make_compound_lookup(compounds)

    rows = [
        compute_sample_burden(
            sample,
            fluorine_fraction_lookup,
            compound_lookup,
        ).to_dict()
        for sample in samples
    ]

    df = pd.DataFrame(rows)

    if args.output_dir is None:
        # Timestamped fallback for exploratory/development runs.
        safe_stem = input_path.name
        safe_stem = safe_stem.replace(".csv.gz", "")
        safe_stem = safe_stem.replace(".xlsx", "")
        safe_stem = safe_stem.replace(".gz", "")
        safe_stem = "".join(
            character
            if character.isalnum() or character in "-_."
            else "_"
            for character in safe_stem
        )

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")

        output_dir = (
            Path("outputs/runs")
            / f"pdh_analysis_{safe_stem}_{timestamp}"
        )

        output_dir.mkdir(parents=True, exist_ok=False)
        working_dir = output_dir
        replace_existing = False

    else:
        # Stable output for the reproducible paper workflow.
        output_dir = args.output_dir
        output_dir.parent.mkdir(parents=True, exist_ok=True)

        working_dir = Path(
            tempfile.mkdtemp(
                prefix=f".{output_dir.name}_",
                dir=output_dir.parent,
            )
        )

        replace_existing = True

    df.to_csv(
        working_dir / "sample_burdens.csv",
        index=False,
    )

    valid = df[df["conversion_factor"].notna()].copy()

    summary = {
        "n_samples": len(df),
        "n_valid_cf": len(valid),
        "median_total_pfas_ng_l":
            valid["total_pfas_ng_l"].median(),
        "mean_total_pfas_ng_l":
            valid["total_pfas_ng_l"].mean(),
        "median_fluorine_burden_ng_l":
            valid["fluorine_burden_ng_l"].median(),
        "mean_fluorine_burden_ng_l":
            valid["fluorine_burden_ng_l"].mean(),
        "median_effective_f_fraction":
            valid["effective_f_fraction"].median(),
        "mean_effective_f_fraction":
            valid["effective_f_fraction"].mean(),
        "median_conversion_factor":
            valid["conversion_factor"].median(),
        "mean_conversion_factor":
            valid["conversion_factor"].mean(),
        "p05_conversion_factor":
            valid["conversion_factor"].quantile(0.05),
        "p95_conversion_factor":
            valid["conversion_factor"].quantile(0.95),
        "median_effective_k":
            valid["effective_k"].median(),
        "mean_effective_k":
            valid["effective_k"].mean(),
    }

    pd.Series(summary).to_json(
        working_dir / "summary.json",
        indent=2,
    )

    # Replace the previous deterministic result only after every output
    # has been generated successfully.
    if replace_existing:
        if output_dir.exists():
            shutil.rmtree(output_dir)

        working_dir.rename(output_dir)

    print(f"Output written to: {output_dir}")
    print()

    for key, value in summary.items():
        print(f"{key}: {value}")

    print()
    print("Dominant PFAS counts:")
    print(valid["dominant_pfas"].value_counts().head(20))


if __name__ == "__main__":
    main()