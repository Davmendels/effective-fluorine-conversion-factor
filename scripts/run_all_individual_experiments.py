from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


DATA_DIR = Path("datasets/raw")
RUNS_DIR = Path("outputs/runs")


DATASETS = [
    ("pdh_export_Veneti.gz", "pdh", "Veneto"),
    (
        "pdh_export_IT_ARPA_Friuli_Venezia.csv.gz",
        "pdh",
        "Italy_Friuli_Venezia",
    ),
    ("pdh_export_UK_EA.csv.gz", "pdh", "UK_EA"),
    ("pdh_export_France_ADES.csv.gz", "pdh", "France_ADES"),
    (
        "pdh_export_France_Naiades.csv.gz",
        "pdh",
        "France_Naiades",
    ),
    (
        "pdh_export_France_ANSES.csv.gz",
        "pdh",
        "France_ANSES",
    ),
    (
        "pdh_export_France_EauRob.csv.gz",
        "pdh",
        "France_EauRob",
    ),
    (
        "pdh_export_France_ICPE.csv.gz",
        "pdh",
        "France_ICPE",
    ),
    ("pdh_export_Flanders.csv.gz", "pdh", "Flanders"),
    (
        "pdh_export_DE_ELWAS.csv.gz",
        "pdh",
        "Germany_ELWAS",
    ),
    ("pdh_export_RIVM.csv.gz", "pdh", "RIVM"),
    (
        "pdh_export_Grand_Lyon.csv.gz",
        "pdh",
        "Grand_Lyon",
    ),
    ("pdh_export_Muir.csv.gz", "pdh", "Muir"),
    (
        "pdh_export_DE_TFA.csv.gz",
        "pdh",
        "Germany_TFA",
    ),
    (
        "0609b147-293e-4913-86ed-c5ce0bbd20cb.xlsx",
        "ucmr5",
        "UCMR5_extract",
    ),
    (
        "c81d295f-dfc7-4c90-a88d-f7d8f50e2c6b.xlsx",
        "wqp",
        "WQP_multimedia",
    ),
]


def main() -> None:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)

    environment = os.environ.copy()

    existing_pythonpath = environment.get("PYTHONPATH")

    if existing_pythonpath:
        environment["PYTHONPATH"] = (
            f"src{os.pathsep}{existing_pythonpath}"
        )
    else:
        environment["PYTHONPATH"] = "src"

    for filename, dataset_format, output_name in DATASETS:
        input_path = DATA_DIR / filename
        output_dir = RUNS_DIR / output_name

        if not input_path.exists():
            raise FileNotFoundError(
                f"Required dataset not found: {input_path}"
            )

        print()
        print("=" * 80)
        print(
            f"Running {filename} [{dataset_format}] "
            f"-> {output_dir}"
        )
        print("=" * 80)

        subprocess.run(
            [
                sys.executable,
                "scripts/analyze_pdh_dataset.py",
                str(input_path),
                "--format",
                dataset_format,
                "--output-dir",
                str(output_dir),
            ],
            check=True,
            env=environment,
        )

    print()
    print("All individual analyses completed successfully.")


if __name__ == "__main__":
    main()