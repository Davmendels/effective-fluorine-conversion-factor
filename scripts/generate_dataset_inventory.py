from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from pfb_simulator.config import DATASET_REGISTRY


# Adjust this date if any files were downloaded on a different day.
DEFAULT_RETRIEVAL_DATE = "2026-07-08"


# Dataset name:
#     local filename,
#     input format,
#     source name,
#     source URL,
#     provenance/preprocessing note
DATASET_INPUTS: dict[
    str,
    tuple[str, str, str, str, str],
] = {
    "Flanders": (
        "pdh_export_Flanders.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Harmonized PDH export.",
    ),
    "France_ADES": (
        "pdh_export_France_ADES.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Groundwater monitoring; harmonized PDH export.",
    ),
    "France_ANSES": (
        "pdh_export_France_ANSES.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Drinking-water monitoring; harmonized PDH export.",
    ),
    "France_EauRob": (
        "pdh_export_France_EauRob.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Tap-water monitoring; harmonized PDH export.",
    ),
    "France_ICPE": (
        "pdh_export_France_ICPE.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Industrial-discharge monitoring; harmonized PDH export.",
    ),
    "France_Naiades": (
        "pdh_export_France_Naiades.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Surface-water monitoring; harmonized PDH export.",
    ),
    "Grand_Lyon": (
        "pdh_export_Grand_Lyon.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Contaminated soil extracts; harmonized PDH export.",
    ),
    "Germany_ELWAS": (
        "pdh_export_DE_ELWAS.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "German water-monitoring data; harmonized PDH export.",
    ),
    "Germany_TFA": (
        "pdh_export_DE_TFA.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Dedicated monitoring dataset containing TFA.",
    ),
    "Italy_Friuli_Venezia": (
        "pdh_export_IT_ARPA_Friuli_Venezia.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Regional Italian monitoring; harmonized PDH export.",
    ),
    "Muir": (
        "pdh_export_Muir.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Marine-water compilation distributed through PDH.",
    ),
    "RIVM": (
        "pdh_export_RIVM.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Dutch monitoring data; harmonized PDH export.",
    ),
    "UK_EA": (
        "pdh_export_UK_EA.csv.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "UK Environment Agency data; harmonized PDH export.",
    ),
    "Veneto": (
        "pdh_export_Veneti.gz",
        "pdh",
        "PFAS Data Hub",
        "https://pdh.cnrs.fr/en/datasets/",
        "Regional Italian monitoring; harmonized PDH export.",
    ),
    "UCMR5_extract": (
        "0609b147-293e-4913-86ed-c5ce0bbd20cb.xlsx",
        "ucmr5",
        "U.S. EPA UCMR 5 occurrence data",
        "https://www.epa.gov/dwucmr/occurrence-data-unregulated-contaminant-monitoring-rule",
        "PFAS-containing records extracted and converted to the common internal format.",
    ),
    "WQP_multimedia": (
        "c81d295f-dfc7-4c90-a88d-f7d8f50e2c6b.xlsx",
        "wqp",
        "Water Quality Portal",
        "https://www.waterqualitydata.us/",
        "PFAS multimedia records extracted and converted to the common internal format.",
    ),
}


FIELDNAMES = [
    "dataset",
    "country",
    "compartment",
    "source",
    "source_url",
    "local_filename",
    "input_format",
    "retrieved_on",
    "file_size_bytes",
    "sha256",
    "n_samples",
    "n_valid_efcf",
    "notes",
]


def sha256_file(path: Path) -> str:
    """
    Calculate the SHA-256 checksum without loading the complete file
    into memory.
    """
    digest = hashlib.sha256()

    with path.open("rb") as file_handle:
        for block in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def read_summary(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing analysis summary:\n  {path}\n"
            "Run scripts/run_all_individual_experiments.py first."
        )

    with path.open("r", encoding="utf-8") as file_handle:
        return json.load(file_handle)


def validate_configuration() -> None:
    registered = set(DATASET_REGISTRY)
    configured = set(DATASET_INPUTS)

    missing_inputs = registered - configured
    unknown_inputs = configured - registered

    if missing_inputs or unknown_inputs:
        messages: list[str] = []

        if missing_inputs:
            messages.append(
                "Missing DATASET_INPUTS entries: "
                + ", ".join(sorted(missing_inputs))
            )

        if unknown_inputs:
            messages.append(
                "Unknown DATASET_INPUTS entries: "
                + ", ".join(sorted(unknown_inputs))
            )

        raise ValueError("\n".join(messages))


def build_inventory(
    raw_dir: Path,
    runs_dir: Path,
    retrieval_date: str,
) -> list[dict[str, Any]]:
    validate_configuration()

    rows: list[dict[str, Any]] = []

    for run_name, metadata in DATASET_REGISTRY.items():
        dataset, country, compartment = metadata

        (
            filename,
            input_format,
            source,
            source_url,
            notes,
        ) = DATASET_INPUTS[run_name]

        raw_path = raw_dir / filename
        summary_path = runs_dir / run_name / "summary.json"

        if not raw_path.is_file():
            raise FileNotFoundError(
                f"Missing raw dataset for {dataset}:\n  {raw_path}"
            )

        summary = read_summary(summary_path)

        n_samples = summary.get("n_samples")
        n_valid = summary.get("n_valid_cf")

        if n_samples is None or n_valid is None:
            raise ValueError(
                f"Required sample counts are absent from:\n"
                f"  {summary_path}"
            )

        row = {
            "dataset": dataset,
            "country": country,
            "compartment": compartment,
            "source": source,
            "source_url": source_url,
            "local_filename": filename,
            "input_format": input_format,
            "retrieved_on": retrieval_date,
            "file_size_bytes": raw_path.stat().st_size,
            "sha256": sha256_file(raw_path),
            "n_samples": int(n_samples),
            "n_valid_efcf": int(n_valid),
            "notes": notes,
        }

        rows.append(row)

        print(
            f"Indexed {dataset:24s}: "
            f"{int(n_samples):>8,} samples, "
            f"{int(n_valid):>8,} valid EFCF"
        )

    return rows


def write_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=FIELDNAMES,
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the environmental dataset inventory used in "
            "the EFCF manuscript."
        )
    )

    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("datasets/raw"),
        help="Directory containing the original input files.",
    )

    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path("outputs/runs"),
        help="Directory containing deterministic dataset analyses.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("supplementary/dataset_inventory.csv"),
        help="Output inventory CSV.",
    )

    parser.add_argument(
        "--retrieved-on",
        default=DEFAULT_RETRIEVAL_DATE,
        help="Dataset retrieval date in YYYY-MM-DD format.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    rows = build_inventory(
        raw_dir=args.raw_dir,
        runs_dir=args.runs_dir,
        retrieval_date=args.retrieved_on,
    )

    write_csv(rows, args.output)

    total_samples = sum(
        int(row["n_samples"])
        for row in rows
    )
    total_valid = sum(
        int(row["n_valid_efcf"])
        for row in rows
    )

    print()
    print(f"Written: {args.output}")
    print(f"Datasets: {len(rows)}")
    print(f"Total source samples: {total_samples:,}")
    print(f"Samples with valid EFCF: {total_valid:,}")


if __name__ == "__main__":
    main()
