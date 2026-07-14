from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class PDHSample:
    sample_id: str
    date: str | None = None
    latitude: float | None = None
    longitude: float | None = None
    matrix: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)

    # Detected concentrations, normally ng/L.
    compounds: dict[str, float] = field(default_factory=dict)

    # Censored/non-detect values, e.g. {"PFOA": 10.0} for <10 ng/L.
    less_than: dict[str, float] = field(default_factory=dict)

    def detected_count(self) -> int:
        return len(self.compounds)

    def censored_count(self) -> int:
        return len(self.less_than)

    def total_detected_concentration(self) -> float:
        return sum(self.compounds.values())


def _read_text(path: str | Path) -> str:
    path = Path(path)

    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as f:
            return f.read()

    with path.open("r", encoding="utf-8") as f:
        return f.read()


def _load_raw_records(path: str | Path) -> list[dict[str, Any]]:
    """
    Load a PDH export.

    Supports:
    - gzip-compressed JSON
    - plain JSON
    - gzip-compressed CSV
    - plain CSV

    The exact PDH format can vary, so this loader is intentionally tolerant.
    """

    path = Path(path)
    text = _read_text(path)
    stripped = text.lstrip()

    if stripped.startswith("["):
        data = json.loads(text)
        if not isinstance(data, list):
            raise ValueError("Expected top-level JSON list.")
        return data

    if stripped.startswith("{"):
        data = json.loads(text)

        if isinstance(data, dict):
            for key in ["data", "records", "rows", "features"]:
                if key in data and isinstance(data[key], list):
                    return data[key]

        raise ValueError("Could not find record list in JSON object.")

    # Fallback: CSV
    compression = "gzip" if path.suffix == ".gz" else None
    df = pd.read_csv(path, compression=compression)
    return df.to_dict(orient="records")


def _coerce_float(value: Any) -> float | None:
    if value is None:
        return None

    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None
        value = value.replace(",", ".")

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _first_present(record: dict[str, Any], keys: list[str]) -> Any:
    for key in keys:
        if key in record and record[key] not in [None, ""]:
            return record[key]
    return None


def _parse_pfas_entries(value: Any) -> tuple[dict[str, float], dict[str, float]]:
    """
    Parse PFAS compound entries.

    Expected PDH-like form:
    [
      {"substance": "PFOA", "value": 12.3},
      {"substance": "PFOS", "less_than": "10.0"}
    ]

    Returns:
    - detected compounds
    - censored / less-than compounds
    """

    detected: dict[str, float] = {}
    less_than: dict[str, float] = {}

    if value is None:
        return detected, less_than

    if isinstance(value, str):
        value = value.strip()
        if not value:
            return detected, less_than

        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            return detected, less_than

    if not isinstance(value, list):
        return detected, less_than

    for entry in value:
        if not isinstance(entry, dict):
            continue

        substance = (
            entry.get("substance")
            or entry.get("name")
            or entry.get("abbr")
            or entry.get("parameter")
            or entry.get("compound")
        )

        if not substance:
            continue

        substance = str(substance).strip()

        measured_value = _coerce_float(entry.get("value"))

        if measured_value is not None and measured_value > 0:
            detected[substance] = measured_value
            continue

        lt_value = (
            entry.get("less_than")
            or entry.get("lessThan")
            or entry.get("loq")
            or entry.get("LOQ")
            or entry.get("limit")
        )

        lt_value = _coerce_float(lt_value)

        if lt_value is not None:
            less_than[substance] = lt_value

    return detected, less_than


def record_to_sample(record: dict[str, Any], index: int) -> PDHSample:
    sample_id = _first_present(
        record,
        ["sample_id", "id", "ID", "station_id", "code", "uuid"],
    )

    if sample_id is None:
        sample_id = f"sample_{index}"

    date = _first_present(
        record,
        ["date", "sampling_date", "sample_date", "datetime", "timestamp"],
    )

    latitude = _coerce_float(
        _first_present(record, ["latitude", "lat", "Latitude", "LAT"])
    )

    longitude = _coerce_float(
        _first_present(record, ["longitude", "lon", "lng", "Longitude", "LON"])
    )

    matrix = _first_present(
        record,
        ["matrix", "medium", "sample_matrix", "water_type"],
    )

    pfas_value = _first_present(
        record,
        [
            "pfas_values",
            "pfas",
            "PFAS",
            "measurements",
            "substances",
            "values",
            "data",
        ],
    )

    compounds, less_than = _parse_pfas_entries(pfas_value)

    return PDHSample(
        sample_id=str(sample_id),
        date=str(date) if date is not None else None,
        latitude=latitude,
        longitude=longitude,
        matrix=str(matrix) if matrix is not None else None,
        raw=record,
        compounds=compounds,
        less_than=less_than,
    )


def load_pdh_samples(path: str | Path) -> list[PDHSample]:
    records = _load_raw_records(path)

    return [
        record_to_sample(record, index=i)
        for i, record in enumerate(records)
        if isinstance(record, dict)
    ]


def summarize_samples(samples: list[PDHSample]) -> dict[str, Any]:
    detected_counts = [s.detected_count() for s in samples]
    censored_counts = [s.censored_count() for s in samples]
    total_detected = [s.total_detected_concentration() for s in samples]

    return {
        "n_samples": len(samples),
        "samples_with_detected_pfas": sum(1 for x in detected_counts if x > 0),
        "samples_with_censored_pfas": sum(1 for x in censored_counts if x > 0),
        "mean_detected_count": sum(detected_counts) / len(samples) if samples else 0,
        "mean_censored_count": sum(censored_counts) / len(samples) if samples else 0,
        "max_detected_count": max(detected_counts) if detected_counts else 0,
        "max_censored_count": max(censored_counts) if censored_counts else 0,
        "mean_total_detected_concentration": (
            sum(total_detected) / len(samples) if samples else 0
        ),
        "max_total_detected_concentration": max(total_detected) if total_detected else 0,
    }
