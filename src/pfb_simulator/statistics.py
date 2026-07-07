from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from pfb_simulator.compounds import PFASCompound


@dataclass(frozen=True)
class DistributionSummary:
    mean: float
    median: float
    std: float
    min: float
    p01: float
    p05: float
    p95: float
    p99: float
    max: float

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class MixtureBounds:
    min_compound_abbr: str
    min_compound_name: str
    min_fluorine_fraction: float
    min_conversion_factor: float

    max_compound_abbr: str
    max_compound_name: str
    max_fluorine_fraction: float
    max_conversion_factor: float

    def to_dict(self) -> dict:
        return asdict(self)


def summarize_distribution(values: np.ndarray) -> DistributionSummary:
    values = np.asarray(values, dtype=float)

    if values.size == 0:
        raise ValueError("Cannot summarize an empty array.")

    return DistributionSummary(
        mean=float(np.mean(values)),
        median=float(np.median(values)),
        std=float(np.std(values)),
        min=float(np.min(values)),
        p01=float(np.quantile(values, 0.01)),
        p05=float(np.quantile(values, 0.05)),
        p95=float(np.quantile(values, 0.95)),
        p99=float(np.quantile(values, 0.99)),
        max=float(np.max(values)),
    )


def compute_mixture_bounds(compounds: list[PFASCompound]) -> MixtureBounds:
    if not compounds:
        raise ValueError("Compound list cannot be empty.")

    min_compound = min(compounds, key=lambda c: c.fluorine_fraction)
    max_compound = max(compounds, key=lambda c: c.fluorine_fraction)

    return MixtureBounds(
        min_compound_abbr=min_compound.abbr,
        min_compound_name=min_compound.name,
        min_fluorine_fraction=min_compound.fluorine_fraction,
        min_conversion_factor=1.0 / min_compound.fluorine_fraction,
        max_compound_abbr=max_compound.abbr,
        max_compound_name=max_compound.name,
        max_fluorine_fraction=max_compound.fluorine_fraction,
        max_conversion_factor=1.0 / max_compound.fluorine_fraction,
    )
