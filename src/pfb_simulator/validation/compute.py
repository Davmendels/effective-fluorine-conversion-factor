
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

from pfb_simulator.compounds import PFASCompound
from pfb_simulator.validation.pdh import PDHSample


SYNONYMS = {
    "PFU(n)DA": "PFUnDA",
    "PFUnA": "PFUnDA",
}


@dataclass(frozen=True)
class SampleBurden:
    sample_id: str
    date: str | None
    latitude: float | None
    longitude: float | None
    matrix: str | None

    detected_count: int
    matched_count: int
    unmatched_count: int

    total_pfas_ng_l: float
    fluorine_burden_ng_l: float
    effective_mw: float | None
    effective_f_fraction: float | None
    conversion_factor: float | None

    effective_k: float | None
    dominant_pfas: str | None
    dominant_fraction: float | None

    unmatched_pfas: list[str]

    def to_dict(self) -> dict:
        return asdict(self)


def normalize_name(name: str) -> str:
    name = name.strip()
    return SYNONYMS.get(name, name)


def make_compound_lookup(compounds: list[PFASCompound]) -> dict[str, PFASCompound]:
    return {
        compound.abbr: compound
        for compound in compounds
    }


def make_fraction_lookup(compounds: list[PFASCompound]) -> dict[str, float]:
    return {
        compound.abbr: compound.fluorine_fraction
        for compound in compounds
    }

def shannon_effective_k(values: list[float]) -> float | None:
    total = sum(values)

    if total <= 0:
        return None

    entropy = 0.0

    for value in values:
        if value <= 0:
            continue

        p = value / total
        entropy -= p * math.log(p)

    return math.exp(entropy)


def compute_sample_burden(
    sample: PDHSample,
    fluorine_fraction_lookup: dict[str, float],
    compound_lookup: dict[str, PFASCompound] | None = None,
) -> SampleBurden:
    total_pfas = 0.0
    fluorine_burden = 0.0
    weighted_mw_sum = 0.0

    matched_values = []
    unmatched = []

    for raw_name, concentration in sample.compounds.items():
        name = normalize_name(raw_name)

        if name not in fluorine_fraction_lookup:
            unmatched.append(raw_name)
            continue

        f_fraction = fluorine_fraction_lookup[name]

        total_pfas += concentration
        fluorine_burden += concentration * f_fraction
        matched_values.append(concentration)

        if compound_lookup is not None and name in compound_lookup:
            weighted_mw_sum += concentration * compound_lookup[name].mw

    effective_mw = (
        weighted_mw_sum / total_pfas
        if total_pfas > 0 and weighted_mw_sum > 0
        else None
    )

    if total_pfas > 0 and fluorine_burden > 0:
        effective_f = fluorine_burden / total_pfas
        cf = total_pfas / fluorine_burden
    else:
        effective_f = None
        cf = None

    if matched_values and total_pfas > 0:
        dominant_value = max(
            (sample.compounds[raw_name], normalize_name(raw_name))
            for raw_name in sample.compounds
            if normalize_name(raw_name) in fluorine_fraction_lookup
        )
        dominant_fraction = dominant_value[0] / total_pfas
        dominant_pfas = dominant_value[1]
    else:
        dominant_fraction = None
        dominant_pfas = None

    return SampleBurden(
        sample_id=sample.sample_id,
        date=sample.date,
        latitude=sample.latitude,
        longitude=sample.longitude,
        matrix=sample.matrix,
        detected_count=sample.detected_count(),
        matched_count=len(matched_values),
        unmatched_count=len(unmatched),
        total_pfas_ng_l=total_pfas,
        fluorine_burden_ng_l=fluorine_burden,
        effective_mw=effective_mw,
        effective_f_fraction=effective_f,
        conversion_factor=cf,
        effective_k=shannon_effective_k(matched_values),
        dominant_pfas=dominant_pfas,
        dominant_fraction=dominant_fraction,
        unmatched_pfas=sorted(set(unmatched)),
    )
