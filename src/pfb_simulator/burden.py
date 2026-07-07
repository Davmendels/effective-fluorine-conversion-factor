from __future__ import annotations

import numpy as np

from pfb_simulator.compounds import PFASCompound


def weighted_fluorine_fraction(
    compounds: list[PFASCompound],
    weights: np.ndarray,
) -> float:
    if len(compounds) != len(weights):
        raise ValueError("Number of compounds and weights must match.")

    weights = np.asarray(weights, dtype=float)

    if np.any(weights < 0):
        raise ValueError("Weights cannot be negative.")

    total = weights.sum()

    if total <= 0:
        raise ValueError("Weight sum must be positive.")

    weights = weights / total

    fractions = np.array([c.fluorine_fraction for c in compounds], dtype=float)

    return float(np.sum(weights * fractions))


def fluorine_burden_from_pfas_mass(
    total_pfas_ng_l: float,
    fluorine_fraction: float,
) -> float:
    return total_pfas_ng_l * fluorine_fraction


def pfas_mass_from_fluorine_burden(
    fluorine_burden_ng_l: float,
    fluorine_fraction: float,
) -> float:
    if fluorine_fraction <= 0:
        raise ValueError("Fluorine fraction must be positive.")

    return fluorine_burden_ng_l / fluorine_fraction


def conversion_factor_from_fluorine_fraction(fluorine_fraction: float) -> float:
    if fluorine_fraction <= 0:
        raise ValueError("Fluorine fraction must be positive.")

    return 1.0 / fluorine_fraction
