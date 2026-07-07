from __future__ import annotations

from collections import Counter

from pfb_simulator.compounds import PFASCompound


def validate_compounds(compounds: list[PFASCompound]) -> None:
    if not compounds:
        raise ValueError("Compound list is empty.")

    abbr_counts = Counter(c.abbr for c in compounds)
    duplicates = [abbr for abbr, count in abbr_counts.items() if count > 1]

    if duplicates:
        raise ValueError(f"Duplicate compound abbreviations found: {duplicates}")

    for c in compounds:
        if c.mw <= 0:
            raise ValueError(f"{c.abbr}: MW must be positive.")

        if c.fluorine_atoms < 1:
            raise ValueError(f"{c.abbr}: fluorine_atoms must be >= 1.")

        if not (0.20 <= c.fluorine_fraction <= 0.85):
            raise ValueError(
                f"{c.abbr}: suspicious fluorine fraction "
                f"{c.fluorine_fraction:.3f}"
            )
