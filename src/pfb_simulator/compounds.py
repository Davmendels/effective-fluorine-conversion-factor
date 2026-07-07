from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


F_ATOMIC_MASS = 18.998403163


@dataclass(frozen=True)
class PFASCompound:
    name: str
    abbr: str
    mw: float
    fluorine_atoms: int
    category: str

    @property
    def fluorine_mass(self) -> float:
        return self.fluorine_atoms * F_ATOMIC_MASS

    @property
    def fluorine_fraction(self) -> float:
        return self.fluorine_mass / self.mw

    @property
    def fluorine_percent(self) -> float:
        return 100.0 * self.fluorine_fraction


def load_compounds(path: str | Path) -> list[PFASCompound]:
    path = Path(path)

    with path.open("r", encoding="utf-8") as f:
        records = json.load(f)

    return [
        PFASCompound(
            name=record["name"],
            abbr=record["abbr"],
            mw=float(record["mw"]),
            fluorine_atoms=int(record["fluorine_atoms"]),
            category=record["category"],
        )
        for record in records
    ]


def filter_compounds(
    compounds: list[PFASCompound],
    *,
    exclude_abbrs: set[str] | None = None,
    include_categories: set[str] | None = None,
    exclude_categories: set[str] | None = None,
) -> list[PFASCompound]:
    exclude_abbrs = exclude_abbrs or set()
    exclude_categories = exclude_categories or set()

    filtered = []

    for compound in compounds:
        if compound.abbr in exclude_abbrs:
            continue

        if compound.category in exclude_categories:
            continue

        if include_categories is not None and compound.category not in include_categories:
            continue

        filtered.append(compound)

    return filtered


def compounds_to_rows(compounds: list[PFASCompound]) -> list[dict]:
    return [
        {
            "abbr": c.abbr,
            "name": c.name,
            "category": c.category,
            "mw": c.mw,
            "fluorine_atoms": c.fluorine_atoms,
            "fluorine_mass": c.fluorine_mass,
            "fluorine_fraction": c.fluorine_fraction,
            "fluorine_percent": c.fluorine_percent,
        }
        for c in compounds
    ]


def print_compound_table(compounds: list[PFASCompound]) -> None:
    header = (
        f"{'Abbr':<14}"
        f"{'Category':<14}"
        f"{'MW':>10}"
        f"{'F atoms':>10}"
        f"{'F %':>10}"
    )
    print(header)
    print("-" * len(header))

    for c in compounds:
        print(
            f"{c.abbr:<14}"
            f"{c.category:<14}"
            f"{c.mw:>10.2f}"
            f"{c.fluorine_atoms:>10d}"
            f"{c.fluorine_percent:>10.2f}"
        )


if __name__ == "__main__":
    compounds = load_compounds("data/pfas_reference.json")
    print_compound_table(compounds)
