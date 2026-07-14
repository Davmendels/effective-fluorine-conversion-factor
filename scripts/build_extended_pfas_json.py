from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


INPUT_CSV = Path("datasets/raw/Chemical_List_PFASSTRUCTV6-2026-07-07.csv")
OUTPUT_JSON = Path("data/pfas_extended.json")
ENV_OUTPUT_JSON = Path("data/pfas_extended_environmental.json")

SOURCE = "EPA CompTox PFASSTRUCTV6 2026-07-07"


def count_f_atoms(formula: str) -> int:
    if not isinstance(formula, str) or not formula.strip():
        return 0

    # Matches F, F2, F15, etc.
    match = re.search(r"F(\d*)", formula)

    if not match:
        return 0

    count = match.group(1)
    return int(count) if count else 1


def clean_string(value) -> str | None:
    if pd.isna(value):
        return None

    value = str(value).strip()
    return value if value else None

def looks_environmental(record: dict) -> bool:
    name = (record.get("name") or "").lower()
    formula = record.get("formula") or ""
    mw = float(record["mw"])
    f_fraction = float(record["fluorine_fraction"])

    acid_terms = [
        "carboxylic acid",
        "sulfonic acid",
        "sulfonate",
        "phosphonic acid",
        "phosphonate",
    ]

    known_terms = [
        "fluorotelomer",
        "perfluoro",
        "hfpo",
        "genx",
        "adona",
        "pfesa",
        "chloroperfluoro",
    ]

    if not any(term in name for term in acid_terms + known_terms):
        return False

    # Keep the water/environmental molecular window.
    if not (150.0 <= mw <= 800.0):
        return False

    # Exclude very weakly fluorinated functionalized species.
    if not (0.45 <= f_fraction <= 0.80):
        return False

    # Avoid obvious polymers/mixtures/salts in names.
    bad_terms = [
        "polymer",
        "copolymer",
        "homopolymer",
        "resin",
        "mixture",
        "reaction products",
        "salt with",
    ]

    if any(term in name for term in bad_terms):
        return False

    return True
    
def main() -> None:
    df = pd.read_csv(INPUT_CSV)

    records = []

    for _, row in df.iterrows():
        name = clean_string(row.get("PREFERRED NAME"))
        casrn = clean_string(row.get("CASRN"))
        formula = clean_string(row.get("MOLECULAR FORMULA"))
        smiles = clean_string(row.get("SMILES"))
        qc_level = clean_string(row.get("QC Level"))

        mw = row.get("AVERAGE MASS")

        if pd.isna(mw):
            continue

        mw = float(mw)
        f_atoms = count_f_atoms(formula)

        if not name:
            continue

        if mw <= 0:
            continue

        if f_atoms <= 0:
            continue

        f_fraction = f_atoms * 18.998403163 / mw

        # Conservative sanity filter.
        # Removes salts, mixtures, polymers, odd records, and non-discrete artefacts.
        if not (150.0 <= mw <= 900.0):
            continue

        if f_atoms < 3:
            continue

        if not (0.35 <= f_fraction <= 0.80):
            continue

        records.append(
            {
                "name": name,
                "abbr": name,
                "mw": round(mw, 6),
                "fluorine_atoms": int(f_atoms),
                "category": "PFASSTRUCTV6",
                "source": SOURCE,
                "formula": formula,
                "casrn": casrn,
                "smiles": smiles,
                "qc_level": qc_level,
                "fluorine_fraction": round(f_fraction, 8),
                "fluorine_percent": round(100.0 * f_fraction, 4),
            }
        )

    # Deduplicate by name + formula + MW
    seen = set()
    deduped = []

    for record in records:
        key = (
            record["name"],
            record["formula"],
            round(record["mw"], 4),
        )

        if key in seen:
            continue

        seen.add(key)
        deduped.append(record)

    deduped.sort(key=lambda r: (r["mw"], r["name"]))

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(deduped, f, indent=2)

    print(f"Input rows      : {len(df):,}")
    print(f"Valid PFAS rows : {len(records):,}")
    print(f"Deduplicated    : {len(deduped):,}")
    print(f"Output          : {OUTPUT_JSON}")

    if deduped:
        f_values = [r["fluorine_fraction"] for r in deduped]
        print()
        print(f"F fraction min  : {min(f_values):.4f}")
        print(f"F fraction max  : {max(f_values):.4f}")
        print(f"F fraction mean : {sum(f_values) / len(f_values):.4f}")

    environmental = [record for record in deduped if looks_environmental(record)]

    with ENV_OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(environmental, f, indent=2)

    print(f"Environmental  : {len(environmental):,}")
    print(f"Env output     : {ENV_OUTPUT_JSON}")

    if environmental:
        f_env = [r["fluorine_fraction"] for r in environmental]
        print()
        print("Environmental subset:")
        print(f"F fraction min  : {min(f_env):.4f}")
        print(f"F fraction max  : {max(f_env):.4f}")
        print(f"F fraction mean : {sum(f_env) / len(f_env):.4f}")


if __name__ == "__main__":
    main()