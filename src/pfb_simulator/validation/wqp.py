from __future__ import annotations

from pathlib import Path

import pandas as pd

from pfb_simulator.validation.pdh import PDHSample

WQP_NAME_MAP = {
    # PFCAs
    "Perfluorobutanoic acid": "PFBA",
    "Perfluorobutanoate": "PFBA",

    "Perfluorovaleric acid": "PFPeA",
    "Perfluoropentanoic acid": "PFPeA",
    "Perfluoropentanoate": "PFPeA",

    "Perfluorohexanoic acid": "PFHxA",
    "Perfluorohexanoate": "PFHxA",

    "Perfluoroheptanoic acid": "PFHpA",
    "Perfluoroheptanoate": "PFHpA",

    "Perfluorooctanoic acid": "PFOA",
    "Perfluorooctanoate": "PFOA",
    "PFOA ion": "PFOA",

    "Perfluorononanoic acid": "PFNA",
    "Perfluorononanoate": "PFNA",

    "Perfluorodecanoic acid": "PFDA",
    "Perfluorodecanoate": "PFDA",

    "Perfluoroundecanoic acid": "PFUnDA",
    "Perfluoroundecanoate": "PFUnDA",

    "Perfluorododecanoic acid": "PFDoDA",
    "Perfluorododecanoate": "PFDoDA",

    "Perfluorotridecanoate": "PFTrDA",
    "Perfluorotetradecanoate": "PFTeDA",

    # PFSAs
    "Perfluorobutanesulfonic acid": "PFBS",
    "Perfluorobutanesulfonate": "PFBS",

    "Perfluoropentanesulfonate": "PFPeS",

    "Perfluorohexanesulfonic acid": "PFHxS",
    "Perfluorohexanesulfonate": "PFHxS",

    "Perfluoroheptanesulfonate": "PFHpS",

    "Perfluorooctane sulfonic acid": "PFOS",
    "Perfluorooctanesulfonate": "PFOS",

    "Perfluorodecanesulfonate": "PFDS",

    # FTS
    "Fluorotelomer sulfonate 8:2": "8:2 FTS",

    "Perfluorooctanesulfonate (PFOS)": "PFOS",
    "Perfluorooctanoate (anionic form)": "PFOA",
    "Perfluoropentanesulfonic acid": "PFPeS",
    "Perfluorononanesulfonate": "PFNS",
    "Perfluorododecanesulfonate": "PFDoDS",
    "6:2 Fluorotelomer sulfonate acid": "6:2 FTS",
    "FtS 6:2 ion": "6:2 FTS",
    "Heptafluorobutyric acid***retired***use Perfluorobutanoic acid": "PFBA",
}

def _first_existing(df: pd.DataFrame, candidates: list[str]) -> str:
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(f"None of these columns found: {candidates}")


def load_wqp_samples(path: str | Path) -> list[PDHSample]:
    path = Path(path)
    df = pd.read_excel(path)

    if "Environmental Media Name" in df.columns:
        df = df[
            df["Environmental Media Name"]
            .astype(str)
            .str.lower()
            .eq("water")
        ].copy()

    sample_col = _first_existing(
        df,
        [
            "Activity Identifier",
            "ActivityIdentifier",
            "Monitoring Location Identifier",
            "MonitoringLocationIdentifier",
            "Sample Result",
        ],
    )

    compound_col = _first_existing(
        df,
        [
            "PFAS Chemical Name",
            "CharacteristicName",
            "Characteristic Name",
            "Contaminant",
            "Substance",
        ],
    )

    value_col = _first_existing(
        df,
        [
            "Result Measure Value (ppt)",
            "Result Measure Value",
            "ResultMeasureValue",
            "Concentration (ng/L)",
            "ResultValue",
        ],
    )

    unit_col = None
    for candidate in [
        "Result Unit of Measure",
        "ResultMeasure/MeasureUnitCode",
        "Result Unit",
        "Unit",
    ]:
        if candidate in df.columns:
            unit_col = candidate
            break
            
    lat_col = next((c for c in ["LatitudeMeasure", "Latitude"] if c in df.columns), None)
    lon_col = next((c for c in ["LongitudeMeasure", "Longitude"] if c in df.columns), None)
    date_col = next((c for c in ["ActivityStartDate", "Date"] if c in df.columns), None)

    samples: list[PDHSample] = []

    for sample_id, group in df.groupby(sample_col):
        compounds: dict[str, float] = {}

        latitude = None
        longitude = None
        date = None

        if lat_col:
            latitude = pd.to_numeric(group[lat_col], errors="coerce").dropna()
            latitude = float(latitude.iloc[0]) if len(latitude) else None

        if lon_col:
            longitude = pd.to_numeric(group[lon_col], errors="coerce").dropna()
            longitude = float(longitude.iloc[0]) if len(longitude) else None

        if date_col:
            dates = group[date_col].dropna()
            date = str(dates.iloc[0]) if len(dates) else None

        for _, row in group.iterrows():
            compound = str(row[compound_col]).strip()

            # Skip TOP assay derived values for now
            if "plus total oxidizable precursors" in compound.lower():
                continue

            compound = WQP_NAME_MAP.get(compound, compound)

            value = pd.to_numeric(row[value_col], errors="coerce")

            if pd.isna(value) or value <= 0:
                continue

            # If using the WQP normalized ppt column, treat it directly as ng/L.
            if value_col == "Result Measure Value (ppt)":
                value = float(value)
            else:
                unit = str(row[unit_col]).strip().lower() if unit_col else "ng/l"

                if unit in ["ug/l", "µg/l", "micrograms per liter", "microgram per liter"]:
                    value = float(value) * 1000.0
                elif unit in ["ng/l", "ppt"]:
                    value = float(value)
                else:
                    continue

            compounds[compound] = value

        samples.append(
            PDHSample(
                sample_id=str(sample_id),
                date=date,
                latitude=latitude,
                longitude=longitude,
                compounds=compounds,
                raw={"source": "WQP", "rows": len(group)},
            )
        )

    # Report unmapped compound names
    unknown = set()

    for sample in samples:
        for compound in sample.compounds:
            if compound not in {
                "PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA",
                "PFUnDA","PFDoDA","PFTrDA","PFTeDA",
                "PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFDS",
                "6:2 FTS","8:2 FTS","HFPO-DA","ADONA"
            }:
                unknown.add(compound)

    if unknown:
        print("\nUnknown PFAS names:")
        for name in sorted(unknown):
            print("  ", name)

    return samples
