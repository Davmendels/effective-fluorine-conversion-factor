from __future__ import annotations

from pathlib import Path

import pandas as pd

from pfb_simulator.validation.pdh import PDHSample


def load_ucmr5_samples(path: str | Path) -> list[PDHSample]:
    path = Path(path)

    df = pd.read_excel(path)

    required = {
        "Sample Result",
        "Contaminant",
        "Concentration (ng/L)",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    samples: list[PDHSample] = []

    for sample_id, group in df.groupby("Sample Result"):
        compounds: dict[str, float] = {}

        for _, row in group.iterrows():
            contaminant = str(row["Contaminant"]).strip()
            value = row["Concentration (ng/L)"]

            if pd.isna(value):
                continue

            value = float(value)

            if value > 0:
                compounds[contaminant] = value

        samples.append(
            PDHSample(
                sample_id=str(sample_id),
                compounds=compounds,
                raw={"source": "UCMR5", "rows": len(group)},
            )
        )

    return samples
