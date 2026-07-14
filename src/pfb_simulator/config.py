from __future__ import annotations


# Stable run-directory name:
#     (display name, country, environmental compartment)
DATASET_REGISTRY: dict[str, tuple[str, str, str]] = {
    "Flanders": (
        "Flanders",
        "Belgium",
        "mixed_water",
    ),
    "France_ADES": (
        "France_ADES",
        "France",
        "groundwater",
    ),
    "France_ANSES": (
        "France_ANSES",
        "France",
        "drinking_water",
    ),
    "France_EauRob": (
        "France_EauRob",
        "France",
        "tap_water",
    ),
    "France_ICPE": (
        "France_ICPE",
        "France",
        "industrial_discharge",
    ),
    "France_Naiades": (
        "France_Naiades",
        "France",
        "surface_water",
    ),
    "Grand_Lyon": (
        "Grand_Lyon",
        "France",
        "soil_extract",
    ),
    "Germany_ELWAS": (
        "Germany_ELWAS",
        "Germany",
        "mixed_water",
    ),
    "Germany_TFA": (
        "Germany_TFA",
        "Germany",
        "surface_water_with_TFA",
    ),
    "Muir": (
        "Muir",
        "Global",
        "sea_water",
    ),
    "Italy_Friuli_Venezia": (
        "Italy_Friuli_Venezia",
        "Italy",
        "mixed_water",
    ),
    "Veneto": (
        "Veneto",
        "Italy",
        "mixed_water",
    ),
    "RIVM": (
        "RIVM",
        "Netherlands",
        "mixed_water",
    ),
    "UCMR5_extract": (
        "UCMR5_extract",
        "USA",
        "drinking_water",
    ),
    "WQP_multimedia": (
        "WQP_multimedia",
        "USA",
        "environmental_hotspot_water",
    ),
    "UK_EA": (
        "UK_EA",
        "United Kingdom",
        "mixed_water",
    ),
}


# Order used in publication figures and dataset-level comparisons.
DATASET_ORDER: list[str] = [
    "Grand_Lyon",
    "Muir",
    "Germany_ELWAS",
    "UK_EA",
    "RIVM",
    "France_Naiades",
    "France_ICPE",
    "Flanders",
    "France_ADES",
    "France_ANSES",
    "France_EauRob",
    "Germany_TFA",
    "Italy_Friuli_Venezia",
    "Veneto",
    "UCMR5_extract",
    "WQP_multimedia",
]


# Display labels for environmental compartments.
COMPARTMENT_LABELS: dict[str, str] = {
    "soil_extract": "Soil extracts",
    "sea_water": "Marine water",
    "surface_water": "Surface water",
    "industrial_discharge": "Industrial discharges",
    "groundwater": "Groundwater",
    "mixed_water": "Mixed water",
    "surface_water_with_TFA": "Surface water (TFA)",
    "tap_water": "Tap water",
    "drinking_water": "Drinking water",
    "environmental_hotspot_water": "Highly contaminated sites",
}


# Muted publication palette shared across figures.
COMPARTMENT_COLORS: dict[str, str] = {
    "soil_extract": "#b79ac8",
    "sea_water": "#7298c9",
    "surface_water": "#89a9cf",
    "industrial_discharge": "#e89a58",
    "groundwater": "#a8bdd8",
    "mixed_water": "#b8c9dd",
    "surface_water_with_TFA": "#86a8d5",
    "tap_water": "#9fc69b",
    "drinking_water": "#83b980",
    "environmental_hotspot_water": "#dc8581",
}


def dataset_metadata(dataset_name: str) -> tuple[str, str, str]:
    """
    Return display name, country, and compartment for one dataset.
    """
    try:
        return DATASET_REGISTRY[dataset_name]
    except KeyError as error:
        raise KeyError(
            f"Unknown dataset name: {dataset_name!r}"
        ) from error


def dataset_compartment(dataset_name: str) -> str:
    """
    Return the environmental compartment associated with one dataset.
    """
    return dataset_metadata(dataset_name)[2]
