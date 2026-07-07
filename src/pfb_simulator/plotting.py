from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_histogram(
    values: np.ndarray,
    output_path: str | Path,
    *,
    title: str,
    xlabel: str,
    bins: int = 80,
) -> None:
    values = np.asarray(values, dtype=float)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.hist(values, bins=bins)
    ax.axvline(np.mean(values), linestyle="--", linewidth=1.5, label="Mean")
    ax.axvline(np.quantile(values, 0.05), linestyle=":", linewidth=1.5, label="5–95% interval")
    ax.axvline(np.quantile(values, 0.95), linestyle=":", linewidth=1.5)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

