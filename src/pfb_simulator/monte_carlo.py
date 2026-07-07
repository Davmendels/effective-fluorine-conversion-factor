from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pfb_simulator.compounds import PFASCompound


@dataclass(frozen=True)
class MonteCarloResult:
    fluorine_fractions: np.ndarray
    conversion_factors: np.ndarray

    @property
    def n_iterations(self) -> int:
        return len(self.fluorine_fractions)


def run_dirichlet_monte_carlo(
    compounds: list[PFASCompound],
    *,
    n_iterations: int = 100_000,
    alpha: float = 1.0,
    seed: int | None = 42,
    batch_size: int = 100_000,
) -> MonteCarloResult:
    """
    Run a Monte Carlo simulation over random PFAS mixtures.

    Each mixture is drawn from a symmetric Dirichlet distribution.

    alpha = 1.0 gives uniform sampling over the simplex.
    alpha < 1.0 gives sparse mixtures dominated by a few compounds.
    alpha > 1.0 gives more evenly mixed compositions.
    """

    if not compounds:
        raise ValueError("Compound list cannot be empty.")

    if n_iterations <= 0:
        raise ValueError("n_iterations must be positive.")

    if alpha <= 0:
        raise ValueError("alpha must be positive.")

    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")

    rng = np.random.default_rng(seed)

    fluorine_vector = np.array(
        [compound.fluorine_fraction for compound in compounds],
        dtype=float,
    )

    n_compounds = len(compounds)

    fluorine_fractions = np.empty(n_iterations, dtype=float)

    start = 0

    while start < n_iterations:
        stop = min(start + batch_size, n_iterations)
        size = stop - start

        weights = rng.dirichlet(
            alpha=np.full(n_compounds, alpha, dtype=float),
            size=size,
        )

        fluorine_fractions[start:stop] = weights @ fluorine_vector

        start = stop

    conversion_factors = 1.0 / fluorine_fractions

    return MonteCarloResult(
        fluorine_fractions=fluorine_fractions,
        conversion_factors=conversion_factors,
    )
