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

@dataclass(frozen=True)
class SubsetMonteCarloResult:
    k_values: list[int]
    summaries: list[dict]


def run_subset_size_experiment(
    compounds: list[PFASCompound],
    *,
    k_values: list[int],
    n_iterations_per_k: int = 100_000,
    alpha: float = 1.0,
    seed: int | None = 42,
) -> list[dict]:
    """
    For each k, repeatedly select k compounds from the full library,
    draw one Dirichlet mixture over those k compounds, and compute the
    weighted fluorine fraction.

    This estimates uncertainty as a function of effective mixture complexity.
    """

    if not compounds:
        raise ValueError("Compound list cannot be empty.")

    rng = np.random.default_rng(seed)

    f_library = np.array(
        [compound.fluorine_fraction for compound in compounds],
        dtype=float,
    )

    n_library = len(f_library)
    rows = []

    for k in k_values:
        if k <= 0:
            raise ValueError("k values must be positive.")

        if k > n_library:
            raise ValueError(f"k={k} exceeds library size {n_library}.")

        fluorine_fractions = np.empty(n_iterations_per_k, dtype=float)

        for i in range(n_iterations_per_k):
            idx = rng.choice(n_library, size=k, replace=False)
            f_subset = f_library[idx]

            weights = rng.dirichlet(np.full(k, alpha, dtype=float))
            fluorine_fractions[i] = float(weights @ f_subset)

        conversion_factors = 1.0 / fluorine_fractions

        rows.append(
            {
                "k": k,
                "n_iterations": n_iterations_per_k,
                "alpha": alpha,
                "mean_f": float(np.mean(fluorine_fractions)),
                "std_f": float(np.std(fluorine_fractions)),
                "p05_f": float(np.quantile(fluorine_fractions, 0.05)),
                "p50_f": float(np.quantile(fluorine_fractions, 0.50)),
                "p95_f": float(np.quantile(fluorine_fractions, 0.95)),
                "mean_cf": float(np.mean(conversion_factors)),
                "p05_cf": float(np.quantile(conversion_factors, 0.05)),
                "p50_cf": float(np.quantile(conversion_factors, 0.50)),
                "p95_cf": float(np.quantile(conversion_factors, 0.95)),
            }
        )

    return rows
