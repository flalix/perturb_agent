"""Shared loader. Input: samples x genes table of log fold change of CPM."""
import numpy as np
import pandas as pd

def load_lfc(path=None, n_samples=80, n_genes=300, seed=0):
    """Load samples x genes lfc(CPM) table from CSV, or build a synthetic demo."""
    if path:
        return pd.read_csv(path, index_col=0)
    rng = np.random.default_rng(seed)
    cols, data = [], []
    for g in range(n_genes):
        kind = rng.choice(["uni", "bi", "tri"], p=[0.6, 0.3, 0.1])
        if kind == "uni":
            x = rng.normal(0, 0.4, n_samples)
        elif kind == "bi":
            shift = rng.choice([-1, 1]) * rng.uniform(2.5, 4)
            comp = rng.random(n_samples) < 0.35
            x = np.where(comp, rng.normal(shift, 0.4, n_samples), rng.normal(0, 0.4, n_samples))
        else:
            comp = rng.choice([-1, 0, 1], n_samples, p=[0.25, 0.5, 0.25])
            x = rng.normal(comp * rng.uniform(3, 4), 0.4, n_samples)
        cols.append(f"GENE_{g:04d}_{kind}")
        data.append(x)
    return pd.DataFrame(np.array(data).T, columns=cols,
                        index=[f"S{i:03d}" for i in range(n_samples)])
