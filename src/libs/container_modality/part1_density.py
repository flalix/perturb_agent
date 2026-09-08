"""Part 1: 5 random genes -> density histogram + uni/multimodal call (KDE peak counting)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from common import load_lfc

def count_modes(x, grid=512, min_prominence=0.05):
    """Number of KDE peaks with normalized prominence >= min_prominence."""
    x = x[~np.isnan(x)]
    kde = gaussian_kde(x)
    xs = np.linspace(x.min(), x.max(), grid)
    dens = kde(xs)
    dens /= dens.max()
    peaks, _ = find_peaks(dens, prominence=min_prominence)
    return len(peaks), xs, dens, peaks

def plot_random_genes(lfc, n=5, seed=0, out="part1_density.png"):
    rng = np.random.default_rng(seed)
    genes = rng.choice(lfc.columns, size=n, replace=False)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 3.5))
    calls = {}
    for ax, g in zip(axes, genes):
        x = lfc[g].to_numpy()
        n_modes, xs, dens, peaks = count_modes(x)
        calls[g] = "unimodal" if n_modes == 1 else f"multimodal ({n_modes} modes)"
        ax.hist(x, bins=30, density=True, alpha=0.5)
        ax.plot(xs, dens * ax.get_ylim()[1], lw=2)      # KDE rescaled to hist height
        ax.plot(xs[peaks], dens[peaks] * ax.get_ylim()[1], "rv")
        ax.set_title(f"{g}\n{calls[g]}", fontsize=9)
        ax.set_xlabel("lfc(CPM)")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    return calls

if __name__ == "__main__":
    lfc = load_lfc()          # load_lfc("cluster_lfc.csv") for real data
    for g, c in plot_random_genes(lfc).items():
        print(f"{g}: {c}")
