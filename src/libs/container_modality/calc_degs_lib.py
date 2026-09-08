"""Part 2: per-gene stats + KMeans k=2..5, silhouette-based modality, per-mode stats and samples.

Outputs
  gene_modes.csv     long table: one row per (gene, mode) with mean/median/std/sem/samples
  gene_summary.csv   one row per gene: k, silhouette, overall stats, n samples in up/down modes
  offset_matrix.csv  samples x genes integer offset of each sample's mode relative to central
                     (0 = central, <0 = down, >0 = up); consumed by part 3
"""
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from common import load_lfc

def basic_stats(x):
    x = np.asarray(x, float)
    return dict(mean=x.mean(), median=np.median(x), std=x.std(ddof=1) if len(x) > 1 else 0.0,
                sem=(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0, n=len(x))

def best_kmeans(x, ks=(2, 3, 4, 5), sil_threshold=0.65, seed=0):
    """Return (k, silhouette, labels). k=1 if no k reaches sil_threshold."""
    X = x.reshape(-1, 1)
    best_k, best_s, best_lab = 1, np.nan, np.zeros(len(x), int)
    for k in ks:
        if k >= len(x):
            continue
        lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(X)
        if len(np.unique(lab)) < 2:
            continue
        s = silhouette_score(X, lab)
        if np.isnan(best_s) or s > best_s:
            best_k, best_s, best_lab = k, s, lab
    if best_k > 1 and best_s < sil_threshold:
        best_k, best_lab = 1, np.zeros(len(x), int)
    return best_k, best_s, best_lab

def analyse_gene(name, series, central_rule="zero", **kw):
    """central_rule: 'zero' -> mode whose mean is closest to 0; 'largest' -> most populated mode."""
    s = series.dropna()
    x, samples = s.to_numpy(float), s.index.to_numpy()
    k, sil, lab = best_kmeans(x, **kw)
    modes = []
    for m in np.unique(lab):
        st = basic_stats(x[lab == m])
        st["label"], st["samples"] = m, samples[lab == m]
        modes.append(st)
    modes.sort(key=lambda d: d["mean"])                       # left -> right
    if central_rule == "largest":
        c_idx = int(np.argmax([d["n"] for d in modes]))
    else:
        c_idx = int(np.argmin([abs(d["mean"]) for d in modes]))
    rows, offsets = [], pd.Series(0, index=samples, dtype=int)
    for i, d in enumerate(modes):
        off = i - c_idx
        mtype = "central" if off == 0 else ("down" if off < 0 else "up")
        offsets[d["samples"]] = off
        rows.append(dict(gene=name, k=k, silhouette=sil, mode_rank=i, offset=off, mode_type=mtype,
                         mean=d["mean"], median=d["median"], std=d["std"], sem=d["sem"],
                         n_samples=d["n"], samples=";".join(map(str, d["samples"]))))
    return rows, offsets

def run(lfc, **kw):
    mode_rows, offsets, summ = [], {}, []
    for g in lfc.columns:
        rows, off = analyse_gene(g, lfc[g], **kw)
        mode_rows += rows
        offsets[g] = off
        ov = basic_stats(lfc[g].dropna())
        summ.append(dict(gene=g, k=rows[0]["k"], silhouette=rows[0]["silhouette"],
                         mean=ov["mean"], median=ov["median"], std=ov["std"], sem=ov["sem"],
                         n_up_modes=sum(r["mode_type"] == "up" for r in rows),
                         n_down_modes=sum(r["mode_type"] == "down" for r in rows),
                         n_up_samples=int((off > 0).sum()), n_down_samples=int((off < 0).sum())))
    gene_modes = pd.DataFrame(mode_rows)
    gene_summary = pd.DataFrame(summ).set_index("gene")
    offset_matrix = pd.DataFrame(offsets).reindex(lfc.index)
    return gene_modes, gene_summary, offset_matrix

if __name__ == "__main__":
    lfc = load_lfc()
    gene_modes, gene_summary, offset_matrix = run(lfc, sil_threshold=0.65, central_rule="zero")
    gene_modes.to_csv("gene_modes.csv", index=False)
    gene_summary.to_csv("gene_summary.csv")
    offset_matrix.to_csv("offset_matrix.csv")
    print(gene_summary["k"].value_counts().sort_index().rename("genes per k"))
    print(gene_modes.head(8).drop(columns="samples"))
