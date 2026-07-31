"""
perturbation_axis.py
====================
Extends `tri_evidence.py` with the PERTURBATIONAL TRANSCRIPTOMIC axis.

The distinction that matters
----------------------------
DepMap measures a *requirement* (viability under knockout). One number per gene per line.
CMap and Tahoe measure a *response* (the whole transcriptome after perturbation).

    requirement  -> "is this gene needed to survive?"        [DepMap, Chronos]
    response     -> "what does the cell DO when you hit it?" [CMap L1000, Tahoe-100M]

A gene can be dispensable for viability yet be the master regulator of a cell state.
DepMap calls that a negative. The response axis calls it a state regulator. That
disambiguation is the single biggest thing these two resources buy you.

    CMap / LINCS 2020    ~1.2M replicate-collapsed signatures, 33.6k compounds AND
                         ~9.3k genes perturbed by shRNA / CRISPR / overexpression,
                         240 cell contexts. Bulk, 978 landmarks + inference, 6h/24h.
                         -> breadth, and the only large GENETIC perturbation compendium.

    Tahoe-100M           100M single cells, 1,100 compounds x 50 cancer lines, ~60k
                         conditions, 24h, full 3' transcriptome, all 50 lines pooled
                         in the same well (Mosaic).
                         -> depth, single-cell response heterogeneity, internally
                            controlled cross-context comparison, lncRNA visibility.

They are complementary, not redundant. CMap has genetic perturbagens and 30x the
compounds; Tahoe has single-cell resolution, a full transcriptome, and a shared-well
design that removes the batch confound from cross-cell-line comparison.

Join key across all five resources: the cell line. DepMap ModelID <-> CCLE name <->
Cellosaurus <-> CMap cell_iname <-> Tahoe cell_name. Harmonise this first or nothing
downstream is valid.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

try:
    from tri_evidence import DepMap, TriEvidence, _cached
except ImportError:  # allow standalone use
    DepMap = TriEvidence = None


# --------------------------------------------------------------------------- #
# Connectivity scoring (CMap WTCS -> NCS)
# --------------------------------------------------------------------------- #
def _weighted_ks(ranked_genes: Sequence[str], gene_set: set[str], weights: np.ndarray | None = None) -> float:
    """GSEA-style weighted enrichment score of `gene_set` in a rank-ordered list."""
    n = len(ranked_genes)
    hits = np.array([g in gene_set for g in ranked_genes], dtype=bool)
    n_hit = hits.sum()
    if n_hit == 0 or n_hit == n:
        return 0.0
    w = np.abs(weights) if weights is not None else np.ones(n)
    p_hit = np.cumsum(np.where(hits, w, 0.0)) / max(np.sum(w[hits]), 1e-12)
    p_miss = np.cumsum(np.where(~hits, 1.0, 0.0)) / (n - n_hit)
    dev = p_hit - p_miss
    return float(dev[np.argmax(np.abs(dev))])


def wtcs(query_up: Iterable[str], query_dn: Iterable[str], signature: pd.Series) -> float:
    """Weighted Connectivity Score of a query gene-set pair against one signature.

    `signature` is a z-score / logFC vector indexed by gene symbol (CMap Level 5 MODZ).
    Positive WTCS = the perturbation RECAPITULATES the query state.
    Negative WTCS = the perturbation REVERSES it.

    Reversal is the drug-repurposing direction. Recapitulation, when the perturbagen is a
    shRNA/CRISPR/OE reagent, is the *upstream regulator inference* direction -- and that
    is the more scientifically interesting one, because it is an experimental answer to
    "what drives my module?" rather than a network-inference guess.
    """
    s = signature.dropna().sort_values(ascending=False)
    ranked, w = s.index.tolist(), s.values
    up, dn = set(query_up), set(query_dn)
    es_up = _weighted_ks(ranked, up, w)
    es_dn = _weighted_ks(ranked, dn, w)
    return 0.0 if np.sign(es_up) == np.sign(es_dn) else (es_up - es_dn) / 2.0


def normalize_cs(scores: pd.Series, groups: pd.Series) -> pd.Series:
    """NCS: divide each WTCS by the mean of same-signed scores within its group
    (cell line x perturbagen type). Required before comparing across contexts."""
    out = scores.copy().astype(float)
    for g, idx in scores.groupby(groups).groups.items():
        v = scores.loc[idx]
        for sign in (1, -1):
            m = np.sign(v) == sign
            if m.sum():
                mu = np.abs(v[m]).mean()
                out.loc[v[m].index] = v[m] / (mu if mu else 1.0)
    return out


# --------------------------------------------------------------------------- #
# CMap / LINCS L1000
# --------------------------------------------------------------------------- #
@dataclass
class CMap:
    """Local Level 5 (MODZ) signature access.

    Recommended acquisition path, in order of preference:
      1. `cmapBQ` -> targeted retrieval from Google BigQuery. Correct choice for a
         pipeline: you never download the 40+ GB GCTx.
      2. CLUE data library GCTx + `cmapPy.pandasGEXpress.parse` (needs a free API key).
      3. GEO: GSE92742 (Phase I), GSE70138 (Phase II), GSE106127 (shRNA/CRISPR subset).

    Quality control that is NOT optional:
      * Filter to exemplar signatures and reasonable TAS (transcriptional activity score).
        A large fraction of L1000 signatures are inert; including them destroys your
        null distribution.
      * Restrict claims to the 978 LANDMARK genes. Of ~11,350 inferred transcripts only
        ~9,196 are 'well inferred', and inference error is structured, not random.
      * `tau` is a percentile against Touchstone, which covers only ~9 core cell lines.
        Outside those contexts, use NCS and say so.
    """
    sig_matrix: pd.DataFrame          # genes x signature_id, Level 5 MODZ z-scores
    sig_meta: pd.DataFrame            # signature_id -> pert_id, pert_iname, pert_type, cell_iname, dose, time

    LANDMARK_ONLY: bool = True

    @classmethod
    def from_parquet(cls, mat_path: str | Path, meta_path: str | Path) -> "CMap":
        return cls(pd.read_parquet(mat_path), pd.read_parquet(meta_path).set_index("sig_id"))

    def _subset(self, pert_types: Sequence[str] | None, cell_lines: Sequence[str] | None) -> pd.Index:
        m = self.sig_meta
        keep = pd.Series(True, index=m.index)
        if pert_types:
            keep &= m["pert_type"].isin(pert_types)
        if cell_lines:
            keep &= m["cell_iname"].isin(cell_lines)
        return m.index[keep].intersection(self.sig_matrix.columns)

    def connectivity(
        self,
        query_up: Sequence[str],
        query_dn: Sequence[str],
        pert_types: Sequence[str] | None = None,
        cell_lines: Sequence[str] | None = None,
        top_n: int = 50,
    ) -> pd.DataFrame:
        """Score a query signature against the compendium.

        pert_types of interest:
          'trt_cp'   compound
          'trt_sh'   shRNA knockdown        -> upstream regulator inference
          'trt_xpr'  CRISPR knockout        -> upstream regulator inference
          'trt_oe'   overexpression         -> sufficiency, not just necessity
        """
        sids = self._subset(pert_types, cell_lines)
        scores = pd.Series(
            {sid: wtcs(query_up, query_dn, self.sig_matrix[sid]) for sid in sids}, name="wtcs"
        )
        meta = self.sig_meta.loc[scores.index]
        grp = meta["cell_iname"].astype(str) + "|" + meta["pert_type"].astype(str)
        out = meta.assign(wtcs=scores, ncs=normalize_cs(scores, grp))
        out = out.sort_values("ncs")
        return pd.concat([out.head(top_n), out.tail(top_n)])

    def upstream_regulators(self, query_up, query_dn, top_n: int = 30) -> pd.DataFrame:
        """Which GENETIC perturbations reproduce my transcriptional module?

        This is the operation people forget CMap can do. It replaces TF-network inference
        (DoRothEA/VIPER/IPA) with an experimental answer: these are genes whose actual
        knockdown/knockout/overexpression produced your state in a real cell. Far stronger
        evidence than a curated regulon, and it nominates causally testable drivers.
        """
        df = self.connectivity(query_up, query_dn, pert_types=["trt_sh", "trt_xpr", "trt_oe"], top_n=200)
        agg = (
            df.groupby(["pert_iname", "pert_type"])["ncs"]
            .agg(["mean", "median", "count"])
            .query("count >= 2")
            .sort_values("median", ascending=False)
        )
        return agg.head(top_n)

    def reversal_compounds(self, query_up, query_dn, top_n: int = 30) -> pd.DataFrame:
        """Compounds that push the transcriptome AWAY from the query state."""
        df = self.connectivity(query_up, query_dn, pert_types=["trt_cp"], top_n=300)
        agg = (
            df.groupby(["pert_iname"])["ncs"]
            .agg(["mean", "median", "count"])
            .query("count >= 3")            # reproducible across cell lines/doses
            .sort_values("median")
        )
        return agg.head(top_n)


# --------------------------------------------------------------------------- #
# Tahoe-100M
# --------------------------------------------------------------------------- #
@dataclass
class Tahoe:
    """Tahoe-100M access.

    Do NOT start from 100M raw cells. The right first move for a systems-biology pipeline
    is to pseudobulk to (cell_line x drug x dose), which collapses the atlas to ~60k
    profiles that behave like a very deep, internally-controlled L1000 -- tractable on a
    workstation, and directly comparable to CMap signatures.

    Then go back to single-cell ONLY for the questions that require it: response
    heterogeneity, subpopulation selection, and program-level bimodality.

    Sources:
      datasets:  tahoebio/Tahoe-100M                (HF, parquet, streamable)
      embeddings: tahoebio/Tahoe-x1-embeddings      (precomputed 3B-model cell embeddings)
      model:     tahoebio/Tahoe-x1                  (70m / 1b / 3b, Apache-2.0)
      scvi:      tahoebio/Tahoe-100M-SCVI-v1

    Design advantages worth exploiting explicitly:
      * All 50 lines are perturbed in the SAME WELL. Drug concentration, timing and batch
        are shared, so cross-cell-line response differences are real biology rather than
        plate effects. This is the cleanest available substrate for asking "which genomic
        context determines response?"
      * Full 3' transcriptome, so lncRNAs are visible (subject to dropout) -- unlike L1000,
        whose 978-landmark design excludes essentially all of them.
    """
    pseudobulk: pd.DataFrame          # (cell_line, drug, dose) x genes, log-CPM
    meta: pd.DataFrame                # matching index metadata incl. control flags
    counts_h5ad: str | Path | None = None   # optional path for single-cell drill-down

    @staticmethod
    def build_pseudobulk(h5ad_or_parquet, out_path: str | Path, min_cells: int = 30) -> Path:
        """Collapse Tahoe to (cell_line, drug, dose) pseudobulk. Run once, cache forever.

        Sum raw counts within group -> CPM -> log1p. Summing counts (not averaging
        normalised values) is the statistically correct collapse for scRNA-seq.
        """
        import anndata as ad
        import scanpy as sc

        adata = ad.read_h5ad(h5ad_or_parquet) if str(h5ad_or_parquet).endswith(".h5ad") else None
        if adata is None:
            raise NotImplementedError("Point this at an .h5ad shard, or stream parquet in chunks.")
        keys = ["cell_name", "drug", "dose"]
        grp = adata.obs.groupby(keys, observed=True).indices
        rows, index = [], []
        for k, idx in grp.items():
            if len(idx) < min_cells:
                continue
            rows.append(np.asarray(adata[idx].X.sum(axis=0)).ravel())
            index.append(k)
        M = pd.DataFrame(rows, index=pd.MultiIndex.from_tuples(index, names=keys), columns=adata.var_names)
        M = np.log1p(M.div(M.sum(axis=1), axis=0) * 1e6)
        M.to_parquet(out_path)
        return Path(out_path)

    def signature(self, cell_line: str, drug: str, dose: str | None = None,
                  control: str = "DMSO_TF") -> pd.Series:
        """logFC vs vehicle within the SAME cell line. Comparable to a CMap Level 5 column."""
        pb = self.pseudobulk
        try:
            ctrl = pb.loc[(cell_line, control)].mean(axis=0)
        except KeyError:
            raise KeyError(f"no vehicle control for {cell_line!r}")
        trt = pb.loc[(cell_line, drug)] if dose is None else pb.loc[(cell_line, drug, dose)]
        trt = trt.mean(axis=0) if isinstance(trt, pd.DataFrame) else trt
        return (trt - ctrl).rename(f"{cell_line}|{drug}|{dose or 'agg'}")

    def context_dependence(self, drug: str, gene_set: Sequence[str]) -> pd.DataFrame:
        """How does response to one drug vary across the 50 lines, for a given program?

        Because all lines share the well, this is an unusually clean readout of
        context-dependent drug response. Join the output to DepMap mutation/CN/expression
        to ask which genomic feature predicts responsiveness.
        """
        recs = []
        for cl in self.pseudobulk.index.get_level_values(0).unique():
            try:
                s = self.signature(cl, drug)
            except KeyError:
                continue
            g = [x for x in gene_set if x in s.index]
            if len(g) < 5:
                continue
            recs.append({"cell_line": cl, "n_genes": len(g),
                         "program_logfc": float(s[g].mean()),
                         "global_logfc_sd": float(s.std())})
        return pd.DataFrame(recs).sort_values("program_logfc")

    def response_heterogeneity(self, cell_line: str, drug: str, gene_set: Sequence[str]) -> dict:
        """Single-cell drill-down: is the response uniform or does a subpopulation escape?

        Bimodality coefficient > ~0.555 suggests a mixture. A bimodal response with a
        clearly unresponsive mode is the transcriptional signature of pre-existing
        resistance -- invisible in bulk L1000 and invisible in DepMap, and a directly
        testable hypothesis about which subclone survives.
        """
        import anndata as ad
        if self.counts_h5ad is None:
            raise ValueError("set counts_h5ad for single-cell queries")
        a = ad.read_h5ad(self.counts_h5ad, backed="r")
        m = (a.obs["cell_name"] == cell_line) & (a.obs["drug"] == drug)
        sub = a[m].to_memory()
        g = [x for x in gene_set if x in sub.var_names]
        score = np.asarray(sub[:, g].X.mean(axis=1)).ravel()
        n = score.size
        from scipy import stats
        skew, kurt = stats.skew(score), stats.kurtosis(score, fisher=True)
        bc = (skew ** 2 + 1) / (kurt + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
        return {"cell_line": cell_line, "drug": drug, "n_cells": int(n),
                "mean_score": float(score.mean()), "bimodality_coefficient": float(bc),
                "interpretation": "likely_mixed_response" if bc > 0.555 else "uniform_response"}


# --------------------------------------------------------------------------- #
# Five-axis integration
# --------------------------------------------------------------------------- #
def state_vs_survival(gene: str, dm: "DepMap", cm: CMap, min_sigs: int = 3) -> dict:
    """THE disambiguation this axis exists for.

    Compares magnitude of transcriptional response to knockdown/knockout (CMap) against
    fitness cost of the same knockout (DepMap). Four quadrants, four different biologies:

      response HIGH, fitness cost HIGH -> essential hub
      response HIGH, fitness cost LOW  -> STATE REGULATOR (subtype identity, plasticity;
                                          e.g. the classic lineage-TF profile). DepMap
                                          alone would have discarded this gene.
      response LOW,  fitness cost HIGH -> structural/metabolic requirement, little
                                          transcriptional feedback
      response LOW,  fitness cost LOW  -> genuinely inert in this system
    """
    sel = dm.selectivity(gene)
    sids = cm.sig_meta.index[
        (cm.sig_meta["pert_iname"] == gene) & (cm.sig_meta["pert_type"].isin(["trt_sh", "trt_xpr"]))
    ].intersection(cm.sig_matrix.columns)
    if len(sids) < min_sigs:
        return {"gene": gene, "status": "insufficient_cmap_signatures", "n_sigs": len(sids)}
    mag = float(cm.sig_matrix[sids].abs().mean().mean())          # mean |z| across signatures
    n_strong = int((cm.sig_matrix[sids].abs() > 2).sum().mean())  # genes moved per signature
    fitness = sel.get("mean_chronos", 0.0)
    hi_resp, hi_fit = mag > 0.6, fitness < -0.4
    quadrant = {
        (True, True): "essential_hub",
        (True, False): "state_regulator",
        (False, True): "structural_requirement",
        (False, False): "inert",
    }[(hi_resp, hi_fit)]
    return {"gene": gene, "n_cmap_signatures": len(sids), "mean_abs_z": mag,
            "genes_moved_per_sig": n_strong, "mean_chronos": fitness, "quadrant": quadrant}


def perturbation_flags(ev: dict) -> list[str]:
    """Contradiction patterns unlocked by adding the response axis."""
    f: list[str] = []
    q = ev.get("state_vs_survival", {}).get("quadrant")
    if q == "state_regulator":
        f.append(
            "STATE_REGULATOR: large transcriptional response to knockdown, no fitness cost. "
            "DepMap-only screening would score this a negative. Prime candidate for a "
            "subtype-identity / plasticity driver -- test by asking whether its CMap "
            "signature matches your cluster-defining program."
        )
    if q == "structural_requirement":
        f.append("STRUCTURAL: fitness cost without transcriptional feedback; unlikely to be "
                 "a tractable state-modulating target")

    rev = ev.get("reversal", {})
    if rev.get("best_ncs", 0) < -1.5 and not rev.get("target_is_dependency", True):
        f.append(
            "CONNECTIVITY_WITHOUT_DEPENDENCY: a compound reverses the signature but its "
            "nominal target is not a DepMap dependency. Suspect polypharmacology or an "
            "off-target driver. Use Tahoe context-dependence to find which lines respond, "
            "then regress response on DepMap mutation/CN to nominate the real target."
        )

    het = ev.get("heterogeneity", {})
    if het.get("interpretation") == "likely_mixed_response":
        f.append(
            "HETEROGENEOUS_RESPONSE: bimodal single-cell response implies a pre-existing "
            "non-responding subpopulation. Bulk assays (L1000, DepMap) average this away. "
            "Hypothesis: intrinsic resistance -- characterise the escaping mode."
        )

    if ev.get("cmap_landmark_coverage", 1.0) < 0.15:
        f.append(
            "POOR_L1000_COVERAGE: your program is mostly outside the 978 landmark genes "
            "(typical for lncRNA-driven programs). CMap connectivity is unreliable here; "
            "use Tahoe, which measures the full transcriptome, as the primary response axis."
        )
    return f


def five_axis_dossier(
    gene: str,
    program_up: Sequence[str],
    program_dn: Sequence[str],
    dm: "DepMap",
    cm: CMap,
    th: Tahoe | None = None,
    landmarks: set[str] | None = None,
) -> dict:
    """Assemble requirement + response + genetic + phenotype evidence in one object.

    Compose with `tri_evidence.triangulate()` for the Open Targets / Monarch axes.
    """
    out: dict = {"gene": gene}
    out["state_vs_survival"] = state_vs_survival(gene, dm, cm)
    out["upstream_regulators"] = cm.upstream_regulators(program_up, program_dn).to_dict("index")
    rev = cm.reversal_compounds(program_up, program_dn)
    out["reversal"] = {"top": rev.head(10).to_dict("index"),
                       "best_ncs": float(rev["median"].min()) if len(rev) else 0.0}
    if landmarks:
        allg = set(program_up) | set(program_dn)
        out["cmap_landmark_coverage"] = len(allg & landmarks) / max(len(allg), 1)
    if th is not None and len(rev):
        top_drug = rev.index[0]
        out["tahoe_context"] = th.context_dependence(top_drug, list(program_up)).head(10).to_dict("records")
    out["flags"] = perturbation_flags(out)
    return out


# --------------------------------------------------------------------------- #
# Additional MCP tools
# --------------------------------------------------------------------------- #
def register_perturbation_tools(mcp, dm, cm: CMap, th: Tahoe | None = None):
    """Attach to the FastMCP server from tri_evidence.build_mcp_server()."""

    @mcp.tool()
    def cmap_upstream_regulators(up_genes: list[str], dn_genes: list[str], top_n: int = 20) -> dict:
        """Which shRNA/CRISPR/OE perturbations REPRODUCE this transcriptional program?
        Experimental upstream-regulator inference. Returns NCS, not tau, unless the query
        is restricted to Touchstone cell lines."""
        return cm.upstream_regulators(up_genes, dn_genes, top_n).reset_index().to_dict("records")

    @mcp.tool()
    def cmap_reversal_compounds(up_genes: list[str], dn_genes: list[str], top_n: int = 20) -> dict:
        """Compounds whose signature is anti-correlated with this program (repurposing
        direction). Requires >=3 independent signatures per compound."""
        return cm.reversal_compounds(up_genes, dn_genes, top_n).reset_index().to_dict("records")

    @mcp.tool()
    def state_regulator_test(gene: str) -> dict:
        """Distinguish state regulator from survival factor by crossing CMap knockdown
        response magnitude against DepMap fitness cost. Use whenever DepMap returns a
        negative for a gene with strong disease evidence."""
        return state_vs_survival(gene, dm, cm)

    if th is not None:
        @mcp.tool()
        def tahoe_context_dependence(drug: str, program_genes: list[str]) -> dict:
            """Cross-cell-line response to one drug for one program, measured in a shared
            well (no batch confound). Join to DepMap genomics to nominate response markers."""
            return th.context_dependence(drug, program_genes).to_dict("records")

        @mcp.tool()
        def tahoe_response_heterogeneity(cell_line: str, drug: str, program_genes: list[str]) -> dict:
            """Single-cell: is the response uniform, or does a subpopulation escape?
            Bimodality implies pre-existing resistance."""
            return th.response_heterogeneity(cell_line, drug, program_genes)

    return mcp
