"""
tri_evidence.py
===============
Triangulating three orthogonal evidence axes for target / mechanism prioritisation:

    DepMap        -> cell-intrinsic CAUSAL requirement (perturbational, context-conditioned)
    Open Targets  -> human GENETIC + clinical + tractability evidence (population-scale, correlative)
    Monarch       -> PHENOTYPE semantics + cross-species model evidence (ontology-grounded)

Design principles
-----------------
1. Bulk beats API. Genome-wide work uses local snapshots (DepMap CSVs, Open Targets
   parquet dumps / BigQuery, Monarch KGX dumps). The REST/GraphQL clients here are for
   *single-entity lookup* -- i.e. what an agent does interactively, not what a pipeline
   does over 20k genes.
2. Every returned object carries provenance (source, release, evidence IDs) so an agent
   can be forced to cite rather than confabulate.
3. Disk cache everything. These resources are versioned quarterly/monthly; reproducibility
   requires pinning, not live queries.

Tested against: Open Targets Platform GraphQL v4 (data 26.03), Monarch API v3,
DepMap 26Q1 file layout. Verify endpoints before production use.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import requests

CACHE = Path(os.environ.get("TRI_CACHE", "~/.cache/tri_evidence")).expanduser()
CACHE.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- #
# plumbing
# --------------------------------------------------------------------------- #
def _cached(namespace: str, key: str, fn, ttl_days: int = 30):
    """Filesystem JSON cache. Keeps agent loops from hammering public endpoints."""
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    p = CACHE / namespace / f"{h}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists() and (time.time() - p.stat().st_mtime) < ttl_days * 86400:
        return json.loads(p.read_text())
    val = fn()
    p.write_text(json.dumps(val))
    return val


def _get(url: str, params: dict | None = None, retries: int = 3) -> dict:
    for i in range(retries):
        r = requests.get(url, params=params, timeout=60)
        if r.status_code == 200:
            return r.json()
        if r.status_code in (429, 502, 503):
            time.sleep(2 ** i)
            continue
        r.raise_for_status()
    raise RuntimeError(f"GET failed after {retries} retries: {url}")


# --------------------------------------------------------------------------- #
# Open Targets Platform (GraphQL v4)
# --------------------------------------------------------------------------- #
class OpenTargets:
    URL = "https://api.platform.opentargets.org/api/v4/graphql"

    def query(self, gql: str, variables: dict) -> dict:
        key = gql + json.dumps(variables, sort_keys=True)

        def _run():
            r = requests.post(self.URL, json={"query": gql, "variables": variables}, timeout=90)
            r.raise_for_status()
            payload = r.json()
            if "errors" in payload:
                raise RuntimeError(payload["errors"])
            return payload["data"]

        return _cached("opentargets", key, _run)

    # -- entity resolution --------------------------------------------------
    def resolve(self, name: str, entity: str = "target") -> dict | None:
        """Symbol / free text -> stable ID. Agents must call this before anything else."""
        gql = """
        query s($q: String!, $e: [String!]) {
          search(queryString: $q, entityNames: $e, page: {index: 0, size: 5}) {
            hits { id name entity description }
          }
        }"""
        hits = self.query(gql, {"q": name, "e": [entity]})["search"]["hits"]
        return hits[0] if hits else None

    # -- target annotation --------------------------------------------------
    def target_profile(self, ensembl_id: str) -> dict:
        """Biotype, constraint (gnomAD o/e), tractability, safety liabilities.

        Genetic constraint is the single most underused field here: a gene with
        LOEUF < 0.35 that shows a strong disease association is a very different
        proposition (likely essential, narrow therapeutic window) from a tolerant one.
        """
        gql = """
        query t($id: String!) {
          target(ensemblId: $id) {
            id approvedSymbol biotype
            geneticConstraint { constraintType score oe oeLower oeUpper }
            tractability { label modality value }
            safetyLiabilities { event eventId datasource literature }
            subcellularLocations { location termSL }
            pathways { pathway pathwayId topLevelTerm }
          }
        }"""
        return self.query(gql, {"id": ensembl_id})["target"]

    # -- association with datatype decomposition ----------------------------
    def target_disease_evidence(self, ensembl_id: str, efo_id: str) -> dict:
        """Overall score PLUS the per-datatype breakdown.

        The breakdown is what matters scientifically. An association score of 0.6 driven
        entirely by `literature` (text mining) is a hypothesis; the same score driven by
        `genetic_association` + `somatic_mutation` is evidence. Never let an agent
        consume the scalar without the decomposition.
        """
        gql = """
        query td($ensemblId: String!, $efoId: String!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            associatedDiseases(efoIds: [$efoId]) {
              rows {
                disease { id name }
                score
                datatypeScores { id score }
              }
            }
          }
        }"""
        d = self.query(gql, {"ensemblId": ensembl_id, "efoId": efo_id})["target"]
        rows = d["associatedDiseases"]["rows"]
        return rows[0] if rows else {"score": 0.0, "datatypeScores": []}

    def disease_targets(self, efo_id: str, size: int = 200) -> pd.DataFrame:
        gql = """
        query d($efoId: String!, $size: Int!) {
          disease(efoId: $efoId) {
            id name
            associatedTargets(page: {index: 0, size: $size}) {
              count
              rows { target { id approvedSymbol biotype } score datatypeScores { id score } }
            }
          }
        }"""
        rows = self.query(gql, {"efoId": efo_id, "size": size})["disease"]["associatedTargets"]["rows"]
        recs = []
        for r in rows:
            rec = {
                "ensembl_id": r["target"]["id"],
                "symbol": r["target"]["approvedSymbol"],
                "biotype": r["target"]["biotype"],
                "ot_score": r["score"],
            }
            rec.update({f"ot_{d['id']}": d["score"] for d in r["datatypeScores"]})
            recs.append(rec)
        return pd.DataFrame(recs)


# --------------------------------------------------------------------------- #
# Monarch Initiative (API v3, Biolink-modelled KG)
# --------------------------------------------------------------------------- #
class Monarch:
    URL = "https://api-v3.monarchinitiative.org/v3/api"

    def search(self, q: str, category: str | None = None, limit: int = 5) -> list[dict]:
        params: dict[str, Any] = {"q": q, "limit": limit}
        if category:
            params["category"] = category
        return _get(f"{self.URL}/search", params).get("items", [])

    def associations(
        self,
        subject: str | None = None,
        object_: str | None = None,
        category: str | None = None,
        limit: int = 500,
    ) -> pd.DataFrame:
        """Generic Biolink association traversal.

        Useful categories:
          biolink:GeneToPhenotypicFeatureAssociation
          biolink:DiseaseToPhenotypicFeatureAssociation
          biolink:CausalGeneToDiseaseAssociation
          biolink:GeneToGeneHomologyAssociation
        """
        params: dict[str, Any] = {"limit": limit}
        if subject:
            params["subject"] = subject
        if object_:
            params["object"] = object_
        if category:
            params["category"] = category
        items = _get(f"{self.URL}/association", params).get("items", [])
        return pd.DataFrame(
            [
                {
                    "subject": a.get("subject"),
                    "subject_label": a.get("subject_label"),
                    "predicate": a.get("predicate"),
                    "object": a.get("object"),
                    "object_label": a.get("object_label"),
                    "species": a.get("subject_taxon_label"),
                    "provenance": a.get("primary_knowledge_source"),
                    "publications": a.get("publications"),
                }
                for a in items
            ]
        )

    def gene_phenotypes(self, gene_curie: str) -> pd.DataFrame:
        """HGNC:xxxx or NCBIGene:xxxx -> HPO/MP terms, human + model organism.

        The cross-species part is the payoff: a mouse-knockout phenotype is orthogonal
        evidence to anything in DepMap (whole-organism, includes non-cell-autonomous
        effects) and to Open Targets genetics (no human ascertainment bias).
        """
        return self.associations(
            subject=gene_curie, category="biolink:GeneToPhenotypicFeatureAssociation"
        )

    def disease_phenotypes(self, mondo_curie: str) -> pd.DataFrame:
        return self.associations(
            subject=mondo_curie, category="biolink:DiseaseToPhenotypicFeatureAssociation"
        )

    def semsim(self, subjects: Sequence[str], objects: Sequence[str], metric: str = "ancestor_information_content") -> dict:
        """Phenotype-profile semantic similarity (PhenIO-backed).

        This is Monarch's genuinely unique contribution: a *quantitative* distance between
        phenotype sets. Lets you ask "does my transcriptional module's phenotype signature
        resemble a known Mendelian disease?" without any keyword matching.
        """
        body = {"subjects": list(subjects), "objects": list(objects), "metric": metric}
        r = requests.post(f"{self.URL}/semsim/compare", json=body, timeout=120)
        r.raise_for_status()
        return r.json()


# --------------------------------------------------------------------------- #
# DepMap (local snapshot -- the portal is behind a captcha as of mid-2026,
# so download once by hand and version-pin the directory)
# --------------------------------------------------------------------------- #
@dataclass
class DepMap:
    root: Path
    release: str = "26Q1"
    _ge: pd.DataFrame | None = field(default=None, repr=False)
    _expr: pd.DataFrame | None = field(default=None, repr=False)
    _model: pd.DataFrame | None = field(default=None, repr=False)

    def __post_init__(self):
        self.root = Path(self.root).expanduser()

    def _load(self, name: str, parquet_cache: bool = True) -> pd.DataFrame:
        pq = self.root / f"{name}.parquet"
        if parquet_cache and pq.exists():
            return pd.read_parquet(pq)
        df = pd.read_csv(self.root / f"{name}.csv", index_col=0)
        # "SYMBOL (ENTREZ)" -> "SYMBOL"
        if not df.columns.str.contains(r"\(").all():
            pass
        else:
            df.columns = df.columns.str.replace(r"\s*\(\d+\)$", "", regex=True)
        if parquet_cache:
            df.to_parquet(pq)
        return df

    @property
    def gene_effect(self) -> pd.DataFrame:
        """Chronos gene effect: rows = ModelID, cols = gene symbol. ~0 = neutral, -1 = median essential."""
        if self._ge is None:
            self._ge = self._load("CRISPRGeneEffect")
        return self._ge

    @property
    def expression(self) -> pd.DataFrame:
        if self._expr is None:
            self._expr = self._load("OmicsExpressionProteinCodingGenesTPMLogp1")
        return self._expr

    @property
    def model(self) -> pd.DataFrame:
        if self._model is None:
            self._model = pd.read_csv(self.root / "Model.csv", index_col=0)
        return self._model

    # -- the three operations that actually matter -------------------------
    def selectivity(self, gene: str) -> dict:
        """Distinguish pan-essential (bad target) from selectively essential (good target).

        Uses the likelihood-ratio-ish heuristic: a gene is *selective* if a minority of
        lines are strongly dependent while the bulk are neutral. Pan-essentials have low
        variance and a mean well below zero -- they will kill everything, including the
        patient's gut epithelium.
        """
        if gene not in self.gene_effect.columns:
            return {"gene": gene, "in_library": False}
        v = self.gene_effect[gene].dropna()
        frac_dep = float((v < -0.5).mean())
        return {
            "gene": gene,
            "in_library": True,
            "n_lines": int(v.size),
            "mean_chronos": float(v.mean()),
            "sd_chronos": float(v.std()),
            "frac_dependent": frac_dep,
            "class": (
                "pan_essential" if v.mean() < -0.6 and frac_dep > 0.8
                else "selective" if 0.02 < frac_dep < 0.5
                else "non_essential"
            ),
            "release": self.release,
        }

    def lineage_enrichment(self, gene: str, lineage: str) -> dict:
        """Is the dependency enriched in a lineage of interest (e.g. 'Pancreas')?"""
        ge = self.gene_effect
        if gene not in ge.columns:
            return {"gene": gene, "in_library": False}
        ids = self.model.index[self.model["OncotreeLineage"] == lineage]
        inl = ge.loc[ge.index.intersection(ids), gene].dropna()
        out = ge.loc[ge.index.difference(ids), gene].dropna()
        from scipy import stats
        t, p = stats.mannwhitneyu(inl, out, alternative="less")
        return {
            "gene": gene,
            "lineage": lineage,
            "n_lineage": int(inl.size),
            "mean_in": float(inl.mean()),
            "mean_out": float(out.mean()),
            "delta": float(inl.mean() - out.mean()),
            "p_mwu": float(p),
        }

    def coessentiality(self, gene: str, top_n: int = 25, min_lines: int = 300) -> pd.DataFrame:
        """Genes whose dependency profile correlates with `gene` across cell lines.

        This is the highest-value systems-biology operation in DepMap and it is
        *functionally orthogonal to co-expression*: co-essentiality recovers complex
        membership and pathway epistasis that co-expression misses entirely, because it
        is measured under perturbation rather than at steady state. Use it to test
        whether a co-expression module you derived is a real functional unit.
        """
        ge = self.gene_effect
        if gene not in ge.columns:
            return pd.DataFrame()
        sub = ge.dropna(axis=1, thresh=min_lines)
        x = sub[gene]
        # Pearson across lines, vectorised
        xc = (x - x.mean()) / x.std()
        Y = (sub - sub.mean()) / sub.std()
        r = Y.mul(xc, axis=0).mean() * (len(sub) / (len(sub) - 1))
        r = r.drop(gene).sort_values(ascending=False)
        return pd.concat([r.head(top_n), r.tail(top_n)]).rename("pearson_r").to_frame()

    def expression_dependency_decoupling(self, gene: str) -> dict:
        """Correlate a gene's own expression with its own dependency.

        Strong negative r => 'expression-addicted' lineage-survival factor (the GATA6 /
        SOX2 / MITF pattern). Near-zero r with high mean expression => the gene is
        expressed but dispensable in vitro, which is exactly the signature you get for
        genes whose function is non-cell-autonomous (secreted, immune, stromal crosstalk).
        """
        ge, ex = self.gene_effect, self.expression
        if gene not in ge.columns or gene not in ex.columns:
            return {"gene": gene, "available": False}
        common = ge.index.intersection(ex.index)
        a, b = ge.loc[common, gene], ex.loc[common, gene]
        m = a.notna() & b.notna()
        r = float(np.corrcoef(a[m], b[m])[0, 1])
        return {
            "gene": gene,
            "available": True,
            "r_expr_vs_dependency": r,
            "mean_expression_log2tpm1": float(b[m].mean()),
            "interpretation": (
                "expression_addicted" if r < -0.3
                else "expressed_but_dispensable_in_vitro" if b[m].mean() > 3 and abs(r) < 0.15
                else "uninformative"
            ),
        }


# --------------------------------------------------------------------------- #
# Triangulation + hypothesis flags
# --------------------------------------------------------------------------- #
@dataclass
class TriEvidence:
    symbol: str
    disease: str
    depmap: dict = field(default_factory=dict)
    opentargets: dict = field(default_factory=dict)
    monarch: dict = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, default=str)


def triangulate(
    symbol: str,
    disease_query: str,
    dm: DepMap,
    ot: OpenTargets | None = None,
    mk: Monarch | None = None,
    lineage: str | None = None,
) -> TriEvidence:
    ot = ot or OpenTargets()
    mk = mk or Monarch()

    tgt = ot.resolve(symbol, "target")
    dis = ot.resolve(disease_query, "disease")
    ev = TriEvidence(symbol=symbol, disease=dis["name"] if dis else disease_query)

    if tgt:
        ev.opentargets["profile"] = ot.target_profile(tgt["id"])
        if dis:
            ev.opentargets["association"] = ot.target_disease_evidence(tgt["id"], dis["id"])

    ev.depmap["selectivity"] = dm.selectivity(symbol)
    ev.depmap["expr_dep"] = dm.expression_dependency_decoupling(symbol)
    if lineage:
        ev.depmap["lineage"] = dm.lineage_enrichment(symbol, lineage)

    hits = mk.search(symbol, category="biolink:Gene")
    if hits:
        curie = hits[0]["id"]
        ph = mk.gene_phenotypes(curie)
        ev.monarch = {
            "curie": curie,
            "n_phenotype_assocs": len(ph),
            "human_phenotypes": ph.loc[ph.species.eq("Homo sapiens"), "object_label"].head(20).tolist() if len(ph) else [],
            "model_phenotypes": ph.loc[~ph.species.eq("Homo sapiens"), "object_label"].head(20).tolist() if len(ph) else [],
        }

    ev.flags = hypothesis_flags(ev)
    return ev


def hypothesis_flags(ev: TriEvidence) -> list[str]:
    """Encode the *disagreement* patterns. Agreement is reassuring; disagreement is where
    the novel hypotheses live, because it means one assay is blind to a real mechanism."""
    f: list[str] = []
    dm_sel = ev.depmap.get("selectivity", {})
    ot_assoc = ev.opentargets.get("association", {})
    dts = {d["id"]: d["score"] for d in ot_assoc.get("datatypeScores", [])}
    genetic = dts.get("genetic_association", 0.0)
    lit = dts.get("literature", 0.0)

    if dm_sel.get("class") == "selective" and genetic > 0.3:
        f.append("CONVERGENT: selective in vitro dependency + human genetic support -> strongest tier")

    if dm_sel.get("class") == "selective" and genetic < 0.05 and ot_assoc.get("score", 0) < 0.2:
        f.append(
            "DEPMAP_ORPHAN: real cellular dependency with no human disease evidence. "
            "Either genuinely novel, or the dependency is a cell-culture artefact. "
            "Check lineage enrichment and co-essentiality neighbours before believing it."
        )

    if ot_assoc.get("score", 0) > 0.4 and dm_sel.get("class") == "non_essential":
        f.append(
            "NON_CELL_AUTONOMOUS: disease-associated but dispensable in monoculture. "
            "Strong prior for a microenvironment / immune / secreted mechanism that "
            "cell-line CRISPR cannot see. Escalate to co-culture or in vivo screens."
        )

    if lit > 0.5 and genetic < 0.05 and ot_assoc.get("score", 0) > 0.3:
        f.append("LITERATURE_INFLATED: association driven by text mining, not orthogonal evidence")

    if dm_sel.get("class") == "pan_essential":
        f.append("PAN_ESSENTIAL: poor therapeutic index; treat as biology, not as a target")

    if ev.depmap.get("expr_dep", {}).get("interpretation") == "expressed_but_dispensable_in_vitro":
        f.append("EXPRESSED_NOT_REQUIRED: consistent with a non-cell-autonomous or redundant role")

    mo = ev.monarch.get("model_phenotypes", [])
    if mo and not ev.monarch.get("human_phenotypes"):
        f.append("MODEL_ONLY_PHENOTYPE: mouse/zebrafish phenotype with no human counterpart yet -- "
                 "candidate for a novel human phenotype hypothesis")

    if not dm_sel.get("in_library", True):
        f.append("NOT_IN_CRISPR_LIBRARY: non-coding / lncRNA / poorly annotated. DepMap is silent by "
                 "construction, not by evidence. Use a CRISPRi/CRISPRa screen or guilt-by-association "
                 "via cis-neighbours instead of concluding 'no dependency'.")

    return f


# --------------------------------------------------------------------------- #
# MCP tool layer (agent interface)
# --------------------------------------------------------------------------- #
def build_mcp_server(depmap_root: str | Path):
    """Expose the above as MCP tools. Requires `pip install fastmcp`.

    Two deliberate design choices:
      * Tools return *structured, provenance-bearing* dicts, never prose. The model does
        the reasoning; the tool does the retrieval. This is what keeps hypotheses falsifiable.
      * `resolve_entity` is a separate mandatory tool. Free-text gene/disease names are the
        single largest source of silent agent error; forcing explicit CURIE resolution makes
        the failure loud instead of quiet.
    """
    from fastmcp import FastMCP

    mcp = FastMCP("tri-evidence")
    dm, ot, mk = DepMap(depmap_root), OpenTargets(), Monarch()

    @mcp.tool()
    def resolve_entity(name: str, kind: str = "target") -> dict:
        """Resolve a gene symbol or disease name to a stable ID (Ensembl / EFO / MONDO).
        ALWAYS call this before any other tool."""
        return ot.resolve(name, kind) or {"error": f"no {kind} match for {name!r}"}

    @mcp.tool()
    def depmap_dependency(gene: str, lineage: str | None = None) -> dict:
        """Chronos dependency profile: pan-essential vs selective, lineage enrichment,
        expression-dependency coupling. Cell-intrinsic evidence only."""
        out = {"selectivity": dm.selectivity(gene),
               "expression_coupling": dm.expression_dependency_decoupling(gene)}
        if lineage:
            out["lineage_enrichment"] = dm.lineage_enrichment(gene, lineage)
        return out

    @mcp.tool()
    def depmap_coessential_partners(gene: str, top_n: int = 20) -> dict:
        """Genes co-essential with `gene` across ~1100 cell lines. Reveals functional
        complex membership orthogonal to co-expression."""
        df = dm.coessentiality(gene, top_n=top_n)
        return {"gene": gene, "release": dm.release, "partners": df.to_dict()["pearson_r"]}

    @mcp.tool()
    def opentargets_evidence(ensembl_id: str, efo_id: str) -> dict:
        """Target-disease association WITH datatype decomposition and tractability.
        Report the decomposition, never the bare score."""
        return {"association": ot.target_disease_evidence(ensembl_id, efo_id),
                "target": ot.target_profile(ensembl_id)}

    @mcp.tool()
    def monarch_phenotypes(gene_or_disease_curie: str) -> dict:
        """Ontology-grounded phenotypes (HPO/MP) for a gene or disease, human and model organism."""
        df = mk.associations(subject=gene_or_disease_curie)
        return {"curie": gene_or_disease_curie, "n": len(df),
                "associations": df.head(50).to_dict(orient="records")}

    @mcp.tool()
    def monarch_phenotype_similarity(profile_a: list[str], profile_b: list[str]) -> dict:
        """Semantic similarity between two HPO/MP term sets. Use to test whether an
        experimentally derived phenotype signature resembles a known disease."""
        return mk.semsim(profile_a, profile_b)

    @mcp.tool()
    def triangulate_gene(symbol: str, disease: str, lineage: str | None = None) -> dict:
        """Full three-axis evidence dossier plus contradiction flags. Start here for
        hypothesis generation; the FLAGS are the interesting part, not the scores."""
        return asdict(triangulate(symbol, disease, dm, ot, mk, lineage))

    return mcp


if __name__ == "__main__":
    dm = DepMap(os.environ.get("DEPMAP_ROOT", "~/data/depmap/26Q1"))
    ev = triangulate("GATA6", "pancreatic adenocarcinoma", dm, lineage="Pancreas")
    print(ev.to_json())
