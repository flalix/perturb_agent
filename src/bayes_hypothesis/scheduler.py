"""
L3 — Scheduler bayesiano de evidência.

Responde: "com orçamento B, qual a PRÓXIMA evidência a adquirir?"

Núcleo formal — ganho esperado de informação (EIG) para um ensaio binário
caracterizado por sensibilidade `se` e especificidade `sp`, sobre uma
hipótese com posterior corrente p:

    P(+)   = p·se + (1-p)·(1-sp)
    p|+    = p·se / P(+)
    p|−    = p·(1-se) / P(−)
    EIG    = H(p) − [ P(+)·H(p|+) + P(−)·H(p|−) ]

com H a entropia binária em bits. Score = EIG / custo.

Duas propriedades que caem de graça e são metodologicamente importantes:

  * EIG é MÁXIMO perto de p = 0,5 e colapsa nos extremos. O scheduler
    naturalmente para de gastar em hipóteses já decididas — o oposto do
    comportamento de torneio Elo, que continua polindo o líder.

  * Um ensaio pouco informativo (se+sp ≈ 1) recebe EIG ≈ 0 por mais barato
    que seja. Não há como "compensar" desenho ruim com volume.

RESERVA DE DIVERSIDADE: uma fração do orçamento é alocada por rodízio sobre
descritores de mecanismo/escala distintos, contra o colapso de diversidade
documentado em buscas evolutivas guiadas por LLM.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .evidence import EvidenceKind
from .hypothesis import Claim, Hypothesis, State
from .provenance import Action, Actor, Ledger


def binary_entropy(p: float) -> float:
    p = min(max(p, 1e-12), 1 - 1e-12)
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


@dataclass(frozen=True)
class Assay:
    """Um ensaio candidato (in silico ou de bancada) e seu perfil operacional."""

    aid: str
    label: str
    kind: EvidenceKind
    sensitivity: float
    specificity: float
    cost_usd: float
    days: float
    applicable_to: tuple[str, ...] = ()      # hids; vazio = qualquer
    identifies_causal: bool = False          # satisfaz G4?
    feasible: bool = True
    blocker: str = ""

    def informativeness(self) -> float:
        """Índice de Youden. ≤0 significa ensaio sem poder discriminante."""
        return self.sensitivity + self.specificity - 1.0


def expected_information_gain(p: float, assay: Assay) -> float:
    """EIG em bits. Ver docstring do módulo."""
    se, sp = assay.sensitivity, assay.specificity
    p_pos = p * se + (1 - p) * (1 - sp)
    p_neg = 1.0 - p_pos
    if p_pos <= 1e-12 or p_neg <= 1e-12:
        return 0.0
    post_pos = p * se / p_pos
    post_neg = p * (1 - se) / p_neg
    expected_post_h = p_pos * binary_entropy(post_pos) + p_neg * binary_entropy(post_neg)
    return max(binary_entropy(p) - expected_post_h, 0.0)


@dataclass(frozen=True)
class Plan:
    hypothesis_id: str
    assay: Assay
    eig_bits: float
    cost_usd: float
    score: float                 # bits por 1k USD
    rationale: str
    slot: str                    # "optimal" | "diversity_reserve"


@dataclass
class EvidenceScheduler:
    budget_usd: float
    diversity_reserve: float = 0.25      # fração do orçamento fora do ranking guloso
    min_eig_bits: float = 0.02           # piso: abaixo disto, o ensaio é teatro
    _spent: float = field(default=0.0, init=False)

    # ------------------------------------------------------------------

    def _candidates(self, hyps: list[Hypothesis], assays: list[Assay]) -> list[Plan]:
        out: list[Plan] = []
        for h in hyps:
            if h.state not in (State.UNDER_TEST, State.REGISTERED):
                continue
            p = h.posterior
            for a in assays:
                if a.applicable_to and h.hid not in a.applicable_to:
                    continue
                if not a.feasible:
                    continue
                eig = expected_information_gain(p, a)

                # bônus de guard: se a hipótese é causal e carece de evidência
                # intervencional, um ensaio que satisfaz G4 vale mais do que o
                # EIG sozinho indica — ele destrava a transição de estado.
                unlock = 1.0
                notes = []
                if h.claim is Claim.CAUSAL and not h.has_interventional_evidence() and a.identifies_causal:
                    unlock = 2.0
                    notes.append("destrava G4 (causal)")
                if not h.has_non_literature_evidence() and a.kind.value not in ("literature", "kg_structure"):
                    unlock *= 1.5
                    notes.append("destrava G3 (não-circularidade)")

                if eig < self.min_eig_bits:
                    continue
                cost = max(a.cost_usd, 1.0)
                score = (eig * unlock) / (cost / 1000.0)
                out.append(
                    Plan(
                        hypothesis_id=h.hid,
                        assay=a,
                        eig_bits=eig,
                        cost_usd=a.cost_usd,
                        score=score,
                        rationale="; ".join(notes) or f"EIG={eig:.3f} bits @ p={p:.2f}",
                        slot="optimal",
                    )
                )
        return sorted(out, key=lambda pl: pl.score, reverse=True)

    def plan(self, hyps: list[Hypothesis], assays: list[Assay], ledger: Ledger, actor: Actor) -> list[Plan]:
        """Seleciona ações até esgotar o orçamento, com reserva de diversidade."""
        by_hid = {h.hid: h for h in hyps}
        ranked = self._candidates(hyps, assays)

        greedy_budget = self.budget_usd * (1 - self.diversity_reserve)
        reserve_budget = self.budget_usd * self.diversity_reserve

        chosen: list[Plan] = []
        used_pairs: set[tuple[str, str]] = set()

        # --- slot 1: guloso por bits/USD ---
        spent = 0.0
        for pl in ranked:
            if spent + pl.cost_usd > greedy_budget:
                continue
            key = (pl.hypothesis_id, pl.assay.aid)
            if key in used_pairs:
                continue
            chosen.append(pl)
            used_pairs.add(key)
            spent += pl.cost_usd

        # --- slot 2: reserva de diversidade ---
        # rodízio sobre (mecanismo, escala) ainda não contemplados no slot 1.
        covered = {
            (by_hid[pl.hypothesis_id].mechanism_descriptor, by_hid[pl.hypothesis_id].scale_descriptor)
            for pl in chosen
        }
        r_spent = 0.0
        for pl in ranked:
            h = by_hid[pl.hypothesis_id]
            desc = (h.mechanism_descriptor, h.scale_descriptor)
            if desc in covered:
                continue
            if r_spent + pl.cost_usd > reserve_budget:
                continue
            key = (pl.hypothesis_id, pl.assay.aid)
            if key in used_pairs:
                continue
            chosen.append(
                Plan(**{**pl.__dict__, "slot": "diversity_reserve",
                        "rationale": pl.rationale + " | reserva de diversidade"})
            )
            used_pairs.add(key)
            covered.add(desc)
            r_spent += pl.cost_usd

        self._spent = spent + r_spent
        ledger.append(
            actor,
            Action.EXPERIMENT_PROPOSED,
            "SCHEDULER",
            {
                "n_plans": len(chosen),
                "spent_usd": self._spent,
                "budget_usd": self.budget_usd,
                "greedy": sum(1 for c in chosen if c.slot == "optimal"),
                "reserve": sum(1 for c in chosen if c.slot == "diversity_reserve"),
                "total_bits": round(sum(c.eig_bits for c in chosen), 4),
            },
        )
        return chosen

    @property
    def spent(self) -> float:
        return self._spent


def unidentifiable_report(h: Hypothesis, assays: list[Assay]) -> str | None:
    """Se nenhum ensaio disponível é informativo para a hipótese, devolve a
    MEDIDA MÍNIMA FALTANTE. Este é o output mais útil do sistema: em vez de
    fabricar uma análise causal, o sistema diz o que precisa ser medido.
    """
    usable = [a for a in assays if a.feasible and expected_information_gain(h.posterior, a) > 0.02]
    if usable:
        return None
    blockers = sorted({a.blocker for a in assays if a.blocker})
    need = []
    if h.claim is Claim.CAUSAL and not h.has_interventional_evidence():
        need.append("uma intervenção (CRISPR/CRISPRi no tipo celular do escopo, "
                    "ou instrumento genético cis-eQTL com colocalização)")
    if not h.has_non_literature_evidence():
        need.append("um dado primário independente da literatura")
    if not need:
        need.append("um ensaio com índice de Youden > 0 para este estimando")
    return (
        f"NÃO IDENTIFICÁVEL com os ensaios disponíveis. Medida mínima faltante: "
        + "; ".join(need)
        + (f". Bloqueios registrados: {', '.join(blockers)}." if blockers else ".")
    )
