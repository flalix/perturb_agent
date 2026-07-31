"""
L2 — Hipótese como objeto tipado de primeira classe, com ciclo de vida.

Os guards da máquina de estados são o coração do desenho. Eles transformam
princípios metodológicos em restrições executáveis, em vez de instruções de
prompt que o modelo pode ignorar:

  G1  Não se registra hipótese sem predição QUANTITATIVA e FALSIFICÁVEL.
  G2  Evidência confirmatória só pode ser anexada após o commit criptográfico
      da predição  ................................................. anti-HARKing.
  G3  Nenhuma hipótese chega a SURVIVING só com literatura/KG  .... anti-circularidade.
  G4  Afirmação CAUSAL exige ≥1 evidência intervencional (perturbação,
      instrumento genético ou bancada)  ............ associação ≠ causalidade.
  G5  SURVIVING exige um número mínimo de tentativas de refutação
      INDEPENDENTES já executadas  ............... o crítico não é decorativo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable

from .evidence import (
    INTERVENTIONAL_KINDS,
    NON_LITERATURE_KINDS,
    DependenceModel,
    EvidenceItem,
    aggregate_log_lr,
)
from .ontology import Scope, Triple
from .provenance import Action, Actor, Ledger, PreRegistration, content_hash


class State(str, Enum):
    PROPOSED = "proposed"
    REGISTERED = "registered"      # predições seladas
    UNDER_TEST = "under_test"
    SURVIVING = "surviving"        # NÃO é "verdadeira"; é "não refutada sob pressão"
    REFUTED = "refuted"
    PARKED = "parked"              # não identificável / sem dado / custo proibitivo
    SUPERSEDED = "superseded"


class Claim(str, Enum):
    ASSOCIATIVE = "associative"
    CAUSAL = "causal"


@dataclass(frozen=True)
class Prediction:
    """Predição quantitativa. `metric`, `operator` e `threshold` tornam a
    avaliação mecânica — sem espaço para reinterpretação retroativa."""

    pid: str
    description: str
    metric: str                 # e.g. "log2FC(HOXB-AS3, P1 vs P2)"
    operator: str               # ">=", "<=", ">", "<", "=="
    threshold: float
    assay: str                  # como se mede
    critical: bool = False      # se falhar, refuta a hipótese (não só reduz posterior)

    def evaluate(self, observed: float) -> bool:
        ops: dict[str, Callable[[float, float], bool]] = {
            ">=": lambda a, b: a >= b,
            "<=": lambda a, b: a <= b,
            ">": lambda a, b: a > b,
            "<": lambda a, b: a < b,
            "==": lambda a, b: math.isclose(a, b, rel_tol=1e-9),
        }
        return ops[self.operator](observed, self.threshold)

    def commit_payload(self) -> dict:
        return {
            "pid": self.pid,
            "metric": self.metric,
            "operator": self.operator,
            "threshold": self.threshold,
            "assay": self.assay,
            "critical": self.critical,
        }


@dataclass(frozen=True)
class RefutationAttempt:
    """Tentativa de derrubar a hipótese. `independent_of` registra o gerador
    para que possamos detectar erros correlacionados gerador/crítico."""

    rid: str
    actor: Actor
    strategy: str               # "alternative_explanation" | "placebo_refuter" |
                                # "negative_control" | "confounder_probe" | "replication"
    independent_of: str         # agent_id do gerador
    succeeded: bool             # True = conseguiu derrubar
    notes: str = ""


@dataclass
class Cost:
    money_usd: float = 0.0
    days: float = 0.0
    feasible: bool = True
    blocker: str = ""


class GuardError(RuntimeError):
    """Violação de guard da máquina de estados. Nunca capture e ignore."""


@dataclass
class Hypothesis:
    hid: str
    statement_nl: str
    formal: Triple
    scope: Scope
    claim: Claim
    log_prior_odds: float
    predictions: list[Prediction] = field(default_factory=list)
    evidence: list[EvidenceItem] = field(default_factory=list)
    refutations: list[RefutationAttempt] = field(default_factory=list)
    state: State = State.PROPOSED
    prereg: PreRegistration | None = None
    cost_to_test: Cost = field(default_factory=Cost)
    mechanism_descriptor: str = "unspecified"   # eixo de diversidade p/ o scheduler
    scale_descriptor: str = "molecular"         # molecular|cellular|tissue|organ|systemic
    parked_reason: str = ""
    _dependence: DependenceModel = field(default_factory=DependenceModel)

    # --- inferência -------------------------------------------------------

    @property
    def log_posterior_odds(self) -> float:
        total, _ = aggregate_log_lr(self.evidence, self._dependence)
        return self.log_prior_odds + total

    @property
    def posterior(self) -> float:
        lo = self.log_posterior_odds
        return 1.0 / (1.0 + math.exp(-lo))

    @property
    def prior(self) -> float:
        return 1.0 / (1.0 + math.exp(-self.log_prior_odds))

    def evidence_breakdown(self) -> dict[str, float]:
        _, bd = aggregate_log_lr(self.evidence, self._dependence)
        return bd

    # --- predicados dos guards -------------------------------------------

    def has_non_literature_evidence(self) -> bool:
        return any(e.kind in NON_LITERATURE_KINDS for e in self.evidence)

    def has_interventional_evidence(self) -> bool:
        return any(e.kind in INTERVENTIONAL_KINDS for e in self.evidence)

    def independent_refutation_count(self) -> int:
        seen = {r.actor.agent_id for r in self.refutations if r.actor.agent_id != r.independent_of}
        return len(seen)

    def failed_critical_prediction(self) -> bool:
        return any(e.decisive for e in self.evidence)

    # --- transições -------------------------------------------------------

    def register(self, ledger: Ledger, actor: Actor) -> PreRegistration:
        """PROPOSED -> REGISTERED. Sela as predições (G1, base de G2)."""
        if self.state is not State.PROPOSED:
            raise GuardError(f"register() exige PROPOSED, estado atual = {self.state}")
        if not self.predictions:
            raise GuardError("G1: hipótese sem predição quantitativa não é registrável")

        payload = {
            "hid": self.hid,
            "statement": self.statement_nl,
            "triple": str(self.formal),
            "scope": self.scope.key(),
            "claim": self.claim.value,
            "predictions": [p.commit_payload() for p in self.predictions],
        }
        ev = ledger.append(actor, Action.PREDICTION_REGISTERED, self.hid, payload)
        self.prereg = PreRegistration(
            hypothesis_id=self.hid,
            commit_hash=content_hash(payload),
            committed_at=ev.timestamp,
            ledger_seq=ev.seq,
            payload_preview={"n_predictions": len(self.predictions)},
        )
        self._transition(ledger, actor, State.REGISTERED, "predições seladas")
        return self.prereg

    def attach_evidence(self, item: EvidenceItem, ledger: Ledger, actor: Actor) -> None:
        """Anexa evidência. G2 barra confirmatória sem pré-registro."""
        if item.confirmatory and self.prereg is None:
            raise GuardError(
                f"G2 (anti-HARKing): evidência confirmatória {item.eid!r} exige pré-registro. "
                "Chame register() antes de tocar no conjunto de confirmação."
            )
        if self.state in (State.REFUTED, State.SUPERSEDED):
            raise GuardError(f"não se anexa evidência a hipótese em {self.state}")

        before = self.log_posterior_odds
        self.evidence.append(item)
        after = self.log_posterior_odds

        ledger.append(
            actor,
            Action.EVIDENCE_ATTACHED,
            self.hid,
            {
                "eid": item.eid,
                "kind": item.kind.value,
                "source": item.source_id,
                "raw_lr": item.raw_lr,
                "effective_lr": item.effective_lr,
                "reliability": item.effective_reliability,
                "independence_class": item.independence_class,
                "confirmatory": item.confirmatory,
            },
        )
        ledger.append(
            actor,
            Action.POSTERIOR_UPDATED,
            self.hid,
            {
                "log_odds_before": before,
                "log_odds_after": after,
                "delta": after - before,
                "posterior": self.posterior,
            },
        )
        if self.state is State.REGISTERED:
            self._transition(ledger, actor, State.UNDER_TEST, "primeira evidência anexada")

    def record_refutation(self, attempt: RefutationAttempt, ledger: Ledger) -> None:
        self.refutations.append(attempt)
        ledger.append(
            attempt.actor,
            Action.REFUTATION_ATTEMPTED,
            self.hid,
            {
                "rid": attempt.rid,
                "strategy": attempt.strategy,
                "succeeded": attempt.succeeded,
                "independent_of": attempt.independent_of,
                "notes": attempt.notes,
            },
        )

    def adjudicate(
        self,
        ledger: Ledger,
        actor: Actor,
        p_survive: float = 0.90,
        p_refute: float = 0.05,
        min_independent_refutations: int = 2,
    ) -> State:
        """Aplica G3, G4, G5 e decide o estado. Idempotente.

        Nota importante: SURVIVING nunca significa "verdadeira". Significa
        "não refutada sob pressão adversarial, com evidência não-circular".
        """
        if self.state not in (State.UNDER_TEST, State.SURVIVING):
            return self.state

        p = self.posterior
        reasons: list[str] = []

        # refutação
        if self.failed_critical_prediction():
            return self._transition(ledger, actor, State.REFUTED, "predição crítica pré-registrada falhou")
        if any(r.succeeded for r in self.refutations):
            return self._transition(ledger, actor, State.REFUTED, "refutação adversarial bem-sucedida")
        if p <= p_refute:
            return self._transition(ledger, actor, State.REFUTED, f"posterior {p:.3f} ≤ {p_refute}")

        # sobrevivência — todos os guards precisam passar
        if p < p_survive:
            reasons.append(f"posterior {p:.3f} < {p_survive}")
        if not self.has_non_literature_evidence():
            reasons.append("G3: sem evidência não-literária (risco de circularidade)")
        if self.claim is Claim.CAUSAL and not self.has_interventional_evidence():
            reasons.append("G4: afirmação causal sem evidência intervencional")
        n_ind = self.independent_refutation_count()
        if n_ind < min_independent_refutations:
            reasons.append(f"G5: {n_ind}/{min_independent_refutations} refutações independentes")

        if reasons:
            if self.state is State.SURVIVING:  # rebaixamento por evidência nova
                self._transition(ledger, actor, State.UNDER_TEST, "; ".join(reasons))
            else:
                ledger.append(actor, Action.STATE_TRANSITION, self.hid,
                              {"from": self.state.value, "to": self.state.value,
                               "blocked_by": reasons})
            return self.state

        return self._transition(ledger, actor, State.SURVIVING, f"posterior {p:.3f}, guards OK")

    def park(self, ledger: Ledger, actor: Actor, reason: str, missing_measurement: str = "") -> State:
        """Saída honesta: não identificável / sem dado. O campo
        `missing_measurement` é o output mais valioso do sistema — é a ponte
        para a bancada."""
        self.parked_reason = reason
        ledger.append(actor, Action.EXPERIMENT_PROPOSED, self.hid,
                      {"reason": reason, "missing_measurement": missing_measurement})
        return self._transition(ledger, actor, State.PARKED, reason)

    def _transition(self, ledger: Ledger, actor: Actor, new: State, reason: str) -> State:
        old, self.state = self.state, new
        ledger.append(actor, Action.STATE_TRANSITION, self.hid,
                      {"from": old.value, "to": new.value, "reason": reason,
                       "posterior": self.posterior})
        return new

    # --- relatório --------------------------------------------------------

    def summary(self) -> dict:
        return {
            "hid": self.hid,
            "state": self.state.value,
            "claim": self.claim.value,
            "statement": self.statement_nl,
            "triple": str(self.formal),
            "scope": self.scope.key(),
            "prior": round(self.prior, 4),
            "posterior": round(self.posterior, 4),
            "n_evidence": len(self.evidence),
            "non_literature": self.has_non_literature_evidence(),
            "interventional": self.has_interventional_evidence(),
            "independent_refutations": self.independent_refutation_count(),
            "prereg_seq": self.prereg.ledger_seq if self.prereg else None,
            "parked_reason": self.parked_reason or None,
        }
