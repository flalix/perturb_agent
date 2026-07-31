"""
L6 — Proveniência e auditoria (transversal).

Princípio de desenho: auditoria NÃO é um agente que lê o output de outros
agentes (isso audita prosa). É um ledger append-only com encadeamento de hash,
onde cada evento aponta para seus pais, formando um DAG de proveniência.

Duas propriedades que isso compra:
  1. Detecção de adulteração: qualquer reescrita retroativa quebra a cadeia.
  2. Selo temporal criptográfico: permite provar que uma predição foi
     registrada ANTES do acesso ao dado de validação (anti-HARKing).

O relatório final deve ser gerado por consulta a este ledger, nunca escrito
"de memória" por um LLM.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Iterator


class Action(str, Enum):
    HYPOTHESIS_PROPOSED = "hypothesis.proposed"
    PREDICTION_REGISTERED = "prediction.registered"  # o commit anti-HARKing
    EVIDENCE_ATTACHED = "evidence.attached"
    POSTERIOR_UPDATED = "posterior.updated"
    REFUTATION_ATTEMPTED = "refutation.attempted"
    STATE_TRANSITION = "state.transition"
    EXPERIMENT_PROPOSED = "experiment.proposed"
    DATA_ACCESSED = "data.accessed"
    HUMAN_DECISION = "human.decision"


@dataclass(frozen=True)
class Actor:
    """Quem agiu. Modelo e versão importam para reprodutibilidade e para
    diagnosticar erros correlacionados entre gerador e crítico."""

    agent_id: str
    role: str                      # "generator" | "critic" | "statistician" | "human" ...
    model: str = "n/a"
    model_version: str = "n/a"

    def fingerprint(self) -> str:
        return f"{self.agent_id}/{self.role}/{self.model}@{self.model_version}"


def _canonical(payload: Any) -> str:
    """Serialização determinística — hash estável entre execuções."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)


def content_hash(payload: Any) -> str:
    return hashlib.sha256(_canonical(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class Event:
    seq: int
    timestamp: str
    actor: Actor
    action: Action
    subject_id: str                # id da hipótese / evidência
    payload_hash: str
    payload: dict[str, Any]
    parents: tuple[int, ...] = ()
    prev_hash: str = ""
    this_hash: str = ""

    def digest(self) -> str:
        return content_hash(
            {
                "seq": self.seq,
                "timestamp": self.timestamp,
                "actor": self.actor.fingerprint(),
                "action": self.action.value,
                "subject_id": self.subject_id,
                "payload_hash": self.payload_hash,
                "parents": list(self.parents),
                "prev_hash": self.prev_hash,
            }
        )


class Ledger:
    """Ledger append-only. Não há método de remoção ou edição — por desenho."""

    def __init__(self) -> None:
        self._events: list[Event] = []

    def append(
        self,
        actor: Actor,
        action: Action,
        subject_id: str,
        payload: dict[str, Any],
        parents: tuple[int, ...] = (),
    ) -> Event:
        prev = self._events[-1].this_hash if self._events else ""
        draft = Event(
            seq=len(self._events),
            timestamp=datetime.now(timezone.utc).isoformat(),
            actor=actor,
            action=action,
            subject_id=subject_id,
            payload_hash=content_hash(payload),
            payload=payload,
            parents=parents,
            prev_hash=prev,
        )
        sealed = Event(**{**asdict(draft), "actor": actor, "action": action, "this_hash": draft.digest()})
        self._events.append(sealed)
        return sealed

    # ---- consultas -------------------------------------------------------

    def __len__(self) -> int:
        return len(self._events)

    def __iter__(self) -> Iterator[Event]:
        return iter(self._events)

    def for_subject(self, subject_id: str) -> list[Event]:
        return [e for e in self._events if e.subject_id == subject_id]

    def first_of(self, subject_id: str, action: Action) -> Event | None:
        for e in self._events:
            if e.subject_id == subject_id and e.action is action:
                return e
        return None

    def verify(self) -> tuple[bool, str]:
        """Revalida a cadeia inteira. Chame antes de emitir qualquer relatório."""
        prev = ""
        for e in self._events:
            if e.prev_hash != prev:
                return False, f"cadeia quebrada em seq={e.seq} (prev_hash divergente)"
            if content_hash(e.payload) != e.payload_hash:
                return False, f"payload adulterado em seq={e.seq}"
            recomputed = Event(**{**asdict(e), "actor": e.actor, "action": e.action, "this_hash": ""}).digest()
            if recomputed != e.this_hash:
                return False, f"hash de evento inválido em seq={e.seq}"
            prev = e.this_hash
        return True, "cadeia íntegra"

    def export(self) -> list[dict[str, Any]]:
        """Exporta em forma serializável (aproxima W3C PROV; mapeie
        actor->prov:Agent, event->prov:Activity, payload->prov:Entity)."""
        out = []
        for e in self._events:
            d = asdict(e)
            d["actor"] = e.actor.fingerprint()
            d["action"] = e.action.value
            out.append(d)
        return out


@dataclass
class PreRegistration:
    """Selo anti-HARKing.

    A predição quantitativa é hasheada e carimbada ANTES de qualquer acesso
    ao conjunto de confirmação. Depois, `matches()` prova que a predição
    avaliada é bit-a-bit a mesma que foi registrada.
    """

    hypothesis_id: str
    commit_hash: str
    committed_at: str
    ledger_seq: int
    payload_preview: dict[str, Any] = field(default_factory=dict)

    def matches(self, payload: Any) -> bool:
        return content_hash(payload) == self.commit_hash
