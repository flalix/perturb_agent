"""
L0 — Camada de identidade e ontologia.

Regra do sistema: NADA entra nas camadas superiores como string livre.
Toda entidade vira um CURIE normalizado. A falha silenciosa mais cara
em pipelines biomédicos é aliasing de símbolo gênico (e.g. MARCH1 -> MARCHF1,
SEPT2 -> SEPTIN2, e a corrupção clássica por Excel).

Em produção, substitua `_SEED_ALIASES` por resolução via BioCypher/BioLink,
HGNC, MONDO, EFO, ChEBI, Uberon, CL. A interface abaixo foi desenhada para
que essa troca seja um drop-in (implemente `Resolver`).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Protocol


class EntityType(str, Enum):
    GENE = "gene"
    TRANSCRIPT = "transcript"
    PROTEIN = "protein"
    LNCRNA = "lncrna"
    PATHWAY = "pathway"
    CELL_TYPE = "cell_type"
    TISSUE = "tissue"
    DISEASE = "disease"
    PHENOTYPE = "phenotype"
    DRUG = "drug"
    PROGRAM = "program"  # p.ex. um cluster/programa transcricional definido pelo projeto


CURIE_RE = re.compile(r"^[A-Za-z0-9._]+:[A-Za-z0-9._\-]+$")


@dataclass(frozen=True, order=True)
class Entity:
    """Entidade normalizada. Imutável e hashable — usada como chave em grafos."""

    curie: str
    etype: EntityType
    label: str = ""

    def __post_init__(self) -> None:
        if not CURIE_RE.match(self.curie):
            raise ValueError(
                f"CURIE malformado: {self.curie!r}. "
                "Esperado prefixo:identificador (e.g. HGNC:4617, MONDO:0009831)."
            )

    def __str__(self) -> str:
        return f"{self.label or self.curie} [{self.curie}]"


class Resolver(Protocol):
    """Contrato de resolução. Implemente sobre BioCypher/HGNC/MONDO em produção."""

    def resolve(self, raw: str, etype: EntityType) -> Entity: ...


# Semente mínima para o demo. Em produção isto vem do KG.
_SEED_ALIASES: dict[tuple[str, EntityType], tuple[str, str]] = {
    ("FAM83A-AS1", EntityType.LNCRNA): ("HGNC:40876", "FAM83A-AS1"),
    ("HOXA10-AS", EntityType.LNCRNA): ("HGNC:53933", "HOXA10-AS"),
    ("HOXB-AS3", EntityType.LNCRNA): ("HGNC:48619", "HOXB-AS3"),
    ("MIR7-3HG", EntityType.LNCRNA): ("HGNC:33584", "MIR7-3HG"),
    ("NEAT1", EntityType.LNCRNA): ("HGNC:30815", "NEAT1"),
    ("H19", EntityType.LNCRNA): ("HGNC:4713", "H19"),
    ("GATA6", EntityType.GENE): ("HGNC:4174", "GATA6"),
    ("KRT17", EntityType.GENE): ("HGNC:6427", "KRT17"),
    ("PRF1", EntityType.GENE): ("HGNC:9360", "PRF1"),
    ("GZMB", EntityType.GENE): ("HGNC:4709", "GZMB"),
    ("PAAD", EntityType.DISEASE): ("MONDO:0009831", "pancreatic adenocarcinoma"),
    ("CD8+ T cell", EntityType.CELL_TYPE): ("CL:0000625", "CD8-positive T cell"),
    ("pancreas", EntityType.TISSUE): ("UBERON:0001264", "pancreas"),
    ("Program-1", EntityType.PROGRAM): ("PROJ:PAAD-P1", "PAAD Program-1 (lncRNA/HOX-AS, immune-cold)"),
    ("Program-2", EntityType.PROGRAM): ("PROJ:PAAD-P2", "PAAD Program-2 (desmoplastic, immune-hot)"),
}


@dataclass
class SeedResolver:
    """Resolver de demonstração. Falha alto e cedo em termos desconhecidos —
    de propósito: entidade não resolvida NUNCA deve virar evidência silenciosa."""

    strict: bool = True
    _unresolved: list[str] = field(default_factory=list)

    def resolve(self, raw: str, etype: EntityType) -> Entity:
        key = (raw.strip(), etype)
        if key in _SEED_ALIASES:
            curie, label = _SEED_ALIASES[key]
            return Entity(curie=curie, etype=etype, label=label)
        self._unresolved.append(raw)
        if self.strict:
            raise KeyError(
                f"Entidade não resolvida: {raw!r} ({etype}). "
                "Registre no KG antes de usar como evidência."
            )
        # modo permissivo: CURIE local, marcado como não-canônico
        slug = re.sub(r"[^A-Za-z0-9._\-]", "_", raw)
        return Entity(curie=f"LOCAL:{slug}", etype=etype, label=raw)


class Direction(str, Enum):
    UP = "up"
    DOWN = "down"
    UNCHANGED = "unchanged"
    UNSPECIFIED = "unspecified"


@dataclass(frozen=True)
class Triple:
    """Forma canônica e falsificável de uma afirmação.

    Sem direção e sinal, uma aresta 'gene A associado a doença B' não gera
    predição testável. Por isso `direction` é obrigatória no construtor.
    """

    subject: Entity
    predicate: str  # idealmente de um vocabulário fechado (BioLink)
    object: Entity
    direction: Direction = Direction.UNSPECIFIED

    def __str__(self) -> str:
        arrow = {"up": "↑", "down": "↓", "unchanged": "=", "unspecified": "?"}[self.direction]
        return f"{self.subject.label} --{self.predicate}({arrow})--> {self.object.label}"


@dataclass(frozen=True)
class Scope:
    """Contexto sem o qual a hipótese não é testável.

    'X regula Y' é vago; 'X regula Y em células ductais pancreáticas humanas
    no contexto de PAAD Program-1' é testável.
    """

    organism: str = "Homo sapiens"
    tissue: Entity | None = None
    cell_type: Entity | None = None
    disease: Entity | None = None
    program: Entity | None = None
    notes: str = ""

    def key(self) -> str:
        parts = [self.organism]
        for e in (self.tissue, self.cell_type, self.disease, self.program):
            parts.append(e.curie if e else "-")
        return "|".join(parts)
