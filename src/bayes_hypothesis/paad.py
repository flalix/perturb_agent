"""
Adaptador L1 — plugue dos seus pipelines de PAAD como fonte de evidência.

Converte resultados de DEG por cluster/programa em `EvidenceItem` com LR
calculada, em vez de "gene significativo" solto.

Três decisões de desenho que carregam metodologia para dentro do código:

1. PODER É OBRIGATÓRIO. `lr_from_test` exige poder e alfa efetivo. Um DEG sem
   poder declarado não tem razão de verossimilhança definida. Se você não
   estimou poder, o adaptador levanta erro — de propósito.

2. AJUSTE COMPOSICIONAL É UM FLAG DE PRIMEIRA CLASSE. Se
   `composition_adjusted=False`, a confiabilidade da evidência é penalizada
   pesadamente, porque "expressão sobe no Program-1" pode ser inteiramente
   um artefato de proporção celular / pureza tumoral. É exatamente o caveat
   classical vs. basal, agora executável em vez de nota de rodapé.

3. CLASSE DE INDEPENDÊNCIA = DATASET + CONTRASTE. Cem genes do mesmo
   contraste no mesmo dataset NÃO são cem evidências independentes; o
   desconto por efeito de desenho cuida disso automaticamente.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..evidence import EvidenceItem, EvidenceKind, lr_from_test
from ..ontology import Direction, Entity, EntityType, Resolver


#: Penalidade de campo para literatura de RNA não-codificante, que está entre
#: as áreas com maior densidade de artigos problemáticos. Aplique só a
#: evidência de tipo LITERATURE/KG_STRUCTURE — não aos seus próprios dados.
NCRNA_FIELD_PENALTY = 0.75


@dataclass
class DEGRow:
    """Uma linha do seu output de DEG (o formato do workbook corrigido)."""

    gene: str
    entity_type: EntityType
    contrast: str                  # e.g. "Program-1_vs_Program-2"
    log2fc: float
    padj: float
    n_group_a: int
    n_group_b: int
    dataset_id: str = "PAAD_local"
    composition_adjusted: bool = False
    purity_adjusted: bool = False


def approx_power(log2fc: float, n_a: int, n_b: int, sd: float = 1.0, alpha: float = 0.05) -> float:
    """Poder aproximado para teste bilateral de duas amostras (aproximação normal).

    Deliberadamente simples e transparente. Substitua por estimativa de poder
    do seu framework (limma/DESeq2 + simulação) quando tiver.
    """
    if n_a < 2 or n_b < 2:
        return 0.05
    ncp = abs(log2fc) / (sd * math.sqrt(1.0 / n_a + 1.0 / n_b))
    z_alpha = 1.959963985  # bilateral 5%
    # Φ(ncp - z) via aproximação com erf
    z = ncp - z_alpha
    return min(max(0.5 * (1.0 + math.erf(z / math.sqrt(2.0))), 0.01), 0.999)


def deg_to_evidence(
    row: DEGRow,
    resolver: Resolver,
    expected_direction: Direction,
    fdr_target: float = 0.05,
    confirmatory: bool = False,
    sd: float = 1.0,
) -> EvidenceItem:
    """Converte uma linha de DEG em evidência com LR explícita."""
    entity: Entity = resolver.resolve(row.gene, row.entity_type)

    observed_dir = Direction.UP if row.log2fc > 0 else Direction.DOWN
    detected = (row.padj < fdr_target) and (observed_dir is expected_direction)

    power = approx_power(row.log2fc, row.n_group_a, row.n_group_b, sd=sd)
    # alfa efetivo: FDR alvo dividido por 2 (direção precisa bater)
    lr = lr_from_test(power=power, alpha_eff=fdr_target / 2.0, detected=detected)

    # Penalidade por não-ajuste — o caveat de deconvolução, executável.
    reliability = 0.95
    caveats: list[str] = []
    if not row.composition_adjusted:
        reliability *= 0.55
        caveats.append("sem ajuste de composição celular (deconvolução)")
    if not row.purity_adjusted:
        reliability *= 0.75
        caveats.append("sem ajuste de pureza tumoral")

    return EvidenceItem(
        eid=f"DEG:{row.dataset_id}:{row.contrast}:{row.gene}",
        kind=EvidenceKind.OWN_OMICS,
        statement=(
            f"{row.gene} {observed_dir.value} em {row.contrast} "
            f"(log2FC={row.log2fc:+.2f}, padj={row.padj:.2e})"
            + (f" — CAVEATS: {'; '.join(caveats)}" if caveats else "")
        ),
        raw_lr=lr,
        source_id=f"{row.dataset_id}#{row.contrast}",
        # todos os genes do mesmo contraste/dataset compartilham a classe:
        # não são evidências independentes entre si.
        independence_class=f"{row.dataset_id}:{row.contrast}",
        reliability=reliability,
        confirmatory=confirmatory,
        entities=(entity,),
        meta={
            "log2fc": row.log2fc,
            "padj": row.padj,
            "power": round(power, 3),
            "composition_adjusted": row.composition_adjusted,
            "purity_adjusted": row.purity_adjusted,
            "caveats": caveats,
        },
    )


def literature_evidence(
    statement: str,
    doi: str,
    raw_lr: float,
    independence_class: str,
    ncrna_field: bool = False,
    retracted: bool = False,
) -> EvidenceItem:
    """Evidência de literatura com penalidade de campo e checagem de retratação.

    No ingest real, cheque `retracted` contra Retraction Watch / CrossRef antes
    de qualquer uso — um artigo retratado deve ir a reliability ≈ 0, não ser
    silenciosamente ponderado.
    """
    penalty = NCRNA_FIELD_PENALTY if ncrna_field else 1.0
    return EvidenceItem(
        eid=f"LIT:{doi}",
        kind=EvidenceKind.LITERATURE,
        statement=statement,
        raw_lr=raw_lr,
        source_id=doi,
        independence_class=independence_class,
        reliability=0.02 if retracted else None,
        field_penalty=penalty,
        meta={"retracted": retracted, "ncrna_field": ncrna_field},
    )
