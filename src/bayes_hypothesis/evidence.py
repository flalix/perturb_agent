"""
L1 — Camada de evidência.

Duas correções que quase nenhum sistema agêntico implementa e que são a
diferença entre um posterior calibrado e um número inventado:

(1) ATENUAÇÃO POR CONFIABILIDADE DA FONTE — com TETO DE INFORMAÇÃO.

    Seja r = P(fonte genuína). Misture as verossimilhanças, não as LRs:

        P(E|H)  = r·λ + (1-r)·q
        P(E|¬H) = r·1 + (1-r)·q

    onde q é a chance de a fonte NÃO-genuína emitir esse mesmo relato. A
    escolha de q é a decisão metodológica central, e a resposta honesta é
    q = λ: fraude e erro produzem preferencialmente resultados marcantes —
    justamente o tipo de achado que um efeito real produziria. Então:

        LR_eff = λ / ( r + (1-r)·λ )

    Consequência (o resultado importante):

        lim(λ→∞) LR_eff = 1 / (1 - r)

    UMA FONTE COM TAXA DE FABRICAÇÃO (1-r) NUNCA PODE ENTREGAR UM FATOR DE
    BAYES MAIOR QUE 1/(1-r), por mais espetacular que seja o resultado.
    Com r = 0,90, o teto é 10 — nenhum paper isolado, com p = 1e-30, move o
    posterior mais que 10:1. Isto captura formalmente por que "altamente
    significativo" e "altamente informativo" não são a mesma coisa quando a
    fonte pode estar errada.

    Aplicamos o teto simetricamente (piso = 1-r) para evidência desconfirmatória.

(2) DESCONTO POR DEPENDÊNCIA.
    Somar log-LR de 30 papers que citam o mesmo experimento original é a
    principal via de superconfiança. Agrupamos evidência em classes de
    independência e aplicamos o efeito de desenho (design effect) de
    amostragem complexa:

        n_eff = n / (1 + ρ(n - 1))

    e escalamos a contribuição do grupo por n_eff/n. Com ρ = 1 (totalmente
    redundante) o grupo inteiro vale por um único item.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from .ontology import Entity


class EvidenceKind(str, Enum):
    """A distinção operacional é: isto é dependente da literatura?"""

    LITERATURE = "literature"
    EXPERT_PRIOR = "expert_prior"
    KG_STRUCTURE = "kg_structure"          # derivado de literatura -> dependente
    OWN_OMICS = "own_omics"                # observacional próprio
    PUBLIC_OMICS = "public_omics"          # observacional público
    PERTURBATION = "perturbation"          # do() já executado (CRISPR, LINCS, Perturb-seq)
    GENETIC_INSTRUMENT = "genetic_instrument"   # MR / coloc
    CLINICAL = "clinical"
    SIMULATION = "simulation"              # ODE, booleano, FBA, modelo preditivo
    WETLAB = "wetlab"


#: Tipos que NÃO dependem da literatura publicada. Uma hipótese não pode
#: alcançar o estado SURVIVING sem pelo menos um destes — é a trava
#: anti-circularidade (literatura -> KG -> hipótese -> "confirmada" pela literatura).
NON_LITERATURE_KINDS: frozenset[EvidenceKind] = frozenset(
    {
        EvidenceKind.OWN_OMICS,
        EvidenceKind.PUBLIC_OMICS,
        EvidenceKind.PERTURBATION,
        EvidenceKind.GENETIC_INSTRUMENT,
        EvidenceKind.CLINICAL,
        EvidenceKind.WETLAB,
    }
)

#: Tipos que constituem intervenção (real ou quase-experimental). Necessários
#: para promover uma afirmação de associação a afirmação causal.
INTERVENTIONAL_KINDS: frozenset[EvidenceKind] = frozenset(
    {EvidenceKind.PERTURBATION, EvidenceKind.GENETIC_INSTRUMENT, EvidenceKind.WETLAB}
)


# --- confiabilidade de base por tipo de fonte --------------------------------
# Ancoragem para LITERATURE: estimativas publicadas de artigos biomédicos
# efetivamente falsos ficam na casa de ~6%; some erro honesto não retratado e
# 0,90 é generoso, não pessimista. Ajuste por campo: literatura de ncRNA
# está entre as mais contaminadas — use `field_penalty`.
DEFAULT_RELIABILITY: dict[EvidenceKind, float] = {
    EvidenceKind.LITERATURE: 0.90,
    EvidenceKind.EXPERT_PRIOR: 0.75,
    EvidenceKind.KG_STRUCTURE: 0.80,
    EvidenceKind.OWN_OMICS: 0.95,
    EvidenceKind.PUBLIC_OMICS: 0.90,
    EvidenceKind.PERTURBATION: 0.93,
    EvidenceKind.GENETIC_INSTRUMENT: 0.92,
    EvidenceKind.CLINICAL: 0.90,
    EvidenceKind.SIMULATION: 0.60,
    EvidenceKind.WETLAB: 0.95,
}


def information_ceiling(reliability: float) -> float:
    """Fator de Bayes máximo que uma fonte com confiabilidade r pode entregar."""
    r = min(max(reliability, 0.0), 1.0 - 1e-9)
    return 1.0 / (1.0 - r)


def attenuate_lr(lr: float, reliability: float) -> float:
    """LR_eff = λ / (r + (1-r)·λ), aplicada simetricamente. Ver docstring do módulo.

    Teto = 1/(1-r); piso = 1-r.
    """
    if lr <= 0:
        raise ValueError("LR deve ser positiva (use decisive=True para refutação lógica)")
    r = min(max(reliability, 0.0), 1.0 - 1e-9)
    if lr >= 1.0:
        return lr / (r + (1.0 - r) * lr)
    inv = 1.0 / lr
    return 1.0 / (inv / (r + (1.0 - r) * inv))


def lr_from_test(power: float, alpha_eff: float, detected: bool) -> float:
    """Converte um resultado de teste estatístico em razão de verossimilhança.

    Detectado:      LR = poder / alfa_efetivo
    Não detectado:  LR = (1 - poder) / (1 - alfa_efetivo)

    Obriga o pipeline a declarar poder. Um DEG "significativo" com poder
    desconhecido não tem LR definida — e isso é uma feature, não um bug.
    """
    power = min(max(power, 1e-6), 1 - 1e-6)
    alpha_eff = min(max(alpha_eff, 1e-6), 1 - 1e-6)
    return power / alpha_eff if detected else (1 - power) / (1 - alpha_eff)


@dataclass
class EvidenceItem:
    """Um item de evidência com LR explícita, proveniência e classe de dependência."""

    eid: str
    kind: EvidenceKind
    statement: str
    raw_lr: float                             # λ antes de atenuar
    source_id: str                            # DOI, accession, caminho do dataset
    independence_class: str                   # itens com mesma classe são tratados como redundantes
    reliability: float | None = None          # None -> usa DEFAULT_RELIABILITY
    field_penalty: float = 1.0                # multiplicador (<1) p/ campos contaminados
    confirmatory: bool = False                # True exige pré-registro prévio
    decisive: bool = False                    # falha de predição pré-registrada e crítica
    entities: tuple[Entity, ...] = ()
    meta: dict = field(default_factory=dict)

    @property
    def effective_reliability(self) -> float:
        base = self.reliability if self.reliability is not None else DEFAULT_RELIABILITY[self.kind]
        return min(max(base * self.field_penalty, 0.0), 1.0)

    @property
    def effective_lr(self) -> float:
        return attenuate_lr(self.raw_lr, self.effective_reliability)

    @property
    def log_lr(self) -> float:
        return math.log(self.effective_lr)

    def __str__(self) -> str:
        return (
            f"[{self.kind.value}] {self.statement} "
            f"(λ={self.raw_lr:.2f} → LR_eff={self.effective_lr:.2f}, r={self.effective_reliability:.2f})"
        )


@dataclass
class DependenceModel:
    """Correlação intra-classe por tipo de evidência.

    rho=1.0 → itens da mesma classe são perfeitamente redundantes.
    Literatura tem rho alto por padrão: papers do mesmo grupo/mesma origem
    experimental são quase o mesmo evento informacional.
    """

    rho: dict[EvidenceKind, float] = field(
        default_factory=lambda: {
            EvidenceKind.LITERATURE: 0.80,
            EvidenceKind.KG_STRUCTURE: 0.90,
            EvidenceKind.EXPERT_PRIOR: 0.70,
            EvidenceKind.OWN_OMICS: 0.60,
            EvidenceKind.PUBLIC_OMICS: 0.50,
            EvidenceKind.PERTURBATION: 0.30,
            EvidenceKind.GENETIC_INSTRUMENT: 0.20,
            EvidenceKind.CLINICAL: 0.40,
            EvidenceKind.SIMULATION: 0.85,
            EvidenceKind.WETLAB: 0.35,
        }
    )

    def deff_scale(self, kind: EvidenceKind, n: int) -> float:
        """Fator n_eff/n para um grupo de n itens redundantes."""
        if n <= 1:
            return 1.0
        r = self.rho.get(kind, 0.5)
        n_eff = n / (1.0 + r * (n - 1))
        return n_eff / n


def aggregate_log_lr(
    items: list[EvidenceItem], dependence: DependenceModel | None = None
) -> tuple[float, dict[str, float]]:
    """Soma log-LR com desconto de dependência por classe de independência.

    Retorna (log_lr_total, contribuição_por_classe).
    """
    dep = dependence or DependenceModel()
    groups: dict[tuple[EvidenceKind, str], list[EvidenceItem]] = {}
    for it in items:
        groups.setdefault((it.kind, it.independence_class), []).append(it)

    total = 0.0
    breakdown: dict[str, float] = {}
    for (kind, cls), grp in groups.items():
        scale = dep.deff_scale(kind, len(grp))
        contrib = sum(it.log_lr for it in grp) * scale
        total += contrib
        breakdown[f"{kind.value}:{cls}"] = contrib
    return total, breakdown
