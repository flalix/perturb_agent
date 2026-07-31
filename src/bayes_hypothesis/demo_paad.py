"""
Demonstração end-to-end com o contexto PAAD.

Mostra os guards fazendo o trabalho: bloqueando HARKing, bloqueando
circularidade, e recusando promover associação a causalidade sem intervenção.
"""

from __future__ import annotations

import json

from sciagent_core import (
    Actor, Assay, Claim, Cost, Direction, EntityType, EvidenceKind,
    EvidenceScheduler, GuardError, Hypothesis, Ledger, Prediction,
    RefutationAttempt, Scope, SeedResolver, State, Triple,
    expected_information_gain, unidentifiable_report,
)
from sciagent_core.adapters.paad import DEGRow, deg_to_evidence, literature_evidence

R = SeedResolver()
LEDGER = Ledger()

GEN = Actor("gen-01", "generator", "model-A", "2026.07")
CRIT_A = Actor("crit-01", "critic", "model-B", "2026.07")     # família diferente
CRIT_B = Actor("crit-02", "critic", "model-C", "2026.07")
STAT = Actor("stat-01", "statistician", "deterministic", "dowhy-0.12")
HUMAN = Actor("flalix", "human", "n/a", "n/a")

BAR = "=" * 78


def hline(title: str) -> None:
    print(f"\n{BAR}\n{title}\n{BAR}")


# ---------------------------------------------------------------- entidades
hoxb_as3 = R.resolve("HOXB-AS3", EntityType.LNCRNA)
prf1 = R.resolve("PRF1", EntityType.GENE)
paad = R.resolve("PAAD", EntityType.DISEASE)
prog1 = R.resolve("Program-1", EntityType.PROGRAM)
cd8 = R.resolve("CD8+ T cell", EntityType.CELL_TYPE)
pancreas = R.resolve("pancreas", EntityType.TISSUE)

scope_p1 = Scope(tissue=pancreas, disease=paad, program=prog1,
                 notes="coorte local + TCGA-PAAD")

# ---------------------------------------------------------------- hipótese 1
H1 = Hypothesis(
    hid="H-001",
    statement_nl=(
        "A superexpressão de HOXB-AS3 em PAAD Program-1 causa redução da "
        "infiltração citotóxica CD8+ (fenótipo immune-cold), e não apenas "
        "acompanha a depleção de parênquima."
    ),
    formal=Triple(hoxb_as3, "negatively_regulates_infiltration_of", cd8, Direction.DOWN),
    scope=scope_p1,
    claim=Claim.CAUSAL,
    log_prior_odds=-2.2,          # prior ~10%: plausível, não estabelecida
    mechanism_descriptor="lncRNA_regulatory",
    scale_descriptor="cellular",
    cost_to_test=Cost(money_usd=48_000, days=90),
)

hline("1. Guard G2 — anti-HARKing")
premature = deg_to_evidence(
    DEGRow("HOXB-AS3", EntityType.LNCRNA, "Program-1_vs_Program-2",
           log2fc=1.9, padj=1e-6, n_group_a=41, n_group_b=57),
    R, Direction.UP, confirmatory=True,
)
try:
    H1.attach_evidence(premature, LEDGER, GEN)
except GuardError as exc:
    print(f"BLOQUEADO -> {exc}")

hline("2. Pré-registro (selo criptográfico das predições)")
H1.predictions = [
    Prediction("P1", "log2FC de HOXB-AS3 em P1 vs P2", "log2FC(HOXB-AS3)", ">=", 1.0,
               "RNA-seq bulk, deconvoluído"),
    Prediction("P2", "Densidade de CD8+ por mm2 menor em P1", "CD8_density_ratio_P1_P2",
               "<=", 0.6, "IHC multiplex / deconvolução"),
    Prediction("P3", "Knockdown de HOXB-AS3 aumenta quimiotaxia de CD8+ >=30%",
               "delta_migration_pct", ">=", 30.0, "transwell + CRISPRi", critical=True),
]
prereg = H1.register(LEDGER, GEN)
print(f"estado={H1.state.value}  commit={prereg.commit_hash[:16]}...  seq={prereg.ledger_seq}")
print(f"verificação do commit: {prereg.matches({'tampered': True})} (deve ser False)")

hline("3. Evidência L1 — literatura (com penalidade de campo ncRNA)")
for st, doi, lr, cls, retr in [
    ("HOX-AS lncRNAs associados a evasão imune em adenocarcinomas", "10.1000/a1", 6.0, "lit:hox-immune", False),
    ("Revisão citando o mesmo achado primário", "10.1000/a2", 5.0, "lit:hox-immune", False),
    ("Segunda revisão, mesma origem experimental", "10.1000/a3", 5.0, "lit:hox-immune", False),
    ("Artigo retratado sobre HOXB-AS3 e proliferação", "10.1000/a4", 12.0, "lit:hox-prolif", True),
]:
    ev = literature_evidence(st, doi, lr, cls, ncrna_field=True, retracted=retr)
    H1.attach_evidence(ev, LEDGER, GEN)
    print(f"  {ev}")

from sciagent_core.evidence import information_ceiling  # noqa: E402
print(f"\n  posterior após literatura: {H1.posterior:.4f}")
print("  DOIS efeitos operando ao mesmo tempo:")
print(f"   (i) TETO DE INFORMAÇÃO: com r=0.68 (literatura x penalidade ncRNA), nenhum")
print(f"       paper isolado pode entregar fator de Bayes > {information_ceiling(0.675):.2f}, por mais")
print(f"       espetacular que seja o p-valor. λ=6 e λ=1000 dão praticamente o mesmo.")
print("   (ii) DESCONTO DE DEPENDÊNCIA: 3 papers na MESMA classe de independência não")
print("        somam; o efeito de desenho os trata quase como um único evento.")
print(f"   O retratado entra com r=0.02 -> LR_eff = {literature_evidence('', 'x', 12.0, 'c', retracted=True).effective_lr:.3f} (≈1, inerte).")

hline("4. Evidência L1 — seus pipelines de PAAD (adaptador DEG)")
rows = [
    DEGRow("HOXB-AS3", EntityType.LNCRNA, "Program-1_vs_Program-2", 1.90, 1.0e-6, 41, 57),
    DEGRow("FAM83A-AS1", EntityType.LNCRNA, "Program-1_vs_Program-2", 1.42, 3.0e-5, 41, 57),
    DEGRow("PRF1", EntityType.GENE, "Program-1_vs_Program-2", -1.35, 8.0e-5, 41, 57),
]
for row in rows:
    ev = deg_to_evidence(row, R, Direction.UP if row.log2fc > 0 else Direction.DOWN,
                         confirmatory=True)
    H1.attach_evidence(ev, LEDGER, STAT)
    print(f"  {ev.statement}")
    print(f"      poder={ev.meta['power']:.2f}  LR_eff={ev.effective_lr:.2f}  r={ev.effective_reliability:.2f}")

print(f"\n  posterior após dados próprios: {H1.posterior:.4f}")
print(f"  Poder ~1.0 e padj=1e-6, mas r=0.39 (sem deconvolução, sem ajuste de pureza)")
print(f"  impõe teto de {information_ceiling(0.392):.2f}. É o resultado correto: se o confundimento")
print("  composicional não foi descartado, o dado NÃO consegue distinguir a hipótese")
print("  célula-intrínseca da alternativa composicional — por mais 'significativo' que seja.")

hline("5. Adjudicação — os guards recusam promover a hipótese")
state = H1.adjudicate(LEDGER, CRIT_A)
print(f"  estado = {state.value}   posterior = {H1.posterior:.4f}")
print("  motivos do bloqueio (últimos eventos do ledger):")
for e in LEDGER.for_subject("H-001")[-1:]:
    for r in e.payload.get("blocked_by", []):
        print(f"    - {r}")

hline("6. Scheduler bayesiano — qual a próxima evidência a adquirir?")
assays = [
    Assay("A-lit", "Mais revisão de literatura", EvidenceKind.LITERATURE,
          0.62, 0.55, cost_usd=200, days=1),
    Assay("A-decon", "Re-deconvolução com BayesPrism + ajuste de pureza", EvidenceKind.OWN_OMICS,
          0.80, 0.75, cost_usd=1_500, days=7),
    Assay("A-mr", "MR cis-eQTL + colocalização (GTEx/eQTLGen)", EvidenceKind.GENETIC_INSTRUMENT,
          0.72, 0.88, cost_usd=3_000, days=14, identifies_causal=True),
    Assay("A-depmap", "Reanálise DepMap/Perturb-seq (do() já executado)", EvidenceKind.PERTURBATION,
          0.78, 0.85, cost_usd=2_500, days=10, identifies_causal=True),
    Assay("A-crispri", "CRISPRi em organoide + co-cultura CD8+", EvidenceKind.WETLAB,
          0.90, 0.92, cost_usd=42_000, days=90, identifies_causal=True),
]
sched = EvidenceScheduler(budget_usd=12_000, diversity_reserve=0.25)
plans = sched.plan([H1], assays, LEDGER, STAT)
print(f"  {'ensaio':<50} {'bits':>7} {'USD':>8} {'bits/1k$':>9}  slot")
for pl in plans:
    print(f"  {pl.assay.label:<50} {pl.eig_bits:>7.3f} {pl.cost_usd:>8.0f} {pl.score:>9.3f}  {pl.slot}")
print(f"\n  gasto: ${sched.spent:,.0f} de ${sched.budget_usd:,.0f}")
print(f"  EIG de 'mais literatura' a p={H1.posterior:.2f}: "
      f"{expected_information_gain(H1.posterior, assays[0]):.4f} bits — quase nada, e é barato.")
print("  O scheduler recusa comprar mais literatura porque bits/USD é o critério, não custo.")

hline("7. Evidência intervencional + refutação adversarial independente")
from sciagent_core import EvidenceItem  # noqa: E402

H1.attach_evidence(
    EvidenceItem(
        eid="MR:HOXB-AS3-cis",
        kind=EvidenceKind.GENETIC_INSTRUMENT,
        statement="MR cis-eQTL: efeito consistente sobre escore de citotoxicidade; PP4=0.86 (coloc)",
        raw_lr=7.5, source_id="GTEx_v10+eQTLGen", independence_class="mr:hoxb-as3",
        confirmatory=True, meta={"PP4": 0.86, "F_stat": 41},
    ), LEDGER, STAT,
)
H1.attach_evidence(
    EvidenceItem(
        eid="PERT:DepMap-CRISPR",
        kind=EvidenceKind.PERTURBATION,
        statement="Knockout em linhagens PDAC: aumento de quimiocinas de recrutamento CD8+",
        raw_lr=5.0, source_id="DepMap_24Q4", independence_class="pert:depmap",
        confirmatory=True,
    ), LEDGER, STAT,
)

for rid, actor, strat in [
    ("R-01", CRIT_A, "confounder_probe"),
    ("R-02", CRIT_B, "negative_control"),
    ("R-03", CRIT_A, "placebo_refuter"),
]:
    H1.record_refutation(
        RefutationAttempt(rid, actor, strat, independent_of="gen-01", succeeded=False,
                          notes="não conseguiu derrubar"), LEDGER,
    )

state = H1.adjudicate(LEDGER, CRIT_A)
print(f"  posterior = {H1.posterior:.4f}   estado = {state.value}")
print(f"  não-literária={H1.has_non_literature_evidence()}  "
      f"intervencional={H1.has_interventional_evidence()}  "
      f"refutações independentes={H1.independent_refutation_count()}")
print("  SURVIVING não significa 'verdadeira' — significa 'não refutada sob pressão'.")

hline("8. Hipótese 2 — refutada por falha de predição crítica pré-registrada")
H2 = Hypothesis(
    hid="H-002",
    statement_nl="A queda de PRF1 em Program-1 é célula-intrínseca, não composicional.",
    formal=Triple(prf1, "downregulated_in", prog1, Direction.DOWN),
    scope=scope_p1, claim=Claim.CAUSAL, log_prior_odds=-1.4,
    mechanism_descriptor="effector_program", scale_descriptor="cellular",
)
H2.predictions = [
    Prediction("Q1", "Queda persiste após deconvolução por BayesPrism",
               "log2FC_PRF1_deconvoluted", "<=", -0.8, "BayesPrism", critical=True),
]
H2.register(LEDGER, GEN)
observed = -0.12
passed = H2.predictions[0].evaluate(observed)
print(f"  predição crítica Q1: log2FC deconvoluído <= -0.8 | observado = {observed:+.2f} -> {'PASSOU' if passed else 'FALHOU'}")
H2.attach_evidence(
    EvidenceItem(
        eid="DECON:PRF1", kind=EvidenceKind.OWN_OMICS,
        statement="Após deconvolução, o efeito desaparece: era composição de linfócitos.",
        raw_lr=0.12, source_id="PAAD_local#deconv", independence_class="own:deconv",
        confirmatory=True, decisive=not passed,
    ), LEDGER, STAT,
)
print(f"  estado = {H2.adjudicate(LEDGER, CRIT_B).value}   posterior = {H2.posterior:.4f}")
print("  Este é o caso que o pré-registro salva: sem o commit prévio, o sistema")
print("  reescreveria a hipótese para caber no resultado.")

hline("9. Saída honesta — não identificável")
H3 = Hypothesis(
    hid="H-003",
    statement_nl="MIR7-3HG medeia crosstalk estroma-tumor no Program-1 in vivo humano.",
    formal=Triple(R.resolve("MIR7-3HG", EntityType.LNCRNA), "mediates_crosstalk_with", cd8, Direction.UNSPECIFIED),
    scope=scope_p1, claim=Claim.CAUSAL, log_prior_odds=-3.0,
    mechanism_descriptor="stromal_crosstalk", scale_descriptor="tissue",
)
H3.predictions = [Prediction("S1", "efeito mediado NIE > 0.2", "NIE", ">=", 0.2, "mediação")]
H3.register(LEDGER, GEN)
blocked = [Assay("A-x", "único ensaio disponível", EvidenceKind.SIMULATION, 0.5, 0.5,
                 cost_usd=100, days=1, feasible=False, blocker="sem coorte com medida temporal pareada")]
msg = unidentifiable_report(H3, blocked)
print(f"  {msg}")
H3.park(LEDGER, HUMAN, "não identificável", "coorte longitudinal pareada pré/pós-tratamento")
print(f"  estado = {H3.state.value}")

hline("10. Auditoria e relatório")
ok, why = LEDGER.verify()
print(f"  integridade do ledger: {ok} ({why})   eventos = {len(LEDGER)}")
print("\n  Relatório final (gerado por consulta ao ledger, não escrito por LLM):")
print(json.dumps([H1.summary(), H2.summary(), H3.summary()], indent=2, ensure_ascii=False))

print("\n  Decomposição de evidência de H-001 (log-LR por classe de independência):")
for cls, contrib in sorted(H1.evidence_breakdown().items(), key=lambda kv: -abs(kv[1])):
    print(f"    {cls:<38} {contrib:+.3f}")
