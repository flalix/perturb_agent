# sciagent_core

Núcleo mínimo e auditável para um sistema multiagente científico. Stdlib apenas.

```
sciagent_core/
├── ontology.py        L0  Entity, Triple, Scope, Resolver (troque por BioCypher)
├── evidence.py        L1  EvidenceItem, LR atenuada, desconto de dependência
├── hypothesis.py      L2  Hypothesis + máquina de estados (guards G1..G5)
├── scheduler.py       L3  EIG bayesiano / custo + reserva de diversidade
├── provenance.py      L6  Ledger append-only com encadeamento de hash
└── adapters/paad.py   L1  plugue dos seus pipelines de PAAD

demo_paad.py           demonstração end-to-end
test_guards.py         14 testes que asseguram que os guards não afrouxam
```

```bash
python3 demo_paad.py
python3 -m unittest test_guards -v
```

## As três ideias que carregam o peso

**1. Teto de informação por confiabilidade da fonte.**
`LR_eff = λ / (r + (1−r)·λ)`, com `lim(λ→∞) = 1/(1−r)`.
Uma fonte com 10% de chance de estar errada nunca entrega fator de Bayes acima
de 10, por mais espetacular que seja o p-valor. Isto formaliza a diferença entre
*altamente significativo* e *altamente informativo*, e resolve automaticamente o
problema de literatura contaminada — sem regras ad hoc.

**2. Desconto de dependência.**
`n_eff = n / (1 + ρ(n−1))`. Trinta papers que citam o mesmo experimento original
valem por aproximadamente um. É a principal via de superconfiança em sistemas
agênticos, e é eliminada por construção.

**3. Guards como código, não como prompt.**
| Guard | Impede |
|---|---|
| G1 | hipótese sem predição quantitativa falsificável |
| G2 | evidência confirmatória sem commit criptográfico prévio → **HARKing** |
| G3 | promoção com evidência só de literatura/KG → **circularidade** |
| G4 | afirmação causal sem evidência intervencional → **associação ≠ causa** |
| G5 | promoção sem refutações independentes → **crítico decorativo** |

Um LLM pode ignorar uma instrução de prompt. Não pode ignorar um `GuardError`.

## Integração com os pipelines reais

**Ontologia → BioCypher.** Implemente o `Protocol` `Resolver`:

```python
class BioCypherResolver:
    def __init__(self, driver): self.driver = driver
    def resolve(self, raw: str, etype: EntityType) -> Entity:
        rec = self.driver.query(
            "MATCH (n) WHERE $raw IN n.synonyms AND n.type = $t RETURN n.id, n.name",
            raw=raw, t=etype.value,
        )
        if not rec:
            raise KeyError(f"não resolvido: {raw}")
        return Entity(curie=rec[0]["n.id"], etype=etype, label=rec[0]["n.name"])
```

**DEG → evidência.** Substitua `approx_power` por estimativa real (simulação
sobre o modelo do DESeq2/limma) e alimente a partir do workbook:

```python
import pandas as pd
from sciagent_core.adapters.paad import DEGRow, deg_to_evidence

df = pd.read_excel("PAAD_clustering_analysis_corrected.xlsx", sheet_name="DEG_N3")
for _, r in df.iterrows():
    row = DEGRow(
        gene=r.gene, entity_type=EntityType.LNCRNA if r.biotype == "lncRNA" else EntityType.GENE,
        contrast="Program-1_vs_Program-2", log2fc=r.log2FoldChange, padj=r.padj,
        n_group_a=n1, n_group_b=n2,
        composition_adjusted=False,   # ← vire True depois do BayesPrism/CIBERSORTx
        purity_adjusted=False,
    )
    h.attach_evidence(deg_to_evidence(row, resolver, Direction.UP, confirmatory=True), ledger, stat)
```

O flag `composition_adjusted` não é cosmético: com ele em `False`, a
confiabilidade cai para ~0,39 e o teto de informação vira **1,64**. Traduzindo:
enquanto a deconvolução não rodar, um DEG com `padj=1e-30` não consegue
distinguir "efeito célula-intrínseca" de "mudança de composição celular" — e o
sistema recusa fingir que consegue. Rodar BayesPrism e virar o flag para `True`
é literalmente o que destrava o posterior.

**Ingest de literatura.** Antes de qualquer uso, cheque `retracted` contra
Retraction Watch/CrossRef. Para lncRNA/miRNA, mantenha `ncrna_field=True` — é
uma das literaturas com maior densidade de artigos problemáticos, e a penalidade
de campo evita que o sistema construa uma hipótese sobre areia.

## O que ainda falta (próximas peças)

1. **Motor de identificabilidade causal** — DAG → estimando → `dowhy` →
   busca de dataset compatível → refutadores automáticos. Emite
   `unidentifiable_report()` quando nada identifica.
2. **Buscador quality-diversity** (MAP-Elites sobre `mechanism_descriptor` ×
   `scale_descriptor`) para o gerador de hipóteses, contra colapso de diversidade.
3. **Harness de holdout temporal** — gerar com literatura até T, medir taxa de
   novidade genuína contra T+1..hoje. É o único jeito honesto de calibrar o
   sistema; não pergunte a novidade ao LLM.
4. **Curvas de calibração** — o posterior declarado prediz a taxa real de
   confirmação? Sem isso, os números são decorativos.

## Nota

`SURVIVING` nunca significa "verdadeira". Significa "não refutada sob pressão
adversarial, com evidência não-circular e ao menos uma intervenção". O relatório
final deve ser gerado por consulta ao ledger (`Ledger.export()`), nunca escrito
de memória por um modelo de linguagem.
