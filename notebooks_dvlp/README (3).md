# 🧪 perturb_agent

**Perturb Agent** is a computational framework to identify **patient-specific pathway perturbations**
from TCGA/cBioPortal transcriptomic data and map them to **potential therapeutic targets**.

Current case study: **pancreatic adenocarcinoma (PAAD)**.

**Status**: 🚧 Under development

**Live demo**: [perturb-agent.onrender.com](https://perturb-agent.onrender.com/)

---

## 📑 Table of contents

- [1. Pipeline overview](#1-pipeline-overview)
- [2. Data sources](#2-data-sources)
- [3. Analysis history and design decisions](#3-analysis-history-and-design-decisions)
  - [3.1 Mutation clustering — ❌ negative result](#31-mutation-clustering---negative-result)
  - [3.2 Bulk DE + ComBat + clustering — ❌ negative result](#32-bulk-de--combat--clustering---negative-result)
  - [3.3 The lncRNA artifact — ❌ root cause found](#33-the-lncrna-artifact---root-cause-found)
  - [3.4 Decision: stop transforming, start deconvolving](#34-decision-stop-transforming-start-deconvolving)
- [4. Cell-type deconvolution with BayesPrism](#4-cell-type-deconvolution-with-bayesprism)
- [5. Program-level analysis: challenges and statistical power](#5-program-level-analysis-challenges-and-statistical-power)
- [6. Compartment-specific factorial DE](#6-compartment-specific-factorial-de-methods-check)
- [7. Results: survival (Kaplan–Meier)](#7-results-survival-kaplanmeier)
- [8. Downstream modules](#8-downstream-modules)
- [9. Roadmap](#9-roadmap)
- [10. Reproducibility](#10-reproducibility)
- [11. References](#11-references)
- [12. Contact](#12-contact)

---

## 1. Pipeline overview

| Layer | Tool | Purpose |
|---|---|---|
| UI | [Streamlit](https://streamlit.io/) | Interactive exploration and visualization |
| Containers | [Docker](https://docker-curriculum.com/) | Reproducibility and R/Bioconductor interface |
| Orchestration | [Nextflow](https://www.nextflow.io/) | Scalable, reproducible data processing |
| ML/AI | Python | Pathway scoring, feature attribution, target prioritization |
| Env | [uv](https://docs.astral.sh/uv/) | Project and dependency management |
| Lint/format | [ruff](https://docs.astral.sh/ruff/) | Code linting and formatting |

Pinned UI dependencies:

```
streamlit==1.55.0
protobuf==3.20.3
click==8.0.4
```

Linting:

```bash
uv run ruff check src/libs/*.py --fix > ruff.txt
uv run ruff format src/libs/*.py
```

---

## 2. Data sources

### 2.1 GDC / TCGA

First data layer. Query any [GDC/TCGA](https://portal.gdc.cancer.gov/analysis_page?app=Projects)
cancer type through the chatbot orchestration layer.

Query flow:

```
program      → project_id            gdc.list_gdc_programs()
primary_site → pid, disease_type     gdc.get_primary_sites(program=program)
cases        → case_id (UUID)        gdc.build_cases(pid=pid, subtype=subtype, stage=stage)
                 ├── subtypes → cancer / tissue subtypes
                 └── stages   → stage_id (AJCC)
samples      → sample_type [tumor | normal], file access
barcodes     → patients
```

Coverage retrieved so far:

| Entity | Count |
|---|---:|
| Primary sites | 57 |
| Cases | 11,428 |
| Samples | 245,657 |
| Annotated mutations | 480,826 |
| Distinct genes | 18,961 |

> ⚠️ **To verify**: the sample/case ratio (~21) is high. Confirm that `samples` counts
> biospecimens and not files.

### 2.2 cBioPortal

GDC restricts controlled-access mutation data. [cBioPortal](https://www.cbioportal.org/)
provides anonymized patient demographics plus mutation annotation, so the two sources are
used together.

Demographics retrieval:

```python
cbio.check_clinical_attributes()
cbio.loop_program_psi_get_cases_samples_mut()
```

**Notebooks**
- `tcga_gdc_and_cBioPortal_mutations_loop.ipynb`
- `mtp_cBio_60_PSI_all_programs_CasesSamplesMut_Demographics`

---

## 3. Analysis history and design decisions

This section documents what **did not** work and why. The negative results drove the current
design, so they are kept on purpose.

### 3.1 Mutation clustering — ❌ negative result

**Approach**

1. Select a disease.
2. Retrieve all annotated mutations.
3. Build a binary pivot table: `barcodes × gene symbols`.
4. Cluster:
   - Pairwise Jaccard distance — `pairwise_distances(X, metric="jaccard")`
   - Hierarchical clustering + dendrogram (seaborn)
   - UMAP embedding
   - Group assignment via **kNN**, k = 8

**Jaccard distance**

Jaccard similarity measures how much two sets overlap; Jaccard distance is its complement
and measures how different they are:

$$J(A,B) = \frac{|A \cap B|}{|A \cup B|} \qquad d_J(A,B) = 1 - J(A,B)$$

**Outcome — did not work.** No stable or reproducible patient structure emerged. The
mutation matrix is extremely sparse: most patients share few or no mutated genes beyond a
handful of drivers, so most pairwise Jaccard distances collapse toward 1 and the embedding
is dominated by noise. Clusters were unstable across seeds and k.

**Interpretation.** Most likely **insufficient data** — too few cases per subtype/stage
stratum, and per-patient mutation counts too low to support a similarity-based partition.
Mutation profiles alone are not a sufficient basis for patient stratification here.

Exploratory figures (disease = Esophagus):

![mutation frequency](./figures/most_mutated_genes_esophagus.png)
*Most mutated genes.*

![heatmap](./figures/esophagus_mutation_plot.png)
*Mutation heatmap.*

![UMAP](./figures/esophagus_UMAP_k=8.png)
*UMAP with kNN, k = 8.*

---

### 3.2 Bulk DE + ComBat + clustering — ❌ negative result

**Approach**

1. Retrieve all patient cases (barcodes).
2. Obtain gene expression (raw counts) per patient.
3. Compute DEGs per patient:
   - Control: TCGA solid-tissue normal samples
   - Method: [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
   - Thresholds: |log2FC| ≥ 1 and FDR < 0.05
4. Batch-correct with **ComBat**, then cluster patients.

**Notebooks**
- `mtp_cBio_63_PAAD_prepare_expression` — differential expression
- `mtp_cBio_62_PAAD_Combat_DE_Cluster` — ComBat + DE + clustering

**Outcome — information loss.** The ComBat + clustering combination removed biological
signal along with batch effect. Bulk expression mixes tumor, stromal and immune
compartments into one average, so per-patient DEGs cannot be attributed to a cell type;
batch correction on top of that averaged away what remained. The resulting clusters were
not interpretable.

---

### 3.3 The lncRNA artifact — ❌ root cause found

**Approach.** HCA + UMAP clustering on the corrected expression matrix. The partition
looked clean and split into apparently meaningful groups.

**Outcome — the clusters were an artifact.** The leading cluster was driven almost entirely
by **lncRNAs**, which is not biologically plausible for PAAD stratification. Tracing it back,
the cause was a **read-counting error: stranded vs. unstranded counts were mixed**, which
inflates antisense and lncRNA signal and creates a spurious dominant axis of variation.

**Lesson.** A visually clean UMAP is not evidence of biology. Every cluster now requires a
plausibility check on its top-loading feature class before it is interpreted.

**Notebook**: `mtp_cBio_64_PAAD_claude_critic_lncRNA_not_bio`

---

### 3.4 Decision: stop transforming, start deconvolving

The three failures above share one cause: **the transformations were doing the damage.**
Jaccard binarization discarded mutation context; ComBat removed real biology; the
count-strand mismatch survived every downstream normalization.

**New rule: do not transform the data.**

- Keep **raw counts**; no ComBat, no ad-hoc batch correction, no rescaling before modeling.
- Fix counting at the source (consistent stranded/unstranded protocol).
- Instead of correcting bulk data, **decompose it** — apply **BayesPrism** to recover
  cell-type composition and cell-type-specific expression, then analyze within cell types.

---

## 4. Cell-type deconvolution with BayesPrism

[BayesPrism](https://github.com/Danko-Lab/BayesPrism) performs Bayesian deconvolution of bulk
transcriptomes against a single-cell reference, returning per-sample **cell-type fractions**
and a per-sample, per-cell-type **expression tensor (Z)**.

- **Input**: raw bulk counts (untransformed)
- **Reference**: Peng *et al.* 2019 PDAC single-cell atlas
- **Output**: cell-type fractions θ + cell-type-specific expression Z
- **Notebook**: `mtp_cBio_65_PAAD_run_bayesprism`

This is the current backbone of the pipeline: everything downstream is computed **within a
cell type or program**, not on the bulk average.

---

## 5. Program-level analysis: challenges and statistical power

Deconvolution yields expression *programs* (cell states) per patient. Testing all of them is
where the analysis becomes fragile.

**Challenges**

1. **Low abundance → no power.** A program present at a small fraction of cells contributes
   few reads, so its deconvolved expression estimate has wide posterior uncertainty. Effect
   sizes are not distinguishable from zero at any usable FDR.
2. **Composition dependence.** Program-level estimates are conditional on θ. When a program's
   fraction is near zero in many patients, its Z values are poorly identified and the
   across-patient comparison is driven by composition rather than by expression.
3. **Multiple testing.** Genes × programs multiplies the hypothesis space; naive FDR control
   across all programs erases the few real signals.
4. **Reference dependence.** Programs absent or underrepresented in the Peng 2019 reference
   cannot be recovered reliably from bulk.

**Filtering strategy applied**

- Require a minimum program fraction per patient and a minimum number of patients above it.
- Require a minimum effective read support for the program's Z estimate.
- Test only surviving programs; apply FDR **within** the retained set.

**Result.** Most programs do not pass. Only two are statistically reliable in this cohort:

| Program | Status | Note |
|---|---|---|
| **CAF** (cancer-associated fibroblasts) | ✅ retained | Sufficient abundance across patients |
| **Basal** | ✅ retained | Sufficient abundance across patients |
| All other programs | ❌ dropped | Insufficient statistical power |

All results below are therefore reported for **CAF** and **Basal** only.

---
## 6. Compartment-specific factorial DE: methods check

### 6.1 Design

Two axis scores, discretized 0/1, fitted independently within each compartment:

```
A = basal_minus_classical   (derived from the malignant compartment)
B = myCAF_minus_iCAF        (derived from the fibroblast compartment)

expr ~ 1 + A + B + A:B
```

```python
R_mal = mc.factorial_state_de(X_mal2, disc)   # malignant compartment
R_fib = mc.factorial_state_de(X_fib2, disc)   # fibroblast compartment

contrasts_of(R_fib)
# ['A_basal_minus_classical', 'B_myCAF_minus_iCAF', 'interaction']
```

n = 108 samples. Gene universe differs by compartment after filtering:
**8,187 genes** in `X_mal2`, **11,464 genes** in `X_fib2`.

> ⚠️ **Coefficient interpretation.** Under 0/1 coding these are *simple effects*, not
> marginal ones: `beta_A` is the A effect at B = 0 (iCAF-high samples), `beta_B` is the
> B effect at A = 0 (classical samples). Sum-to-zero coding (−0.5/+0.5) is required
> before "averaged over the other axis" language is correct. `beta_AB` is invariant
> to the coding choice.

> ⚠️ **Unequal denominators.** The two compartments are FDR-corrected over different
> gene sets (8,187 vs. 11,464). Null results are therefore not directly comparable in
> strength across compartments; the intersection should be used for any cross-compartment
> claim.

### 6.2 Result: the diagonal is circular

| | tested on `X_mal2` (8,187 genes) | tested on `X_fib2` (11,464 genes) |
|---|---|---|
| **A** — malignant-derived | **143 hits** — circular | not significant¹ |
| **B** — fibroblast-derived | not significant¹ | **3 hits** — circular |

<sub>¹ Observed only conditional on the other contrast's selection; the unconditional
genome-wide test has not yet been run (see §6.5).</sub>

Each contrast reaches significance only on the matrix it was derived from. This is
guaranteed by construction and carries no biological information — it confirms that
`factorial_state_de` recovers the program it was handed, and nothing beyond that.

Effect sizes support the reading. The malignant diagonal reaches p = 4.0e−08 with
|coef| up to 3.79; the fibroblast diagonal tops out at p = 2.5e−06, |coef| ≈ 0.57.
The fibroblast axis is weaker because it is the weaker program, not because it is
less circular.

**Neither list below is a validation candidate set.**

### 6.3 Diagonal cell — malignant contrast A on `X_mal2` (143 hits)

Top 10 by p-value:

| # | Gene | Ensembl ID | A coef / FDR | raw *p* |
|---|---|---|---|---|
| 1 | PLA2G10 | ENSG00000069764 | −1.642 / 3.3e−04 | 4.0e−08 |
| 2 | LYZ | ENSG00000090382 | −1.833 / 8.8e−04 | 2.1e−07 |
| 3 | MCU | ENSG00000156026 | −0.646 / 1.9e−03 | 8.5e−07 |
| 4 | ANG | ENSG00000214274 | −0.987 / 1.9e−03 | 1.0e−06 |
| 5 | CREB3L1 | ENSG00000157613 | −1.039 / 1.9e−03 | 1.4e−06 |
| 6 | FMO5 | ENSG00000131781 | −1.179 / 1.9e−03 | 1.3e−06 |
| 7 | BCAS1 | ENSG00000064787 | −1.669 / 2.0e−03 | 1.7e−06 |
| 8 | REG4 | ENSG00000134193 | −3.788 / 8.0e−03 | 8.5e−06 |
| 9 | GPX2 | ENSG00000176153 | −1.446 / 8.0e−03 | 9.7e−06 |
| 10 | CAPN5 | ENSG00000149260 | −1.018 / 8.0e−03 | 9.7e−06 |

The full 143-gene list reproduces the Moffitt basal/classical signature at both poles:

- **Classical**: GATA6, AGR2, AGR3, TFF1, TFF2, TFF3, CLDN18, LGALS4, REG4, TSPAN8,
  HNF4G, CTSE, PGC, MYO1A, GPA33, A4GNT, CDHR5, SLC44A4
- **Basal**: KRT5, KRT6A, KRT81, FAM83A, GSDMC, SERPINB4

Recovering the classical program from a contrast defined as basal-minus-classical on
the same matrix is circular by construction.

**Flags for later.**
- **LYZ** (−1.833) is a myeloid gene. Given the low NNLS cross-validation reliability of
  the macrophage compartment (ρ = 0.197), compartment leakage is more likely than
  malignant-intrinsic expression.
- **SHH** and **IHH** both appear in the list. Hedgehog ligands are epithelium-derived
  and signal outward to stroma — one of the few well-documented tumour → stroma axes in
  PDAC. Circular here, but the first place to look if the off-diagonal test in §6.5
  returns anything.

### 6.4 Diagonal cell — fibroblast contrast B on `X_fib2` (3 hits)

Top 10 by B p-value:

| # | Gene | Ensembl ID | A coef / FDR | B coef / FDR | B raw *p* | Sig. |
|---|---|---|---|---|---|:--:|
| 1 | ADAMTS12 | ENSG00000151388 | +0.112 / 0.845 | +0.573 / 0.029 | 2.5e−06 | ✓ |
| 2 | NTM | ENSG00000182667 | −0.056 / 0.936 | +0.576 / 0.040 | 7.0e−06 | ✓ |
| 3 | TTYH2 | ENSG00000141540 | +0.063 / 0.927 | −0.547 / 0.049 | 1.3e−05 | ✓ |
| — | *FDR 0.05 cutoff* | | | | | |
| 4 | SERPINE2 | ENSG00000135919 | +0.061 / 0.967 | −1.077 / 0.051 | 1.8e−05 | |
| 5 | HTR2A | ENSG00000102468 | −0.226 / 0.806 | −0.836 / 0.101 | 5.3e−05 | |
| 6 | PLPP4 | ENSG00000203805 | +0.061 / 0.943 | +0.668 / 0.101 | 5.3e−05 | |
| 7 | STAT4 | ENSG00000138378 | +0.031 / 0.973 | −0.607 / 0.101 | 7.0e−05 | |
| 8 | ITGA11 | ENSG00000137809 | −0.068 / 0.943 | +0.721 / 0.101 | 7.0e−05 | |
| 9 | COL11A1 | ENSG00000060718 | +0.091 / 0.962 | +1.264 / 0.107 | 8.4e−05 | |
| 10 | SOX11 | ENSG00000176887 | +0.069 / 0.960 | +0.919 / 0.112 | 9.7e−05 | |

ITGA11, COL11A1 and SERPINE2 are canonical myCAF markers — again recovered by
construction.

The FDR 0.05 line falls between ranks 3 and 4 on a difference of 0.002 (TTYH2 0.049 vs.
SERPINE2 0.051), despite SERPINE2 carrying nearly twice the effect size. That boundary
is a threshold artifact, not a biological distinction. Rank order within the shortlist
is not stable at n = 108.

**STAT4** is a T-cell transcription factor appearing in a fibroblast contrast; with the
T-cell compartment at ρ = 0.136, leakage is the more parsimonious explanation.

### 6.5 Required next step: the off-diagonal

The informative comparison is each contrast tested on the compartment it was **not**
derived from, unconditionally:

```python
R_mal[R_mal.fdr_B_myCAF_minus_iCAF < 0.05]       # stroma axis → malignant expression
R_fib[R_fib.fdr_A_basal_minus_classical < 0.05]  # tumour axis → fibroblast expression
```

Restrict both to the 8,187 ∩ 11,464 gene intersection so the two FDR corrections share
a denominator.

Given the replicated malignant.prolif ↔ iCAF composition coupling (meta r = −0.476,
I² = 0), `B → X_mal2` is the more likely direction of any real effect. If both cells are
empty genome-wide, the gene-level cross-compartment null is established, and the θ-level
coupling stands as a purely compositional relationship with no within-compartment
expression correlate — a coherent and reportable negative result.

The `interaction` contrast is returned by `contrasts_of` but is not reported: at n = 108
across four factorial cells it is underpowered genome-wide, and on the diagonal it
inherits the circularity of its parent contrast.
## 7. Results: survival (Kaplan–Meier)

Patients were stratified by program-level expression of each candidate gene (high vs. low)
and by CAF / Basal program fraction.

- **Method**: Kaplan–Meier estimator, log-rank test; Cox proportional hazards for adjusted HR
- **Endpoint**: overall survival (TCGA-PAAD clinical follow-up)
- **Covariates**: stage (AJCC), age, sex

<!-- TODO: fill in n, HR, CI and p per stratification -->

| Stratification | n (high / low) | HR (95% CI) | log-rank *p* |
|---|---|---|---|
| `GENE_1`, CAF program | — | — | — |
| `GENE_2`, Basal program | — | — | — |
| `GENE_3`, Basal program | — | — | — |
| CAF fraction (high vs low) | — | — | — |

![Kaplan-Meier curves](./figures/paad_kaplan_meier.png)
*Overall survival by program-level gene expression.*

---

## 8. Downstream modules

### 8.1 Pathway perturbation modeling

1. Retrieve [Reactome](https://reactome.org/) pathways and gene sets.
2. Map DEGs onto pathways.
3. For each pathway:
   - Identify DEGs present in the pathway.
   - Expand DEG signal over the Reactome functional interaction graph, including first-order
     (1-hop) neighbors, restricted to pathway-local topology.
4. Pathway selection: find the minimum N by hypergeometric statistics; keep pathways with
   n ≥ N genes.
5. Build a perturbation profile per selected pathway, highlighting DEGs vs.
   network-propagated genes.

### 8.2 Patient representation and clustering

1. Represent each patient as a pathway perturbation vector.
2. Cluster patients on pathway-level features.

> Note: this now runs on **program-level** DEGs, not bulk DEGs (see §3.2).

### 8.3 Biological and therapeutic annotation

For each patient cluster and pathway:

1. **Gene–phenotype**: [dbGaP](https://dbgap.ncbi.nlm.nih.gov/home/) ·
   [PhenoScanner](https://github.com/phenoscanner/phenoscanner) ·
   [PheGenI](https://www.ncbi.nlm.nih.gov/gap/phegeni)
2. **Gene–disease**: [MalaCards](https://www.malacards.org/) · [DisGeNET](https://disgenet.com/)
3. **Drug associations**: [LINCS](https://clue.io/lincs) ·
   [Allosteric Database](https://mdl.shsmu.edu.cn/ASD/) ·
   [DGIdb](https://dgidb.org/about/overview/introduction) ·
   [DrugBank](https://go.drugbank.com/) · [ChEMBL](https://www.ebi.ac.uk/chembl/)
4. **Mechanism of action**: inferred from LINCS perturbation profiles

### 8.4 Visualization and reporting

1. **Pathway-level views**: perturbed genes, network structure, DEG vs. propagated genes.
2. **Cluster-level summaries**: shared pathways/genes, candidate drugs and MOA.
3. **LLM-generated summaries** using [Tahoe](https://www.tahoebio.ai/) perturbation data —
   [Tahoe Bio](https://huggingface.co/tahoebio), a gigascale single-cell perturbational atlas
   (May 2025) — for biological interpretation and therapeutic hypotheses.

---

## 9. Roadmap

### 9.1 Next — cell-type-specific differential expression

Compute DEGs directly on the BayesPrism cell-type-specific expression tensor rather than on
bulk counts, so that each DEG is attributable to a defined cell type/program. This addresses
the core limitation identified in §3.2, where bulk DE could not separate tumor from stromal
signal.

### 9.2 Final — validation against pure PDAC single-cell data

Compare deconvolution-derived programs, DEGs and candidate targets against **pure pancreatic
cancer single-cell experiments**. The single-cell data acts as ground truth: a program or
DEG recovered from bulk should reproduce in real single-cell PDAC measurements. Agreement
validates the deconvolution approach; disagreement bounds what bulk deconvolution can
legitimately claim.

### 9.3 Open questions

1. Given a primary site, subtype and stage — are all samples similar? What structure does
   clustering actually recover?
2. Do the recovered clusters resemble **EXCEPTIONAL_RESPONDERS** cohorts, or **organoid**
   models?

---

## 10. Reproducibility

<!-- TODO: complete with actual build/run commands -->

```bash
# environment
uv sync

# run the Streamlit app
uv run streamlit run src/app.py

# containerized run (R / Bioconductor: DESeq2, BayesPrism)
docker build -t perturb_agent .
docker run -p 8501:8501 perturb_agent
```

**Docker image contents**: R + Bioconductor (DESeq2, BayesPrism), Python environment,
Nextflow runtime.

---

## 11. References

**BayesPrism** — Chu T., Wang Z., Pe'er D., Danko C.G. *Cell type and gene expression
deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and
single-cell RNA sequencing in oncology.* Nat Cancer. 2022.
[Paper](https://www.nature.com/articles/s43018-022-00356-3) ·
[GitHub](https://github.com/Danko-Lab/BayesPrism)

**PDAC single-cell reference** — Peng J., Sun B.F., Chen C.Y., *et al.* *Single-cell RNA-seq
highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal
adenocarcinoma.* Cell Res. 2019 Sep;29(9):725–738. doi: 10.1038/s41422-019-0195-y.
PMID: 31273297; PMCID: PMC6796938.

---

## 12. Contact

**Flavio Lichtenstein, PhD**
📧 flalix@gmail.com
📍 São Paulo, Brazil
