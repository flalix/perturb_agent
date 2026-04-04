# 🧪  perturb_agent

Perturb Agent is a computational framework to identify **patient-specific pathway perturbations** from TCGA transcriptomic data and map them to **potential therapeutic targets**.

**Status**: 🚧 Under development


## ⚙️ Pipeline Overview

The pipeline integrates:

1. **[Streamlit](https://streamlit.io/)** for interactive exploration and visualization  
2. **[Nextflow](https://www.nextflow.io/)** for scalable, reproducible data processing  
3. **Python (ML/AI layer)** — pathway scoring, feature attribution, and target prioritization
4. **[uv](https://docs.astral.sh/uv/)** Python project and dependency management

---

## ⚙️ First results


### Interfacing GDC TCGA data, one gathered:
- 11428 cases.
- 245657 samples.
- 480826 annotated mutations.
- 18961 different genes.

---

## 💡 GDC flow

> project → project_id - gdc.list_gdc_progams()  
> primary_sites → pid and disease_type - gdc.get_primary_sites(program=program)  
> cases → case_id (UUID) - gdc.build_cases(pid=pid, subtype=subtype, stage=stage)  
  - subtypes → cancer subtypes, tissue subtypes
  - stages →  stage_id (AJCC)

> samples → sample type: [tumor, normal] and file access  
> barcodes → patients  
> annotated mutations (from [cBioPortal](https://www.cbioportal.org/))  

---

## 💡 Core components

2. **Chatbot**
   * Query any [GDC/TCGA](https://portal.gdc.cancer.gov/analysis_page?app=Projects) cancer type
   * Acts as the orchestration layer for pipeline execution

---

3. **tool1 — Mutation clusterization**
   1. Given a disease
   2. Given a subtype or state
      * Or nogthing
   3. Clusterize muation profiles
      * Calculate DEGs
      * Enriched Pathways
      * Query chatbots

---

4. **tool1 — Differential Expression (per patient)**
   1. Retrieve all patient cases (barcodes)
   2. For each patient:
      * Obtain gene expression (raw counts)
   3. Compute DEGs per patient:
      * Control: TCGA solid tissue normal samples
      * Method: [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
      * Thresholds:
        * |log2FC| ≥ 1
        * FDR < 0.05

---


5. **tool2 — Pathway Perturbation Modeling**
    1. Retrieve Reactome pathways and gene sets
    2. Map DEGs onto [Reactome](https://reactome.org/) pathways
    3. For each pathway:
        * Identify DEGs present in the pathway
        * Expand DEG signal using the Reactome functional interaction graph:
           * Include first-order neighbors (1-hop) of DEGs within the pathway graph
           * Expansion is restricted to pathway-local topology
    4. Pathway selection
        * Find the minimum N according to the hypergeometric statistics
        * Select pathway havein n >= N genes
    5. For each selected pathway:
        * Construct a perturbation profile including:
           * Highligh DEGs
           * network-propagated genes (neighbors)

---


6. **tool3 - Patient Representation & Clustering**
    1. Represent each patient as a pathway perturbation vector
    2. Cluster patients based on pathway-level features

---


7. **tool4 — Biological and Therapeutic Annotation**

  - For each patient cluster and pathway:

    1. Gene–phenotype associations:
        * [dbGaP](https://dbgap.ncbi.nlm.nih.gov/home/)
        * [PhenoScanner](https://github.com/phenoscanner/phenoscanner)
        * [PheGenI](https://www.ncbi.nlm.nih.gov/gap/phegeni)
    2. Gene–disease associations:
       * [MalaCards](https://www.malacards.org/)
       * [DisGeNET](https://disgenet.com/)
    3. Drug associations:
       * [LINCS](https://clue.io/lincs) (perturbation signatures)
       * [Allosteric Database](https://mdl.shsmu.edu.cn/ASD/)
       * [Drug-Gene Interaction Database (DGIdb)](https://dgidb.org/about/overview/introduction)
       * [DrugBank](https://go.drugbank.com/)
       * [ChEMBL](https://www.ebi.ac.uk/chembl/)
    4. Mechanism of action (MOA):
       * inferred from LINCS perturbation profiles

---

8. tool5 — Visualization & Reporting

  - Dashboard includes:

    1. Pathway-level views:
        * perturbed genes
           * network structure
           * DEG vs propagated genes
    2. Cluster-level summaries:
        * shared pathways/genes
        * candidate drugs and MOA
    3. LLM-generated summaries using [TAHOE](https://www.tahoebio.ai/) perturbed datasets
        * [Tahoe Bio LLMs](https://huggingface.co/tahoebio) a gigascale single cell perturbational atlas (May 2025)
           * biological interpretation
           * therapeutic hypotheses

---

9. Possible **'new'** questions:

    1. Given a primary site, a subtype, and stage
        * are all samples similar?
        * what kind of info returns a clusterization?
    2. For each clusterization, are they similar to:
        * EXCEPTIONAL_RESPONDERS?
        * Organoids?

---

## ⚙️  Under development

> PhD Flavio Lichtenstein  
> email: flalix@gmail.com  
> phone: +55-11-96560-1960  
> local: Brazil/Sao Paulo  
