# 🧪  perturb_agent

Perturb Agent is a computational framework to identify patient-specific pathway perturbations from TCGA transcriptomic data and map them to potential therapeutic targets.

**Status**: 🚧 Under development


## ⚙️ Pipeline Overview

The pipeline integrates:

1. **Streamlit** for interactive exploration and visualization  
2. **Nextflow** for scalable, reproducible data processing  
3. **Python (ML/AI layer)** — pathway scoring, feature attribution, and target prioritization

## 💡 Core components

1. Chatbot: ask for any TCGA cancer/type
2. tool1 — Differential Expression (per patient)
   1. Retrieve all patient cases (barcodes)
   2. For each patient:
      * Obtain gene expression (raw counts)
   3. Compute DEGs per patient:
      * Control: TCGA solid tissue normal samples
      * Method: DESeq2
      * Thresholds:
        * |log2FC| ≥ 1
        * FDR < 0.05
3. tool2 — Pathway Perturbation Modeling
    1. Retrieve Reactome pathways and gene sets
    2. Map DEGs onto Reactome pathways
    3. For each pathway:
        * Identify DEGs present in the pathway
        * Expand DEG signal using the Reactome functional interaction graph:
           * Include first-order neighbors (1-hop) of DEGs within the pathway graph
           * Expansion is restricted to pathway-local topology
  4. Pathway selection:
     1. Find the minimum N according to the hypergeometric statistics
     2. Select pathway havein n >= N genes
  5. For each selected pathway:
     1. Construct a perturbation profile including:
        * Highligh DEGs
       * network-propagated genes (neighbors)
4. tool3 - Patient Representation & Clustering
    1. Represent each patient as a pathway perturbation vector
    2. Cluster patients based on pathway-level features
5. tool4 — Biological and Therapeutic Annotation
    1. For each patient cluster and pathway:
       * Gene–phenotype associations:
       * dbGaP, PhenoScanner, PheGenI
    2. Gene–disease associations:
       * MalaCards, DisGeNET
    3. Drug associations:
       * LINCS (perturbation signatures)
       * Allosteric Database
       * EpiDBase
    4. Mechanism of action (MOA):
       * inferred from LINCS perturbation profiles
6. tool5 — Visualization & Reporting
    1. Dashboard includes:
       * Pathway-level views:
          * perturbed genes
          * network structure
          * DEG vs propagated genes
       * Cluster-level summaries:
          * shared pathways/genes
          * candidate drugs and MOA
       * LLM-generated summaries:
          * biological interpretation
          * therapeutic hypotheses


## ⚙️  Under development

> PhD Flavio Lichtenstein  
> email: flalix@gmail.com  
> phone: +55-11-96560-1960  
> local: Brazil/Sao Paulo  
