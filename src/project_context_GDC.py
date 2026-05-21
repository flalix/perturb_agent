#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/05/08
# Udated  on 2026/05/08
# @author: Flavio Lichtenstein
# @local: Home sweet home

from dataclasses import dataclass
from pathlib import Path
import os
import numpy as np

@dataclass
class ProjectContext:

    # subtype / histology dictionaries
    SUBTYPE_GENES: dict
    HISTOLOGY_GENES: dict
    TUMOR_CLASS: dict
    GLOBAL_SUBTYPE: dict
    HISTOLOGY: dict
    SITE_MAP: dict

    # plotting
    colors: list
    nc_biotype_list: list

def load_project_context() -> ProjectContext:
    """
    Load project-level configuration from a YAML dictionary and return
    a structured ProjectContext object.

    Returns
    -------
    ProjectContext
    """


    # -------------------------
    # Mutation subtype dictionaries
    # -------------------------
    SUBTYPE_GENES = {
        "TCGA-BRCA": {
            "Luminal_A": {
                "PIK3CA",
                "GATA3",
                "MAP3K1",
                "FOXA1",
                "TBX3",
                "RUNX1",
                "CBFB",
            },
            "Luminal_B": {
                "TP53",
                "RB1",
                "CCND1",
                "MYC",
                "ERBB2",
                "PTEN",
            },
            "HER2": {
                "ERBB2",
                "GRB7",
                "PIK3CA",
                "PTEN",
                "TP53",
            },
            "Basal": {
                "TP53",
                "BRCA1",
                "RB1",
                "PTEN",
                "NF1",
            },
        },
        "TCGA-LUAD": {
            "EGFR_driven": {"EGFR", "ERBB2", "RBM10"},
            "KRAS_STK11_KEAP1_like": {"KRAS", "STK11", "KEAP1", "SMARCA4"},
            "TP53_proximal_inflammatory_like": {"TP53", "NF1", "BRAF", "RIT1"},
            "RTK_fusion_or_MAPK": {"ALK", "ROS1", "RET", "MET", "BRAF"},
        },
        "TCGA-COAD": {
            "Canonical_CIN_like": {"APC", "TP53", "KRAS", "SMAD4", "PIK3CA"},
            "MSI_like": {"BRAF", "RNF43", "ARID1A", "PIK3CA", "ACVR2A"},
            "WNT_or_TGFbeta_like": {"APC", "CTNNB1", "FBXW7", "SMAD4", "TGFBR2"},
        },
        "TCGA-READ": {
            "Canonical_CIN_like": {"APC", "TP53", "KRAS", "SMAD4", "PIK3CA"},
            "MSI_like": {"BRAF", "RNF43", "ARID1A", "PIK3CA", "ACVR2A"},
            "WNT_or_TGFbeta_like": {"APC", "CTNNB1", "FBXW7", "SMAD4", "TGFBR2"},
        },
        "TCGA-LIHC": {
            "CTNNB1_WNT_like": {"CTNNB1", "AXIN1", "APC"},
            "TP53_proliferative_like": {"TP53", "RB1", "CCNE1"},
            "Chromatin_remodeling_like": {"ARID1A", "ARID2", "BAP1"},
            "TERT_telomere_like": {"TERT", "TERF2"},
        },
    }

    HISTOLOGY_GENES = {
        "TCGA-BRCA": {
            "Lobular": {"CDH1", "PIK3CA", "FOXA1", "TBX3", "GATA3", "MAP3K1"}
        },
        "TCGA-ESCA": {
            "Adenocarcinoma": {"TP53", "CDKN2A", "SMAD4", "ERBB2", "KRAS", "ARID1A"},
            "Squamous": {"TP53", "NFE2L2", "NOTCH1", "PIK3CA", "FAT1", "KMT2D", "ZNF750"},
        },
        "TCGA-TGCT": {
            "Seminoma": {"KIT", "KRAS", "NRAS", "RAC1"},
            "Nonseminoma": {"TP53", "MDM2"},
        },
    }

    TUMOR_CLASS = {
        "adenocarcinoma": ["adenocarcinoma"],
        "squamous_cell_carcinoma": ["squamous"],
        "urothelial_carcinoma": ["urothelial"],
        "hepatocellular_carcinoma": ["hepatocellular"],
        "renal_cell_carcinoma": ["renal cell"],
        "thyroid_carcinoma": ["thyroid carcinoma"],
        "adrenal_cortical_carcinoma": ["adrenal cortical carcinoma"],
        "glioma": ["glioma", "astrocytoma", "glioblastoma"],
        "sarcoma": ["sarcoma"],
        "osteosarcoma": ["osteosarcoma"],
        "leukemia": ["leukemia"],
        "lymphoma": ["lymphoma"],
        "myeloma": ["myeloma"],
        "neuroendocrine_tumor": ["neuroendocrine"],
        "melanoma": ["melanoma"],
        "germ_cell_tumor": ["germ cell"],
        "meningioma": ["meningioma"],
        "other": [],
    }

    GLOBAL_SUBTYPE = {
        "endometrioid": ["endometrioid"],
        "serous": ["serous"],
        "clear_cell": ["clear cell"],
        "ductal": ["ductal"],
        "lobular": ["lobular"],
        "squamous": ["squamous"],
        "neuroendocrine": ["neuroendocrine"],
        "sarcoma": ["sarcoma"],
        "adenocarcinoma_generic": ["adenocarcinoma"],
    }

    HISTOLOGY = {
        "epithelial": [
            "adenocarcinoma",
            "serous",
            "endometrioid",
            "ductal",
            "lobular",
            "clear_cell",
        ],
        "squamous": ["squamous"],
        "mesenchymal": ["sarcoma"],
        "neuroendocrine": ["neuroendocrine"],
    }

    SITE_MAP = {
        "TCGA-UCEC": {
            "endometrioid": "uterine_endometrioid",
            "serous": "uterine_serous",
            "clear_cell": "uterine_clear_cell",
        },
        "TCGA-BRCA": {
            "ductal": "breast_ductal",
            "lobular": "breast_lobular",
        },
        "TCGA-LUAD": {
            "adenocarcinoma_generic": "lung_adenocarcinoma",
        },
        "TCGA-LUSC": {
            "squamous": "lung_squamous",
        },
    }

    colors = ["red", "green", "blue", "orange", "pink", "purple", "black", "cyan", "tomato", "lime", 
              "magenta", "yellow", "gray", "brown", "olive", "navy", "teal", "maroon", "silver",]

    nc_biotype_list = ['lncRNA', 'transcribed_unitary_pseudogene',
                        'unprocessed_pseudogene', 
                        'transcribed_unprocessed_pseudogene', 'processed_pseudogene',
                        'snoRNA', 'miRNA',
                        'polymorphic_pseudogene', 'transcribed_processed_pseudogene',
                        'snRNA', 'misc_RNA', 'IG_V_pseudogene',
                        'unitary_pseudogene', 'rRNA_pseudogene',
                        'translated_unprocessed_pseudogene']

    return ProjectContext(
        SUBTYPE_GENES=SUBTYPE_GENES,
        HISTOLOGY_GENES=HISTOLOGY_GENES,
        TUMOR_CLASS=TUMOR_CLASS,
        GLOBAL_SUBTYPE=GLOBAL_SUBTYPE,
        HISTOLOGY=HISTOLOGY,
        SITE_MAP=SITE_MAP,
        colors=colors,
        nc_biotype_list=nc_biotype_list,
    )