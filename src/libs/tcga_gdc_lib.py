#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

import glob
from langchain_anthropic import data
import os, requests, json, re
import warnings
from tabnanny import verbose
import pandas as pd
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import plotly.express as px

# import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances

from scipy.stats import hypergeom

import hdbscan
import umap

from setuptools import glob
from typing import List, Tuple, Any, Iterable, Optional

from libs.Basic import *
from libs.stat_lib import *

class GDC(object):
	def __init__(self, root0:Path=Path('../data/')):
		
		self.url_gdc_project = "https://api.self.cancer.gov/projects"
		self.url_gdc_cases = "https://api.self.cancer.gov/cases"
		self.url_gdc_files = "https://api.self.cancer.gov/files"
		self.url_gdc_data  = "https://api.self.cancer.gov/data/%s"
		self.url_data      = "https://api.gdc.cancer.gov/data"

		self.url_cbioportal = "https://www.cbioportal.org/api"

		self.prog_id, self.psi_id = '', ''

		self.root0 = Path(root0)

		# root_data will be: ../data/TCGA
		self.root_data = Path()
		self.root_summary = Path()
		self.root_psi = Path()

		self.clean_gdc_files()

		self.fname_all_cases = '%s_summ_cases.tsv'
		self.fname_all_samples = '%s_summ_samples.tsv'
		self.fname_all_mutations = '%s_summ_mutations.tsv'

		self.colors = ['red', 'green', 'blue', 'orange', 'pink', 'purple', 'black', 'cyan', 
				 	   'tomato', 'lime', 'magenta', 'yellow', 'gray', 'brown', 'olive',
					   'navy', 'teal', 'maroon', 'silver']
		

		self.SUBTYPE_GENES2 = {'TCGA-BRCA':
				{
				"Luminal_A": {
					"PIK3CA","GATA3","MAP3K1","CDH1","TBX3","RUNX1","FOXA1"
				},
				"Luminal_B": {
					"PIK3CA","TP53","GATA3","RB1","CCND1","ERBB2"
				},
				"HER2": {
					"ERBB2","PIK3CA","TP53","PTEN"
				},
				"TNBC_Basal": {
					"TP53","BRCA1","RB1","PTEN","NF1"
				},
				"Lobular": {
					"CDH1","PIK3CA","FOXA1","TBX3","GATA3","MAP3K1"
				}
			}
		}


		# Curated starter dictionaries for mutation-based enrichment / labeling
		# Not official TCGA subtype definitions.
		# Use as a practical rule-based layer on top of your clustering.

		self.SUBTYPE_GENES = {
			"TCGA-ACC": {
				"WNT_beta_catenin": {"CTNNB1", "ZNRF3", "APC"},
				"TP53_cell_cycle": {"TP53", "RB1", "CDKN2A"},
				"Chromatin_remodeling": {"MEN1", "DAXX", "ATRX", "TERT"}
			},

			"TCGA-PCPG": {
				"Pseudohypoxia": {"VHL", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "EPAS1"},
				"Kinase_signaling": {"RET", "NF1", "HRAS", "MAX", "TMEM127"},
				"WNT_or_other": {"CSDE1", "ATRX", "SETD2", "KMT2D"}
			},

			"TCGA-BLCA": {
				"Luminal_papillary_like": {"FGFR3", "KDM6A", "STAG2", "PIK3CA", "ELF3"},
				"Basal_squamous_like": {"TP53", "RB1", "NFE2L2", "KEAP1", "FAT1", "KMT2D"},
				"Neuronal_like": {"TP53", "RB1", "ERCC2"}
			},

			"TCGA-LGG": {
				"IDH_mut_astrocytoma_like": {"IDH1", "IDH2", "TP53", "ATRX"},
				"IDH_mut_1p19q_oligodendroglioma_like": {"IDH1", "IDH2", "CIC", "FUBP1", "TERT"},
				"IDH_wildtype_progressive": {"EGFR", "PTEN", "NF1", "TERT", "PDGFRA"}
			},

			"TCGA-GBM": {
				"RTK_EGFR_like": {"EGFR", "PTEN", "TERT", "MDM4"},
				"PDGFRA_proneural_like": {"PDGFRA", "IDH1", "TP53", "ATRX"},
				"Mesenchymal_NF1_like": {"NF1", "PTEN", "RB1", "TP53"}
			},

			"TCGA-BRCA": {
				"Luminal_A": {"PIK3CA", "GATA3", "MAP3K1", "FOXA1", "TBX3", "RUNX1", "CBFB"},  # "PIK3CA","GATA3","MAP3K1","CDH1","TBX3","RUNX1","FOXA1"
				"Luminal_B": {"TP53", "RB1", "CCND1", "MYC", "ERBB2", "PTEN"}, #  {"TP53", "RB1", "CCND1", "MYC", "ERBB2", "PTEN"},
				"HER2": {"ERBB2", "GRB7", "PIK3CA", "PTEN", "TP53"}, # "ERBB2","PIK3CA","TP53","PTEN"
				"Basal": {"TP53", "BRCA1", "RB1", "PTEN", "NF1"} # "TP53","BRCA1","RB1","PTEN","NF1"
			},

			"TCGA-LUAD": {
				"EGFR_driven": {"EGFR", "ERBB2", "RBM10"},
				"KRAS_STK11_KEAP1_like": {"KRAS", "STK11", "KEAP1", "SMARCA4"},
				"TP53_proximal_inflammatory_like": {"TP53", "NF1", "BRAF", "RIT1"},
				"RTK_fusion_or_MAPK": {"ALK", "ROS1", "RET", "MET", "BRAF"}
			},

			"TCGA-LUSC": {
				"Oxidative_stress_like": {"NFE2L2", "KEAP1", "CUL3"},
				"PI3K_squamous_like": {"PIK3CA", "PTEN", "SOX2", "TP63"},
				"Cell_cycle_like": {"TP53", "CDKN2A", "RB1", "FBXW7"}
			},

			"TCGA-MESO": {
				"BAP1_chromatin_like": {"BAP1", "SETD2", "PBRM1"},
				"Hippo_pathway_like": {"NF2", "LATS1", "LATS2"},
				"DNA_damage_or_other": {"TP53", "SETDB1", "DDX3X"}
			},

			"TCGA-CESC": {
				"PI3K_squamous_like": {"PIK3CA", "EP300", "FBXW7", "KMT2C"},
				"TGF_beta_or_EMT_like": {"FAT1", "PTEN", "ARID1A"},
				"Adenocarcinoma_like": {"KRAS", "ERBB2", "ELF3"}
			},

			"TCGA-COAD": {
				"Canonical_CIN_like": {"APC", "TP53", "KRAS", "SMAD4", "PIK3CA"},
				"MSI_like": {"BRAF", "RNF43", "ARID1A", "PIK3CA", "ACVR2A"},
				"WNT_or_TGFbeta_like": {"APC", "CTNNB1", "FBXW7", "SMAD4", "TGFBR2"}
			},

			"TCGA-READ": {
				"Canonical_CIN_like": {"APC", "TP53", "KRAS", "SMAD4", "PIK3CA"},
				"MSI_like": {"BRAF", "RNF43", "ARID1A", "PIK3CA", "ACVR2A"},
				"WNT_or_TGFbeta_like": {"APC", "CTNNB1", "FBXW7", "SMAD4", "TGFBR2"}
			},

			"TCGA-DLBC": {
				"ABC_like": {"MYD88", "CD79B", "CARD11", "PIM1", "PRDM1"},
				"GCB_like": {"EZH2", "BCL2", "CREBBP", "KMT2D", "MEF2B", "TNFRSF14"},
				"Other_BCR_or_NFkB_like": {"TNFAIP3", "NFKBIE", "B2M", "SOCS1"}
			},

			"TCGA-SKCM": {
				"BRAF_like": {"BRAF", "PTEN", "MAP2K1"},
				"NRAS_like": {"NRAS", "PPP6C", "RAC1"},
				"NF1_like": {"NF1", "RASA2"},
				"Triple_wildtype_like": {"KIT", "GNAQ", "GNA11", "SF3B1"}
			},

			"TCGA-ESCA": {
				"EAC_like": {"TP53", "CDKN2A", "SMAD4", "ERBB2", "KRAS", "ARID1A"},
				"ESCC_like": {"TP53", "NFE2L2", "NOTCH1", "PIK3CA", "FAT1", "KMT2D"}
			},

			"TCGA-UVM": {
				"Galphaq_pathway": {"GNAQ", "GNA11", "CYSLTR2", "PLCB4"},
				"BAP1_high_risk_like": {"BAP1"},
				"SF3B1_intermediate_like": {"SF3B1"},
				"EIF1AX_low_risk_like": {"EIF1AX"}
			},

			"TCGA-LAML": {
				"NPM1_FLT3_like": {"NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2"},
				"RUNX1_spliceosome_like": {"RUNX1", "ASXL1", "SRSF2", "U2AF1", "STAG2"},
				"TP53_complex_karyotype_like": {"TP53", "PPM1D"},
				"CEBPA_like": {"CEBPA", "GATA2"}
			},

			"TCGA-KICH": {
				"TP53_PTEN_like": {"TP53", "PTEN"},
				"TERT_mitochondrial_or_other": {"TERT", "MT-ND5", "MT-ND1"},
				"Chromatin_like": {"SETD2", "KMT2C", "ARID1A"}
			},

			"TCGA-KIRP": {
				"Type1_MET_like": {"MET", "KRAS", "EGFR"},
				"Type2_CDKN2A_SETD2_like": {"CDKN2A", "SETD2", "BAP1", "TFE3", "FH"},
				"Chromatin_or_mTOR_like": {"KDM6A", "PBRM1", "MTOR"}
			},

			"TCGA-KIRC": {
				"VHL_PBRM1_like": {"VHL", "PBRM1", "SETD2", "KDM5C"},
				"BAP1_aggressive_like": {"BAP1", "TP53"},
				"mTOR_pathway_like": {"MTOR", "TSC1", "TSC2", "PTEN", "PIK3CA"}
			},

			"TCGA-HNSC": {
				"HPV_positive_like": {"PIK3CA", "TRAF3", "CYLD", "E2F1"},
				"HPV_negative_classic": {"TP53", "CDKN2A", "FAT1", "NOTCH1", "CASP8"},
				"NSD1_differentiated_like": {"NSD1", "HRAS", "KMT2D"}
			},

			"TCGA-LIHC": {
				"CTNNB1_WNT_like": {"CTNNB1", "AXIN1", "APC"},
				"TP53_proliferative_like": {"TP53", "RB1", "CCNE1"},
				"Chromatin_remodeling_like": {"ARID1A", "ARID2", "BAP1"},
				"TERT_telomere_like": {"TERT", "TERF2"}
			},

			"TCGA-CHOL": {
				"IDH_FGFR_like": {"IDH1", "IDH2", "FGFR2", "BAP1"},
				"KRAS_TP53_like": {"KRAS", "TP53", "SMAD4"},
				"Chromatin_or_MAPK_like": {"ARID1A", "BRAF", "PBRM1"}
			},

			"TCGA-PAAD": {
				"KRAS_core": {"KRAS", "TP53", "CDKN2A", "SMAD4"},
				"DNA_repair_deficient_like": {"BRCA1", "BRCA2", "PALB2", "ATM", "ATR"},
				"Chromatin_or_other": {"ARID1A", "KDM6A", "RNF43", "GNAS"}
			},

			"TCGA-PRAD": {
				"ETS_fusion_like": {"ERG", "ETV1", "ETV4", "PTEN"},
				"SPOP_FOXA1_like": {"SPOP", "FOXA1", "CHD1"},
				"DNA_repair_like": {"BRCA2", "ATM", "CDK12", "PALB2"},
				"IDH1_mutant_like": {"IDH1"}
			},

			"TCGA-SARC": {
				"Leiomyosarcoma_like": {"TP53", "RB1", "ATRX"},
				"Dedifferentiated_liposarcoma_like": {"MDM2", "CDK4", "HMGA2"},
				"Undifferentiated_pleomorphic_like": {"TP53", "NF1", "RB1"}
			},

			"TCGA-OV": {
				"HR_deficient_like": {"BRCA1", "BRCA2", "RAD51C", "RAD51D", "PALB2"},
				"CCNE1_amplified_like": {"CCNE1", "TP53"},
				"RB1_NF1_like": {"RB1", "NF1", "TP53"},
				"Other_HGSOC_core": {"TP53", "CDK12", "PTEN"}
			},

			"TCGA-STAD": {
				"EBV_like": {"PIK3CA", "ARID1A", "BCOR"},
				"MSI_like": {"PIK3CA", "ARID1A", "KRAS", "RNF43"},
				"Genomically_stable_diffuse_like": {"CDH1", "RHOA", "CLDN18"},
				"CIN_like": {"TP53", "ERBB2", "EGFR", "MET", "FGFR2"}
			},

			"TCGA-TGCT": {
				"Seminoma_like": {"KIT", "KRAS", "NRAS"},
				"Nonseminoma_like": {"TP53", "MDM2", "KRAS"},
				"Fusion_or_chromatin_like": {"NANOG", "SOX17", "BCOR"}
			},

			"TCGA-THYM": {
				"GTF2I_indolent_like": {"GTF2I"},
				"Aggressive_TP53_like": {"TP53", "CYLD"},
				"RAS_or_epigenetic_like": {"HRAS", "NRAS", "BAP1"}
			},

			"TCGA-THCA": {
				"BRAF_like": {"BRAF", "RET", "NTRK1", "NTRK3"},
				"RAS_like": {"NRAS", "HRAS", "KRAS", "EIF1AX", "PPARG"},
				"Advanced_dedifferentiated_like": {"TERT", "TP53", "PIK3CA", "AKT1"}
			},

			"TCGA-UCS": {
				"Serous_like": {"TP53", "PIK3CA", "PPP2R1A", "FBXW7"},
				"PI3K_pathway_like": {"PTEN", "PIK3CA", "PIK3R1"},
				"Chromatin_like": {"ARID1A", "ARID1B", "KMT2D"}
			},

			"TCGA-UCEC": {
				"POLE_ultramutated": {"POLE"},
				"MSI_hypermutated": {"MMR", "MLH1", "MSH2", "MSH6", "PMS2", "KRAS", "ARID1A"},
				"Copy_number_low_endometrioid_like": {"PTEN", "CTNNB1", "PIK3CA", "ARID1A", "KRAS"},
				"Copy_number_high_serous_like": {"TP53", "PPP2R1A", "FBXW7", "ERBB2"}
			}
		}

		self.HISTOLOGY_GENES = {
			"TCGA-BRCA": {
				# Histology axis, not molecular subtype
				"Lobular": {"CDH1", "PIK3CA", "FOXA1", "TBX3", "GATA3", "MAP3K1"}
				# Ductal_like = default / absence of lobular signal
			},

			"TCGA-ESCA": {
				"Adenocarcinoma": {"TP53", "CDKN2A", "SMAD4", "ERBB2", "KRAS", "ARID1A"},
				"Squamous": {"TP53", "NFE2L2", "NOTCH1", "PIK3CA", "FAT1", "KMT2D", "ZNF750"}
			},

			"TCGA-TGCT": {
				"Seminoma": {"KIT", "KRAS", "NRAS", "RAC1"},
				"Nonseminoma": {"TP53", "MDM2"}
			}
		}




	def clean_gdc_files(self):

		self.gdc_file_name = ''
		self.gdc_file_id = ''
		self.gdc_data_type = ''

		self.s_case = ''

		self.fname_programs = 'gdc_programs.txt'

		# primary_site
		self.fname_prim_site  = 'primary_site_program_%s.tsv'
		self.fname_cases0     = 'cases_for_%s.tsv'
		self.fname_subtype0   = 'subtype_for_%s.tsv'
		self.fname_samples0   = 'samples_for_%s.tsv'
		self.fname_vcf_files0 = 'vcf_files_for_%s.tsv'
		self.fname_rnaseq_exp_files = 'rnaseq_exp_files_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.fname_cases_deprecated = 'cases_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		# fname = self.fname_fileid%(file_type, self.psi_id, case_id, file_id)
		self.fname_fileid   = '%s_for_%s_case_%s_file_%s.%s'
		self.fname_mut_anal0 = 'mutations_anal_for_study_%s.tsv'
		self.fname_mut_summ0 = 'mutations_summ_for_study_%s.tsv'

		self.gdc_fname = ''
		self.gdc_filename = ''
		self.gdc_ouptut_fname = ''
		self.gdc_ouptut_filename = ''

		self.exp_unit = ""
		self.value_col = ""

		# program, primary site, subtype, stage, case_id, samples
		self.df_psi, self.df_subt, self.df_cases = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
		self.df_stage, self.df_samples, self.df_files = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

		'''
		primary_diagnosis
			↓
		tumor_class        ← NEW (ACC lives here)
			↓
		subtype_global     ← histology (may be "other")
			↓
		subtype_tissue

		'''

		self.TUMOR_CLASS = {
			# epithelial cancers
			"adenocarcinoma": ["adenocarcinoma"],
			"squamous_cell_carcinoma": ["squamous"],
			"urothelial_carcinoma": ["urothelial"],
			"hepatocellular_carcinoma": ["hepatocellular"],
			"renal_cell_carcinoma": ["renal cell"],
			"thyroid_carcinoma": ["thyroid carcinoma"],
			
			# organ-specific epithelial entities
			"adrenal_cortical_carcinoma": ["adrenal cortical carcinoma"],
			"glioma": ["glioma", "astrocytoma", "glioblastoma"],
			
			# mesenchymal
			"sarcoma": ["sarcoma"],
			"osteosarcoma": ["osteosarcoma"],
			
			# hematologic
			"leukemia": ["leukemia"],
			"lymphoma": ["lymphoma"],
			"myeloma": ["myeloma"],
			
			# neuroendocrine
			"neuroendocrine_tumor": ["neuroendocrine"],
			
			# melanoma
			"melanoma": ["melanoma"],
			
			# germ cell
			"germ_cell_tumor": ["germ cell"],
			
			# CNS specific
			"meningioma": ["meningioma"],
			
			# fallback
			"other": []
		}

		self.GLOBAL_SUBTYPE = {
			"endometrioid": ["endometrioid"],
			"serous": ["serous"],
			"clear_cell": ["clear cell"],
			"ductal": ["ductal"],
			"lobular": ["lobular"],
			"squamous": ["squamous"],
			"neuroendocrine": ["neuroendocrine"],
			"sarcoma": ["sarcoma"],
			"adenocarcinoma_generic": ["adenocarcinoma"]
		}

		self.HISTOLOGY = {
			"epithelial": [
				"adenocarcinoma", "serous", "endometrioid",
				"ductal", "lobular", "clear_cell"
			],
			"squamous": ["squamous"],
			"mesenchymal": ["sarcoma"],
			"neuroendocrine": ["neuroendocrine"]
		}

		self.SITE_MAP = {
			"TCGA-UCEC": {
				"endometrioid": "uterine_endometrioid",
				"serous": "uterine_serous",
				"clear_cell": "uterine_clear_cell"
			},
			"TCGA-BRCA": {
				"ductal": "breast_ductal",
				"lobular": "breast_lobular"
			},
			"TCGA-LUAD": {
				"adenocarcinoma_generic": "lung_adenocarcinoma"
			},
			"TCGA-LUSC": {
				"squamous": "lung_squamous"
			}
		}

	def text_normalization(self, x):
		if pd.isna(x):
			return ""
	
		x = x.lower().strip()
		x = re.sub(r"\bnos\b", "", x)
		x = re.sub(r"[,;()]", " ", x)
		x = re.sub(r"\s+", " ", x)
		return x.strip()


	def map_global_subtype(self, text:str) -> str:
		for k, patterns in self.GLOBAL_SUBTYPE.items():
			if any(p in text for p in patterns):
				return k
		return "other"
	
	def map_tumor_class(self, text:str) -> str:
		for k, patterns in self.TUMOR_CLASS.items():
			if any(p in text for p in patterns):
				return k
		return "other"
	
	
	def map_histology(self, subtype:str) -> str:
		for h, patterns in self.HISTOLOGY.items():
			if any(p in subtype for p in patterns):
				return h
		return "other"

	def map_tissue_subtype(self, global_subtype:str) -> str:
		if self.psi_id in self.SITE_MAP:
			return self.SITE_MAP[self.psi_id].get(global_subtype, global_subtype)
	
		return global_subtype

	def validate_consistency(self, global_subtype:str, disease_type:str) -> str:
		disease_type = self.text_normalization(disease_type)

		if "squamous" in disease_type and global_subtype != "squamous":
			return "conflict"

		return "ok"

	def simplify_stage(self, stage:str) -> Any:
		if not isinstance(stage, str):
			return None
		
		stage = stage.strip().upper()
		
		if stage.startswith("STAGE X"):
			return "missing"
		elif stage.startswith("STAGE IV"):
			return "IV"
		elif stage.startswith("STAGE III"):
			return "III"
		elif stage.startswith("STAGE II"):
			return "II"
		elif stage.startswith("STAGE I"):
			return "I"
		elif stage.startswith("UNKNOWN"):
			return "missing"
		
		return None

	def set_program(self, prog_id:str):
		self.prog_id = prog_id

		self.root_data    = create_dir(self.root0, prog_id)
		self.root_summary = create_dir(self.root_data, 'summary')
		
		self.clean_gdc_files()


	def get_primary_sites(self, prog_id:str='TCGA',force:bool=False, verbose:bool=False) -> pd.DataFrame:

		self.df_psi = pd.DataFrame()

		self.set_program(prog_id)

		fname = self.fname_prim_site%(prog_id)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_psi = pdreadcsv(fname, self.root_data, verbose=verbose)
			self.df_psi = df_psi

			return df_psi

		filters = {
					"op": "in",
					"content": { "field": "program.name", "value": [self.prog_id] }
				}

		params = {
			"filters": json.dumps(filters),
			"fields": "project_id,name,primary_site,disease_type",
			"format": "JSON",
			"size": 1000
		}

		response = None
		try:
			res = requests.get(self.url_gdc_project, params=params)
			response = res.json()
			
			if 'data' not in response.keys():
				print(f"No data found while searching for '{self.prog_id}'")
				print(">>> response", response)
				self.df_psi = pd.DataFrame()
				return self.df_psi

			hits = response["data"]["hits"]
			print(">>> hits", len(hits))

			df_psi = pd.DataFrame(hits)
			# fix list columns
			for col in df_psi.columns:
				df_psi[col] = df_psi[col].apply(
					lambda x: ", ".join(x) if isinstance(x, list) else x  )

			df_psi = df_psi.rename(columns={"id": "psi_id"})

			df_psi = df_psi.sort_values(["primary_site", "disease_type"])

			_ = pdwritecsv(df_psi, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"Error searching for '{self.prog_id}': {e}")
			print(">>> response", response)
			self.df_psi = pd.DataFrame()
			return self.df_psi

		self.df_psi = df_psi
	
		return df_psi
	

	def set_primary_site(self, psi_id:Any=None, primary_site:Any=None, verbose:bool=False) -> bool:

		self.psi_id = ''
		self.primary_site, self.disease_type, self.disease_name = '', '', ''

		if isinstance(psi_id, str) and psi_id != '':
			dfa = self.df_psi[self.df_psi.psi_id == psi_id]
			if dfa.empty:
				print("No primary site information found for:", psi_id)
				return False
		elif isinstance(primary_site, str) and primary_site != '':
			dfa = self.df_psi[self.df_psi.primary_site == primary_site]
			if dfa.empty:
				print("No primary site information found for:", primary_site)
				return False
		else:
			if verbose: print("No primary site information provided.")
			return False


		row = dfa.iloc[0]

		self.psi_id = row.psi_id
		self.primary_site = row.primary_site
		self.disease_type = row.disease_type
		self.disease_name = row.name

		self.root_psi = self.root_data / self.psi_id
		os.makedirs(self.root_psi, exist_ok=True)

		self.set_filenames()

		return True

	def get_gdc_progams(self, force:bool=False, verbose:bool=False) -> List:

		filename = os.path.join(self.root0, self.fname_programs)

		if os.path.exists(filename) and not force:
			txt = read_txt(filename, verbose=verbose)
			prog_list = eval(txt)
			return prog_list

		params = {
			"facets": "program.name",
			"size": 0
		}

		try:
			res = requests.get(self.url_gdc_project, params=params)
			buckets = res.json()["data"]["aggregations"]["program.name"]["buckets"]

			prog_list = [b["key"] for b in buckets]

			write_txt(str(prog_list), filename, verbose=verbose)

		except Exception as e:
			print(f"No programs found. Error: {e}")
			return []

		return prog_list


	def list_disease_types(self, psi_id:str) -> List:

		self.psi_id = psi_id
	
		try:
			row = self.df_psi[self.df_psi.psi_id == psi_id].iloc[0]
			deas_type_list = row.disease_type

			if isinstance(deas_type_list, str):
				deas_type_list = eval(deas_type_list)
		except:
			print("No disease types were found.")
			deas_type_list = []

		self.deas_type_list = deas_type_list

		return deas_type_list



	def set_filenames(self):
		self.fname_cases = self.fname_cases0%(self.psi_id)
		self.filename_cases = os.path.join(self.root_psi, self.fname_cases)

		self.fname_subt = self.fname_subtype0%(self.psi_id)
		self.filename_subt = os.path.join(self.root_psi, self.fname_subt)

	def apply_filter_cases(self, df_cases: pd.DataFrame) -> pd.DataFrame:
		df_cases = df_cases[df_cases.validity == 'valid'].copy()
		df_cases = df_cases[df_cases["consistency"] == "ok"]
		df_cases.reset_index(drop=True, inplace=True)
		
		# frac_threshold:float=0.01,
		# df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()
		# df_cases = df_cases[df_cases["frac"] > frac_threshold]
		# df_cases.reset_index(drop=True, inplace=True)'
		# df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()

		return df_cases
		
	def get_cases_and_subtypes(self, batch_size:int=200, 
							   do_filter:bool=True, debug:bool=False, 
				    		   force:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
		'''
		calc all subtypes, given and psi_id --> df_cases
		group by ["subtype_global", "subtype_tissue", "stage"] --> df_subt

		filter: NOS (Not Otherwise Specified) → the pathologist could not (or did not) assign a more specific subtype.
		e.g.: "Yes, it's an adenocarcinoma — but we don’t have finer classification"

		input: psi_id = primary site ID
		output: df_cases, df_subt, df_prof
		'''

		self.set_filenames()

		if os.path.exists(self.filename_cases) and os.path.exists(self.filename_subt) and not force:
			df_cases = pdreadcsv(self.fname_cases, self.root_psi, verbose=verbose)
			if 'pid' in df_cases.columns:
				df_cases = df_cases.rename(columns={'pid': 'psi_id'})
				pdwritecsv(df_cases, self.fname_cases, self.root_psi)

			self.df_cases = df_cases

			if do_filter:
				df_cases = self.apply_filter_cases(df_cases)

			df_subt = self.groupby_case_by_subtypes(df_cases)
			df_prof = self.build_profile(df_cases)

			self.df_cases = df_cases
			self.df_subt = df_subt
			self.df_prof = df_prof

			return df_cases, df_subt, df_prof
				

		def build_tcga_ontology(df):

			"""
			'id', 'primary_site', 'disease_type', 'case_id', 'diagnoses',
			'project.project_id', 'subtype_global', 'stage_ajcc', 'tumor_grade',
			'stage_clin', 'figo_stage', 'tumor_stage', 'stage'],


			print("-------- build -------------")
			print(df.columns)
			print("---------------------------")
			"""

			df["primary_site_norm"] = df["primary_site"].apply(self.text_normalization)
			df["disease_type_norm"] = df["disease_type"].apply(self.text_normalization)
			df["diagnosis_norm"]    = df["primary_diagnosis"].apply(self.text_normalization)

			df["tumor_class"]    = df["diagnosis_norm"].apply(self.map_tumor_class)
			df["subtype_global"] = df["diagnosis_norm"].apply(self.map_global_subtype)

			"""
			for psi_id=='TCGA-ACC' if subtype_global = other --> change to adrenal_cortical_carcinoma
			"""
			df["subtype_global"] = [df.iloc[i]["tumor_class"] 
						            if  (df.iloc[i]["psi_id"] =='TCGA-ACC' and df.iloc[i]["subtype_global"] =='other')
						            else df.iloc[i]["subtype_global"] for i in range(len(df))]

			# histology
			df["histology"] = df["subtype_global"].apply(self.map_histology)

			# tissue-specific subtype
			df["subtype_tissue"] = df.apply(
				lambda r: self.map_tissue_subtype(r["subtype_global"], r["psi_id"]),
				axis=1
			)

			# consistency check
			df["consistency"] = df.apply(
				lambda r: self.validate_consistency(r["subtype_global"], r["disease_type"]),
				axis=1
			)

			return df

		# nos -> removes valid dominant classes, like "Endometrioid adenocarcinoma, NOS"
		def classify_validity(row, debug:bool=False) -> str:

			if debug:
				print("---------------------")
				print(row)
				print("---------------------")

			diag = row["diagnosis_norm"]

			if diag in ["", "unknown", "not reported"]:
				return "invalid"

			if row["subtype_global"] == "other":
				return "ambiguous"

			return "valid"

		def extract_any(x, key, debug:bool=False):
			if debug: print("#", x)

			if isinstance(x, list):
				for d in x:
					if isinstance(d, dict) and key in d and d[key] is not None:
						return d[key]
			elif isinstance(x, dict):
				return x.get(key)
			return None

		def extract_all(x, key, main_diag):
			# print(">>> extract_all")
			# print(f">>> main_diag '{main_diag}' - key '{key}' - '{x}'")

			dicf = {"diagnosis": 'unknown',
					"ajcc": 'unknown'
					}
			
			i=0
			if isinstance(x, list):
				for d in x:
					# print("int loop", i)
					if isinstance(d, dict):
						diag_aux = d.get("primary_diagnosis")

						if diag_aux == main_diag:
							if key in d.keys():
								dicf = {"diagnosis": main_diag,
										"ajcc": d.get(key)
										}
								# print(">> exit1:", i, dicf, '\n\n')
								return dicf
							else:
								dicf = {"diagnosis": main_diag,
										"ajcc": 'unknown'
										}                        
					i+=1

			elif isinstance(x, dict):
				dicf = {"diagnosis": x.get("primary_diagnosis"),
						"ajcc": x.get(key)
						}

			# print(">> exit2:", i, dicf, '\n\n')
			return dicf
				
		def calc_main_diagnosis(df_cases: pd.DataFrame) -> str:
			diag_list = df_cases["diagnoses"].map(lambda x: extract_any(x, "primary_diagnosis"))
			diag_list = [x for x in diag_list if isinstance(x, str) and x.strip()]
			dic = Counter(diag_list)

			dfa = pd.DataFrame({
				'diag': list(dic.keys()),
				'n': list(dic.values())
			})

			dfa = dfa.sort_values('n', ascending=False)
			
			return dfa.iloc[0].diag

		def unpack_diagnoses(df, main_diag):

			# series = df["diagnoses"].map(lambda x: extract_all(x, "ajcc_pathologic_stage", main_diag))
			series = [extract_all(diags, "ajcc_pathologic_stage", main_diag) for diags in df_cases["diagnoses"] ]

			df["subtype_global"] = [dic["diagnosis"] for dic in series]
			df["stage_ajcc"]     = [dic["ajcc"]      for dic in series]

			df["primary_diagnosis"] = df["diagnoses"].map(lambda x: extract_any(x, "primary_diagnosis"))
			df["tumor_grade"]       = df["diagnoses"].map(lambda x: extract_any(x, "tumor_grade"))
			df["stage_clin"]        = df["diagnoses"].map(lambda x: extract_any(x, "ajcc_clinical_stage"))
			df["figo_stage"]        = df["diagnoses"].map(lambda x: extract_any(x, "figo_stage"))
			df["tumor_stage"]       = df["diagnoses"].map(lambda x: extract_any(x, "tumor_stage"))

			df["stage"] = df["stage_ajcc"] \
							.fillna(df["stage_clin"]) \
							.fillna(df["figo_stage"]) \
							.fillna(df["tumor_stage"])

			# df["stage"] = df["stage"].fillna('unknown')
			return df

		#-------------------------- batch loop ---------------------------
		filters = {
			"op": "in",
			"content": {
				"field": "cases.project.project_id",
				"value": [psi_id]
			}
		}

		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None
		df_cases = pd.DataFrame()

		try:
			while True:
				print(".", end='')

				params = {
					"filters": json.dumps(filters),
					"fields": ",".join([
						"case_id",
						"project.project_id",
						"primary_site",
						"disease_type",
						"diagnoses.primary_diagnosis",
						"diagnoses.tumor_grade",
						"diagnoses.ajcc_pathologic_stage",
						"diagnoses.ajcc_clinical_stage"
						"diagnoses.figo_stage",
						"diagnoses.tumor_stage",
					]),
					"format": "JSON",
					"size": size_,
					"from": from_
				}

				res = requests.get(self.url_gdc_cases, params=params)
				response = res.json()

				if 'data' not in response.keys():
					print(f"No data found while searching for '{psi_id}'")
					print(">>> response", response)
					self.df_cases = pd.DataFrame()
					self.df_subt  = pd.DataFrame()
					self.df_prof  = pd.DataFrame()
					return self.df_cases, self.df_subt, self.df_prof
			

				hits = response.get("data", {}).get("hits", [])

				if total is None:
					total = response["data"]["pagination"]["total"]
				
				if not hits:
					break

				all_hits.extend(hits)
				from_ += size_

			print("\n")

			if all_hits == []:
				print(f"No subtypes found for {psi_id} ")
				self.df_cases = response
				self.df_subt  = pd.DataFrame()
				self.df_prof  = pd.DataFrame()
				return self.df_cases, self.df_subt, self.df_prof
	
			#------------ lost data? ------------------
			N = len(all_hits)

			if N < total:
				print(f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} ")
			else:
				print(f"👉 Returned {N} / Total paginated {total}")

			#------------ having all hits -------------

			df_cases = pd.json_normalize(all_hits)
			self.df_cases = df_cases

			"""
			print("> 1")
			print("----------- 1 ---------------")
			print('rows', len(df_cases), '\ncolumns', df_cases.columns)
			print("---------------------------")


			"fields": ",".join([
				"case_id",
				"project.project_id",
				"primary_site",
				"disease_type",
				"diagnoses.primary_diagnosis",
				"diagnoses.tumor_grade",
				"diagnoses.ajcc_pathologic_stage",
				"diagnoses.ajcc_clinical_stage"
				"diagnoses.figo_stage",
				"diagnoses.tumor_stage",
			]),

			['id', 'primary_site', 'disease_type', 'case_id', 'diagnoses', 'project.project_id']
			"""

			# rename for sanity
			self.df_cases = df_cases

			#------------------- main_diag -------------------------------------------------------
			main_diag = calc_main_diagnosis(df_cases)
			self.main_diag = main_diag

			df_cases = unpack_diagnoses(df_cases, main_diag)
			"""
			'id', 'primary_site', 'disease_type', 'case_id', 'diagnoses',
			'project.project_id', 'subtype_global', 'stage_ajcc', 'tumor_grade',
			'stage_clin', 'figo_stage', 'tumor_stage', 'stage'],
			"""

			df_cases = df_cases.rename(columns={"project.project_id": "psi_id"})

			if debug:
				print("----------- 2 ---------------")
				print(df_cases.head(3).T)
				print("---------------------------")

			self.df_cases2 = df_cases

			if debug:
				print("----------- 3 ---------------")
				print(df_cases.head(3).T)
				print("---------------------------")

			df_cases = build_tcga_ontology(df_cases)
			
			df_cases["validity"] = df_cases.apply(classify_validity, axis=1)

			df_cases["n"] = 1
			df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()
			df_cases = df_cases.sort_values("n", ascending=False).reset_index(drop=True)
			df_cases.reset_index(drop=True, inplace=True)

			df_cases = df_cases.drop(columns=['id'])
			df_subt = self.groupby_case_by_subtypes(df_cases)

			_ = pdwritecsv(df_cases, self.fname_cases, self.root_psi, verbose=verbose)
			_ = pdwritecsv(df_subt,  self.fname_subt, self.root_psi, verbose=verbose)


		except Exception as e:
			print(f"Error for searching diags for '{psi_id}'. error: {e}")
			self.df_cases = df_cases
			self.df_subt  = pd.DataFrame()
			self.df_prof  = pd.DataFrame()
			return self.df_cases, self.df_subt, self.df_prof

		if do_filter:
			df_cases = self.apply_filter_cases(df_cases)

		df_prof = self.build_profile(df_cases)

		self.df_cases = df_cases
		self.df_subt = df_subt
		self.df_prof = df_prof
	
		return df_cases, df_subt, df_prof
	

	def group_file_types(self, df_samples:pd.DataFrame) -> pd.DataFrame:
		dic =  Counter(df_samples.data_type)
		dfu = pd.DataFrame(dic.items(), columns=["data_type", "n"])
		dfu = dfu.sort_values('n', ascending=False).reset_index(drop=True)
		return dfu

	def groupby_sstate(self, df_cases:pd.DataFrame):

		df_cases["sstage"] = df_cases["stage"].map(lambda x: self.simplify_stage(x))
		df_subt = df_cases.groupby(["psi_id", "subtype_global", "tumor_class", "subtype_tissue", "sstage"], dropna=False).size().reset_index(name="n")
		df_subt = df_subt.sort_values("n", ascending=False).reset_index(drop=True)

		self.df_subt = df_subt

		return df_subt


	def groupby_case_by_subtypes(self, df_cases:pd.DataFrame):

		# df_subt = df_cases[cols].copy().drop_duplicates()
		# df_subt = df_subt.sort_values(cols).reset_index(drop=True)		

		cols = ['psi_id', 'subtype_global', 'tumor_class', 'subtype_tissue', 'stage']
		df_subt = df_cases.groupby(cols, dropna=False).size().reset_index(name="n")
		df_subt = df_subt.sort_values("n", ascending=False).reset_index(drop=True)

		self.df_subt = df_subt

		return df_subt
	

	def build_profile(self, df_cases: pd.DataFrame) -> pd.DataFrame:

		def clean_case_profile(diagnoses) -> dict:

			result = {
				"primary_diagnosis": None,
				"stage": None,
				"tumor_grade": None,
				"diagnosis_conflict": False,
				"n_diagnoses": 0
			}

			# ---------- validate ----------
			if not isinstance(diagnoses, list) or len(diagnoses) == 0:
				return result

			result["n_diagnoses"] = len(diagnoses)

			# ---------- extract diagnoses ----------
			diag_list = [
				d.get("primary_diagnosis")
				for d in diagnoses
				if isinstance(d, dict) and d.get("primary_diagnosis") is not None
			]

			if not diag_list:
				return result

			# ---------- main diagnosis (mode) ----------
			main_diag = Counter(diag_list).most_common(1)[0][0]
			result["primary_diagnosis"] = main_diag

			# ---------- conflict detection ----------
			unique_diags = set(diag_list)
			if len(unique_diags) > 1:
				result["diagnosis_conflict"] = True

			# ---------- extract attributes ONLY for main diagnosis ----------
			for d in diagnoses:
				if not isinstance(d, dict):
					continue

				if d.get("primary_diagnosis") == main_diag:

					# stage (pathologic > clinical fallback)
					stage = (
						d.get("ajcc_pathologic_stage") or
						d.get("ajcc_clinical_stage") or
						d.get("figo_stage") or
						d.get("tumor_stage")
					)

					result["stage"] = stage

					# tumor grade
					grade = d.get("tumor_grade")
					if grade and result["tumor_grade"] is None:
						result["tumor_grade"] = grade

			return result

		profiles = df_cases["diagnoses"].apply(clean_case_profile)
		df_prof = pd.DataFrame(profiles.tolist())
		self.df_prof = df_prof

		return df_prof

	def get_stage_from_cases(self, case_id_list: list):

		df2 = self.df_cases[ ~pd.isnull(self.df_cases.stage)]
		lista = []
		for case_id in case_id_list:
			dfa = df2[ df2.case_id == case_id]

			if dfa.empty:
				lista.append(None)
			else:
				lista.append(dfa.iloc[0].stage)
		return lista

	def set_s_case(self, subtype_global:str, tumor_class:str, subtype_tissue:str):

		self.s_case = f"{self.psi_id}_{self.primary_site}_subtype_{subtype_global}_tumor_{tumor_class}_subtype_tissue_{subtype_tissue}"

		if len(self.s_case) > 180:
			self.s_case = f"{self.psi_id}_{self.primary_site[:40]}_subtype_{subtype_global[:40]}_tumor_{tumor_class[:40]}_tissue_{subtype_tissue[:40]}"

		self.s_case = title_replace(self.s_case)

	def get_samples_for_subtypes(self, subtype_global:str, tumor_class:str, subtype_tissue:str, 
								 batch_cases:int=50, batch_size:int=200, 
								 force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		return all samples given a list of cases
		for psi_id, subtype_global, tumor_class, subtype_tissue

		input: psi_id, subtype_global, tumor_class, subtype_tissue
		output: df_samples
		'''
		self.df_samples = pd.DataFrame()

		self.set_s_case(subtype_global, tumor_class, subtype_tissue)

		df_cases, _, _ = self.get_cases_and_subtypes(batch_size=batch_size, do_filter=False, 
											         force=False, verbose=verbose)
	
		self.df_cases = df_cases

		if df_cases is None or df_cases.empty:
			print(f"No cases found while searching for '{self.s_case}'")
			return self.df_samples
		
		"""
		lista=[]
		if isinstance(sstage, str):
			if sstage.startswith('I'):
				stage = 'Stage ' + sstage
			elif sstage ==  'missing':
				lista = ['unknown', 'X']
		stage = sstage

		if len(lista) > 0:
			df_cases = df_cases[(df_cases.subtype_global == subtype_global) & 
								(df_cases.tumor_class == tumor_class) &
								(df_cases.subtype_tissue == subtype_tissue) &
								(df_cases.stage.isin(lista))]
		else:
			df_cases = df_cases[(df_cases.subtype_global == subtype_global) & 
								(df_cases.tumor_class == tumor_class) &
								(df_cases.subtype_tissue == subtype_tissue) &
								(df_cases.stage == stage)]
		"""

		df_cases = df_cases[(df_cases.subtype_global == subtype_global) & 
							(df_cases.tumor_class == tumor_class) &
							(df_cases.subtype_tissue == subtype_tissue) ]

		df_cases = df_cases.copy().reset_index(drop=True)
		self.df_cases = df_cases

		if df_cases.empty:
			print(f"No cases found for {self.s_case}")
			return self.df_samples

		fname = self.fname_samples0%(self.s_case)
		fname = title_replace(fname)
		filename = os.path.join(self.root_psi, fname)

		if os.path.exists(filename) and not force:
			df_samples = pdreadcsv(fname, self.root_psi, verbose=verbose)
			self.df_samples = df_samples

			return df_samples

		case_id_list = list(df_cases.case_id)
		case_id_list.sort()

		# s_case_id_list3 = f"[{','.join(case_id_list[:3])}]"

		N_cases = len(case_id_list)
		print(f">>> {N_cases} cases")


		#-------------------------- batch loop ---------------------------
		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None
		df_samples = pd.DataFrame()

		# print("Searching: ", end='')

		ini = -batch_cases
		end = 0

		while(True):
			ini += batch_cases
			end += batch_cases

			if ini >= N_cases:
				break
			
			if end > N_cases:
				end = N_cases

			print(f"{ini}-{end} ", end='')

			lista = case_id_list[ini: end]
			# print("\n>>>", len(lista), lista)

			filters = {
				"op": "in",
				"content": {
					"field": "cases.case_id",
					"value": lista
				}
			}			

			try:
				while True:
					print(".", end='')

					params = {
						"filters": json.dumps(filters),
						"fields": ",".join([
							"file_id",
							"file_name",
							"data_type",
							"data_format",
							"experimental_strategy",
							"cases.case_id",
							"cases.submitter_id",
							"cases.samples.sample_id",
							"cases.samples.submitter_id",
							"cases.samples.sample_type",
						]),
						"format": "JSON",
						"size": size_,
						"from": from_
					}
					
					res = requests.get(self.url_gdc_files, params=params)
					response = res.json()

					if 'data' not in response.keys():
						print(f"No data found while searching for '{self.psi_id}' cases {case_id_list}")
						print(">>> response", response)
						self.df_samples = pd.DataFrame()
						return self.df_samples
				

					hits = response.get("data", {}).get("hits", [])

					if total is None:
						total = response["data"]["pagination"]["total"]
					
					if not hits:
						break

					all_hits.extend(hits)
					from_ += size_

				# print("\n")

				if all_hits == []:
					print(f"No files were found for {psi_id} cases {case_id_list}")
					self.df_samples = pd.DataFrame()
					return self.df_samples
		
				#------------ lost data? ------------------
				N = len(all_hits)

				if N < total:
					print(f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} ")
				else:
					if verbose: print(f"👉 Returned {N} / Total paginated {total}")

				#------------ having all hits -------------
				records = []

				for hit in all_hits:
					for case in hit.get("cases", []):
						for sample in case.get("samples", []):
							records.append({
								"case_id": case["case_id"],
								"submitter_id": case["submitter_id"],
								"sample_id":   sample["sample_id"],
								"sample_type": sample["sample_type"],
								"barcode_sample":  sample["submitter_id"],
								"file_id":    hit["file_id"],
								"file_name":  hit["file_name"],
								"data_type":  hit["data_type"],
								"data_format": hit["data_format"],
							})

				df_samples = pd.DataFrame(records)
				self.df_samples = df_samples
				cols = list(df_samples.columns)

				# 🔹 Metadata 
				df_samples['psi_id'] = psi_id
				df_samples['subtype_global'] = subtype_global
				df_samples['tumor_class'] = tumor_class
				df_samples['subtype_tissue'] = subtype_tissue
				df_samples['stage'] = self.get_stage_from_cases(df_samples.case_id.tolist())

				cols = ["psi_id", "subtype_global", "tumor_class", "subtype_tissue", "stage"] + cols

				df_samples = df_samples.sort_values(["case_id", "sample_type"], ascending=[False,False]).reset_index(drop=True)
				df_samples.reset_index(drop=True, inplace=True)

				_ = pdwritecsv(df_samples, fname, self.root_psi, verbose=verbose)

			except Exception as e:
				print(f"Error for searching files for {self.s_case}'. error: {e}")
				self.df_samples = df_samples
				return df_samples

		self.df_samples = df_samples

		return df_samples

	def get_table_given_fileID(self, file_type:str, case_id:str, file_id:str,
							   timeout:int=120, force:bool=False, verbose:bool=False) -> Any:
		"""
		Retrieve any kind of table like: RNA or Proteomic expression
		input: case_id and file_id
		output: the desired file
		"""

		if not file_id and not isinstance(file_id, str):
			print(f"No file_id defined.")

			self.df_table = pd.DataFrame()
			return self.df_table
		
		is_expression = False
		file_type = file_type.strip()

		if file_type == 'Gene Expression Quantification':
			is_expression = True
			type_of_file = 'tsv'
		elif file_type == 'Raw Simple Somatic Mutation':
			is_expression = True
			type_of_file = 'tar.gz'
		else:
			print(f"Develope the method for this file type {file_type}")
			raise Exception("\n------------ stop ---------------\n")

		fname = self.fname_fileid%(file_type, self.psi_id, case_id, file_id, type_of_file)
		fname = title_replace(fname)
		filename = os.path.join(self.root_psi, fname)

		if os.path.exists(filename) and not force:

			if is_expression:
				df_table = pdreadcsv(fname, self.root_psi, verbose=verbose)
				self.df_table = df_table
				return df_table
			else:
				return filename
	

		if verbose: print("Downloading: ", end='')
		try:
			url_file = self.url_gdc_data%(file_id)

			with requests.get(url_file, stream=True, timeout=timeout) as r:
				if r.status_code != 200:
					print("Error:", r.status_code)
					try:
						print(r.text[:500])
					except Exception:
						pass
					return None

				with open(filename, "wb") as f:
					for chunk in r.iter_content(chunk_size=8192):
						f.write(chunk)

			if is_expression:
				df_table = pd.read_csv(filename, sep="\t", comment="#")
				df_table = self.clean_expression_table(df_table)

				cols = ["gene_id", "symbol", "gene_type", "unstranded", "counts", "stranded_second",
						"tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"]
				df_table = df_table[cols]

				_ = pdwritecsv(df_table, fname, self.root_psi, verbose=verbose)
			else:
				return filename

		except Exception as e:
			print(f"Download error for '{self.psi_id}', '{file_type}', case {case_id} and {file_id}: {e}")
			if is_expression:
				self.df_table = pd.DataFrame()
				return self.df_table
			else:
				return "Error: " + filename

		if is_expression:
			self.df_table = df_table
			return df_table
		
		return filename
	
	def clean_expression_table(self, df:pd.DataFrame) -> pd.DataFrame:
		
		# Remove summary rows (N_*)
		df = df[~df["gene_id"].str.startswith("N_")]

		# Keep only valid Ensembl genes
		df = df[df["gene_id"].str.startswith("ENSG")]

		# Remove version from gene_id (ENSG... -> ENSG...)
		df["gene_id"] = df["gene_id"].str.split(".").str[0]

		df = df.rename(columns={"gene_name": "symbol", 'stranded_first': 'counts'})

		return df


	def get_table_searching_for_fileID(self, psi_id:str, data_type:str, sample_type:str, file_id:str,  verbose:bool=False) -> pd.DataFrame:
			
		data_type2 = title_replace(data_type)
		sample_type2 = title_replace(sample_type)

		files = [x for x in os.listdir(self.root_psi) if file_id in x and data_type2 in x and sample_type2 in x]

		if len(files) == 0:
			print(f"No files found for {file_id}.")
			self.df_table = pd.DataFrame()
			return self.df_table

		if len(files) > 1:
			print(f"Multiple files found for {file_id}. Using the first one.")

		fname = files[0]
		df_table = pdreadcsv(fname, self.root_psi, verbose=verbose)
		self.df_table = df_table

		return df_table

	def get_tumor_normal_tables(self, df_samples:pd.DataFrame, case_id_list:List[str], data_type:str, 
							    verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame]:
		'''
		Retrieve tumor and normal tables for a given case ID and data type.

		input:
			df_samples: DataFrame containing sample information
			case_id_list: List[str], list of case IDs to filter
			data_type: str, data type to filter
			verbose: bool, whether to print verbose messages
		output:
			Tuple[pd.DataFrame, pd.DataFrame]: normal and tumor tables

		'''

		self.df_samples = df_samples
		df2 = df_samples[df_samples.data_type == data_type]
		
		if df2.empty:
			print(f"No data found for data type='{data_type}'.")
			self.df_normal = pd.DataFrame()
			self.df_tumor = pd.DataFrame()
			return self.df_normal, self.df_tumor
		
		df2 = df2[df2.case_id.isin(case_id_list)]

		if df2.empty:
			print(f"No cases found for data type='{data_type}'.")
			self.df_normal = pd.DataFrame()
			self.df_tumor = pd.DataFrame()
			return self.df_normal, self.df_tumor

		case_id_list = np.unique(df2.case_id)
		if verbose: print(f"There are {len(case_id_list)} unique case IDs")

		def return_normal_tumor(x:Any, term:str):
			if not isinstance(x,str):
				return False
			
			x = x.lower()

			if 'blood' in x:
				return False

			if term in x:
				return True

			return False


		df_normal = df2[ [return_normal_tumor(x, 'normal') for x in df2.sample_type] ]

		if df_normal.empty:

			df3 = df_samples[df_samples.data_type == data_type]
			df_normal = df3[ [return_normal_tumor(x, 'normal') for x in df3.sample_type] ]

			if df_normal.empty:
				print(f"No normal samples found for case ID {case_id_list} and data type {data_type}.")
				self.df_normal = pd.DataFrame()
			else:
				self.df_normal = df_normal

		df_tumor  = df2[ [return_normal_tumor(x, 'tumor') for x in df2.sample_type] ]

		if df_tumor.empty:
			print(f"No tumor samples found for case ID {case_id_list} and data type {data_type}.")
			self.df_tumor = pd.DataFrame()
			return self.df_normal, self.df_tumor

		cols = ['file_id', 'data_type', 'sample_type']
		self.df_tumor  = df_tumor[cols].copy().reset_index(drop=True)
		self.df_normal = df_normal[cols].copy().reset_index(drop=True)

		return self.df_normal, self.df_tumor

	def merge_normal_tumor_tables(self, psi_id:str, df_normal:pd.DataFrame, df_tumor:pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

		cols = ["gene_id", "symbol", "gene_type", "counts"]
		common_cols = ["gene_id", "symbol", "gene_type"]


		dfa_normal = pd.DataFrame()

		for i, row in df_normal.iterrows():
			data_type = row.data_type
			sample_type = row.sample_type
			file_id = row.file_id
			
			dft = self.get_table_searching_for_fileID(psi_id=psi_id, data_type=data_type, sample_type=sample_type, file_id=file_id, verbose=False)

			if dft.empty:
				print(f"No data found for file_id: {i} {file_id} {data_type} and {sample_type}")
				continue

			dft = dft[cols]
			dft = dft.rename(columns={"counts": f"counts_{i+1}"})

			if dfa_normal.empty:
				dfa_normal = dft
			else:
				dfa_normal = dfa_normal.merge(dft, on=common_cols, how="outer") 


		dfa_tumor = pd.DataFrame()

		for i, row in df_tumor.iterrows():
			data_type = row.data_type
			sample_type = row.sample_type
			file_id = row.file_id
			
			dft = self.get_table_searching_for_fileID(psi_id=psi_id, data_type=data_type, sample_type=sample_type, file_id=file_id, verbose=False)

			if dft.empty:
				print(f"No data found for file_id: {i} {file_id} {data_type} and {sample_type}")
				continue

			dft = dft[cols]
			dft = dft.rename(columns={"counts": f"counts_{i+1}"})

			if dfa_tumor.empty:
				dfa_tumor = dft
			else:
				dfa_tumor = dfa_tumor.merge(dft, on=common_cols, how="outer") 

		return dfa_normal, dfa_tumor


	def get_case_id(self, barcode:str) -> str:

		self.clean_gdc_files()

		if not isinstance(barcode, str) or len(barcode) < 3:
			print(f"Barcode bad formated {barcode}.")
			return ''
		
		try:
			filters = {
				"op": "in",
				"content": {
					"field": "submitter_id",
					"value": [barcode]
				}
			}

			params = {
				"filters": json.dumps(filters),
				"fields": "case_id,submitter_id",
				"format": "JSON",
				"size": 1
			}

			response = requests.get(self.url_gdc_cases, params=params)
			data = response.json()

			hits = data.get("data", {}).get("hits", [])
			if not hits:
				raise ValueError(f"No case found for barcode {barcode}")
			
		except Exception as e:
			print(f"No data found for {barcode}. error: {e}")
			return ''

		return hits[0]["case_id"]


	#======= dummy - remove soon =================
	def get_files(self, case_id:str, data_type:str="Gene Expression Quantification") -> List:
		lista = self.get_expression_files_by_case(case_id, data_type="Gene Expression Quantification")
		return lista

	def get_expression_files(self, case_id:str, data_type:str="Gene Expression Quantification") -> List:
		lista = self.get_expression_files_by_case(case_id, data_type="Gene Expression Quantification")
		return lista


	def get_expression_files_by_case(self, case_id:str, exp_unit:str='TPM',
								     data_type:str="Gene Expression Quantification") -> List:

		if not isinstance(case_id, str) or len(case_id) < 3:
			print(f"UUID bad formated {case_id}.")
			return []


		if not isinstance(data_type, str) or len(data_type) < 3:
			print(f"Data Type bad formated {data_type}.")
			return []

		filters = {
			"op": "and",
			"content": [
				{
					"op": "in",
					"content": {
						"field": "cases.case_id",
						"value": [case_id]
					}
				},
				{
					"op": "in",
					"content": {
						"field": "data_type",
						"value": [data_type]
					}
				}
			]
		}

		try:
			params = {
				"filters": json.dumps(filters),
				"fields": "file_id,file_name,data_type",
				"format": "JSON",
				"size": 100
			}

			response = requests.get(self.url_gdc_files, params=params)
			data = response.json()

		except Exception as e:
			print(f"No data found for {case_id}. error: {e}")
			return []

		response = data["data"]["hits"]

		try:
			dic = response[0]
			if isinstance(dic, str):
				dic = eval(dic)

			self.gdc_file_name = dic['file_name']
			self.gdc_file_id = dic['file_id']
			self.gdc_data_type = dic['data_type']

			fname = f"{self.gdc_file_id}.dat"
			self.gdc_fname = fname.replace(".dat", ".tsv")

			self.gdc_filename = os.path.join(self.root_psi, fname)

			if exp_unit == 'TPM':
				self.exp_unit = exp_unit
				self.value_col ="tpm_unstranded"
			else:
				self.exp_unit = "???"
				self.value_col ="???"
				raise Exception("Error in which count col, define as TPM.")


			self.gdc_ouptut_fname = f"{fname.replace('.dat', '')}_{exp_unit}.tsv"
			self.gdc_ouptut_filename = os.path.join(self.root_psi, self.gdc_ouptut_fname)


		except Exception as e:
			print(f"No response. error: {e}")

		return response
	




	def resolve_mutation_profile(self, study_id: str) -> str:

		candidates = [
			f"{study_id}_mutations",
			f"{study_id}_mutations_extended",
		]

		for mp in candidates:
			url = f"{self.url_cbioportal}/molecular-profiles/{mp}"
			if requests.get(url, timeout=20).ok:
				return mp

		raise ValueError(f"No mutation profile found for {study_id}")


	def get_mutations_from_samples(self, barcode_sample_list: Iterable[str], study_id: str,
		session: Optional[requests.Session]=None, timeout:int=60) -> pd.DataFrame:
		"""
		Fetch mutation records from cBioPortal for a list of sample IDs.

		Parameters
		----------
		barcode_sample_list : Iterable[str]
			cBioPortal sample IDs, e.g. ["TCGA-GC-A3BM-01", "TCGA-XF-A9SY-01"].
		study_id : str
			cBioPortal study ID, e.g. "blca_tcga".
		molecular_profile_id : str | None
			Mutation profile ID. If None, defaults to "{study_id}_mutations".
		base_url : str
			cBioPortal API base URL.
		session : requests.Session | None
			Optional requests session.
		timeout : int
			Request timeout in seconds.

		Returns
		-------
		pd.DataFrame
			Mutation table. Empty DataFrame if nothing is returned.

		Notes
		-----
		- This function assumes all sample_ids belong to the same study.
		- If your samples come from multiple studies, call the function per study.
		"""

		# molecular_profile_id = f"{study_id}_mutations"
		molecular_profile_id = self.resolve_mutation_profile(study_id)

		http = session or requests.Session()

		url = f"{self.url_cbioportal}/molecular-profiles/{molecular_profile_id}/mutations/fetch"

		payload = {
			"sampleIds": barcode_sample_list
		}

		headers = {
			"Accept": "application/json",
			"Content-Type": "application/json",
		}

		resp = http.post(url, json=payload, headers=headers, timeout=timeout)

		# Helpful error message from cBioPortal
		if not resp.ok:
			msg = ""
			try:
				msg = resp.json()
			except Exception:
				msg = resp.text

			raise RuntimeError(
				f"cBioPortal request failed: HTTP {resp.status_code} | "
				f"profile={molecular_profile_id} | details={msg}"
			)

		data = resp.json()
		if not data:
			print(f"Error: cBioPortal URL: {url}")
			print(f"No mutations found for molecular profile '{molecular_profile_id}' barcodes: {barcode_sample_list}.")
			return pd.DataFrame()

		df = pd.DataFrame(data)

		cols_ori = list(df.columns)

		if "tumorAltCount" in df.columns:
			cols_ori.remove('tumorAltCount')
		else:
			df["tumorAltCount"] = None

		cols = cols_ori + ['tumorAltCount']
		df = df[cols]

		# self.df = df
		# raise Exception('stop3')

		# the selected cols + others not listed
		# cols = [c for c in df.columns if c in df.columns.to_list()] + [c for c in df.columns if c not in cols]
		# self.df = df

		df['keyword'] = [x.split(' ')[0] if isinstance(x, str) else x for x in df['keyword']]

		dic_rename = {'uniqueSampleKey':'unique_sample_key',
				'uniquePatientKey':'unique_patient_key',  'molecularProfileId':'molecular_profile_id',
				'sampleId':'barcode_sample', 'patientId':'barcode', 'entrezGeneId':'entrez_gene_id',
				'studyId':'psi_id', 'center':'center', 'mutationStatus':'mutation_status',
				'validationStatus':'validation_status', 'tumorRefCount':'tumor_ref_count', 
				'normalRefCount':'normal_ref_count', 'startPosition':'start',
				'endPosition':'end', 'referenceAllele':'ref_allele', 
				'proteinChange':'protein_mut', 'mutationType':'mutation_type',
				'ncbiBuild':'ncbi_build', 'variantType':'variant_type', 'keyword':'symbol', 
				'chr':'chr', 'variantAllele':'variant_allele',
				'refseqMrnaId':'refseq_mrna_id', 'proteinPosStart':'protein_pos_start', 
				'proteinPosEnd':'protein_pos_end', 'tumorAltCount':'tumor_alt_count'}

		rename_cols = [dic_rename.get(col, col) for col in df.columns]
		
		df.columns = rename_cols

		df['sample'] = [x.split('-')[-1] for x in df['barcode_sample'] ]

		order_cols = ['psi_id', 'molecular_profile_id', 'barcode', 'sample', 'barcode_sample',
				'symbol', 'refseq_mrna_id', 'entrez_gene_id', 
				'protein_mut', 'mutation_type', 'mutation_status',
				'ref_allele', 'variant_allele', 'variant_type', 
				'chr', 'start', 'end', 
				'validation_status', 'protein_pos_start', 'protein_pos_end', 'tumor_alt_count',
				'ncbi_build', 'center', 'tumor_ref_count', 'unique_sample_key', 'unique_patient_key']

		df = df[order_cols]
		# self.df = df

		return df

	def change_cbioportal_studyid(self, study_id: str) -> str:
		"""
		Normalize TCGA study IDs to cBioPortal PanCancer Atlas studies.

		In cBioPortal:
			COAD = colon adenocarcinoma
			READ = rectum adenocarcinoma

			👉 In PanCancer Atlas they are merged into one cohort:
		"""

		dic = {
			"acc_tcga": "acc_tcga_pan_can_atlas_2018",
			"luad_tcga": "luad_tcga_pan_can_atlas_2018",
			"lusc_tcga": "lusc_tcga_pan_can_atlas_2018",
			"coad_tcga": "coadread_tcga_pan_can_atlas_2018",
			"read_tcga": "coadread_tcga_pan_can_atlas_2018",
			"brca_tcga": "brca_tcga_pan_can_atlas_2018",
			"gbm_tcga":  "gbm_tcga_pan_can_atlas_2018",
			"ov_tcga":   "ov_tcga_pan_can_atlas_2018",
			"skcm_tcga": "skcm_tcga_pan_can_atlas_2018",
			"ucec_tcga": "ucec_tcga_pan_can_atlas_2018",
			"stad_tcga": "stad_tcga_pan_can_atlas_2018",
			"blca_tcga": "blca_tcga_pan_can_atlas_2018",
			"hnsc_tcga": "hnsc_tcga_pan_can_atlas_2018",
			"kirc_tcga": "kirc_tcga_pan_can_atlas_2018",
			"kirp_tcga": "kirp_tcga_pan_can_atlas_2018",
			"lihc_tcga": "lihc_tcga_pan_can_atlas_2018",
			"prad_tcga": "prad_tcga_pan_can_atlas_2018",
			"thca_tcga": "thca_tcga_pan_can_atlas_2018",
			"esca_tcga": "esca_tcga_pan_can_atlas_2018",
			"paad_tcga": "paad_tcga_pan_can_atlas_2018",
			"kich_tcga": "kich_tcga_pan_can_atlas_2018",  # kidney chromophobe
			"sarc_tcga": "sarc_tcga_pan_can_atlas_2018",
			"pcpg_tcga": "pcpg_tcga_pan_can_atlas_2018",
			"tgct_tcga": "tgct_tcga_pan_can_atlas_2018",
			"thym_tcga": "thym_tcga_pan_can_atlas_2018",
			"meso_tcga": "meso_tcga_pan_can_atlas_2018",
			"ucs_tcga":  "ucs_tcga_pan_can_atlas_2018",
			"uvm_tcga":  "uvm_tcga_pan_can_atlas_2018",
			"chol_tcga": "chol_tcga_pan_can_atlas_2018",
			"dlbc_tcga": "dlbc_tcga_pan_can_atlas_2018",
		}

		return dic.get(study_id, study_id)
	
	def set_mutation_filenames(self):
		self.fname_mut_anal = self.fname_mut_anal0%(self.s_case)
		self.fname_mut_anal = title_replace(self.fname_mut_anal)
		self.filename_mutanal = self.root_psi / self.fname_mut_anal

		self.fname_mut_summ = self.fname_mut_summ0%(self.s_case)
		self.fname_mut_summ = title_replace(self.fname_mut_summ)
		self.filename_mutsumm = self.root_psi / self.fname_mut_summ


	def get_df_mut_transform_mutation_table(self, study_id: str, barcode_sample_list: List[str], 
		session: Optional[requests.Session] = None, timeout: int=60,
		force:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame]:

		self.study_id0 = study_id

		'''
		if TCGA remove the last characters if len > 2
		TCGA-OR-A5J2-01A -> TCGA-OR-A5J2-01
		'''
		barcode_list = self.prepare_barcode_sample_list(barcode_sample_list)	

		if study_id[0].isupper():
			mat = study_id.lower().split('-')
			# cBioPortal disease - tcga 
			study_id = mat[1] + '_' + mat[0]

		study_id = self.change_cbioportal_studyid(study_id)
		self.study_id = study_id

		print(f"\n>>> {study_id} --> {self.s_case} len = {len(barcode_list)} - {barcode_list[:5]}...")

		self.set_mutation_filenames()

		if os.path.exists(self.filename_mutanal) and os.path.exists(self.filename_mutsumm) and not force:
			dff    = pdreadcsv(self.fname_mut_summ, self.root_psi, verbose=verbose)
			df_mut = pdreadcsv(self.fname_mut_anal, self.root_psi, verbose=verbose)
			return dff, df_mut
		
		'''
			df_mut cols: ["sample_id", "barcode_sample", "psi_id", "mol_profile_id","gene",
			"entrez_gene_id", "protein_mut", "mutation_type", "mutation_status",
			"variant_type", "chr", "start", "end",
			"ref_allele", "tumor_seq_allele"]		
		'''
		df_mut = self.get_mutations_from_samples(barcode_sample_list=barcode_sample_list, study_id=study_id,
											 	session=session, timeout=timeout)
		
		self.df_mut = df_mut
		
		if df_mut.empty:
			print("No mutations found for these samples.")
			return pd.DataFrame(), pd.DataFrame()

		#--------------- map main cols from df_mut ------------------------
		"""
		order_cols = ['barcode_sample', 'barcode_sample', 'psi_id', 'mol_profile_id',
			'symbol', 'refseq_mrna_id', 'entrez_gene_id', 
			'protein_mut', 'mutation_type', 'mutation_status',
			'ref_allele', 'variant_allele', 'variant_type', 
			'chr', 'start', 'end', 
			'validation_status', 'protein_pos_start', 'protein_pos_end', 'tumor_alt_count',
			'ncbi_build', 'center', 'tumorRefCount', 'unique_sample_key',
			'unique_patient_key']
		"""
		dff = (df_mut
					.groupby(['psi_id', 'barcode', "barcode_sample", 'symbol', 'refseq_mrna_id', "entrez_gene_id", "protein_mut", 'mutation_type', "variant_type", "chr"])
					.size()
					.reset_index(name="n_mutations")
				)
		
		dff = dff[dff.barcode.notna()]
		dff = dff[dff.barcode_sample.notna()]
		dff = dff[dff.entrez_gene_id.notna()]

		dff = dff.sort_values(["barcode", "symbol", "protein_mut"])
		dff = dff.reset_index(drop=True)
		
		self.dff = dff

		_ = pdwritecsv(dff,    self.fname_mut_summ, self.root_psi, verbose=False)
		_ = pdwritecsv(df_mut, self.fname_mut_anal, self.root_psi, verbose=False)

		return dff, df_mut


	def cbioportal_studies(self):
		url = "https://www.cbioportal.org/api/studies"

		res = requests.get(url, headers={"Accept": "application/json"})
		res.raise_for_status()

		studies = res.json()

		study_ids = [s["studyId"] for s in studies]

		return study_ids
	
	def loop_program_psi_samples(self, prog_id:str='TCGA', force:bool=False, 
			verbose:bool=True) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
		
		df_psi = self.get_primary_sites(prog_id=prog_id, force=force, verbose=verbose)

		df_cases, df_subt, df_prof = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

		fname_all_cases = self.fname_all_cases%(self.prog_id)
		filename_cases = os.path.join(self.root_summary, fname_all_cases)

		fname_all_samples = self.fname_all_samples%(self.prog_id)
		filename_samples = os.path.join(self.root_summary, fname_all_samples)

		fname_all_mutations = self.fname_all_mutations%(self.prog_id)
		filename_mutations = os.path.join(self.root_summary, fname_all_mutations)

		if os.path.exists(filename_cases) and os.path.exists(filename_samples) and \
		   os.path.exists(filename_mutations) and not force:

			df_all_cases = pdreadcsv(fname_all_cases, self.root_summary)
			df_all_samples = pdreadcsv(fname_all_samples, self.root_summary)
			df_all_mutations = pdreadcsv(fname_all_mutations, self.root_summary)

			self.df_all_cases = df_all_cases
			self.df_all_samples = df_all_samples
			self.df_all_mutations = df_all_mutations
			
			return df_all_cases, df_all_samples, df_all_mutations


		lista = np.arange(len(df_psi))

		df_list_cases, df_list_samples, df_list_mutations = [], [], []

		for ipsi in lista:
			row = df_psi.iloc[ipsi]
			psi_id = row.psi_id
			primary_site = row.primary_site

			self.set_primary_site(psi_id)

			print(f'{ipsi}) {primary_site}', end=' - ')

			df_cases, df_subt, _ = self.get_cases_and_subtypes(batch_size=200, do_filter=False, force=force, verbose=verbose)

			if df_cases.empty:
				print("No cases found for PSI_ID:", psi_id)
				continue

			if isinstance(df_cases, pd.DataFrame):
				df_list_cases.append(df_cases)
			else:
				print("Unexpected type for df_cases:", type(df_cases))
				raise Exception("Stope: unexpected type for df_cases")


			for isubt, row in df_subt.iterrows():
				subtype_global = row.subtype_global
				tumor_class = row.tumor_class
				subtype_tissue = row.subtype_tissue


				df_samples = self.get_samples_for_subtypes(subtype_global=subtype_global,
														   tumor_class=tumor_class, subtype_tissue=subtype_tissue,
														   batch_size=200, force=force, verbose=verbose)
				print(f'{isubt}) {self.s_case}')
				
				if df_samples.empty:
					print("No samples found for PSI_ID:", psi_id)
					continue

				df_list_samples.append(df_samples)

				df2 = df_samples[~df_samples.sample_type.str.contains('Blood', case=False, na=False)]

				if df2.empty:
					print("No samples having non-blood types for PSI_ID:", psi_id)
					continue

				barcode_sample_list = list(np.unique(df2.barcode_sample))

				print("Getting mutations", end=' ')
				dff, _ = self.get_df_mut_transform_mutation_table(study_id=self.psi_id, barcode_sample_list=barcode_sample_list, force=force, verbose=verbose)

				if dff.empty:
					print("Could not find mutations for :", self.s_case)
					continue

				df_list_mutations.append(dff)

		if len(df_list_cases) > 0:
			df_all_cases = pd.concat(df_list_cases, ignore_index=True)
			df_all_cases = df_all_cases.drop_duplicates()
			df_all_cases = df_all_cases.reset_index(drop=True)
		else:
			df_all_cases = pd.DataFrame()

		if len(df_list_samples) > 0:
			df_all_samples = pd.concat(df_list_samples, ignore_index=True)
			df_all_samples = df_all_samples.drop_duplicates()
			df_all_samples = df_all_samples.reset_index(drop=True)
		else:
			df_all_samples = pd.DataFrame()

		if len(df_list_mutations) > 0:
			df_all_mutations = pd.concat(df_list_mutations, ignore_index=True)
			df_all_mutations = df_all_mutations.drop_duplicates()
			df_all_mutations = df_all_mutations.reset_index(drop=True)
		else:
			df_all_mutations = pd.DataFrame()

		_ = pdwritecsv(df_all_cases, fname_all_cases, self.root_summary)
		_ = pdwritecsv(df_all_samples, fname_all_samples, self.root_summary)
		_ = pdwritecsv(df_all_mutations, fname_all_mutations, self.root_summary)

		self.df_all_cases = df_all_cases
		self.df_all_samples = df_all_samples
		self.df_all_mutations = df_all_mutations
		   
		return df_all_cases, df_all_samples, df_all_mutations


	def get_filtered_tables(self, primary_site: str, 
						    verbose: bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:

		dfa = self.df_psi[self.df_psi.primary_site == primary_site]

		if dfa.empty:
			self.psi_id = ''
			print("No primary site information found for:", primary_site)
			return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), []
		
		row = dfa.iloc[0]
		self.set_primary_site(psi_id = row.psi_id)

		df_cases, df_all_samples, df_all_mut, all_barcode_list = \
			self.get_filtered_tables_subtypes(sample_type='tumor', do_filter=True, verbose=verbose)
		

		if df_cases.empty:
			print("No cases found for primary site:", self.primary_site)
		
		if df_all_samples.empty:
			print("No samples found for primary site:", self.primary_site)

		
		if df_all_mut.empty:
			print("No mutations found for primary site:", self.primary_site)

		if len(all_barcode_list) == 0:
			print("No barcodes found for primary site:", self.primary_site)

		return df_cases, df_all_samples, df_all_mut, all_barcode_list


	def get_filtered_tables_subtypes(self, sample_type:str='tumor', do_filter:bool=True,
								     verbose:bool=True) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:

		self.df_cases   = pd.DataFrame()
		self.df_all_samples = pd.DataFrame()
		self.df_all_mut = pd.DataFrame()
		self.all_barcode_list = []

		if not os.path.exists(self.filename_cases):
			print("Error: could not find cases file:", self.filename_cases)
			return self.df_cases, self.df_all_samples, self.df_all_mut, self.all_barcode_list

		df_cases = pdreadcsv(self.fname_cases, self.root_psi, verbose=verbose)
		if 'pid' in df_cases.columns:
			df_cases = df_cases.rename(columns={'pid': 'psi_id'})
			pdwritecsv(df_cases, self.fname_cases, self.root_psi)


		if do_filter:
			df_cases = self.apply_filter_cases(df_cases)

		self.df_cases = df_cases

		df_subt = self.groupby_case_by_subtypes(df_cases)

		df_list_samples =[]
		df_list_mut = []
		list_all_barcodes = []

		for _, row in df_subt.iterrows():
			subtype_global = row.subtype_global
			tumor_class    = row.tumor_class
			subtype_tissue = row.subtype_tissue

			self.set_s_case(subtype_global, tumor_class, subtype_tissue)

			df3 = df_cases[ (df_cases.subtype_global == subtype_global) & 
							(df_cases.tumor_class == tumor_class) &
							(df_cases.subtype_tissue == subtype_tissue) ]
			
			if df3.empty:
				print(f"No cases found for {subtype_global} {tumor_class} {subtype_tissue}")
				continue
					
			df3 = df3.copy().reset_index(drop=True)

			case_id_list = np.unique(df3.case_id)

			fname = self.fname_samples0%(self.s_case)
			fname = title_replace(fname)
			filename = os.path.join(self.root_psi, fname)

			if not os.path.exists(filename):
				print("Error: could not find samples file:", filename)
				continue
		
			df_samples = pdreadcsv(fname, self.root_psi, verbose=verbose)

			if df_samples.empty:
				print("Error: could not read samples file:", filename)
				continue

			df_samples = df_samples[ (df_samples.case_id.isin(case_id_list)) & (df_samples.sample_type.str.contains(sample_type, case=False)) ]
			self.df_samples = df_samples

			if df_samples.empty:
				print("Error: could not filter df_samples")
				continue
			
			df_samples = df_samples.copy().reset_index(drop=True)
			self.barcode_list = self.prepare_barcode_sample_list(df_samples.barcode_sample.tolist())	

			self.set_mutation_filenames()

			if os.path.exists(self.filename_mutsumm):
				df_mut = pdreadcsv(self.fname_mut_anal, self.root_psi, verbose=verbose)
			else:
				print("No mutation analysis file found for:", self.s_case)
				df_mut = pd.DataFrame()
			
			df_list_samples.append(df_samples)
			if not df_mut.empty:
				df_list_mut.append(df_mut)
			list_all_barcodes += list(self.barcode_list)

		df_all_samples = pd.concat(df_list_samples, ignore_index=True) if df_list_samples else pd.DataFrame()
		df_all_mut = pd.concat(df_list_mut, ignore_index=True) if df_list_mut else pd.DataFrame()
		list_all_barcodes = list(np.unique(list_all_barcodes))

		return self.df_cases, df_all_samples, df_all_mut, list_all_barcodes


	def prepare_barcode_sample_list(self, barcode_sample_list: list[str]) -> list[str]:

		'''
		if TCGA remove the last characters if len > 2
		TCGA-OR-A5J2-01A -> TCGA-OR-A5J2-01
		'''

		barcode_sample_list = [self.to_cbioportal_barcode_sample(x) for x in barcode_sample_list]
		if not barcode_sample_list:
			raise ValueError("barcode_sample_list is empty.")
		
		# 01A, 01B, 01Z → all collapse to 01
		barcode_sample_list = list(np.unique(barcode_sample_list))

		return barcode_sample_list
	

	def to_cbioportal_barcode_sample(self, x: str) -> str:
		parts = x.split("-")

		if len(parts) >= 4 and parts[0] == "TCGA":
			sample_code = parts[3][:2]   # 01A -> 01, 11A -> 11
			return "-".join([parts[0], parts[1], parts[2], sample_code])
		
		return x
	

	def build_pivot_table(self, df_all_mut: pd.DataFrame, min_barcodes:int=2, min_genes:int=2) -> pd.DataFrame:
		"""
		Build a barcode x gene binary mutation matrix (0/1).

		Rows represent barcodes (samples), columns represent gene symbols.
		A value of 1 indicates that at least one mutation was observed for that
		barcode-gene pair.

		Parameters
		----------
		df_all_mut : pd.DataFrame
			Input mutation table. Must contain at least:
			- 'barcode': sample identifier
			- 'symbol': gene symbol

		Returns
		-------
		pd.DataFrame
			Binary mutation matrix with:
			- index   = barcode
			- columns = gene symbol
			- values  = 0 or 1

			Returns an empty DataFrame if the input is empty or required columns
			are missing.
		"""
		if df_all_mut is None or df_all_mut.empty:
			return pd.DataFrame()

		required_cols = {"barcode", "symbol"}
		missing_cols = required_cols - set(df_all_mut.columns)
		if missing_cols:
			raise ValueError(
				f"build_pivot_table requires columns {sorted(required_cols)}, "
				f"but is missing {sorted(missing_cols)}."
			)

		# Keep only the columns needed for the mutation matrix
		dfa = df_all_mut.loc[:, ["barcode", "symbol"]].copy()
		dfa["barcode"] = dfa["barcode"].astype(str).str.strip()
		dfa["symbol"] = dfa["symbol"].astype(str).str.strip()

		dfa = dfa[
			(dfa["barcode"] != "")
			& (dfa["symbol"] != "")
			& (dfa["barcode"].str.lower() != "nan")
			& (dfa["symbol"].str.lower() != "nan")
		]

		if dfa.empty:
			return pd.DataFrame()
		
		# Mark presence of at least one mutation per barcode-gene pair
		dfa["mutated"] = 1

		dfpiv = dfa.pivot_table(
			index="barcode",
			columns="symbol",
			values="mutated",
			aggfunc="max",
			fill_value=0,
		).astype(np.uint8)

		# Remove empty samples and genes
		dfpiv = dfpiv.loc[dfpiv.sum(axis=1) >= min_genes, :]

		dfpiv = dfpiv.loc[:, dfpiv.sum(axis=0) >= min_barcodes]

		if dfpiv.shape[0] < 3:
			print("dfpiv has less than 3 samples.")
			return pd.DataFrame()

		if dfpiv.shape[1] < 3:
			print("dfpiv has less than 3 genes.")
			return pd.DataFrame()

		'''
		It sorts your matrix in a consistent order:

		axis=0 → sort rows (barcodes)
		axis=1 → sort columns (genes)

		So after this:
			barcodes are alphabetically (or lexicographically) ordered
			genes are alphabetically ordered
		'''
		dfpiv = dfpiv.sort_index(axis=0).sort_index(axis=1)

		return dfpiv


	def calc_HDBSCAN(self, dfpiv: pd.DataFrame, min_cluster_size:int=8, 
				     min_samples:int=3) -> tuple[List, List, Any]:
		"""
		Cluster with HDBSCAN, not KMeans
		pairwise_distances with jaccard
		Multidimensional Scaling (MDS) 
		If there are a few dense groups plus many ambiguous samples, HDBSCAN can work better.

		In HDBSCAN:
			min_cluster_size → minimum size of a cluster
			min_samples → minimum local neighborhood density
			A point is considered “core” if it has at least min_samples neighbors nearby.

		input: dfpiv, min_cluster_size (minimum number of samples in a cluster)
		output: embedding and labels
		"""
		
		X = dfpiv.to_numpy(dtype=bool)

		n_samples = X.shape[0]
		n_genes = X.shape[1]

		if n_samples < 3:
			print("Need at least 3 non-empty samples to compute HDBSCAN + clustering.")
			return [], [], None
		
		if min_cluster_size > n_samples:
			print(f"min_cluster_size={min_cluster_size} is larger than number of samples ({X.shape[0]}). Using min_cluster_size={X.shape[0]}.")
			min_cluster_size = max(2, n_samples - 1)

		if n_genes < 3:
			print("Need at least 3 non-empty genes to compute HDBSCAN + clustering.")
			return [], [], None
		
		min_cluster_size = min(min_cluster_size, n_samples)

		print(">>> calc_HDBSCAN MIN_CLUSTER_SIZE", min_cluster_size)

		D = pairwise_distances(X, metric="jaccard")

		embedding = MDS(
			n_components=2,
			dissimilarity="precomputed",
			n_init=8,
			init='classical_mds',
			random_state=42,
		).fit_transform(D)

		if isinstance(embedding, tuple):
			print("embedding return as a tuple")
			embedding = embedding[0]

		embedding = np.asarray(embedding)

		if embedding.shape[0] < 3:
			print("Too few valid embedded samples after filtering.")
			return [], [], None


		clusterer = hdbscan.HDBSCAN(
			min_cluster_size=min_cluster_size,
			min_samples=min_samples if min_samples is not None else min_cluster_size,
			metric="euclidean"
		)
		labels = clusterer.fit_predict(embedding)

		return embedding.tolist(), labels.tolist(), D

	def calc_UMAP(self, dfpiv: pd.DataFrame, k:int=8) -> tuple[List, List]:
		# Binary mutation matrix for Jaccard
		# Force numeric/binary and remove bad values
		
		X = dfpiv.to_numpy(dtype=np.uint8)

		n_samples = X.shape[0]
		n_genes = X.shape[1]

		if n_samples < 3:
			print("Need at least 3 non-empty samples to compute UMAP + clustering.")
			return [], []
		
		if k > n_samples:
			print(f"k={k} is larger than number of samples ({X.shape[0]}). Using k={X.shape[0]}.")
			k = max(2, n_samples - 1)

		if n_genes < 3:
			print("Need at least 3 non-empty genes to compute UMAP + clustering.")
			return [], []
		
		k = min(k, n_samples)

		n_neighbors = min(15, n_samples - 1)
		n_neighbors = max(2, n_neighbors)    

		reducer = umap.UMAP(
			n_neighbors=n_neighbors,
			min_dist=0.1,
			metric="jaccard",
			random_state=42,
			init="random",      # important
			output_dens=False   # avoid tuple output
		)

		with warnings.catch_warnings():
			warnings.filterwarnings("ignore", category=UserWarning, module="umap")
			embedding = reducer.fit_transform(X)

		if isinstance(embedding, tuple):
			print("embedding return as a tuple")
			embedding = embedding[0]

		embedding = np.asarray(embedding)

		'''
		good = np.isfinite(embedding).all(axis=1)
		embedding = embedding[good]
		'''
		if not np.isfinite(embedding).all():
			print("UMAP embedding contains NaN or infinite values.")
			return [], []

		if embedding.shape[0] < 3:
			print("Too few valid embedded samples after filtering.")
			return [], []

		labels = KMeans(
			n_clusters=min(k, embedding.shape[0] - 1),
			random_state=42,
			n_init=10,
		).fit_predict(embedding)

		return embedding.tolist(), labels.tolist()

	def plot_UMAP(self, dfpiv: pd.DataFrame, k:int=8, figsize:tuple=(14, 10)) -> Tuple[Any, Any, Any]:

		n_samples = dfpiv.shape[0]
		n_genes = dfpiv.shape[1]

		embedding, labels = self.calc_UMAP(dfpiv, k)
		embedding = np.array(embedding)

		if len(embedding) == 0 or len(labels) == 0:
			print("No valid UMAP embedding or labels.")
			return None, None, None

		fig, ax = plt.subplots(figsize=figsize)

		# cmap = plt.cm.get_cmap("tab10", k)
		sc = plt.scatter(
			embedding[:, 0],
			embedding[:, 1],
			c=[self.colors[label] for label in labels],
			s=20
		)

		ax.set_title(f"Clustering using UMAP mutation profiles: (k={k})\nPrimary Site: '{self.primary_site}' #{n_samples} samples and #{n_genes} genes")
		ax.set_xlabel("UMAP1")
		ax.set_ylabel("UMAP2")

		counts = Counter(labels)

		legend_handles = []

		for cluster_id in sorted(counts.keys()):
			color = self.colors[cluster_id]
			n = counts[cluster_id]

			patch = mpatches.Patch(
				color=color,
				label=f"Cluster {cluster_id} (n={n})"
			)
			legend_handles.append(patch)

		ax.legend(
			handles=legend_handles,
			title="Groups",
			loc="best"
		)

		plt.show()
		return fig, embedding, labels


	def plot_HDBSCAN(self, dfpiv: pd.DataFrame, min_cluster_size:int=8, min_samples:int=3, figsize:tuple=(14, 10)) -> Tuple[Any, Any, Any, Any]:

		n_samples = dfpiv.shape[0]
		n_genes = dfpiv.shape[1]

		embedding, labels, d = self.calc_HDBSCAN(dfpiv=dfpiv, min_cluster_size=min_cluster_size, min_samples=min_samples)
		embedding = np.array(embedding)

		if len(embedding) == 0 or len(labels) == 0:
			print("No valid HDBSCAN embedding or labels.")
			return None, None, None, None

		fig, ax = plt.subplots(figsize=figsize)

		sc = plt.scatter(
			embedding[:, 0],
			embedding[:, 1],
			c=[self.colors[label] for label in labels],
			s=20
		)

		stri = f"Clustering using HDBSCAN mutation profiles"
		stri += f"\nPrimary Site: {self.psi_id} - '{self.primary_site}' #{n_samples} samples and #{n_genes} genes"
		stri += f"\nmin_cluster_size={min_cluster_size} and min_samples={min_samples}"
		ax.set_title(stri)
		ax.set_xlabel("embedding1")
		ax.set_ylabel("embedding2")

		counts = Counter(labels)

		legend_handles = []

		for cluster_id in sorted(counts.keys()):
			color = self.colors[cluster_id]
			n = counts[cluster_id]

			patch = mpatches.Patch(
				color=color,
				label=f"Cluster {cluster_id} (n={n})"
			)
			legend_handles.append(patch)

		ax.legend(
			handles=legend_handles,
			title="Groups",
			loc="best"
		)
		plt.show()

		return fig, embedding, labels, d
	
	def cluster_mutation_table(self, dfpiv:pd.DataFrame, labels, cluster:int=1,
                           min_barcodes:int=2) -> pd.DataFrame:

		if len(labels) != dfpiv.shape[0]:
			stri  = "Error: build_pivot_table filter empty lines."
			stri += "\nNumber of labels does not match number of samples. "
			stri += "\n----------- stop execution -----------\n"
			raise Exception(stri)

		labels = pd.Series(labels, index=dfpiv.index)
		sel_barcodes = labels[labels == cluster].index
		dff = dfpiv.loc[sel_barcodes].copy()
		dff = dff.loc[:, dff.sum(axis=0) >= min_barcodes]

		return dff
	
	def calc_shannon_entropy_from_dfstat(self, dfstat: pd.DataFrame) -> pd.DataFrame:
		rows = []

		for (k, cluster), dfsub in dfstat.groupby(["k", "cluster"]):
			deg_list = dfsub["degree"].to_list()

			if len(deg_list) == 0:
				H = np.nan
				Hmax = np.nan
				Hnorm = np.nan
				n_genes = 0
			else:
				w = np.array(deg_list, dtype=float)
				p = w / w.sum()
				H = -np.sum(p * np.log2(p))
				n_genes = len(p)
				Hmax = np.log2(n_genes) if n_genes > 1 else 0.0
				Hnorm = H / Hmax if Hmax > 0 else 0.0

			rows.append({
				"k": k,
				"cluster": cluster,
				"n_genes": n_genes,
				"cluster_size": dfsub["cluster_size"].iloc[0],
				"entropy": H,
				"entropy_max": Hmax,
				"entropy_norm": Hnorm
			})

		if rows == []:
			return pd.DataFrame()

		dfh = pd.DataFrame(rows)
		dfh = dfh.sort_values("entropy_norm", ascending=True)

		return dfh
	
	def score_k_from_entropy_table(self, dfh: pd.DataFrame) -> pd.DataFrame:
		rows = []

		for k, sub in dfh.groupby("k"):
			total_n = sub["cluster_size"].sum()

			weighted_mean_entropy = (
				(sub["entropy_norm"] * sub["cluster_size"]).sum() / total_n
				if total_n > 0 else np.nan
			)

			mean_entropy = sub["entropy_norm"].mean()
			std_entropy = sub["entropy_norm"].std()
			min_cluster_size = sub["cluster_size"].min()
			max_cluster_size = sub["cluster_size"].max()
			n_clusters = sub.shape[0]
			n_small_clusters = (sub["cluster_size"] < 3).sum()

			rows.append({
				"k": k,
				"n_clusters": n_clusters,
				"weighted_mean_hnorm": weighted_mean_entropy,
				"mean_hnorm": mean_entropy,
				"std_hnorm": std_entropy,
				"min_cluster_size": min_cluster_size,
				"max_cluster_size": max_cluster_size,
				"n_small_clusters_lt3": n_small_clusters
			})

		return pd.DataFrame(rows).sort_values("weighted_mean_hnorm", ascending=True)	

	def entropy_analysis_for_primary_site(self, cluster_type:str, primary_site:str, Kmin:int=2, Kmax:int=10, 
									      min_barcodes:int=2, min_genes:int=2,
							    		  verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
		
		_, _, df_all_mut, _ = self.get_filtered_tables(primary_site=primary_site, verbose=verbose)

		dfempty = pd.DataFrame()

		if df_all_mut.empty:
			return dfempty, dfempty, dfempty, dfempty, dfempty

		dfpiv = self.build_pivot_table(df_all_mut, min_barcodes=min_barcodes, min_genes=min_genes)
		self.dfpiv = dfpiv

		if dfpiv.shape[0] < 3 or dfpiv.shape[1] < 3:
			return dfempty, dfempty, dfempty, dfpiv, df_all_mut
		
		if Kmax >= dfpiv.shape[0]:
			Kmax = dfpiv.shape[0] - 1

		if Kmax <= 3:
			return dfempty, dfempty, dfempty, dfpiv, df_all_mut

		df_list = []
		for k in range(Kmin, Kmax + 1):
			if cluster_type == 'UMAP':
				_, labels = self.calc_UMAP(dfpiv, k)
			elif cluster_type == 'HDBSCAN':
				_, labels, _ = self.calc_HDBSCAN(dfpiv, k)
			else:
				raise Exception(f"\n---------- Define the cluster_type like UMAP or HDBSCAN, got: {cluster_type}")
			
			if labels is None or len(labels) == 0:
				continue

			for cluster in np.unique(labels):

				dfc = self.cluster_mutation_table(
					dfpiv=dfpiv,
					labels=labels,
					cluster=cluster,
					min_barcodes=min_barcodes
				)

				n_cluster = dfc.shape[0]

				if n_cluster < 3 or dfc.shape[1] < 3:
					continue

				gene_degree = dfc.sum(axis=0).sort_values(ascending=False)
				gene_freq = (gene_degree / n_cluster).sort_values(ascending=False)

				df_cluster_stat = pd.DataFrame(
				{	"k": k,
					"cluster": cluster,
					"gene": gene_degree.index,
					"degree": gene_degree.values,
					"cluster_size": n_cluster,
					"freq": gene_freq.values
				})
				df_list.append(df_cluster_stat)

		if len(df_list) == 0:
			return dfempty, dfempty, dfempty, dfpiv, df_all_mut

		dfstat = pd.concat(df_list, ignore_index=True)

		dfh = self.calc_shannon_entropy_from_dfstat(dfstat)

		dfw = self.score_k_from_entropy_table(dfh)

		cols = dfw.columns.to_list()

		dfw['psi_id'] = self.psi_id
		dfw['primary_site'] = self.primary_site
		dfw['min_barcodes'] = min_barcodes
		dfw['min_genes'] = min_genes

		cols = ['psi_id', 'primary_site', 'min_barcodes', 'min_genes'] + cols
		dfw = dfw[cols]

		return dfw, dfh, dfstat, dfpiv, df_all_mut
	

	def buid_purity_table(self, dfpiv: pd.DataFrame, labels:list) -> pd.DataFrame:
		lab_list = np.unique(labels)
		pu_list = []
		n_list = []
		pairs = []

		for label in lab_list:
			idx = labels == label
			X = dfpiv[idx].to_numpy(dtype=bool)

			n_bardodes = len(X)
			n_list.append(n_bardodes)
			pairs.append(n_bardodes * (n_bardodes - 1) / 2)

			if len(X) > 1:
				dist = pairwise_distances(X, metric="jaccard")
				similarity = 1 - dist
				purity = similarity[np.triu_indices_from(similarity, k=1)].mean()
			else:
				purity = 0

			pu_list.append(np.round(purity, 3))


		dfa = pd.DataFrame({
			"label": lab_list,
			"n_barcodes": n_list,
			"purity": pu_list,
			"n_pairs": pairs
		})

		max_pairs = dfa["n_pairs"].max()
		dfa["purity_norm"] = dfa["purity"] * (dfa["n_pairs"] / max_pairs)

		dfa = dfa.sort_values(by="purity_norm", ascending=False)

		return dfa
	

	def plot_purity(self, dfpur: pd.DataFrame, dfclu:pd.DataFrame, good_clusters: list, min_perc: float = 0.10):
		
		# top_n_genes = 30

		ngood = len(good_clusters)
		nrows = int(np.ceil(ngood / 2))
		height = 6*nrows
			
		fig, axes = plt.subplots(nrows, 2, figsize=(12, height))
		axes = axes.flatten()

		for ax, label in zip(axes, good_clusters):
			if label == -1:
				continue
			
			dfb = dfclu[label]

			dfb = dfb[ dfb.values > min_perc]
			dfb = dfb.sort_values(ascending=False)
			# dfb = dfb.sort_values(ascending=False).head(top_n_genes)
			
			ax.bar(dfb.index, dfb.values)

			stri = f"Label {label} | purity_norm={dfpur.loc[dfpur['label']==label].iloc[0].purity_norm:.3f}"

			ax.set_title(stri)
			ax.set_ylabel("Representative percentage")
			ax.set_xlabel("Genes")
			ax.tick_params(axis="x", rotation=70)

			print(stri, ", ".join(dfb.index.to_list()))

		plt.tight_layout()


		return fig


	def enrichment_test(self, sample_genes, subtype_genes, background_genes):
		N = len(background_genes)
		K = len(subtype_genes)
		n = len(sample_genes)
		overlap_genes = set(sample_genes) & set(subtype_genes)
		overlap = len(overlap_genes)

		# P(X >= k)
		pval = hypergeom.sf(overlap - 1, N, K, n)

		return pval, overlap_genes




	def cluster_analysis(self, cluster_type:str, primary_site:str, k:int=5, Kmin:int=2, Kmax:int=10, 
						 min_barcodes:int=3, min_genes:int=5, pur_threshold:float = 0.05, min_represent_perc = 0.10,
						 verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, 
									                  pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
		
		dfw, dfh, dfstat, dfpiv, df_all_mut = self.entropy_analysis_for_primary_site(cluster_type,  primary_site,
																		 Kmin, Kmax, 
																		 min_barcodes, min_genes, verbose)

		dfempty = pd.DataFrame()

		if dfpiv.empty:
			print("Did not define the pivot table")
			return dfempty, dfempty, dfempty, dfempty, dfempty, dfempty, dfempty, dfempty

		if cluster_type == 'UMAP':
			print(f"Chose {cluster_type} with k={k}")
			_, labels = self.calc_UMAP(dfpiv, k)

		elif cluster_type == 'HDBSCAN':
			min_cluster_size=5
			min_samples=3
			print(f"Chose {cluster_type} with min_cluster_size={min_cluster_size} and min_samples={min_samples}")

			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				_, labels, _ = self.calc_HDBSCAN(dfpiv, min_cluster_size=min_cluster_size, min_samples=min_samples)
		else:
			print("Did not defined the clustering method")
			return dfempty, dfempty, dfempty, dfw, dfh, dfstat, dfpiv, df_all_mut

		#----------- cluster ----------------
		dfpiv2 = dfpiv.copy()
		dfpiv2["cluster"] = labels
		dfclu = dfpiv2.groupby("cluster").mean().T

		dfpur = self.buid_purity_table(dfpiv, labels)
		if dfpur.empty:
			return dfempty, dfempty, dfempty, dfw, dfh, dfstat, dfpiv, df_all_mut
	
		good_clusters = dfpur.loc[dfpur["purity_norm"] >= pur_threshold, "label"]

		#------------ hypergeometric statistics ---------------------
		background_genes = np.unique(df_all_mut.symbol.to_list())

		dic	= self.SUBTYPE_GENES.get(self.psi_id, {})

		if dic=={}:
			print(f"No subtype genes found for PSI ID: {self.psi_id}")
			return dfempty, dfpur, dfclu, dfw, dfh, dfstat, dfpiv, df_all_mut
		

		lista = []
		for label in good_clusters:
			
			purity_norm = dfpur.loc[dfpur['label']==label].iloc[0].purity_norm
					
			dfb = dfclu[label]

			n_barcodes = len(dfpiv2[dfpiv2["cluster"] == label])
			
			dfb = dfb[ dfb.values > min_represent_perc]
			# dfb = dfb.sort_values(ascending=False)

			sample_genes = dfb.index.to_list()
		
			for subtype, annotated_genes in dic.items():
				pval, overlap_genes = self.enrichment_test(sample_genes, annotated_genes, background_genes)
				# print(f"Subtype: {subtype}, overlap: {overlap}, p-value: {pval}")

				overlap = len(overlap_genes)

				if overlap >= 2:
					mat = [cluster_type, k, n_barcodes, label, purity_norm, subtype, overlap, len(sample_genes), len(annotated_genes), len(background_genes), pval, overlap_genes]
					lista.append(mat)


		if lista == []:
			df = pd.DataFrame()
		else:
			df = pd.DataFrame(lista, columns=["cluster_type", "k", "n_barcodes", "label", "purity_norm", "subtype", "overlap", "sample_genes", "annotated_genes", "background_genes", "pval", "overlap_genes"])
			df['fdr'] = fdr(df['pval'])
		
		return df, dfpur, dfclu, dfw, dfh, dfstat, dfpiv, df_all_mut



	def get_VCF_files(self, subtype_global:str, tumor_class:str, subtype_tissue:str, 
					batch_cases:int=20, batch_size:int=20, timeout:int=100,
					force:bool=False, verbose:bool=False) -> pd.DataFrame:

		df_vcf = pd.DataFrame()
		self.df_vcf = df_vcf

		self.set_s_case(subtype_global, tumor_class, subtype_tissue)

		df_cases, _, _ = self.get_cases_and_subtypes(batch_size=batch_size, do_filter=False, 
													force=False, verbose=verbose)

		self.df_cases = df_cases

		if df_cases is None or df_cases.empty:
			print(f"No cases found while searching for '{self.s_case}'")
			return df_vcf
			
		case_id_list = list(df_cases.case_id)
		case_id_list.sort()

		N_cases = len(case_id_list)
		print(f">>> {N_cases} cases")

		self.fname_vcf_files0 = 'vcf_files_for_%s.tsv'
		fname = self.fname_vcf_files0%(self.s_case)
		fname = title_replace(fname)
		filename = os.path.join(self.root_psi, fname)

		if os.path.exists(filename) and not force:
			df_vcf = pdreadcsv(fname, self.root_psi, verbose=verbose)
			self.df_vcf = df_vcf

			return df_vcf
		#-------------------------- batch loop ---------------------------
		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None

		ini = -batch_cases
		end = 0
		res = None

		while(True):
			ini += batch_cases
			end += batch_cases

			if ini >= N_cases:
				break
			
			if end > N_cases:
				end = N_cases

			print(f"{ini}-{end} ", end='')

			lista = case_id_list[ini: end]

			filters = {
				"op": "and",
				"content": [
					{"op": "in",
						"content": {"field": "cases.case_id", "value": lista}
					},
					{"op": "=",
						"content": {"field": "files.data_format", "value": "VCF"}
					},
					{"op": "in",
						"content": {"field": "files.data_type",
									"value": [
										"Raw Simple Somatic Mutation",
										"Masked Somatic Mutation"
									]
						}
					}
				]
			}
		
			try:
				while True:
					print(".", end='')

					params = {
						"filters": json.dumps(filters),
						"fields": ",".join([
							"file_id",
							"file_name",
							"data_type",
							"data_format",
							"analysis.workflow_type",
							"cases.case_id",
							"cases.submitter_id",
							"cases.samples.sample_id",
							"cases.samples.submitter_id",
							"cases.samples.sample_type",
						]),
						"format": "JSON",
						"size": size_,
						"from": from_
					}
					
					res = requests.get(self.url_gdc_files, params=params, timeout=timeout)
					response = res.json()

					if 'data' not in response.keys():
						print(f"No data found while searching for '{self.psi_id}' cases {case_id_list}")
						print(">>> response", response)
						return df_vcf
				
					hits = response.get("data", {}).get("hits", [])

					if total is None:
						total = response["data"]["pagination"]["total"]
					
					if not hits:
						break

					all_hits.extend(hits)
					from_ += size_

				print("\n")

				if all_hits == []:
					print(f"No files were found for {self.psi_id} and {len(case_id_list)} cases")
					return df_vcf

				#------------ lost data? ------------------
				N = len(all_hits)

				if N < total:
					print(f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} ")
				else:
					if verbose: print(f"👉 Returned {N} / Total paginated {total}")

				#------------ having all hits -------------
				records = []

				for hit in all_hits:
					for case in hit.get("cases", []):
						for sample in case.get("samples", []):
							records.append({
								"case_id": case["case_id"],
								"submitter_id": case["submitter_id"],
								"sample_id":   sample["sample_id"],
								"sample_type": sample["sample_type"],
								"barcode_sample":  sample["submitter_id"],
								"file_id":    hit["file_id"],
								"file_name":  hit["file_name"],
								"data_type":  hit["data_type"],
								"data_format": hit["data_format"],
								"workflow_type": hit["analysis"]["workflow_type"],
							})
	
				df_vcf = pd.DataFrame(records)
				cols = list(df_vcf.columns)

				# 🔹 Metadata 
				df_vcf['psi_id'] = self.psi_id
				df_vcf['subtype_global'] = subtype_global
				df_vcf['tumor_class'] = tumor_class
				df_vcf['subtype_tissue'] = subtype_tissue

				cols = ["psi_id", "subtype_global", "tumor_class", "subtype_tissue"] + cols

				df_vcf = df_vcf.sort_values(["case_id", "sample_type"], ascending=[False,False]).reset_index(drop=True)
				df_vcf.reset_index(drop=True, inplace=True)

				_ = pdwritecsv(df_vcf, fname, self.root_psi, verbose=verbose)

			except Exception as e:
				print(f"Error for searching files for {self.s_case}'. error: {e}")
				print("Trye to diminish the batch size - like ~20")

				try:
					res.raise_for_status()   # raises if HTTP 4xx/5xx
					print("status:", res.status_code)
					print("content-type:", res.headers.get("Content-Type"))
					print("text:", res.text[:500])   # inspect response before parsing
				except:
					print("Could not access http response.")
								
				return df_vcf

		self.df_vcf = df_vcf

		return df_vcf
