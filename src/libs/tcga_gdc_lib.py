#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

import glob
import os, requests, json, re
from tabnanny import verbose
import pandas as pd
from collections import Counter
from pathlib import Path

from setuptools import glob
from typing import List, Tuple, Any, Iterable, Optional

from libs.Basic import *

class GDC(object):
	def __init__(self, root0:Path=Path('../data/')):
		
		self.url_gdc_project = "https://api.gdc.cancer.gov/projects"
		self.url_gdc_cases = "https://api.gdc.cancer.gov/cases"
		self.url_gdc_files = "https://api.gdc.cancer.gov/files"
		self.url_gdc_data  = "https://api.gdc.cancer.gov/data/%s"

		self.url_cbioportal = "https://www.cbioportal.org/api"

		self.prog_id, self.psi_id = '', ''

		self.root0 = Path(root0)

		# root_data will be: ../data/TCGA
		self.root_data = ''
		self.root_summary = ''
		self.root_psi = ''

		self.clean_gdc_files()

		self.fname_all_cases = '%s_summ_cases.tsv'
		self.fname_all_samples = '%s_summ_samples.tsv'
		self.fname_all_mutations = '%s_summ_mutations.tsv'

	def clean_gdc_files(self):

		self.gdc_file_name = ''
		self.gdc_file_id = ''
		self.gdc_data_type = ''

		self.s_case = ''

		self.fname_programs = 'gdc_programs.txt'

		# primary_site
		self.fname_prim_site = 'primary_site_program_%s.tsv'
		self.fname_cases0    = 'cases_for_%s.tsv'
		self.fname_subtype0  = 'subtype_for_%s.tsv'
		self.fname_samples0  = 'samples_for_%s.tsv'
		self.fname_rnaseq_exp_files = 'rnaseq_exp_files_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.fname_cases_deprecated = 'cases_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.fname_fileid   = '%s_for_%s_case_%s_sample_type_%s_stage_%s_fileid_%s.tsv'
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

		df_cases, _, _ = self.get_cases_and_subtypes(batch_size=200, do_filter=False, 
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
						print(f"No data found while searching for '{psi_id}' cases {case_id_list}")
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

	def get_table_given_fileID(self, psi_id:str, case_id:str, 
								sample_type:str, stage:str,
								file_type:str, file_id:str,
								force:bool=False, verbose:bool=False) -> pd.DataFrame:
		"""
		Retrieve table like: RNA or Proteomic expression
		input: file_id
		output: dataframe
		"""

		if not file_id and not isinstance(file_id, str):
			print(f"No file_id defined.")

			self.df_table = pd.DataFrame()
			return self.df_table

		fname = self.fname_fileid%(file_type, psi_id, case_id, sample_type, stage, file_id)
		fname = title_replace(fname)
		filename = os.path.join(self.root_psi, fname)

		if os.path.exists(filename) and not force:
			df_table = pdreadcsv(fname, self.root_psi, verbose=verbose)
			self.df_table = df_table

			return df_table

		is_expression = True if file_type == 'Gene Expression Quantification' else False

		if verbose: print("Searching: ", end='')
		try:
			url_file = f"https://api.gdc.cancer.gov/data/{file_id}"

			res = requests.get(url_file)
			data = res.content

			with open(filename, "wb") as f:
				f.write(data)

			if is_expression:
				df_table = pd.read_csv(filename, sep="\t", comment="#")
				df_table = self.clean_expression_table(df_table)

				cols = ["gene_id", "symbol", "gene_type", "unstranded", "counts", "stranded_second",
						"tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"]
				df_table = df_table[cols]

				_ = pdwritecsv(df_table, fname, self.root_psi, verbose=verbose)
			else:
				df_table = pdreadcsv(fname, self.root_psi, verbose=verbose)

			if verbose: print(f"found {len(df_table)} files in selected samples.")

		except Exception as e:
			print(f"Error for '{psi_id}', '{file_type}', and {file_id}: {e}")
			self.df_table = pd.DataFrame()
			return self.df_table

		self.df_table = df_table

		return df_table
	
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


	def get_case_uuid(self, barcode:str) -> str:

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


	def get_files(self, case_uuid:str, data_type:str="Gene Expression Quantification") -> List:
		lista = self.get_expression_files_by_uuid(case_uuid, data_type="Gene Expression Quantification")
		return lista

	def get_expression_files_by_uuid(self, case_uuid:str, data_type:str="Gene Expression Quantification") -> List:

		if not isinstance(case_uuid, str) or len(case_uuid) < 3:
			print(f"UUID bad formated {case_uuid}.")
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
						"value": [case_uuid]
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
			print(f"No data found for {case_uuid}. error: {e}")
			return []

		response = data["data"]["hits"]

		try:
			dic = response[0]
			if isinstance(dic, str):
				dic = eval(dic)

			self.build_gdc_filenames(dic)

		except Exception as e:
			print(f"No response. error: {e}")

		return response
	
	def build_gdc_filenames(self, dic:dict, exp_unit:str='TPM'):

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

		fname_all_cases = self.fname_all_cases%(self.psi_id)
		filename_cases = os.path.join(self.root_summary, fname_all_cases)

		fname_all_samples = self.fname_all_samples%(self.psi_id)
		filename_samples = os.path.join(self.root_summary, fname_all_samples)

		fname_all_mutations = self.fname_all_mutations%(self.psi_id)
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