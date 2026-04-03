#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

import glob
import os, requests, json, re
import pandas as pd
from collections import Counter

from setuptools import glob
from typing import List, Tuple, Any, Iterable, Optional

from libs.Basic import *

class GDC(object):
	def __init__(self, root_data:str='../data/'):
		
		self.url_gdc_project = "https://api.gdc.cancer.gov/projects"
		self.url_gdc_cases = "https://api.gdc.cancer.gov/cases"
		self.url_gdc_files = "https://api.gdc.cancer.gov/files"
		self.url_gdc_data  = f"https://api.gdc.cancer.gov/data/%s"

		self.url_cbioportal = "https://www.cbioportal.org/api"

		self.root_data = root_data

		self.clean_gdc_files()

	def clean_gdc_files(self):
		self.gdc_file_name = ''
		self.gdc_file_id = ''
		self.gdc_data_type = ''

		self.fname_programs = 'programs.txt'

		# primary_site
		self.fname_primary_site = 'primary_site_program_%s.tsv'
		self.fname_cases        = 'cases_for_%s.tsv'
		self.fname_subtype      = 'subtype_for_%s.tsv'
		self.fname_stage        = 'stage_for_%s_subtype_%s.tsv'
		self.fname_pid_samples  = 'samples_for_%s.tsv'
		self.fname_rnaseq_exp_files = 'rnaseq_exp_files_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.fname_cases_deprecated = 'cases_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.fname_fileid   = '%s_for_%s_case_%s_sample_type_%s_stage_%s_fileid_%s.tsv'
		self.fname_mutation = 'mutations_for_study_%s_samples_%s_to_%s.tsv'
		self.fname_extmut   = 'mutations_all_fields_for_study_%s_samples_%s_to_%s.tsv'

		self.gdc_fname = ''
		self.gdc_filename = ''
		self.gdc_ouptut_fname = ''
		self.gdc_ouptut_filename = ''

		self.exp_unit = ""
		self.value_col = ""

		# program, primary site, subtype, stage, case_id, samples
		self.df_ps, self.df_subt, self.df_cases = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
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

	def map_tissue_subtype(self, global_subtype:str, pid:str) -> str:
		if pid in self.SITE_MAP:
			return self.SITE_MAP[pid].get(global_subtype, global_subtype)
	
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

	def set_cancer_type(self, case:str="Breast Cancer"):
		self.case = case

		self.root_case = os.path.join(self.root_data, self.dir_case)
		os.makedirs(self.root_case, exist_ok=True)

		self.dir_case = case.lower().replace(' ', '_')

		self.clean_gdc_files()

	def get_gdc_progams(self, force:bool=False, verbose:bool=False) -> List:

		filename = os.path.join(self.root_data, self.fname_programs)

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


	def get_primary_sites(self, program:str='TCGA', force:bool=False, verbose:bool=False) -> pd.DataFrame:

		fname = self.fname_primary_site%(program)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_ps = pdreadcsv(fname, self.root_data, verbose=verbose)
			self.df_ps = df_ps

			return df_ps

		filters = {
					"op": "in",
					"content": { "field": "program.name", "value": [program] }
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
				print(f"No data found while searching for '{program}'")
				print(">>> response", response)
				self.df_ps = pd.DataFrame()
				return self.df_ps

			hits = response["data"]["hits"]

			df_ps = pd.DataFrame(hits)
			# fix list columns
			for col in df_ps.columns:
				df_ps[col] = df_ps[col].apply(
					lambda x: ", ".join(x) if isinstance(x, list) else x  )

			df_ps = df_ps.rename(columns={"id": "pid"})

			df_ps = df_ps.sort_values(["primary_site", "disease_type"])

			_ = pdwritecsv(df_ps, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"Error searching for '{program}': {e}")
			print(">>> response", response)
			self.df_ps = pd.DataFrame()
			return self.df_ps

		self.df_ps = df_ps
	
		return df_ps
	
	def list_disease_types(self, pid = 'TCGA-BLCA') -> List:
	
		try:
			row = self.df_ps[self.df_ps.project_id == pid].iloc[0]
			deas_type_list = row.disease_type

			if isinstance(deas_type_list, str):
				deas_type_list = eval(deas_type_list)
		except:
			print("No disease types were found.")
			deas_type_list = []

		return deas_type_list

	def get_cases_and_subtypes(self, pid:str, batch_size:int=200, 
							   do_filter:bool=True, debug:bool=False, 
				    		   force:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
		'''
		calc all subtypes, given and pid --> df_cases
		group by ["subtype_global", "subtype_tissue", "stage"] --> df_subt

		filter: NOS (Not Otherwise Specified) → the pathologist could not (or did not) assign a more specific subtype.
		e.g.: "Yes, it's an adenocarcinoma — but we don’t have finer classification"

		input: pid = primary site ID
		output: df_cases, df_subt, df_prof
		'''

		def apply_filter(df_cases: pd.DataFrame) -> pd.DataFrame:
			df_cases = df_cases[df_cases.validity == 'valid'].copy()
			df_cases = df_cases[df_cases["consistency"] == "ok"]
			df_cases.reset_index(drop=True, inplace=True)
			
			# frac_threshold:float=0.01,
			# df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()
			# df_cases = df_cases[df_cases["frac"] > frac_threshold]
			# df_cases.reset_index(drop=True, inplace=True)
			# df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()

			return df_cases

		fname_cases = self.fname_cases%(pid)
		filename_cases = os.path.join(self.root_data, fname_cases)


		fname_subt = self.fname_subtype%(pid)
		filename_subt = os.path.join(self.root_data, fname_subt)

		if os.path.exists(filename_cases) and os.path.exists(filename_subt) and not force:
			df_cases = pdreadcsv(fname_cases, self.root_data, verbose=verbose)
			self.df_cases = df_cases

			if do_filter:
				df_cases = apply_filter(df_cases)

			df_subt = self.groupby_state(df_cases)
			df_prof = self.build_profile(df_cases)

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
			for pid=='TCGA-ACC' if subtype_global = other --> change to adrenal_cortical_carcinoma
			"""
			df["subtype_global"] = [df.iloc[i]["tumor_class"] 
						            if  (df.iloc[i]["pid"] =='TCGA-ACC' and df.iloc[i]["subtype_global"] =='other')
						            else df.iloc[i]["subtype_global"] for i in range(len(df))]

			# histology
			df["histology"] = df["subtype_global"].apply(self.map_histology)

			# tissue-specific subtype
			df["subtype_tissue"] = df.apply(
				lambda r: self.map_tissue_subtype(r["subtype_global"], r["pid"]),
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
				"value": [pid]
			}
		}

		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None
		df_cases = pd.DataFrame()

		# print("Searching: ", end='')
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
					print(f"No data found while searching for '{pid}'")
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
				print(f"No subtypes found for {pid} ")
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

			df_cases = df_cases.rename(columns={"project.project_id": "pid"})

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
			df_subt = self.groupby_state(df_cases)

			_ = pdwritecsv(df_cases, fname_cases, self.root_data, verbose=verbose)
			_ = pdwritecsv(df_subt,  fname_subt, self.root_data, verbose=verbose)


		except Exception as e:
			print(f"Error for searching diags for '{pid}'. error: {e}")
			self.df_cases = df_cases
			self.df_subt  = pd.DataFrame()
			self.df_prof  = pd.DataFrame()
			return self.df_cases, self.df_subt, self.df_prof

		if do_filter:
			df_cases = apply_filter(df_cases)

		self.df_cases = df_cases
		self.df_subt  = df_subt
		
		df_prof = self.build_profile(df_cases)

		return df_cases, df_subt, df_prof
	

	def group_file_types(self, df_samples:pd.DataFrame) -> pd.DataFrame:
		dic =  Counter(df_samples.data_type)
		dfu = pd.DataFrame(dic.items(), columns=["data_type", "n"])
		dfu = dfu.sort_values('n', ascending=False).reset_index(drop=True)
		return dfu

	def groupby_sstate(self, df_cases:pd.DataFrame):

		df_cases["sstage"] = df_cases["stage"].map(lambda x: self.simplify_stage(x))
		df_subt = df_cases.groupby(["pid", "subtype_global", "tumor_class", "subtype_tissue", "sstage"], dropna=False).size().reset_index(name="n")
		df_subt = df_subt.sort_values("n", ascending=False).reset_index(drop=True)

		self.df_subt = df_subt

		return df_subt


	def groupby_state(self, df_cases:pd.DataFrame):

		# df_subt = df_cases[cols].copy().drop_duplicates()
		# df_subt = df_subt.sort_values(cols).reset_index(drop=True)		

		cols = ['pid', 'subtype_global', 'tumor_class', 'subtype_tissue', 'stage']
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

	def get_samples_for_pid_subtypes(self, pid:str, subtype_global:str, tumor_class:str, 
									 subtype_tissue:str, 
									 batch_cases:int=5, batch_size:int=200, 
									 force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		return all samples given a list of cases
		for pid, subtype_global, tumor_class, subtype_tissue

		input: pid, subtype_global, tumor_class, subtype_tissue
		output: df_samples
		'''
		self.df_samples = pd.DataFrame()

		df_cases, _, _ = self.get_cases_and_subtypes(pid=pid, batch_size=200, do_filter=False, 
											         force=False, verbose=verbose)
	
		self.df_cases = df_cases


		s_case = f"{pid} subtype '{subtype_global}' tumor '{tumor_class}' subtype_tissue '{subtype_tissue}'"

		print(">>>", s_case)

		if df_cases is None or df_cases.empty:
			print(f"No cases found while searching for '{s_case}'")
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
			print(f"No cases found for {s_case}")
			return self.df_samples

		fname = self.fname_pid_samples%(s_case)
		fname = title_replace(fname)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_samples = pdreadcsv(fname, self.root_data, verbose=verbose)
			self.df_samples = df_samples

			return df_samples

		case_id_list = list(df_cases.case_id)
		case_id_list.sort()

		s_case_id_list3 = f"[{','.join(case_id_list[:3])}]"

		N_cases = len(case_id_list)
		print(f">>> {N_cases} cases", s_case_id_list3, ".....")


		#-------------------------- batch loop ---------------------------
		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None
		df_samples = pd.DataFrame()

		print("Searching: ", end='')

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
			print("\n>>>", len(lista), lista)

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
						print(f"No data found while searching for '{pid}' cases {case_id_list}")
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

				print("\n")

				if all_hits == []:
					print(f"No files were found for {pid} cases {case_id_list}")
					self.df_samples = pd.DataFrame()
					return self.df_samples
		
				#------------ lost data? ------------------
				N = len(all_hits)

				if N < total:
					print(f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} ")
				else:
					print(f"👉 Returned {N} / Total paginated {total}")

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
								"barcode_id":  sample["submitter_id"],
								"file_id":    hit["file_id"],
								"file_name":  hit["file_name"],
								"data_type":  hit["data_type"],
								"data_format": hit["data_format"],
							})

				df_samples = pd.DataFrame(records)
				self.df_samples = df_samples
				cols = list(df_samples.columns)

				# 🔹 Metadata 
				df_samples['pid'] = pid
				df_samples['subtype_global'] = subtype_global
				df_samples['tumor_class'] = tumor_class
				df_samples['subtype_tissue'] = subtype_tissue
				df_samples['stage'] = self.get_stage_from_cases(df_samples.case_id.tolist())

				cols = ["pid", "subtype_global", "tumor_class", "subtype_tissue", "stage"] + cols

				df_samples = df_samples.sort_values(["case_id", "sample_type"], ascending=[False,False]).reset_index(drop=True)
				df_samples.reset_index(drop=True, inplace=True)

				_ = pdwritecsv(df_samples, fname, self.root_data, verbose=verbose)

			except Exception as e:
				print(f"Error for searching files for {s_case}'. error: {e}")
				self.df_samples = df_samples
				return df_samples

		self.df_samples = df_samples

		return df_samples

	def get_table_given_fileID(self, pid:str, case_id:str, 
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

		fname = self.fname_fileid%(file_type, pid, case_id, sample_type, stage, file_id)
		fname = title_replace(fname)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_table = pdreadcsv(fname, self.root_data, verbose=verbose)
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

				_ = pdwritecsv(df_table, fname, self.root_data, verbose=verbose)
			else:
				df_table = pdreadcsv(fname, self.root_data, verbose=verbose)

			if verbose: print(f"found {len(df_table)} files in selected samples.")

		except Exception as e:
			print(f"Error for '{pid}', '{file_type}', and {file_id}: {e}")
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


	def get_table_searching_for_fileID(self, pid:str, data_type:str, sample_type:str, file_id:str,  verbose:bool=False) -> pd.DataFrame:
			
		data_type2 = title_replace(data_type)
		sample_type2 = title_replace(sample_type)

		files = [x for x in os.listdir(self.root_data) if file_id in x and data_type2 in x and sample_type2 in x]

		if len(files) == 0:
			print(f"No files found for {file_id}.")
			self.df_table = pd.DataFrame()
			return self.df_table

		if len(files) > 1:
			print(f"Multiple files found for {file_id}. Using the first one.")

		fname = files[0]
		df_table = pdreadcsv(fname, self.root_data, verbose=verbose)
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

		case_id_list = np.unique(df2.case_id.to_list())
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

	def merge_normal_tumor_tables(self, pid:str, df_normal:pd.DataFrame, df_tumor:pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

		cols = ["gene_id", "symbol", "gene_type", "counts"]
		common_cols = ["gene_id", "symbol", "gene_type"]


		dfa_normal = pd.DataFrame()

		for i, row in df_normal.iterrows():
			data_type = row.data_type
			sample_type = row.sample_type
			file_id = row.file_id
			
			dft = self.get_table_searching_for_fileID(pid=pid, data_type=data_type, sample_type=sample_type, file_id=file_id, verbose=False)

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
			
			dft = self.get_table_searching_for_fileID(pid=pid, data_type=data_type, sample_type=sample_type, file_id=file_id, verbose=False)

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

		self.clean_gdc_files()

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

		self.gdc_filename = os.path.join(self.root_case, fname)

		if exp_unit == 'TPM':
			self.exp_unit = exp_unit
			self.value_col ="tpm_unstranded"
		else:
			self.exp_unit = "???"
			self.value_col ="???"
			raise Exception("Error in which count col, define as TPM.")


		self.gdc_ouptut_fname = f"{fname.replace('.dat', '')}_{exp_unit}.tsv"
		self.gdc_ouptut_filename = os.path.join(self.root_case, self.gdc_ouptut_fname)


	def get_mutations_from_samples(self, sample_ids: Iterable[str], study_id: str,
		session: Optional[requests.Session] = None, timeout: int = 60,
		verbose: bool = False) -> pd.DataFrame:
		"""
		Fetch mutation records from cBioPortal for a list of sample IDs.

		Parameters
		----------
		sample_ids : Iterable[str]
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

		molecular_profile_id = f"{study_id}_mutations"

		http = session or requests.Session()

		url = f"{self.url_cbioportal}/molecular-profiles/{molecular_profile_id}/mutations/fetch"
		print(">>>", url)

		payload = {
			"sampleIds": sample_ids
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
			return pd.DataFrame()

		df = pd.DataFrame(data)

		cols = ['sampleId', 'patientId', 'studyId', 'molecularProfileId',
				'entrezGeneId', 'keyword', 'proteinChange', 'mutationType',
				'mutationStatus', 'center', 'tumorRefCount', 'normalAltCount',
				'normalRefCount', 'variantType', 'chr', 'startPosition',
				'endPosition', 'referenceAllele', 'uniqueSampleKey',
				'uniquePatientKey', 'validationStatus', 'tumorAltCount',
				'ncbiBuild', 'variantAllele', 'refseqMrnaId', 'proteinPosStart',
				'proteinPosEnd']


		# the selected cols + others not listed
		# cols = [c for c in cols if c in df.columns.to_list()] + [c for c in df.columns if c not in cols]
		df = df[cols]

		df['keyword'] = [x.split(' ')[0] if isinstance(x, str) else x for x in df['keyword']]


		rename_cols = ['sample_id', 'barcode', 'pid', 'mol_profile_id',
					'entrez_gene_id', 'symbol', 'protein_mut', 'mutation_type',
					'mutation_status', 'center', 'tumor_ref_count', 'normal_alt_count',
					'normal_ref_count', 'variant_type', 'chr', 'start',
					'end', 'ref_allele', 'unique_sample_key',
					'unique_patient_key', 'validation_status', 'tumor_alt_count',
					'ncbi_build', 'variant_allele', 'refseq_mrna_id', 'protein_pos_start', 'protein_pos_end']

		df.columns = rename_cols


		order_cols = ['sample_id', 'barcode', 'pid', 'mol_profile_id',
				'symbol', 'refseq_mrna_id', 'entrez_gene_id', 
				'protein_mut', 'mutation_type', 'mutation_status',
				'ref_allele', 'variant_allele', 'variant_type', 
				'chr', 'start', 'end', 
				'validation_status', 'protein_pos_start', 'protein_pos_end', 'tumor_alt_count',
				'ncbi_build',  'center', 'tumor_ref_count', 'normal_alt_count', 'normal_ref_count', 'unique_sample_key',
				'unique_patient_key']

		df = df[order_cols]

		return df

	
	def get_df_mut_transform_mutation_table(self, sample_ids: Iterable[str], study_id: str,
		session: Optional[requests.Session] = None, timeout: int = 60,
		force: bool = False, verbose: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:

		sample_ids = [str(x).strip() for x in sample_ids if str(x).strip()]
		if not sample_ids:
			raise ValueError("sample_ids is empty.")
		
		if len(sample_ids[0].split('-')[-1]) == 3:
			sample_ids = [x[:-1] for x in sample_ids]

		sample_ids = list(sample_ids)
		sample_ids.sort()

		sample_ids_ini = sample_ids[0]
		sample_ids_end = sample_ids[1] if len(sample_ids) > 1 else sample_ids[0]

		if study_id[0].isupper():
			mat = study_id.lower().split('-')
			# cBioPortal disease - tcga 
			study_id = mat[1] + '_' + mat[0]

		if study_id == "luad_tcga":
			study_id = "luad_tcga_pan_can_atlas_2018"

		if verbose:
			print(f">>> {study_id} len = {len(sample_ids)} - {sample_ids[:5]}...")


		fname_mut = self.fname_mutation%(study_id, sample_ids_ini, sample_ids_end)
		fname_mut = title_replace(fname_mut)
		filename_mutation = os.path.join(self.root_data, fname_mut)

		fname_extmut = self.fname_extmut%(study_id, sample_ids_ini, sample_ids_end)
		fname_extmut = title_replace(fname_extmut)
		filename_extmut = os.path.join(self.root_data, fname_extmut)


		if os.path.exists(filename_mutation) and os.path.exists(filename_extmut) and not force:
			dff    = pdreadcsv(fname_mut, self.root_data, verbose=verbose)
			df_mut = pdreadcsv(fname_mut, self.root_data, verbose=verbose)
			return dff, df_mut
		
		'''
			df_mut cols: ["sample_id", "barcode", "pid", "mol_profile_id","gene",
			"entrez_gene_id", "protein_mut", "mutation_type", "mutation_status",
			"variant_type", "chr", "start", "end",
			"ref_allele", "tumor_seq_allele"]		
		'''
		df_mut = self.get_mutations_from_samples(sample_ids=sample_ids, study_id=study_id,
											 	session=session, timeout=timeout, verbose=verbose)
		
		if df_mut.empty:
			print("No mutation data found for these samples..")
			return pd.DataFrame(), pd.DataFrame()

		#--------------- map main cols from df_mut ------------------------
		"""
		order_cols = ['sampleId', 'barcode', 'pid', 'mol_profile_id',
			'symbol', 'refseq_mrna_id', 'entrez_gene_id', 
			'protein_mut', 'mutation_type', 'mutation_status',
			'ref_allele', 'variant_allele', 'variant_type', 
			'chr', 'start', 'end', 
			'validation_status', 'protein_pos_start', 'protein_pos_end', 'tumor_alt_count',
			'ncbi_build', 'center', 'tumorRefCount', 'normalAltCount', 'normalRefCount', 'unique_sample_key',
			'unique_patient_key']
		"""
		dff = (df_mut
					.groupby(['pid', "sample_id", "barcode", 'symbol', 'refseq_mrna_id', "entrez_gene_id", "protein_mut", 'mutation_type', "variant_type", "chr"])
					.size()
					.reset_index(name="n_mutations")
				)
		
		dff = dff[dff.barcode.notna()]
		dff = dff[dff.entrez_gene_id.notna()]

		dff = dff.sort_values(["barcode", "sample_id", "symbol", "protein_mut"])
		dff = dff.reset_index(drop=True)
		
		_ = pdwritecsv(dff, fname_mut, self.root_data, verbose=verbose)
		_ = pdwritecsv(df_mut, fname_extmut, self.root_data, verbose=verbose)

		return dff, df_mut


	def cbioportal_studies(self):
		url = "https://www.cbioportal.org/api/studies"

		res = requests.get(url, headers={"Accept": "application/json"})
		res.raise_for_status()

		studies = res.json()

		study_ids = [s["studyId"] for s in studies]

		return study_ids
	
	"""
	def get_expression_files_given_samples(self, pid:str, subtype:str, stage:str,
										   case_ids: List, batch_size:int=20,
										   force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		Retrieve RNA-seq expression files for given case_ids
		input: case_ids
		output: dataframe
		'''

		fname = self.fname_rnaseq_exp_files%(pid, subtype, stage)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_files = pdreadcsv(fname, self.root_data, verbose=verbose)
			self.df_files = df_files

			return df_files

		if not case_ids:
			print(f"No samples ids")
			return pd.DataFrame()

		filters = {"op": "and",
					"content": [
						{ "op": "in",
							"content": {"field": "cases.case_id", "value": case_ids}
						},
						{ "op": "in",
							"content": {"field": "data_type", "value": ["Gene Expression Quantification"]}
						},
						{ "op": "in",
							"content": {"field": "experimental_strategy", "value": ["RNA-Seq"]}
						},
						{ "op": "in",
							"content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}
						},
						{ "op": "in",
						  "content": { "field": "cases.samples.sample_type", "value": ["Primary Tumor", "Solid Tissue Normal"]	}
						},		
						{ "op": "in",
							"content": {"field": "data_format", "value": ["TSV"]}
						},
						{ "op": "in",
							"content": {"field": "access", "value": ["open"] }
						}
					]
				}
		

		#------------- batch search 
		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None

		print("Searching: ", end='')
		try:
			while True:
				print(".", end='')

				params = {
					"filters": json.dumps(filters),
					"fields": ",".join([
						"file_id", 	"file_name", "file_size", "md5sum",
						"cases.case_id", "cases.submitter_id",
						"cases.samples.sample_id", 	"cases.samples.submitter_id", 
						"analysis.workflow_type"
					]),
					"expand": "cases,cases.samples",
					"format": "JSON",
					"from": from_,
					"size": size_,
				}

				res = requests.get(self.url_gdc_files, params=params)
				data = res.json()
	
				hits = data.get("data", {}).get("hits", [])

				if total is None:
					total = data["data"]["pagination"]["total"]

				if not hits:
					break

				all_hits.extend(hits)
				from_ += size_

			print("\n")

			if all_hits == []:
				print(f"No cases found for {pid} / {subtype} / {stage} / samples {len(case_ids)} ")
				return pd.DataFrame()
			
			#------------ lost data? ------------------
			N = len(all_hits)

			if N < total:
				print(f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} ")
			else:
				print(f"👉 Returned {N} / Total paginated {total}")

			#------------ having all hits -------------
			rows = []

			for hit in all_hits:

				workflow = hit.get("analysis", {}).get("workflow_type")

				for case in hit.get("cases", []):
					for sample in case.get("samples", []):

						if case.get("case_id") not in case_ids:
							continue  # 🔥 critical fix

						rows.append({
							"case_id": case.get("case_id"),
							"sample_id": sample.get("sample_id"),
							"file_id": hit.get("file_id"),
							"file_name": hit.get("file_name"),
							"case_submitter_id": case.get("submitter_id"),
							"sample_submitter_id": sample.get("submitter_id"),
							"workflow": workflow
						})

			df_files = pd.DataFrame(rows)
			cols = list(df_files.columns)

			# 🔹 Metadata 
			df_files["project_id"] = pid 
			df_files["subtype"] = subtype 
			df_files["stage"] = stage 

			cols = ["project_id", "subtype", "stage"] + cols
			df_files = df_files[cols]

			df_files = df_files.drop_duplicates(subset="file_id")
			df_files = df_files.sort_values("sample_id", ascending=False).reset_index(drop=True)
			
			_ = pdwritecsv(df_files, fname, self.root_data, verbose=verbose)
			if verbose: print(f"Found {len(df_files)} files in selected samples.")

		except Exception as e:
			print(f"Error for '{pid}', '{subtype}', and {stage}: {e}")
			self.df_files = pd.DataFrame()
			return self.df_files

		self.df_files = df_files

		return df_files

	
	def download_file(self, force:bool=False, verbose:bool=True) -> bool:

		if os.path.exists(self.gdc_ouptut_filename) and not force:
			if verbose: print(f"File already downloaded and cleaned: {self.gdc_filename}")
			return True

		if not isinstance(self.gdc_file_id, str) or len(self.gdc_file_id) < 3:
			print(f"Filename bad formated {self.gdc_file_id}.")
			return False

		if "." in self.gdc_file_id:
			print(f"Filename without '.' {self.gdc_file_id}.")
			return False

		ret = True

		try:
			url = self.url_gdc_data%(self.gdc_file_id)
			response = requests.get(url, stream=True)

			with open(self.gdc_filename, "wb") as f:
				for chunk in response.iter_content(chunk_size=8192):
					if chunk:
						f.write(chunk)

			if verbose: print(f"File saved at {self.gdc_filename}.")

		except Exception as e:
			print(f"No data found for {self.gdc_file_id}. error: {e}")
			ret = False
			
		return ret

	def clean_gdc_tsv(self, force:bool=False, verbose:bool=True) -> pd.DataFrame:
		'''
		clean comment lines and transform in tsv file
		auto-detect comments

		original table --> self.gdc_filename
		final tsv table --> self.gdc_ouptut_filename

		'''		
		if os.path.exists(self.gdc_ouptut_filename) and not force:
			df2 = pdreadcsv(self.gdc_ouptut_filename)
			if verbose: print(f"Reading table {df2.shape}: {self.gdc_ouptut_filename}")
			return df2
		
		try:
			df = pd.read_csv(self.gdc_filename, sep="\t", comment="#")

			df = self.clean_expression_table_value_col(df)

			_ = pdwritecsv(df, self.gdc_ouptut_filename)

			if verbose: print(f"Read table {df.shape}: {self.gdc_ouptut_filename}")
		except Exception as e:
			print(f"Could not prepare and save '{self.gdc_ouptut_filename}' - error: {e}")
			return pd.DataFrame()
		
		return df2
	
	def show_fnames(self):
		print(f"fname       '{self.gdc_file_name}'")
		print(f"fname id    '{self.gdc_file_id}'")
		print(f"data type   '{self.gdc_data_type}'")
		print(f"downloaded  '{self.gdc_fname}'")
		print(f"final fname '{self.gdc_ouptut_fname}'")

	def clean_expression_table_value_col(self, df:pd.DataFrame) -> pd.DataFrame:
		
		# Remove summary rows (N_*)
		df = df[~df["gene_id"].str.startswith("N_")]

		# Keep only valid Ensembl genes
		df = df[df["gene_id"].str.startswith("ENSG")]

		# Remove version from gene_id (ENSG... -> ENSG...)
		df["gene_id"] = df["gene_id"].str.split(".").str[0]

		# Optional: keep only relevant columns
		cols = ["gene_id", "gene_name", "gene_type", self.value_col]
		df2 = df[cols].copy()

		# Rename for clarity
		df2 = df2.rename(columns={"gene_name": "symbol", self.value_col: self.exp_unit})

		return df2


	self.fname_samples      = 'samples_for_%s_subtype_%s_stage_%s.tsv'
	def get_samples_deprecated(self, pid:str, subtype:str, stage:str, batch_size:int=200,
				    force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		calc all cases, given and pid, subtype, stage, and case_id

		input: pid = primary site ID, subtype, stage, and case_id
		output: dataframe
		'''

		fname = self.fname_samples%(pid, subtype, stage)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_case = pdreadcsv(fname, self.root_data, verbose=verbose)
			self.df_case = df_case

			return df_case
				

		filters = { "op": "and", 
			 		"content": [ 
						{ "op": "in", 
						  "content": { "field": "cases.project.project_id", "value": [pid] } 
						}, 
						{ "op": "in", 
						  "content": { "field": "diagnoses.primary_diagnosis", "value": [subtype] } 
						},
						{ "op": "in", 
						  "content": { "field": "diagnoses.ajcc_pathologic_stage", "value": [stage] } 
						},
					]
				   }

		#------------- batch search 
		all_hits = []
		from_ = 0
		size_ = batch_size
		total = None

		print("Searching: ", end='')
		try:
			while True:
				print(".", end='')
				params = {
					"filters": json.dumps(filters),
					"fields": ",".join([
						"case_id",
						"submitter_id",
						"samples.sample_id",
						"samples.submitter_id",
						"samples.sample_type"
					]),
					"from": from_,
					"size": size_,
				}

				res = requests.get(self.url_gdc_cases, params=params)
				response = res.json()

				if 'data' not in response:
					print(f"No data found for '{pid}', '{subtype}', and {stage}")
					self.df_samples = pd.DataFrame()
					return self.df_samples

				hits = response.get("data", {}).get("hits", [])

				if total is None:
					total = response["data"]["pagination"]["total"]
				
				if not hits:
					break

				all_hits.extend(hits)
				from_ += size_

			print("\n")

			if all_hits == []:
				print(f"No cases found for {pid} / {subtype} / {stage} ")
				self.df_samples = pd.DataFrame()
				return self.df_samples
			
			#------------ lost data? ------------------
			N = len(all_hits)

			if N < total:
				print(f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} ")
			else:
				print(f"👉 Returned {N} / Total paginated {total}")

			#------------ having all hits -------------
			rows = []

			for hit in all_hits:
				for s in hit.get("samples", []):
					rows.append({
						"case_id": hit["case_id"],
						"case_submitter_id": hit["submitter_id"],
						"sample_id": s.get("sample_id"),
						"sample_submitter_id": s.get("submitter_id"),
						"sample_type": s.get("sample_type")
					})

			df_samples = pd.DataFrame(rows)
			cols = list(df_samples.columns)

			# 🔹 Metadata 
			df_samples["project_id"] = pid 
			df_samples["subtype"] = subtype 
			df_samples["stage"] = stage 

			cols = ["project_id", "subtype", "stage"] + cols
			df_samples = df_samples[cols]
			
			df_samples = df_samples.sort_values("sample_id", ascending=False).reset_index(drop=True)
			
			_ = pdwritecsv(df_samples, fname, self.root_data, verbose=verbose)
			if verbose: print(f"Found {len(df_samples)} samples for many cases.")

		except Exception as e:
			print(f"Error while searching with '{pid}', '{subtype}', and {stage}: {e}")
			df_samples = pd.DataFrame()

		self.df_samples = df_samples

		return df_samples


	def get_cases_deprecated(self, pid:str, subtype:str, stage:str, size_:int=10000,
				    force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		calc all cases, given and pid, subtype, and stage

		input: pid = primary site ID, subtype, and stage
		output: dataframe
		'''

		fname = self.fname_cases%(pid, subtype, stage)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_case = pdreadcsv(fname, self.root_data, verbose=verbose)
			self.df_case = df_case

			return df_case
		

		stage_clean = stage.lower().replace("stage", "").strip()
		stage1 = f"stage {stage_clean.upper()}"
		stage2 = f"Stage {stage_clean.upper()}"
		
		# Stage III, or Stage IIIa and Stage IIIb //  stage IIIa and stage IIIb
		filters = { "op": "and", 
			 		"content": [ 
						{ "op": "in", 
						  "content": { "field": "cases.project.project_id", "value": [pid] } }, 
						{ "op": "in", 
						  "content": { "field": "diagnoses.primary_diagnosis", "value": [subtype] } },
						{ "op": "or",
						  "content": [
								{ "op": "like",
									"content": { "field": "diagnoses.ajcc_pathologic_stage", "value": stage1 + "%" }
								},
								{ "op": "like",
								"content": { "field": "diagnoses.ajcc_pathologic_stage", "value": stage2 + "%" }
								},
							]
						}
					]
				   }

		params = {
			"filters": json.dumps(filters), 
			"fields": "case_id",
    		"size": size_
			}

		try:
			res = requests.get(self.url_gdc_cases, params=params)
			data = res.json()

			hits = data.get("data", {}).get("hits", [])

			if not hits:
				print(f"No cases found for {pid} / {subtype} / {stage}")
				return pd.DataFrame()

			case_list = [h["case_id"] for h in hits]

			df_case = pd.DataFrame({"case_id": case_list, "n":1})
			# 🔹 Fraction 
			df_case["frac"] = df_case["n"] / df_case["n"].sum() 
			cols = list(df_case.columns)
			
			# 🔹 Metadata 
			df_case["project_id"] = pid 
			df_case["subtype"] = subtype 
			df_case["stage"] = stage 

			cols = ["project_id", "subtype", "stage"] + cols
			df_case = df_case[cols]
			
			df_case = df_case.sort_values("n", ascending=False).reset_index(drop=True)
			
			_ = pdwritecsv(df_case, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"No data found for '{pid}', '{subtype}', and {stage}. error: {e}")
			df_case = pd.DataFrame()

		self.df_case = df_case

		return df_case
			

	"""
