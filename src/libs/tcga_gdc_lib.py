#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

import os, requests, json, re
import pandas as pd
from collections import Counter
from typing import List, Tuple  # , Any

from libs.Basic import *

class GDC(object):
	def __init__(self, root_data:str='../data/'):
		
		self.url_gdc_project = "https://api.gdc.cancer.gov/projects"
		self.url_gdc_cases = "https://api.gdc.cancer.gov/cases"
		self.url_gdc_files = "https://api.gdc.cancer.gov/files"
		self.url_gdc_data  = f"https://api.gdc.cancer.gov/data/%s"

		self.root_data = root_data

		self.clean_gdc_files()

	def clean_gdc_files(self):
		self.gdc_file_name = ''
		self.gdc_file_id = ''
		self.gdc_data_type = ''

		self.fname_programs = 'programs.txt'

		# primary_site
		self.fname_primary_site = 'primary_site_program_%s.tsv'
		self.fname_cases        = 'cases_for_PS_%s.tsv'
		self.fname_subtype      = 'subtype_for_PS_%s.tsv'
		self.fname_stage        = 'stage_for_PS_%s_Subtype_%s.tsv'
		self.fname_samples      = 'samples_for_PS_%s_Subtype_%s_Stage_%s.tsv'
		self.fname_rnaseq_exp_files = 'rnaseq_exp_files_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.fname_cases_deprecated = 'cases_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.gdc_fname = ''
		self.gdc_filename = ''
		self.gdc_ouptut_fname = ''
		self.gdc_ouptut_filename = ''

		self.exp_unit = ""
		self.value_col = ""

		# program, primary site, subtype, stage, case_id, samples
		self.df_ps, self.df_subt, self.df_cases = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
		self.df_stage, self.df_sample, self.df_files = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

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
		output: df_cases, df_subt
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

			df_subt, df_prof = self.build_profile(df_cases)

			return df_cases, df_subt, df_prof
				

		def build_tcga_ontology(df):

			# ['id', 'primary_site', 'disease_type', 'case_id', 'pid', 'primary_diagnosis', 'tumor_grade', 'stage']

			print("-------- build -------------")
			print(df.columns)
			print("---------------------------")

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
			print(">>> main_diag", main_diag, key)
			dicf = {"diagnosis": 'unknown',
					"ajcc": 'unknown'
					}
			
			if isinstance(x, list):
				print(">>> list", x)
				for d in x:
					if isinstance(d, dict):
						print(">>> found1111?", d)

						if d.get("primary_diagnosis") == main_diag and key in d.keys():
							print(">>> found2222?", d)
							dicf = {"diagnosis": main_diag,
									"ajcc": d.get(key)
									}
							break

			elif isinstance(x, dict):
				print(">>> dict", x)
				dicf = {"diagnosis": x.get("primary_diagnosis"),
						"ajcc": x.get(key)
						}

			return dicf		

		def unpack_diagnoses(df):

			#------------------- main_diag -------------------------------------------------------
			diag_list = df["diagnoses"].map(lambda x: extract_any(x, "primary_diagnosis"))
			diag_list = [x for x in diag_list if isinstance(x, str) and x.strip()]
			dic = Counter(diag_list)

			dfa = pd.DataFrame({
				'diag': list(dic.keys()),
				'n': list(dic.values())
			})

			dfa = dfa.sort_values('n', ascending=False)
			main_diag = dfa.iloc[0].diag

			#------------------- main_diag end ---------------------------------------------------

			dicf = df["diagnoses"].map(lambda x: extract_all(x, "ajcc_clinical_stage", main_diag))

			df["subtype_global"] = dicf["diagnosis"]
			df["stage_ajcc"]     = dicf["ajcc"]

			df["tumor_grade"] = df["diagnoses"].map(lambda x: extract_any(x, "tumor_grade"))
			df["stage_clin"]  = df["diagnoses"].map(lambda x: extract_any(x, "ajcc_clinical_stage"))
			df["figo_stage"]  = df["diagnoses"].map(lambda x: extract_any(x, "figo_stage"))
			df["tumor_stage"] = df["diagnoses"].map(lambda x: extract_any(x, "tumor_stage"))

			df["stage"] = df["stage_ajcc"] \
							.fillna(df["stage_clin"]) \
							.fillna(df["figo_stage"]) \
							.fillna(df["tumor_stage"])

			df["stage"] = df["stage"].fillna('unknown')

			# df_cases = df_cases.drop('diagnoses', axis=1)
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

		print("Searching: ", end='')
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
					self.df_cases = response
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

			print("> 1")
			print("----------- 1 ---------------")
			print('rows', len(df_cases), '\ncolumns', df_cases.columns)
			print("---------------------------")

			'''
			# flatten lists
			for col in df_cases.columns:
				df_cases[col] = df_cases[col].apply(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x)


			print("> 2")
			print("---------------------------")
			print(df_cases)
			print("---------------------------")
			'''

			"""
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

			df_cases = unpack_diagnoses(df_cases)

			if debug:
				print("> 2")
				print("----------- 2 ---------------")
				print(df_cases.head(3).T)
				print("---------------------------")

			self.df_cases2 = df_cases

			# rename for sanity
			df_cases = df_cases.rename(columns={
				"project.project_id": "pid",
			})

			if debug:
				print("> 3")
				print("----------- 3 ---------------")
				print(df_cases.head(3).T)
				print("---------------------------")

			df_cases = build_tcga_ontology(df_cases)

			df_cases["validity"] = df_cases.apply(classify_validity, axis=1)

			df_cases["n"] = 1
			df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()
			df_cases = df_cases.sort_values("n", ascending=False).reset_index(drop=True)
			df_cases.reset_index(drop=True, inplace=True)

			_ = pdwritecsv(df_cases, fname_cases, self.root_data, verbose=verbose)
			df_subt = df_cases.groupby(["subtype_global", "tumor_class", "subtype_tissue", "stage"]).size().reset_index(name="n")
			_ = pdwritecsv(df_subt,  fname_subt, self.root_data, verbose=verbose)


		except Exception as e:
			print(f"Error for searching cases for '{pid}'. error: {e}")
			self.df_cases = df_cases
			self.df_subt  = pd.DataFrame()
			self.df_prof  = pd.DataFrame()
			return self.df_cases, self.df_subt, self.df_prof

		if do_filter:
			df_cases = apply_filter(df_cases)

		self.df_cases = df_cases
		self.df_subt  = df_subt
		
		df_subt, df_prof = self.build_profile(df_cases)

		return df_cases, df_subt, df_prof
	

	def build_profile(self, df_cases: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:

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

					if stage and result["stage"] is None:
						result["stage"] = stage

					# tumor grade
					grade = d.get("tumor_grade")
					if grade and result["tumor_grade"] is None:
						result["tumor_grade"] = grade

			return result

		df_subt = df_cases.groupby(["subtype_global", "tumor_class", "subtype_tissue", "stage"]).size().reset_index(name="n")
		self.df_subt  = df_subt

		profiles = df_cases["diagnoses"].apply(clean_case_profile)
		df_prof = pd.DataFrame(profiles.tolist())
		self.df_prof = df_prof

		return df_subt, df_prof


	def get_stages(self, pid:str, subtype:str, do_filter:bool=True, 
				     force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		calc all stages, given and pid and a subtype

		input: pid = primary site ID and subtype 
		output: dataframe
		'''

		fname = self.fname_stage%(pid, subtype)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_stage = pdreadcsv(fname, self.root_data, verbose=verbose)

			if do_filter:
				df_stage = df_stage[df_stage.is_valid == True].copy()
				df_stage.reset_index(drop=True, inplace=True)
				df_stage["frac"] = df_stage["n"] / df_stage["n"].sum()	

			self.df_stage = df_stage

			return df_stage
				

		filters = { "op": "and", 
			 		"content": [ 
						 { "op": "in", 
						   "content": { "field": "cases.project.project_id", "value": [pid] } }, 
						 { "op": "in",
						   "content": { "field": "diagnoses.primary_diagnosis", "value": [subtype] } } 
					]
				   }

		params = {
			"filters": json.dumps(filters), 
			"facets": "diagnoses.ajcc_pathologic_stage", 
			"size": 0 
			}

		try:
			res = requests.get(self.url_gdc_cases, params=params)
			response = res.json()

			if 'data' not in response.keys():
				print(f"No data found for '{pid}' and '{subtype}'")
				print(">>> response", response)
				self.df_stage = pd.DataFrame()
				return self.df_stage
			

			aggs = response.get("data", {}).get("aggregations", {}) 
			
			buckets = aggs.get("diagnoses.ajcc_pathologic_stage", {}).get("buckets", []) 
			
			if not buckets: 
				if verbose: print(f"No stages found for {pid} / {subtype}") 
				return pd.DataFrame() 
			
			stage_list = [b["key"] for b in buckets] 
			stage_count = [b["doc_count"] for b in buckets] 
			
			df_stage = pd.DataFrame({"stage": stage_list, "n": stage_count })

			df_stage["stage_raw"] = df_stage["stage"]
			df_stage["stage_clean"] = df_stage["stage"].str.strip().str.upper()
			# 🔹 Normalize 
			df_stage["stage"] = df_stage["stage"].str.strip().str.lower()

			# 🔹 Validity flag 
			invalid_terms = ["not reported", "unknown", "_missing"]
			pattern = "|".join(invalid_terms) 
			
			df_stage["is_valid"] = ~df_stage["stage"].str.contains(pattern, case=False, regex=True) 
			
			# 🔹 Fraction 
			df_stage["frac"] = df_stage["n"] / df_stage["n"].sum() 
			
			# 🔹 Metadata 
			df_stage["project_id"] = pid 
			df_stage["subtype"] = subtype 

			cols = list(df_stage.columns)
			cols = ["project_id", "subtype"] + cols[:-2]
			df_stage = df_stage[cols]
			
			df_stage = df_stage.sort_values("n", ascending=False).reset_index(drop=True)
			
			_ = pdwritecsv(df_stage, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"Error searching stages for '{pid}' and '{subtype}': {e}")
			df_stage = pd.DataFrame()

		if do_filter:
			df_stage = df_stage[df_stage.is_valid == True].copy()
			df_stage.reset_index(drop=True, inplace=True)
			df_stage["frac"] = df_stage["n"] / df_stage["n"].sum()

		self.df_stage = df_stage

		return df_stage
	

	def get_samples(self, pid:str, subtype:str, stage:str, batch_size:int=200,
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
					self.df_sample = pd.DataFrame()
					return self.df_sample

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
				self.df_sample = pd.DataFrame()
				return self.df_sample
			
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

			df_sample = pd.DataFrame(rows)
			cols = list(df_sample.columns)

			# 🔹 Metadata 
			df_sample["project_id"] = pid 
			df_sample["subtype"] = subtype 
			df_sample["stage"] = stage 

			cols = ["project_id", "subtype", "stage"] + cols
			df_sample = df_sample[cols]
			
			df_sample = df_sample.sort_values("sample_id", ascending=False).reset_index(drop=True)
			
			_ = pdwritecsv(df_sample, fname, self.root_data, verbose=verbose)
			if verbose: print(f"Found {len(df_sample)} samples for many cases.")

		except Exception as e:
			print(f"Error while searching with '{pid}', '{subtype}', and {stage}: {e}")
			df_sample = pd.DataFrame()

		self.df_sample = df_sample

		return df_sample

	def get_expression_files_given_samples(self, pid:str, subtype:str, stage:str,
										   case_ids: List, batch_size:int=20,
										   force:bool=False, verbose:bool=False) -> pd.DataFrame:
		"""
		Retrieve RNA-seq expression files for given case_ids
		input: case_ids
		output: dataframe
		"""

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
		
		"""
		"""

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
		'''		
		if os.path.exists(self.gdc_ouptut_filename) and not force:
			df2 = pdreadcsv(self.gdc_ouptut_filename)
			if verbose: print(f"Reading table {df2.shape}: {self.gdc_ouptut_filename}")
			return df2
		
		try:
			df = pd.read_csv(self.gdc_filename, sep="\t", comment="#")

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

			_ = pdwritecsv(df2, self.gdc_ouptut_filename)

			if verbose: print(f"Read table {df2.shape}: {self.gdc_ouptut_filename}")
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
			
			# 🔹 Metadata 
			df_case["project_id"] = pid 
			df_case["subtype"] = subtype 
			df_case["stage"] = stage 

			cols = list(df_case.columns)
			cols = ["project_id", "subtype", "stage"] + cols[:-3]
			df_case = df_case[cols]
			
			df_case = df_case.sort_values("n", ascending=False).reset_index(drop=True)
			
			_ = pdwritecsv(df_case, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"No data found for '{pid}', '{subtype}', and {stage}. error: {e}")
			df_case = pd.DataFrame()

		self.df_case = df_case

		return df_case
	
