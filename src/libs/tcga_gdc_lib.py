#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

import os, requests, json
import pandas as pd

from typing import List, Tuple, Any

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
		self.fname_subtype      = 'subtype_for_PS_%s.tsv'
		self.fname_stage        = 'stage_for_PS_%s_Subtype_%s.tsv'
		self.fname_cases        = 'cases_for_PS_%s_Subtype_%s_Stage_%s.tsv'
		self.fname_samples      = 'samples_for_PS_%s_Subtype_%s_Stage_%s.tsv'
		self.fname_rnaseq_exp_files = 'rnaseq_exp_files_for_PS_%s_Subtype_%s_Stage_%s.tsv'

		self.gdc_fname = ''
		self.gdc_filename = ''
		self.gdc_ouptut_fname = ''
		self.gdc_ouptut_filename = ''

		self.exp_unit = ""
		self.value_col = ""

		# program, primary site, subtype, stage, case_id, samples
		self.df_ps, self.df_subt, self.df_stage = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
		self.df_case, self.df_sample, self.df_files = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

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
					"content": { "field": "program.name", "value": [program]
					}
				}

		params = {
			"filters": filters,
			"fields": "project_id,name,primary_site,disease_type",
			"format": "JSON",
			"size": 1000
		}

		try:
			res = requests.get(self.url_gdc_project, params=params)
			projects = res.json()["data"]["hits"]

			df_ps = pd.DataFrame(projects)
			cols = ["project_id", "primary_site", "disease_type"]
			df_ps = df_ps[cols]

			df_ps = df_ps.sort_values(["primary_site", "disease_type"])

			_ = pdwritecsv(df_ps, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"No data found for '{program}'. error: {e}")
			df_ps = pd.DataFrame()

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


	def get_subtypes(self, pid:str, do_filter:bool=True, force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		calc all subtypes, given and pid

		filter: NOS (Not Otherwise Specified) → the pathologist could not (or did not) assign a more specific subtype.
		e.g.: "Yes, it's an adenocarcinoma — but we don’t have finer classification"

		input: pid = primary site ID
		output: dataframe
		'''

		fname = self.fname_subtype%(pid)
		filename = os.path.join(self.root_data, fname)

		if os.path.exists(filename) and not force:
			df_subt = pdreadcsv(fname, self.root_data, verbose=verbose)

			if do_filter:
				df_subt = df_subt[df_subt.is_valid == True].copy()
				df_subt.reset_index(drop=True, inplace=True)
				df_subt["frac"] = df_subt["n"] / df_subt["n"].sum()

			self.df_subt = df_subt

			return df_subt
				

		filters = {
			"op": "in",
			"content": {
				"field": "cases.project.project_id",
				"value": [pid]
			}
		}

		params = {
			"filters": json.dumps(filters),
			"facets": "diagnoses.primary_diagnosis", 
			"size": 0
		}

		try:
			res = requests.get(self.url_gdc_cases, params=params)
			data = res.json()

			buckets = data["data"]["aggregations"]["diagnoses.primary_diagnosis"]["buckets"]

			subtype_list  = [b["key"] for b in buckets]
			subtype_count = [b["doc_count"] for b in buckets]

			df_subt = pd.DataFrame({'subtype': subtype_list, 'n':subtype_count})

			df_subt["subtype"] = (df_subt["subtype"].str.lower().str.strip())
			df_subt["subtype_raw"] = df_subt["subtype"]  # exact GDC value

			invalid_terms = ["not reported", "unknown", "nos"]
			df_subt["is_valid"] = ~df_subt["subtype"].apply( lambda x: any(term in x for term in invalid_terms)
)
			df_subt["frac"] = df_subt["n"] / df_subt["n"].sum()

			df_subt = df_subt.sort_values("n", ascending=False).reset_index(drop=True)
			df_subt.reset_index(drop=True, inplace=True)
			
			# 🔹 Metadata 
			df_subt["project_id"] = pid 

			cols = list(df_subt.columns)
			cols = ["project_id"] + cols[:-1]
			df_subt = df_subt[cols]

			_ = pdwritecsv(df_subt, fname, self.root_data, verbose=verbose)

		except Exception as e:
			print(f"No data found for '{pid}'. error: {e}")
			df_subt = pd.DataFrame()

		if do_filter:
			df_subt = df_subt[df_subt.is_valid == True].copy()
			df_subt.reset_index(drop=True, inplace=True)
			df_subt["frac"] = df_subt["n"] / df_subt["n"].sum()

		self.df_subt = df_subt

		return df_subt
	

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
			data = res.json()

			aggs = data.get("data", {}).get("aggregations", {}) 
			
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
			print(f"No data found for '{pid}' and '{subtype}'. error: {e}")
			df_stage = pd.DataFrame()

		if do_filter:
			df_stage = df_stage[df_stage.is_valid == True].copy()
			df_stage.reset_index(drop=True, inplace=True)
			df_stage["frac"] = df_stage["n"] / df_stage["n"].sum()

		self.df_stage = df_stage

		return df_stage
	

	def get_cases(self, pid:str, subtype:str, stage:str, size_:int=10000,
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
				print(f"No cases found for {pid} / {subtype} / {stage} ")
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
			print(f"No data found for '{pid}', '{subtype}', and {stage}. error: {e}")
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

