#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/19
# @author: Flavio Lichtenstein
# @local: Home sweet home

import os, requests, json
import pandas as pd

from typing import List, Tuple

from libs.Basic import *


class GDC(object):
	def __init__(self, case:str="Breast Cancer", root_data:str='../data/'):
		
		self.url_gdc_cases = "https://api.gdc.cancer.gov/cases"
		self.url_gdc_files = "https://api.gdc.cancer.gov/files"
		self.url_gdc_data  = f"https://api.gdc.cancer.gov/data/%s"


		self.case = case

		self.dir_case = case.lower().replace(' ', '_')
		self.root_case = os.path.join(root_data, self.dir_case)
		os.makedirs(self.root_case, exist_ok=True)

		self.clean_gdc_files()

	def clean_gdc_files(self):
		self.gdc_file_name = ''
		self.gdc_file_id = ''
		self.gdc_data_type = ''

		self.gdc_fname = ''
		self.gdc_filename = ''
		self.gdc_ouptut_fname = ''
		self.gdc_ouptut_filename = ''

		self.exp_unit = ""
		self.value_col = ""		


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