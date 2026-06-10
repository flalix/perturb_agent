#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-

import json
import requests
import io
import time
from pathlib import Path
from os.path import join as osjoin
from os.path import exists as exists
import pandas as pd
import numpy as np
from typing import Tuple, List # Optional, Iterable, Set, Any

from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from scipy.stats import f_oneway

import multiprocessing
from multiprocessing import Pool, freeze_support

from IPython.display import HTML

from libs.Basic import pdwritecsv, pdreadcsv, write_txt, read_txt, title_replace # create_dir
from libs.venn_lib import get_venn_sections, defineClass
from libs.MTD_lib import MTD
from libs.stat_lib import chi2_or_fisher_exact_test


class enricheR(MTD):
	def __init__(self, disease:str, gene_protein:str, s_omics:str, project:str, s_project:str, 
			     root0:Path, root0_data:Path,  prog_id:str, psi_id:str,
				 case_list:List, dic_case_list:dict, has_age:bool=True, has_gender:bool=True, exp_normalization:bool=False, 
				 std_filename:str='', std_filename_list:list=[],
				 geneset_num:int=0, ptw_min_num_of_degs_cut:int=3,
				 tolerance_pPMI:float=.15, s_pathw_enrichm_method:str='enricher',
				 LFC_cut_inf:float=0.40, fdr_ptw_cutoff_list:List=[],
				 num_of_genes_list:List=[3], lfc_list = [], fdr_list = [],
				 min_lfc_modulation:float=0.40, type_sat_ptw_index:str='linear_sat', 
				 saturation_lfc_param:float=5., enr_db_list:List=[], pPMI_normalized:bool=False):
			
		super().__init__(disease=disease, gene_protein=gene_protein, s_omics=s_omics, 
				project=project, s_project=s_project, root0=root0, root0_data=root0_data,
				prog_id=prog_id, psi_id=psi_id, case_list=case_list, dic_case_list=dic_case_list, 
				has_age=has_age, has_gender=has_gender, exp_normalization=exp_normalization,
				std_filename=std_filename, std_filename_list=std_filename_list,
				geneset_num=geneset_num, ptw_min_num_of_degs_cut=ptw_min_num_of_degs_cut,
				tolerance_pPMI=tolerance_pPMI, s_pathw_enrichm_method=s_pathw_enrichm_method,
				LFC_cut_inf=LFC_cut_inf, fdr_ptw_cutoff_list=fdr_ptw_cutoff_list,
				num_of_genes_list=num_of_genes_list, lfc_list=lfc_list, fdr_list=fdr_list, 
				min_lfc_modulation=min_lfc_modulation, type_sat_ptw_index=type_sat_ptw_index,
				saturation_lfc_param=saturation_lfc_param, enr_db_list=enr_db_list, pPMI_normalized=pPMI_normalized)
				

		self.fname_new_proteomics = 'new_taubate_LFC_%s_x_%s_%s.tsv'

		self.pathid = ''
		self.kgml, self.pathway_kgml = '', ''
		self.dfgc = pd.DataFrame()

		self.ENRICHR_URL_add_list  = 'https://maayanlab.cloud/Enrichr/addList'
		self.ENRICHR_URL_user_list = 'https://maayanlab.cloud/Enrichr/view?userListId=%s'
		self.ENRICHR_URL_query	   = 'https://maayanlab.cloud/Enrichr/enrich'
		self.ENRICHR_query_string  = '?userListId=%s&backgroundType=%s'

		self.old_pathway_cols = ['pathway', 'overlap', 'ptw_pval_cut', 'ptw_FDR_cut', 'Old P-value', 'Old Adjusted P-value', 'odds_ratio', 'combined_score', 'genes'],
		self.sel_pathway_cols = ['pathway', 'overlap', 'ptw_pval_cut', 'ptw_FDR_cut','odds_ratio', 'combined_score', 'genes']

		self.parse_gc_fields = ['ENTRY', 'SYMBOL', 'NAME', 'SEQUENCE', 'ORTHOLOGY', 'ORGANISM',
								'TYPE', 'REMARK', 'COMMENT', 'PATHWAY', 'NETWORK', 'BRITE', 'DBLINKS',
								'ELEMENT', 'DISEASE', 'POSITION', 'MOTIF', 'STRUCTURE', 'AASEQ', 'NTSEQ']

		self.dfkegg =  pd.DataFrame()
		self.df_uniprot = pd.DataFrame()

		self.degs_in_pathways_random, self.degs_not_in_pathways_random = [], []
		self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random = 0, 0
		
		self.n_degs_in_pathways_bca, self.n_degs_not_in_pathways_bca = 0, 0
		self.n_degs, self.n_degs_bca = 0, 0


		self.MAX_GENES_ENRICHR_SAFE = 2000

		self.biotype_annot = ['protein_coding','IG_C_gene', 'IG_D_gene', 'IG_J_gene',
			'IG_V_gene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene',
			'lncRNA', 'miRNA', 'misc_RNA', 'rRNA', 'tRNA', 'scaRNA', 'snRNA', 'snoRNA',]


	def set_df_enr_ipath(self, ipath: int) -> Tuple[str, List]:
		try:
			pathway = self.df_enr.iloc[ipath].pathway
			genes   = self.df_enr.iloc[ipath].genes.split(';')
		except:
			pathway, genes = '', []

		self.pathway = pathway
		self.genes = genes

		return pathway, genes

	def get_rest_kegg_pathways(self) -> pd.DataFrame:
		try:
			result = REST.kegg_list("pathway").read()

			dfkegg = self.to_df(result)

			if not dfkegg.empty:
				dfkegg.columns = ['path_id', 'pathway']
				_ = pdwritecsv(dfkegg, self.kegg_fname, self.root_kegg)
		except:
			print("Error: accessing/reading REST.kegg_list('pathway')")
			dfkegg = pd.DataFrame()

		self.dfkegg = dfkegg

		return dfkegg
	

	def get_kegg_pathays(self) -> pd.DataFrame:
		if not self.dfkegg.empty:
			return self.dfkegg

		filefull = osjoin(self.root_kegg, self.kegg_fname)

		if exists(filefull):
			dfkegg = pdreadcsv(self.kegg_fname, self.root_kegg)
		else:
			dfkegg = self.get_rest_kegg_pathways()

		self.dfkegg = dfkegg
		return dfkegg

	def set_kegg_ipath(self, ipath: int) -> Tuple[object, object]:
		if self.dfkegg.empty:
			dfkegg = self.get_kegg_pathays()
			if dfkegg.empty:
				return '???', '???'

		try:
			path_id = self.dfkegg.iloc[ipath].path_id
			pathway = self.dfkegg.iloc[ipath].pathway
		except:
			pathway, path_id = '???', '???'

		self.pathway = pathway
		self.path_id = path_id

		return path_id, pathway

	def set_kegg_path_id(self, path_id: int) -> str:
		if self.dfkegg.empty:
			dfkegg = self.get_kegg_pathays()
			if dfkegg.empty:
				return '???'
					
		dfa = self.dfkegg[self.dfkegg.path_id == path_id]
		if dfa.empty:
			pathway = '???'
		else:
			pathway = dfa.iloc[0].pathway

		self.path_id = path_id
		self.pathway = pathway

		return pathway


	def find_kegg_pathway_by_name(self, pathway, verbose=False) -> str:
		if self.dfkegg.empty:
			dfkegg = self.get_kegg_pathays()
			if dfkegg.empty:
				return '???'

		dfa = self.dfkegg[self.dfkegg.pathway == pathway]
		if dfa.empty:
			if verbose: print("Could not find pathway: '%s'"%(pathway))
			path_id = '???'
		else:
			path_id = dfa.iloc[0].path_id
			path_id = path_id.replace('path:', '')

		self.pathid  = path_id
		self.path_id = path_id
		self.pathway = pathway

		return path_id


	'''--- get_id_list - https://maayanlab.cloud/Enrichr/enrich -----'''
	def open_session_upload_symbols(self, deg_list: list, description:object=None) -> Tuple[str, str]:
		self.deg_list = deg_list
		self.shortId, self.userListId = None, None

		if not isinstance(deg_list, list):
			print("deg_list must be a list.")
			return None, None

		if description is None:
			description = self.s_project + ' ' + self.case

		genes_str = '\n'.join(deg_list)

		payload = {
			'list': (None, genes_str),
			'description': (None, description)
		}

		try:
			response = requests.post(self.ENRICHR_URL_add_list, files=payload)
			if not response.ok:
				print('Error analyzing gene list')
				return '', ''

			data = json.loads(response.text)
			shortId = data['shortId']
			userListId = data['userListId']
		except:
			print("Problems in request.")
			return '', ''

		self.shortId, self.userListId = shortId, userListId

		return shortId, userListId

	def is_ok_symbols(self, deg_list):

		response = requests.get(self.ENRICHR_URL_user_list % self.userListId)
		if not response.ok:
			print('Error getting gene list')
			return False, []

		try:
			data = json.loads(response.text)
			genes = data['genes']
		except:
			print("Problems in request.")
			return False, []

		good_genes = [x for x in genes if x not in deg_list]
		return True, good_genes

	def open_reactome_enrichment_analysis(self, case:str,  species:str='Homo sapiens', 
										  force:bool=False, verbose:bool=False) -> pd.DataFrame:
		fname = 'taubate_%s_reactome_enrichment_results.csv'%(case)
		fname_sig = fname.replace('.csv', '') + '_sig.tsv'

		filefull = osjoin(self.root_enrich, fname_sig)
		if exists(filefull) and not force:
			df_enr = pdreadcsv(fname_sig, self.root_enrich, verbose=verbose)
			self.df_enr = df_enr

			dicgenes = {df_enr.iloc[i].pathway_id: df_enr.iloc[i].submitted_entities_found.split(';') for i in range(len(df_enr)) }
			self.dicgenes = dicgenes
			return df_enr

		filefull = osjoin(self.root_enrich, fname)
		if not exists(filefull):
			print("Could not find: %s"%(filefull))
			self.df_enr = pd.DataFrame()
			self.dicgenes = {}
			return pd.DataFrame()

		df_enr = pdreadcsv(fname, self.root_enrich, sep=',')

		cols = ['pathway_id', 'pathway', 'found_entities', 'total_entities', 'interactors_found', 'interactors_total',
			   'entities_ratio', 'ptw_pval_cut', 'ptw_FDR_cut', 'reactions_found', 'ractions_total', 'reactions_ratio', 'species_identifier',
			   'species_name', 'submitted_entities_found', 'mapped_entities', 'submitted_entities_hit_interactor', 'interacts_with',
			   'found_reaction_identifiers']

		df_enr.columns = cols

		df_enr = df_enr[df_enr.species_name == species]
		df_enr = df_enr[ (df_enr.ptw_FDR_cut < self.ptw_FDR_cut) & (df_enr.num_of_genes >= self.ptw_min_num_of_degs_cut)]
		df_enr = df_enr.sort_values('ptw_FDR_cut', ascending=True)
		df_enr.reset_index(inplace=True, drop=True)

		cols.remove("species_identifier")
		cols.remove("species_name")
		df_enr = df_enr[cols]

		ret = pdwritecsv(df_enr, fname_sig, self.root_enrich)
		self.df_enr = df_enr

		dicgenes = {df_enr.iloc[i].pathway_id: df_enr.iloc[i].submitted_entities_found.split(';') for i in range(len(df_enr)) }
		self.dicgenes = dicgenes
		return df_enr

	def calc_default_enrichment_analysis(self, geneset_num_list:List=[0, 1, 2, 4, 5, 7], force:bool=False, verbose:bool=False):

		LFC_cut = 1.
		lfc_FDR_cut = 0.05
		ptw_FDR_cut = 0.05

		for case in self.case_list:
			# icount=-1

			for geneset_num in geneset_num_list:
				self.set_db(geneset_num)

				if verbose: print(f">>> {case}, geneset = {self.geneset_lib}, lfc cutoff = {self.LFC_cut} and fdr cutoff = {self.lfc_FDR_cut}")

				# ret, degs, degs_ensembl, dflfc
				_, _, degs_ensembl, _ =  self.open_case_params(case, LFC_cut=LFC_cut, 
												lfc_FDR_cut=lfc_FDR_cut, 
												ptw_FDR_cut=ptw_FDR_cut, verbose=verbose)

				self.calc_EA_dataset_symbol(degs_ensembl, default=True, force=force, verbose=verbose)


	def calc_all_enrichment_analysis(self, geneset_num_list:List=[0, 1, 2, 4, 5, 7], 
								     method:str='spearman', iecho:int=-1,
									 force:bool=False, verbose:bool=False):

		dfsim = self.open_simulation_table()
		if dfsim.empty:
			print(f"Error: no previous simulation.")
			raise Exception("Stop: calc_all_enrichment_analysis()")

		for case in self.case_list:
			icount=-1

			'''
				open_case_simple:
					set case
					open dflfc_ori
			'''
			if not self.open_case_simple(case):
				print(f"Problems for case='{case}'")
				continue

			''' get filtered fdr's written in json files '''
			df_fdr = self.open_fdr_lfc_correlation(case, self.LFC_cut_inf)
			if df_fdr is None or df_fdr.empty:
				print(f"Problems with df_fdr json for case='{case}'")
				continue

			''' filter simulation df for case, df_fdr.fdr, LFC_cut_inf, and ptw_min_num_of_degs_cut '''
			dfsim2 = dfsim[ (dfsim.normalization == self.normalization) & (dfsim.case == case) &
							(dfsim.lfc_FDR_cut.isin(df_fdr.fdr) ) &
							(dfsim.LFC_cut >= self.LFC_cut_inf) &
							(dfsim.n_degs >= self.ptw_min_num_of_degs_cut)]

			if dfsim2.empty:
				print(f"No previous simulation for case='{case} with fdr in df_fdr.fdr and LFC_cut >= {self.LFC_cut_inf} and  num of {self.s_deg_dap}s >= {self.ptw_min_num_of_degs_cut}'")
				continue

			N = len(dfsim2)
			print(f"# {case}: {N} lfc-fdr-degs simulations; LFC_cut >= {self.LFC_cut_inf} and {self.s_deg_dap}s >= {self.ptw_min_num_of_degs_cut}")

			for geneset_num in geneset_num_list:
				self.set_db(geneset_num)
				print(f"\ngeneset = {self.geneset_lib}")

				for i in range(N):
					row = dfsim2.iloc[i]
					LFC_cut = row.LFC_cut
					lfc_FDR_cut = row.lfc_FDR_cut

					icount += 1

					'''
						if df_fdr has no correlation, than only calculate enrichment analysis for LFC_cut == 1
					'''
					df_fdr2 = df_fdr[ (df_fdr.case == case) & 
									  (df_fdr.fdr  == lfc_FDR_cut) &
									  (df_fdr.method == method) ]

					if df_fdr2.empty:
						print(f"Could not find df_fdr for {case} lfc_FDR_cut {lfc_FDR_cut:.2f} and method={method}")
						continue

					if pd.isnull(df_fdr2.iloc[0]['corr']):
						if LFC_cut != 1:
							continue

					if verbose: print(f"# {case}, geneset = {self.geneset_lib}, lfc cutoff = {LFC_cut} and fdr cutoff = {lfc_FDR_cut}")

					if isinstance(row.degs_ensembl, list):
						degs_ensembl2 = row.degs_ensembl
					else:
						degs_ensembl2 = eval(row.degs_ensembl)

					# degs2, degs_ensembl, dflfc
					_, degs_ensembl, _ = self.list_of_degs_set_params(LFC_cut, lfc_FDR_cut, verbose=False)


					if len(degs_ensembl) != len(degs_ensembl2):
						print(f"Error: mismatch len for {self.s_deg_dap}s for {case}, lfc cutoff = {LFC_cut} and fdr cutoff = {lfc_FDR_cut} -> {len(degs_ensembl)} pre calc != {len(degs_ensembl2)} calc now")
						continue

					if len(degs_ensembl) < self.ptw_min_num_of_degs_cut:
						if verbose: print(f"Warning:  degs_ensembl {len(degs_ensembl)} >= {self.ptw_min_num_of_degs_cut} {self.s_deg_dap}s - {degs_ensembl}")
						continue
		
					# if i > 10:break
					self.calc_EA_dataset_symbol(degs_ensembl, calc_many_sig=True, default=False, force=force, verbose=verbose)
					if iecho > 9 and icount%iecho==0:
						print("")
						print(f"{icount} / {N}) {self.s_deg_dap}s {len(degs_ensembl)} for {case}, lfc cutoff = {LFC_cut} and fdr cutoff = {lfc_FDR_cut}")
						self.echo_degs()
						print("")
						# self.echo_enriched_pathways()
						# print("")
					else:
						print(".",end="")

			print("")


	def calc_random_EA_dataset_symbol(self, i_sim:int, deg_list:List, verbose:bool=False):
		'''
		Calc Enerichment Analysis given a dataset and a DEG list.
		Saven in enrichment_random dir
		allway cut and save the default cutoffs
		save as many send i_sim - to perform statistics
		'''
		self.df_enr0 = pd.DataFrame()
		self.df_enr = pd.DataFrame()

		if len(deg_list) < self.ptw_min_num_of_degs_cut:
			print(f"Error: Random simulation {i_sim}, please send a deg_list with >= {self.ptw_min_num_of_degs_cut} {self.s_deg_dap}s")
			return

		self.deg_list = deg_list
		df_enr0 = self.calc_get_enriched_pathway(deg_list, return_value=True, verbose=verbose, i_sim=i_sim)
		self.df_enr0 = df_enr0

		if df_enr0 is None or df_enr0.empty:
			return

		df_enr = self.df_enr0[(self.df_enr0.fdr  < self.ptw_FDR_cut) &
							  (self.df_enr0.pval < self.ptw_pval_cut) &
							  (self.df_enr0.num_of_genes >= self.ptw_min_num_of_degs_cut)].copy()

		if verbose: print("<<<", self.ptw_FDR_cut, self.ptw_pval_cut, self.ptw_min_num_of_degs_cut, len(df_enr))

		if df_enr.empty:
			return

		df_enr = df_enr.sort_values(['fdr', 'num_of_genes'], ascending=[True, False])
		df_enr.reset_index(inplace=True, drop=True)
		self.df_enr = df_enr

		_, fname_cutoff = self.set_enrichment_name()
		fname_cutoff = fname_cutoff.replace('.tsv', f'_random_{i_sim}.tsv')
		_ = pdwritecsv(df_enr, fname_cutoff, self.root_enrich_random, verbose=verbose)

		return

	def calc_EA_dataset_symbol(self, deg_list:List, calc_many_sig:bool=True, default:bool=False, 
							  force:bool=False, verbose:bool=False):
		'''
		Calc Enerichment Analysis given a dataset and a DEG list.
		if calc_many_sig: calc the fdr cutoffs = [0.05, 0.75]
		'''

		if len(deg_list) < self.ptw_min_num_of_degs_cut:
			print(f"Warning in calc_EA_dataset_symbol(): for {self.case} - LFC_cut={self.LFC_cut:.2f} and lfc_FDR_cut=={self.lfc_FDR_cut:.2f}, please send a deg_list with >= {self.ptw_min_num_of_degs_cut} {self.s_deg_dap}s")
			return

		deg_list = list(deg_list)
		self.deg_list = deg_list
		df_enr0 = self.calc_get_enriched_pathway(deg_list, return_value=True, force=force, verbose=verbose)
		self.df_enr0 = df_enr0

		if (calc_many_sig or default) and df_enr0 is not None and not df_enr0.empty:
			self.calc_many_sig_enrich_pathways(default=default, force=force, verbose=verbose)
	
		return

	def calc_get_enriched_pathway(self, deg_list:list, return_value:bool=True,
								  force:bool=False, verbose:bool=False, i_sim:int=-1) -> pd.DataFrame:
		'''
		Calc enrichment analysis -> send call to Enrichr Web Service
		if i_sim == -1:
			no simulation
			if alrealdy exists read table if return_value else return pd.DataFrame()
		else:
			force has no meaning
			always send a call to the Web Service
		'''

		self.df_enr = pd.DataFrame()
		self.all_enr_degs, self.enr_found_degs, self.enr_not_found_degs = [], [], []

		fname, _ = self.set_enrichment_name()  # fname_cutoff

		if i_sim == -1:
			root_enrich = self.root_enrich
			filefull = osjoin(root_enrich, fname)

			if exists(filefull) and not force:  #  and not calc_many_sig
				if not return_value:
					return pd.DataFrame()

				return pdreadcsv(fname, self.root_enrich)
		else:
			#---------- is a simulation ------------------------
			root_enrich = self.root_enrich_random
			fname = fname.replace(".tsv", f"_random_{i_sim}.tsv")
			filefull = osjoin(root_enrich, fname)

		shortId, userListId = self.open_session_upload_symbols(deg_list)

		if shortId == '' or userListId == '':
			print(f"Problems in open_session_upload_symbols().")
			return pd.DataFrame()

		# print(f"Enrichr web service: for '{filefull}'")
		try:
			response = requests.get( self.ENRICHR_URL_query + self.ENRICHR_query_string%(userListId, self.geneset_lib) )
		except:
			print(f"Problems in request geneset {self.geneset_lib} and userListId {str(userListId)}.")
			return pd.DataFrame()

		if not response.ok:
			print('Error fetching enrichment results: %s in %s'%(userListId, self.geneset_lib))
			return pd.DataFrame()

		try:
			data = json.loads(response.text)
		except:
			print('Error getting data: loading json response.')
			return pd.DataFrame()

		dic = {}
		try:
			mats = data[self.geneset_lib]

			if len(mats) == 0:
				print("No results")
				return pd.DataFrame()

			for mat in mats:
				# print(mat)
				count = mat[0]
				desc  = mat[1]
				pval2  = mat[2]
				odds_ratio = mat[3]
				combined_score = mat[4]
				genes  = mat[5]
				fdr2  = mat[6]
				# unk0  = mat[7]
				# unk1  = mat[8]

				dic[count] ={}
				dic2 = dic[count]

				if self.geneset_lib == 'Reactome_2022':
					# Neuronal System R-HSA-112316
					mat = desc.split(' ')
					pathway = " ".join(mat[:-1])
					_id  = mat[-1]
					dic2['pathway'] = pathway
					dic2['pathway_id'] = _id
				elif self.geneset_lib.startswith('GO_'):
					# Positive Regulation of Cellular Process (GO:0048522)
					mat = desc.split(' ')
					pathway = " ".join(mat[:-1])
					_id  = mat[-1].strip().replace(')','').replace('(','')
					dic2['pathway'] = pathway
					dic2['pathway_id'] = _id
				elif self.geneset_lib.startswith('WikiPathways_'):
					# Small Cell Lung Cancer WP4658
					mat = desc.split(' ')
					pathway = " ".join(mat[:-1])
					_id  = mat[-1].strip()
					dic2['pathway'] = pathway
					dic2['pathway_id'] = _id
				else:
					dic2['pathway'] = desc

				dic2['pval'] = pval2
				dic2['fdr'] = fdr2
				dic2['odds_ratio'] = odds_ratio
				dic2['combined_score'] = combined_score
				dic2['genes'] = genes
				dic2['num_of_genes'] = len(genes)

		except:
			print('Error Enrichr response')
			return pd.DataFrame()

		if len(dic) == 0:
			return pd.DataFrame()

		df_enr = pd.DataFrame(dic).T

		if self.geneset_lib.startswith('Reactome') and self.geneset_lib != 'Reactome_2022':
			df_reactome = self.open_reactome_pathway()

			df_enr = pd.merge(df_enr, df_reactome[ ['pathway', 'pathway_id'] ], how="outer", on='pathway')
			cols = ['pathway', 'pathway_id', 'pval', 'fdr', 'odds_ratio', 'combined_score', 'genes', 'num_of_genes']
			df_enr = df_enr[cols]
			df_enr = df_enr[ ~pd.isnull(df_enr.fdr) ]

		df_enr = df_enr.sort_values(['fdr', 'num_of_genes'], ascending=[True, False])
		df_enr.reset_index(inplace=True, drop=True)
		_ = pdwritecsv(df_enr, fname, self.root_enrich, verbose=verbose)

		return df_enr


	def calc_many_sig_enrich_pathways(self, ptw_pval_cut:float=0.05, default:bool=False, 
									  force:bool=False, verbose:bool=False):
		'''
			loop on fdr_ptw_cutoff_list: 0.05 ... to ... 0.75
			loop on num_of_genes_list: default [3] one can be more flexible [2] or rigid [3, 4, 5]
			for each filter it filter and save the df_enr table
		'''

		if default:
			fdr_ptw_cutoff_list = [0.05]
		else:
			fdr_ptw_cutoff_list = self.fdr_ptw_cutoff_list

		if self.num_of_genes_list == []:
			self.num_of_genes_list = [3]

		for ptw_FDR_cut in fdr_ptw_cutoff_list:
			for ptw_min_num_of_degs_cut in self.num_of_genes_list:

				''' to get the correct fname_cutoff '''
				self.set_pathway_cutoff_params(ptw_FDR_cut, ptw_pval_cut, ptw_min_num_of_degs_cut)

				_, fname_cutoff = self.set_enrichment_name() # fname
				filefull = osjoin(self.root_enrich, fname_cutoff)

				if exists(filefull) and not force:
					continue

				df_enr = self.df_enr0[(self.df_enr0.fdr < ptw_FDR_cut) &
									  (self.df_enr0.pval < ptw_pval_cut) &
									  (self.df_enr0.num_of_genes >= ptw_min_num_of_degs_cut)].copy()


				if df_enr.empty:
					if verbose: 
						print("<<< Empty:", fname_cutoff, ptw_FDR_cut, ptw_pval_cut, ptw_min_num_of_degs_cut, len(df_enr))
					continue

				df_enr.reset_index(inplace=True, drop=True)
				_ = pdwritecsv(df_enr, fname_cutoff, self.root_enrich, verbose=verbose)


	''' old: set_enriched_pathway_line '''
	def get_enriched_pathway_line(self, i_line) -> object:
		'''
			assign:
				self.genes_in_pathway, self.pathway, self.pathway_id
		'''
		if not isinstance(self.df_enr, pd.DataFrame) or self.df_enr.empty:
			print("No self.df_enr was found or empty")
			return None

		try:
			row = self.df_enr.iloc[i_line]
		except:
			print("df_enr has not this line.")
			return None

		genes = row.genes
		if isinstance(genes, str):
			genes = eval(genes)

		self.genes_in_pathway = genes
		self.pathway = row.pathway
		try:
			self.pathway_id = row.pathway_id
		except:
			self.pathway_id = None

		return row

	def open_db_pathway(self, verbose:bool=False) -> pd.DataFrame:

		if self.geneset_lib.startswith('Reactome'):
			return self.open_reactome_pathway(verbose=verbose)
		
		if self.geneset_lib.startswith('KEGG'):
			return self.open_kegg_pathway(verbose=verbose)
		
		if self.geneset_lib.startswith('BioPlanet'):
			return self.open_bioplanet_pathway(verbose=verbose)

		print(f"DB {self.geneset_lib} is not developed")
		return pd.DataFrame()

	def open_reactome_pathway(self, verbose:bool=False) -> pd.DataFrame:
		fname = 'reactome_pathways_human.tsv'
		filefull = osjoin(self.root_reactome, fname)

		if not exists(filefull):
			print("Could not found: %s"%(filefull))
			return pd.DataFrame()

		df_reactome = pdreadcsv(fname, self.root_reactome, verbose=verbose)
		self.df_reactome = df_reactome
	
		return df_reactome

	def open_bioplanet_pathway(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:
		fname = 'pathway.tsv'
		filefull = osjoin(self.root_bioplanet, fname)

		if exists(filefull) and not force:
			df_biop = pdreadcsv(fname, self.root_bioplanet, verbose=verbose)
			df_biop.columns = ['pathway_id', 'pathway', 'gene_id', 'symbol']
			self.df_biop = df_biop
			return df_biop

		df_biop = pdreadcsv('pathway.csv', self.root_bioplanet, sep=',')
		cols = ['pathway_id', 'pathway', 'gene_id', 'symbol']
		df_biop.columns = cols

		df_biop = df_biop.drop_duplicates()
		df_biop.reset_index(inplace=True, drop=True)
		df_biop.columns = ['pathway_id', 'pathway', 'gene_id', 'symbol']

		_ = pdwritecsv(df_biop, fname, self.root_bioplanet)
		self.df_biop = df_biop

		return df_biop

	def open_bioplanet_category(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:
		cols = ['pathway_id', 'pathway', 'category']

		fname = 'pathway-category.tsv'
		filefull = osjoin(self.root_bioplanet, fname)

		if exists(filefull) and not force:
			dfbiop_cat = pdreadcsv(fname, self.root_bioplanet, verbose=verbose)
			dfbiop_cat.columns = cols
			self.dfbiop_cat = dfbiop_cat
			return dfbiop_cat

		dfbiop_cat = pdreadcsv('pathway-category.csv', self.root_bioplanet, sep=',')

		dfbiop_cat.columns = cols

		dfbiop_cat = dfbiop_cat.drop_duplicates()
		dfbiop_cat.reset_index(inplace=True, drop=True)
		_ = pdwritecsv(dfbiop_cat, fname, self.root_bioplanet, verbose=verbose)

		self.dfbiop_cat = dfbiop_cat
		return dfbiop_cat

	def open_bioplanet_disease(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		fname = 'pathway-disease-mapping.tsv'
		filefull = osjoin(self.root_bioplanet, fname)

		if exists(filefull) and not force:
			dfdisease = pdreadcsv(fname, self.root_bioplanet, verbose=verbose)
			self.dfdisease = dfdisease
			return dfdisease

		dfdisease = pdreadcsv(fname, self.root_bioplanet)
		cols = ['pathway_id', 'pathway', 'geneid', 'symbol', 'mim_id', 'disease', 'symbols', 'chr_location', 'phenotype']
		dfdisease.columns = cols

		dfdisease = dfdisease.drop_duplicates()
		dfdisease.reset_index(inplace=True, drop=True)
		_ = pdwritecsv(dfdisease, fname, self.root_bioplanet, verbose=verbose)

		self.dfdisease = dfdisease
		return dfdisease

	# A bit of code that will help us display the PDF output
	def PDF(self, filename):
		return HTML(f'<iframe src={filename} width=700 height=350></iframe>')

	# Some code to return a Pandas dataframe, given tabular text
	def to_df(self, result:str) -> pd.DataFrame:

		try:
			df = pd.read_table(io.StringIO(result), header=None)
		except:
			print("Could not read:", result)
			df = pd.DataFrame()

		return df

	def open_kegg_pathway(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		filefull = osjoin(self.root_kegg, self.fname_kegg_pathways)

		if exists(filefull) and not force:
			df_enr = pdreadcsv(self.fname_kegg_pathways, self.root_kegg, verbose=verbose)
			self.df_enr = df_enr
			return df_enr

		df_enr = self.get_rest_kegg_pathways()
		self.df_enr = df_enr

		return df_enr

	def get_kegg_pathway_image(self) -> object:
		try:
			return REST.kegg_get(self.pathid, "image").read()
		except:
			print("Error: problems with internet?")
			return None

	def get_kegg_kgml(self, force:bool=False, verbose:bool=False) -> bool:

		pathid_hsa = self.pathid.replace('map', 'hsa')
		self.pathid_hsa = pathid_hsa

		if self.pathid is None or self.pathid == '':
			print("pathid is not defined '%s'"%(str(self.pathid_hsa)))
			self.kgml, self.pathway_kgml = '', ''
			return False

		fname	 = "kgml_%s.xml"%(title_replace(self.pathway))
		filefull  = osjoin(self.root_enrich, fname)

		if exists(filefull) and not force:
			kgml = read_txt(fname, self.root_enrich, verbose=verbose)
			kgml = "\n".join(kgml)

			self.kgml = kgml
			self.pathway_kgml = KGML_parser.read(kgml)

			return True

		try:
			self.kgml = REST.kegg_get(pathid_hsa, "kgml").read()
			self.pathway_kgml = KGML_parser.read(self.kgml)
			ret = True
		except:
			print("Could not find KGML for '%s'"%(pathid_hsa))
			self.kgml, self.pathway_kgml = '', ''
			ret = False


		ret = write_txt(self.kgml, fname, self.root_enrich, verbose=True)

		return ret

	def get_kegg_human_gene_annotation(self, gene_id:str, force:bool=False, verbose:bool=False) -> bool:
		''' like "hsa:5624" '''
		kegg_id = gene_id
		self.kegg_id = kegg_id
		self.dfgc = pd.DataFrame()
		filefull = osjoin(self.root_kegg, self.fname_kegg_gene_comp)

		if exists(filefull) and not force:
			df = pdreadcsv(self.fname_kegg_gene_comp, self.root_kegg, verbose=verbose)
			dfgc = df[df.kegg_id == kegg_id].copy()

			if not dfgc.empty:
				dfgc.reset_index(inplace=True, drop=True)
				self.dfgc = dfgc
				return True

		if 'hsa' not in gene_id:
			self.gene_id = None
			print("This is not a human gene_id (without 'hsa'): '%s'"%(gene_id))
			return False

		try:
			result = REST.kegg_get(gene_id).read()
			self.kegg_id = gene_id
			self.result = result
			ret = self.parse_kegg_gene_compound_annotation(gene_id, result, force=force, verbose=verbose)
		except:
			print("Could not get REST.kegg_get(gene = '%s')"%(gene_id))
			self.kegg_id = ''
			self.result = ''
			ret = False

		return ret

	def get_kegg_compound_annotation(self, compound_id:str, force:bool=False, verbose:bool=False) -> bool:

		kegg_id = compound_id
		self.kegg_id = kegg_id
		self.dfgc = pd.DataFrame()
		filefull = osjoin(self.root_kegg, self.fname_kegg_gene_comp)

		if exists(filefull) and not force:
			df = pdreadcsv(self.fname_kegg_gene_comp, self.root_kegg)
			dfgc = df[df.kegg_id == kegg_id].copy()

			if not dfgc.empty:
				dfgc.reset_index(inplace=True, drop=True)
				self.dfgc = dfgc
				return True

		try:
			result = REST.kegg_get(compound_id).read()
			self.kegg_id = compound_id
			self.result = result
			ret = self.parse_kegg_gene_compound_annotation(compound_id, result, force=force, verbose=verbose)

		except:
			print("Could not get REST.kegg_get(compound = '%s')"%(compound_id))
			self.kegg_id = None
			self.result = None
			ret = False

		return ret

	def parse_kegg_gene_compound_annotation(self, kegg_id, result, force:bool=False, verbose:bool=False) -> bool:
		kegg_id = kegg_id
		result = result
		force = force
		verbose = verbose

		return False

	'''
	def parse_kegg_gene_compound_annotation_TO_REVIEW(self, kegg_id, result, force:bool=False, verbose:bool=False) -> bool:

		if not isinstance(result, str):
			print("Parse result must be a string")
			return False

		if len(result) <= 10:
			print("Parse result must exists: '%s'"%(result))
			self.dfgc = pd.DataFrame()
			return False

		filefull = osjoin(self.root_kegg, self.fname_kegg_gene_comp)

		if exists(filefull) and not force:
			df = pdreadcsv(self.fname_kegg_gene_comp, self.root_kegg, verbose=verbose)
			dfgc = df[df.kegg_id != kegg_id].copy()
		else:
			dfgc = pd.DataFrame()

		comp_list = result.split('\n')

		i = -1; dic={}; stop = False
		for comp in comp_list:
			comp = comp.strip()

			terms = comp.split(' ')
			j = -1
			for term in terms:
				if term == '///':
					stop = True
					break

				j += 1
				if term == '': continue
				if j == 0:
					if term in self.parse_gc_fields:
						i += 1
						dic[i] = {}
						dic2 = dic[i]
						dic2['term'] = term
						dic2['val'] = ''
						continue

					#print("**", i, ')', dic2['term'])
					mat = []
					for term2 in terms[j:]:
						if term2 != '': mat.append(term2)

					stri = " ".join(mat)

					if dic2['val'] == '':
						dic2['val'] = stri
						# print("### start", i, dic2['term'], dic2['val'])
					else:
						dic2['val'] += ' \n' + stri
						# print("### mergin", i, dic2['term'], dic2['val'])

				# print("$$ %d) '%s'"%(i, stri))

				break


			if stop: break

		df = pd.DataFrame(dic).T
		df['kegg_id'] = kegg_id
		fields = ['kegg_id', 'term', 'val']
		df = df[fields]

		if dfgc is None:
			dfgc = df
		else:
			dfgc = pd.concat([dfgc, df])

		dfgc = dfgc.sort_values(['kegg_id', 'term'])
		dfgc.reset_index(inplace=True, drop=True)
		_ = pdwritecsv(dfgc, self.fname_kegg_gene_comp, self.root_kegg)

		dfgc = df[df.kegg_id == kegg_id].copy()
		dfgc.reset_index(inplace=True, drop=True)
		self.dfgc = dfgc

		return True
	'''

	def find_kegg_gene_name_symbol(self, gene_id, force:bool=False, verbose:bool=False) -> Tuple[str, str, str]:

		ret = self.get_kegg_human_gene_annotation(gene_id, force=force, verbose=verbose)

		if ret:
			dfq = self.dfgc[ self.dfgc.kegg_id == gene_id ]
		else:
			dfq = pd.DataFrame()

		if not dfq.empty:
			gene_name	  = "; ".join( dfq[dfq.term == 'NAME'].val)
			gene_synonyms = "; ".join( dfq[dfq.term == 'SYMBOL'].val)
			gene_symbol   = gene_synonyms.split('; ')[0]
		else:
			gene_name, gene_symbol, gene_synonyms = '', '', ''

		return gene_name, gene_symbol, gene_synonyms

	def find_kegg_compound_name_type(self, compound_id, force:bool=False, verbose:bool=False) -> Tuple[str, str]:

		ret = self.get_kegg_compound_annotation(compound_id, force=force, verbose=verbose)

		if ret:
			dfq = self.dfgc[ self.dfgc.kegg_id == compound_id ]
		else:
			dfq = pd.DataFrame()

		if not dfq.empty:
			compound_name = "; ".join( dfq[dfq.term == 'NAME'].val)
			compound_type = "; ".join( dfq[dfq.term == 'TYPE'].val)

			compound_name = compound_name.replace("\n", ", ")
			compound_type = compound_type.replace("\n", ", ")
		else:
			compound_name, compound_type = '', ''

		return compound_name, compound_type


	def find_kegg_gene_diseases(self, gene_id_list, force:bool=False, verbose:bool=False) -> list:

		if not isinstance(gene_id_list, list):
			print("Error: gene_id_list must be a list")
			return []

		diseases = []
		for gene_id in gene_id_list:
			ret = self.get_kegg_human_gene_annotation(gene_id, force=force, verbose=verbose)
			if not ret: continue

			dfd = self.dfgc[ (self.dfgc.kegg_id == gene_id) & (self.dfgc.term == 'DISEASE')]
			if not dfd.empty:
				diseases.append(dfd.val)

		return diseases


	def return_gene_compound(self, node_names:List, node_type:str) -> Tuple[str, str, str, str, str]:

		s_gene_name, s_gene_symbol, s_gene_synonyms, s_compound_name, s_compound_type = '', '', '', '', ''

		# print(node_type, type(node_names), node_names)
		if node_type != 'gene' and node_type != 'compound':
			return s_gene_name, s_gene_symbol, s_gene_synonyms, s_compound_name, s_compound_type

		for kegg_id in node_names:
			gene_name, gene_symbol, gene_synonyms, compound_name, compound_type = None, None, None, None, None

			if node_type == 'compound':
				kegg_id = kegg_id[4:]  # without 'cpd:'
				compound_name, compound_type = self.find_kegg_compound_name_type(kegg_id)

				if s_compound_name == '':
					s_compound_name = compound_name
					s_compound_type = compound_type
				else:
					s_compound_name += '||' + compound_name
					s_compound_type += '||' + compound_type

			else:
				gene_name, gene_symbol, gene_synonyms = self.find_kegg_gene_name_symbol(kegg_id)

				if s_gene_name == '':
					s_gene_name = gene_name
					s_gene_symbol = gene_symbol
					s_gene_synonyms = gene_synonyms
				else:
					s_gene_name += '||' + gene_name
					s_gene_symbol += '||' + gene_symbol
					s_gene_synonyms += '||' + gene_synonyms


		return s_gene_name, s_gene_symbol, s_gene_synonyms, s_compound_name, s_compound_type


	# @enforce.runtime_validation
	def define_KEGG_pathway(self, ipath: int, verbose:bool = False) -> bool:
		pathway, genes = self.set_df_enr_ipath(ipath)

		if pathway is None:
			self.pathway = ''
			self.pathway_genes = []
			self.pathid = ''
			return False

		if verbose:
			print("pathway:", pathway)
			print("genes:",   ", ".join(genes))

		self.pathway = pathway
		self.pathway_genes = genes

		pathid = self.find_kegg_pathway_by_name(pathway)
		self.pathid = pathid

		if verbose:
			print(self.pathid, pathway)

		return True

	def calc_col_name(self, i:int) -> str:
		group = self.dfpa.iloc[i].group.strip()
		num   = int(self.dfpa.iloc[i].group_num)

		if num < 10:
			num = '0'+str(num)
		else:
			num = str(num)

		return group + '_' + num

	def open_proteomics_table_log2(self, fname:str, lim_elder:float=55, lim_obese:float=30,
								   verbose:bool=False) -> pd.DataFrame:
		'''
		fname = 'pacientes_x_amostra_proteomica.tsv'
		'''
		fnamefull = osjoin(self.root_data, fname)
		if not exists(fnamefull):
			print(f"File not found: '{fnamefull}'")
			return pd.DataFrame()

		dfpa = pdreadcsv(fname, self.root_data, verbose=verbose)

		dfpa.columns = ['group', 'group_num', '_id', 'sex', 'age', 'spec_id', 'bmi']
		dfpa['pac_id'] = [x.split('_')[1] for x in dfpa._id]
		dfpa['gender'] = [0 if x == 'F' else 1 for x in dfpa.sex]
		dfpa['elder']  = [0 if x < lim_elder else 1 for x in dfpa.age]
		dfpa['obese']  = [0 if x < lim_obese else 1 for x in dfpa.bmi]

		self.dfpa = dfpa
		dfpa['group_num'] = [ self.calc_col_name(i) for i in range(len(dfpa))]
		self.dfpa = dfpa

		return dfpa


	def open_log2_table_limma(self, fname:str, want_ctrl:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		fname = 'Log2_limma_Ctrl_x_%s.tsv'%(self.case)
		'''

		fnamefull = osjoin(self.root_data, fname)
		if not exists(fnamefull):
			print(f"File not found: '{fnamefull}'")
			return pd.DataFrame()

		dflimma = pdreadcsv(fname, self.root_data, verbose=verbose)

		if want_ctrl:
			cols = [x.replace('C_','ctrl_') if x.startswith('C_') else x for x in dflimma.columns]
			dflimma.columns = cols
			case2 = 'ctrl'
		else:
			case2 = self.case

		cols = [x for x in list(dflimma.columns) if case2 in x]

		dflimma = dflimma[ ['Majority.protein.IDs'] + cols ]
		dflimma.columns =  ['uniprot_id'] + cols
		self.dflimma = dflimma

		return dflimma

	def select_patients(self, _all:bool=True, only_obese:bool=False, obese:int=1,
					    only_elder:bool=False, elder:int=1, want_ctrl:bool=False) -> Tuple[pd.DataFrame, str]:

		case2 = 'ctrl' if want_ctrl else self.case

		dfpa_sel = self.dfpa.copy()

		if _all:
			dfpa_sel = dfpa_sel[(dfpa_sel.group == case2) ]
			s_title = 'all samples'
		else:
			if only_obese:
				s_title = 'obese' if obese == 1 else 'not obese'
				dfpa_sel = dfpa_sel[(dfpa_sel.group == case2) & (dfpa_sel.obese == obese)].copy()
			elif only_elder:
				s_title = 'elder' if elder == 1 else 'not elder'
				dfpa_sel = dfpa_sel[(dfpa_sel.group == case2) & (dfpa_sel.elder == elder)].copy()
			else:
				s_title = 'obese' if obese == 1 else 'not obese'
				s_title += ' and ' + 'elder' if elder == 1 else 'not elder'
				dfpa_sel = dfpa_sel[(dfpa_sel.group == case2) &(dfpa_sel.elder == elder) & (dfpa_sel.obese==obese)].copy()

		dfpa_sel.reset_index(inplace=True, drop=True)
		self.dfpa_sel = dfpa_sel

		return dfpa_sel, s_title


	def clean_select_patients(self):
		cols = list(self.dfpa_sel.group_num)

		cols = [x for x in cols if x in self.dflimma.columns]

		dflimma_sel = self.dflimma[  ['uniprot_id'] + cols ].copy()
		dflimma_sel = dflimma_sel.dropna(axis=0)
		self.dflimma_sel = dflimma_sel

		return dflimma_sel


	def calc_anova_patients(self) -> Tuple[float, float]:
		cols = list(self.dflimma_sel.columns)
		cols = cols[1:]
		df2 = self.dflimma_sel[cols]

		try:
			stat, pval = f_oneway(*df2)
		except:
			print('Error: could not calculate ANOVA')
			stat, pval = -1, 1.

		return stat, pval


	def calc_anova_patients_OLD(self) -> Tuple[float, float]:
		cols = list(self.dflimma_sel.columns)
		cols = cols[1:]

		df2 = self.dflimma_sel[cols]
		ncols = len(cols)

		stat, pval = -1, 1.

		if ncols == 3:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2])
		elif ncols == 4:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3])
		elif ncols == 5:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4])
		elif ncols == 6:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5])
		elif ncols == 7:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5], df2.iloc[:,6])
		elif ncols == 8:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5], df2.iloc[:,6], df2.iloc[:,7])
		elif ncols == 9:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5], df2.iloc[:,6], df2.iloc[:,7], df2.iloc[:,8])
		elif ncols == 10:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5], df2.iloc[:,6], df2.iloc[:,7], df2.iloc[:,8], df2.iloc[:,9])
		elif ncols == 11:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5], df2.iloc[:,6], df2.iloc[:,7], df2.iloc[:,8], df2.iloc[:,9], df2.iloc[:,10])
		elif ncols == 12:
			stat, pval = f_oneway(df2.iloc[:,0], df2.iloc[:,1], df2.iloc[:,2], df2.iloc[:,3], df2.iloc[:,4], df2.iloc[:,5], df2.iloc[:,6], df2.iloc[:,7], df2.iloc[:,8], df2.iloc[:,9], df2.iloc[:,10], df2.iloc[:,11])
		else:
			print(f'Error: wrong number of cols = ncols')

		return stat, pval


	def boxplot_patients(self, s_title:str, pval:float) -> object:

		cols = list(self.dflimma_sel.columns)
		cols = cols[1:]

		df2 = self.dflimma_sel[cols]

		fig = plt.figure(figsize=(12,8))
		ax = plt.axes()
		plt.boxplot(df2, showmeans=True, meanline=True)

		title = 'Boxplot for %s\n%s - ANOVA = %.2e'%(self.case, s_title, pval)

		plt.title(title)
		plt.xlabel('samples')
		plt.xticks(np.arange(1, len(cols)+1))
		ax.set_xticklabels(cols)
		plt.ylabel('log2(protein)')

		return fig

	def define_case(self, fname_limma:str, fname_samples:str, case:str, _all:bool=True, 
				    only_obese:bool=False, obese:int=1, only_elder:bool=False, elder:int=1, 
					lim_elder:float=55, lim_obese:float=3, verbose:bool=False) -> bool:
		'''
		fname_limma   = 'Log2_limma_Ctrl_x_%s.tsv'%(self.case)
		fname_samples = 'pacientes_x_amostra_proteomica.tsv'
		'''
		self.case = case

		dfpa = self.open_proteomics_table_log2(fname_samples, lim_elder=lim_elder, lim_obese=lim_obese, verbose=verbose)
		if dfpa.empty:
			return False

		print("Opened, groups:", dfpa.group.unique())

		pd_empty = pd.DataFrame()
		self.dflimma_ctrl, self.dfpa_sel_ctrl, self.dflimma_sel_ctrl = pd_empty, pd_empty, pd_empty
		self.dflimma_case, self.dfpa_sel_case, self.dflimma_sel_case = pd_empty, pd_empty, pd_empty

		want_ctrl = True
		dflimma = self.open_log2_table_limma(fname=fname_limma, want_ctrl=want_ctrl, verbose=verbose)
		if dflimma is None: 
			return False
		
		dfpa_sel, _ = self.select_patients(_all=_all, only_obese=only_obese, obese=obese, 
									 	   only_elder=only_elder, elder=elder, want_ctrl=want_ctrl)
		if dfpa_sel is None: return False
		dflimma_sel = self.clean_select_patients()
		if dflimma_sel is None or dflimma_sel.empty:
			return False

		self.dflimma_ctrl = dflimma
		self.dfpa_sel_ctrl = dfpa_sel
		self.dflimma_sel_ctrl = dflimma_sel

		want_ctrl = False
		dflimma = self.open_log2_table_limma(fname=fname_limma, want_ctrl=want_ctrl)
		if dflimma is None:
			return False
		
		dfpa_sel, _ = self.select_patients(_all=_all, only_obese=only_obese, obese=obese, 
									 	   only_elder=only_elder, elder=elder, want_ctrl=want_ctrl)
		if dfpa_sel is None: return False
		dflimma_sel = self.clean_select_patients()
		if dflimma_sel is None or dflimma_sel.empty:
			return False

		self.dflimma_case = dflimma
		self.dfpa_sel_case = dfpa_sel
		self.dflimma_sel_case = dflimma_sel

		return True

	def calc_elder_data(self, case:str, elder:int=1, verbose:bool=False) -> Tuple[pd.DataFrame, List]:
		if self.df_uniprot.empty:
			self.df_uniprot = self.gene.open_uniprot(verbose=verbose)

		self.case = case
		self.elder = elder

		s_elder = 'elder' if elder == 1 else 'not_elder'
		self.s_complement = s_elder

		fname = self.fname_new_proteomics%(case, 'ctrl', s_elder)
		filefull = osjoin(self.root_data, fname)

		if not exists(filefull):
			print("File does not exists '%s'"%(filefull))
			return pd.DataFrame(), []

		dfp = pdreadcsv(fname, self.root_data, verbose=verbose)

		dfn = pd.merge(dfp, self.df_uniprot, how='inner', on='uniID')
		print(len(dfp), len(dfn))
		dfn = dfn[(dfn.fdr < 0.05) & (np.abs(dfn.lfc) >= 1)]
		print(len(dfn))
		dfn = dfn.sort_values(['fdr', 'lfc'], ascending=[True, False])

		symbols = list(dfn.symbol)
		symbols.sort()

		return dfn, symbols

	def plot_venn2(self, lista1: List, lista2:List, names:List=['elder', 'adult'],
				   title0:str="Comparing Proteins elder x adult - case %s in dataset %s",
				   verbose:bool=True) -> Tuple[object, str]:

		set1 = set(lista1)
		set2 = set(lista2)

		inter_lista = list(set1.intersection(set2))
		inter_lista.sort()
		set1_only_lista = [x for x in set1 if x not in inter_lista]
		set1_only_lista.sort()
		set2_only_lista = [x for x in set2 if x not in inter_lista]
		set2_only_lista.sort()


		mat = get_venn_sections( (set1, set2) )

		vals = []; geneList =[]; index=[]
		dic = {}
		for i in range(len(mat)):
			genes = mat[i][1]
			index.append(mat[i][0])
			vals.append(len(mat[i][1]))
			geneList.append(mat[i][1])
			dic[mat[i][0]] = len(mat[i][1])

		if verbose:
			print(index)
			print(vals)
			print("")
			for key in dic.keys():
				print(key, dic[key])


		fig = plt.figure(figsize=(14,10))

		v2 = venn2( (set1, set2), set_labels = None, alpha=0.3)

		for i in range(len(mat)):
			# v2.get_patch_by_id(index[i]).set_color(i)
			v2.get_patch_by_id(index[i]).set_edgecolor('none')
			v2.get_label_by_id(index[i]).set_text('%s\n%d'%( defineClass(index[i], names), vals[i])) ### error???

			label = v2.get_label_by_id(index[i])
			label.set_fontsize(18)
			# label.set_family('arial')
			# label.set_x(label.get_position()[0] + 0.1)

		title = title0%(self.case, self.geneset_lib)
		plt.title(title, fontsize=20)

		fnamefig = title_replace(title)
		filefull = osjoin(self.root_figure, fnamefig + '.png')
		plt.savefig(filefull, dpi=300, format='png', facecolor='white')


		stri = "inter_lista (%d): %s\n"%(len(inter_lista), "; ".join(inter_lista))
		#--- elder
		stri += "\nonly elder (%d): %s\n"%(len(set1_only_lista), "; ".join(set1_only_lista))
		#--- adult
		stri += "\nonly adults (%d): %s\n"%(len(set2_only_lista), "; ".join(set2_only_lista))

		write_txt(stri, fnamefig+'.txt', self.root_figure, verbose=verbose)

		return fig, stri


	def calc_enriched_pathways_random_genes(self, i_sim:int, case:str, LFC_cut_default:float=1,
											lfc_FDR_cut_default:float=0.05, ptw_FDR_cut_default:float=0.05,
											prompt_verbose:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		calc DEGs with the default cutoff
		usually less than the best cutoff
		calc diff = n_degs_bca - n_degs_default
		fulfill the list of default degs with random diff not used genes
		calc the enriched analysis with the random fulfill gene list
		filter the enrichment analyis table with the same best cutoff
		calc in and not_in genes in pathway to perform the chi-square table

		'''
		# ret, degs_default, degs_ensembl_default, dflfc_default 
		_, _, _, _ = self.open_case_params(case, LFC_cut=LFC_cut_default,
										   lfc_FDR_cut=lfc_FDR_cut_default, 
										   ptw_FDR_cut=ptw_FDR_cut_default, verbose=verbose)
		self.degs_in_pathways_default = self.degs_in_pathways
		self.degs_not_in_pathways_default = self.degs_not_in_pathways
		self.degs_default = self.degs

		''' get best cutoff parameters'''

		# ret, degs_bca, degs_ensembl_bca, dflfc_bca 
		_, degs_bca, _, _= self.open_case(case)

		degs_in_pathways_bca	 = self.degs_in_pathways
		degs_not_in_pathways_bca = self.degs_not_in_pathways

		lista1 = degs_in_pathways_bca + degs_not_in_pathways_bca
		lista1.sort()

		lista2 = degs_bca
		lista2.sort()

		assert lista1 == lista2

		self.degs_in_pathways_bca	  = degs_in_pathways_bca
		self.degs_not_in_pathways_bca = degs_not_in_pathways_bca
		self.degs_bca = degs_bca

		self.n_degs_in_pathways_bca	    = len(degs_in_pathways_bca)
		self.n_degs_not_in_pathways_bca = len(degs_not_in_pathways_bca)
		self.n_degs_bca = len(degs_bca)

		n_genes = self.n_degs
		n_diff = n_genes - len(self.degs_default)
		if prompt_verbose: print(f"Randomizing {n_diff} = {n_genes} - {len(self.degs_default)} default DEGs")

		if n_diff == 0:
			flag_ok = True
			bca_has_more = True
			print("There is the same number of DEGs BCA and default. No improvement")
		else:
			flag_ok = False

			if n_diff > 0:
				bca_has_more = True
			else:
				bca_has_more = False
				print("DEGs BCA is less than DEGs default! No improvement")

		if flag_ok or not bca_has_more:
			self.degs_in_pathways_random, self.degs_not_in_pathways_random = [], []
			self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random = 0, 0
			return pd.DataFrame()

		""" remove the default and get only symbols with known ensembl_id """
		dfa = self.dflfc_ori.copy()
		dfa = dfa[ (~pd.isnull(dfa.ensembl_id)) & (~pd.isnull(dfa.symbol)) & \
				   (~dfa.symbol.isin(self.degs_in_pathways_default)) &
				   (~dfa.symbol.isin(self.degs_default)) ]
		dfa.reset_index(inplace=True, drop=True)

		rows = np.random.randint(0, len(dfa), n_diff)
		random_genes = list(dfa.iloc[rows].symbol)
		random_genes.sort()

		random_fulfill_genes = self.degs_default + random_genes
		if prompt_verbose: print(f"Calculating {len(random_fulfill_genes)} total default + random genes.")

		self.calc_random_EA_dataset_symbol(i_sim, random_fulfill_genes, verbose=False)

		if self.df_enr is None or self.df_enr.empty:
			if prompt_verbose: print("No pathway was found.")
			df_enr = pd.DataFrame()
			self.degs_in_pathways_random, self.degs_not_in_pathways_random = [], random_fulfill_genes
			self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random  = 0, len(random_fulfill_genes)
		else:
			df_enr = self.df_enr
			if prompt_verbose: print(f"There are {len(df_enr)} enriched pathways.")

			all_enr_degs = []
			for i in range(len(df_enr)):
				genes = df_enr.iloc[i].genes
				if isinstance(genes, str):
					genes = eval(genes)
				all_enr_degs += genes

			all_enr_degs = list(np.unique(all_enr_degs))

			degs_in_pathways_random	 = [x for x in random_fulfill_genes if x	 in all_enr_degs]
			degs_not_in_pathways_random = [x for x in random_fulfill_genes if x not in all_enr_degs]

			self.degs_in_pathways_random = degs_in_pathways_random
			self.degs_not_in_pathways_random = degs_not_in_pathways_random

			self.n_degs_in_pathways_random = len(degs_in_pathways_random)
			self.n_degs_not_in_pathways_random = len(degs_not_in_pathways_random)

		return df_enr


	def build_matrix_calc_chi_square(self, n_degs_in_pathways_bca:int, n_degs_not_in_pathways_bca:int,
									 n_degs_in_pathways_default:int, 
									 n_degs_not_in_pathways_random:int) -> Tuple[pd.DataFrame, bool, int, float, float, str]:

		mat = [ [n_degs_in_pathways_bca, n_degs_not_in_pathways_bca],
				[n_degs_in_pathways_default, n_degs_not_in_pathways_random] ]
	
		dfmat = pd.DataFrame(mat)
		dfmat.rename(index={0: 'BCA', 1: 'random'}, inplace=True)
		dfmat.columns = ['degs_in', 'degs_out']

		# ret_chi, dof, stat, pvalue, stri_stat = chisquare_2by2(dfmat)

		vals0 = dfmat.degs_in.to_list()
		vals1 = dfmat.degs_out.to_list()

		s_stat, statistic, pvalue, dof, _ = chi2_or_fisher_exact_test(vals0, vals1)

		ret = True if pvalue < 0.05  and pvalue >= 0 else False

		return dfmat, ret, dof, statistic, pvalue, s_stat


	def parallel_n_simulations(self, n_sim:int, case:str, LFC_cut_default:float=1,
					lfc_FDR_cut_default:float=0.05, ptw_FDR_cut_default:float=0.05) -> pd.DataFrame:

		dic = {}
		for i_sim in range(n_sim):
			# print(i_sim,end=' ')
			# df_enr
			print(".",end="")    
			_ = self.calc_enriched_pathways_random_genes(i_sim, case, LFC_cut_default, 
														lfc_FDR_cut_default, ptw_FDR_cut_default)

			dfmat, ret_chi, dof, stat, pvalue, stri_stat = self.build_matrix_calc_chi_square( \
								self.n_degs_in_pathways_bca, self.n_degs_not_in_pathways_bca, \
								self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random)
			if pvalue < 0:
				print(f"\n {i_sim}) Error", self.n_degs_in_pathways_bca, self.n_degs_not_in_pathways_bca,
										self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random)
			dic[i_sim] = {}
			dic2 = dic[i_sim]

			dic2['stat_sig'] = ret_chi
			dic2['pvalue'] = f'{pvalue:.3e}'
			dic2['stat'] = stat
			dic2['dof'] = dof
			dic2['stri_stat'] = stri_stat
			dic2['dfmat'] = dfmat

		if len(dic) == 0:
			return pd.DataFrame()

		dff = pd.DataFrame(dic).T
		return dff

	def run_parallel_n_simulations(self, n_sim:int, case:str, LFC_cut_default:float=1,
						  lfc_FDR_cut_default:float=0.05, ptw_FDR_cut_default:float=0.05) -> pd.DataFrame:

		t0 = time.time()

		nprocessors = multiprocessing.cpu_count()
		nprocessors = nprocessors - 2

		nprocessors = 1 if nprocessors <= 1 else nprocessors

		pool = Pool(nprocessors)

		lista_results = []

		n_sim_parallel = int(np.ceil(n_sim/nprocessors))

		n_sim_tot = 0
		for i in range(nprocessors):
			n_sim_tot += n_sim_parallel
			if n_sim_tot > n_sim:
				n_sim2 = n_sim_tot - n_sim
			else:
				n_sim2 = n_sim_parallel

			# print(i, n_sim2)
			print(".",end="")

			result = pool.apply_async(self.parallel_n_simulations, (n_sim2, case, LFC_cut_default, lfc_FDR_cut_default, ptw_FDR_cut_default) )
			lista_results.append(result)

			if n_sim_tot >= n_sim:
				break

		freeze_support()

		df_list = []
		for i in range(len(lista_results)):
			result = lista_results[i]
			dff = result.get(timeout=10000)
			df_list.append(dff)

		if df_list == []:
			dff = pd.DataFrame()
			print("\n\n---- dff Error - parallel returned empty -----------")
		else:
			dff = pd.concat(df_list)
			dff.reset_index(inplace=True, drop=True)

			seconds = time.time()-t0
			print(f"\n\n---- End: {seconds:.1f} sec -----------\n")

		return dff

	def run_n_simulations(self, n_sim:int, case:str, LFC_cut_default:float=1,
						  lfc_FDR_cut_default:float=0.05, ptw_FDR_cut_default:float=0.05,
						  force:bool=False, verbose:bool=False) -> pd.DataFrame:
		'''
		run n chi-square simulations
		todo: paralelize
		'''

		self.degs_in_pathways_bca,   self.degs_not_in_pathways_bca,   self.degs_bca = [], [], []
		self.n_degs_in_pathways_bca, self.n_degs_not_in_pathways_bca, self.n_degs_bca = -1, -1, -1
		self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random = -1, -1
		self.n_degs = -1

		fname = self.fname_enr_simulation%(case, self.normalization)
		filefull = osjoin(self.root_ressum, fname)

		if exists(filefull) and not force:
			return pdreadcsv(fname, self.root_ressum, verbose=verbose)


		dff = self.run_parallel_n_simulations(n_sim=n_sim, case=case, LFC_cut_default=LFC_cut_default,
						  lfc_FDR_cut_default=lfc_FDR_cut_default, ptw_FDR_cut_default=ptw_FDR_cut_default)
		'''
		dic = {}
		for i_sim in range(n_sim):
			print(i_sim,end=' ')
			# df_enr
			_ = self.calc_enriched_pathways_random_genes(i_sim, case, LFC_cut_default, 
														 lfc_FDR_cut_default, 
														 ptw_FDR_cut_default, 
														 prompt_verbose=prompt_verbose, verbose=verbose)

			dfmat, ret_chi, dof, stat, pvalue, stri_stat = self.build_matrix_calc_chi_square( \
								self.n_degs_in_pathways_bca, self.n_degs_not_in_pathways_bca, \
								self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random)
			if pvalue < 0:
				print(f"\n {i_sim}) Error", self.n_degs_in_pathways_bca, self.n_degs_not_in_pathways_bca,
										self.n_degs_in_pathways_random, self.n_degs_not_in_pathways_random)
				continue

			dic[i_sim] = {}
			dic2 = dic[i_sim]

			dic2['stat_sig'] = ret_chi
			dic2['pvalue'] = f'{pvalue:.3e}'
			dic2['stat'] = stat
			dic2['dof'] = dof
			dic2['stri_stat'] = stri_stat
			dic2['dfmat'] = dfmat

		dff = pd.DataFrame(dic).T
		'''

		if dff is not None and not dff.empty:
			_ = pdwritecsv(dff, fname, self.root_ressum, verbose=verbose)

		return dff


	def multiple_best_points(self, selected_toi_col:str='toi4_median', max_points:int=8,
							 with_params:bool=False, verbose:bool=False) -> pd.DataFrame:

		dic = {}
		for ipoint in range(1,max_points+1):
			dfconfig = self.calc_best_cutoffs_params(selected_toi_col=selected_toi_col, n_best_sample=ipoint, save_config=False, verbose=verbose)
			dfconfig = dfconfig[dfconfig.med_max_ptw == 'median']

			for j in range(len(dfconfig)):
				row = dfconfig.iloc[j]
				
				if ipoint == 1:
					dic[j] = {}
					dic2 = dic[j]
					dic2['case'] = row.case
				else:
					dic2 = dic[j]

				dic2[f'nptw_{ipoint}'] = row.n_pathways
				dic2[f'ndegs_{ipoint}'] = row.n_degs_in_pathways
				dic2[f'nptw_per_ndegs_{ipoint}'] = round(row.n_pathways / row.n_degs_in_pathways, 3)

				if with_params:
					dic2[f'fdrc_{ipoint}'] = row.lfc_FDR_cut
					dic2[f'lfcc_{ipoint}'] = row.LFC_cut
					dic2[f'fdrp_{ipoint}'] = row.ptw_FDR_cut
				
		dfpoints = pd.DataFrame(dic).T
		return dfpoints
	
	def calc_best_points(self, dfpoints:pd.DataFrame) -> Tuple[List, pd.DataFrame]:

		all_cols  = list(dfpoints.columns)
		cols_per  = [x for x in all_cols if x.startswith('nptw_per_ndegs')]
		cols_degs = [x for x in all_cols if x.startswith('ndegs_')]
		cols_ptws = [x for x in all_cols if x.startswith('nptw_') and len(x) < 10]

		dic={}
		n_best_sample_list=[]

		for i in range(len(dfpoints)):
			maxi = 0
			for j in range(len(cols_per)):
				col = cols_per[j]
				val = dfpoints.iloc[i][col]
				if val > maxi:
					cmax = col
					maxi = val
					dic[i] = {}
					dic[i]['case'] = dfpoints.iloc[i]['case']
					dic[i]['col'] = col
					dic[i]['sample'] = int(col.split('_')[-1])
					dic[i]['ptw_per_degs'] = val
					# print(j, col, cols_degs[j], cols_ptws[j] )
					dic[i]['ndegs'] = dfpoints.iloc[i][cols_degs[j]]
					dic[i]['nptw'] = dfpoints.iloc[i][cols_ptws[j]]

		dff = pd.DataFrame(dic).T
		n_best_sample_list = list(dff['sample'])

		return n_best_sample_list, dff


