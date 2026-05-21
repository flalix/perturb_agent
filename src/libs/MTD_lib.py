#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2023/06/25
# Udated  on 2024/10/11; 2024/05/06; 2024/03/29; 2023/08/28; 2023/08/16
# @author: Flavio Lichtenstein
# @local: Bioinformatics: CENTD/Molecular Biology; Instituto Butatan

import numpy as np
import os
import shutil
from pathlib import Path
from os.path import join as osjoin
from os.path import exists as exists
import pandas as pd
from typing import  Tuple, List  # Optional, Iterable, Set, Any
from datetime import datetime
import psutil

import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colors
from   matplotlib_venn import venn2 # venn2_circles

from pycirclize import Circos

import plotly.express as px
import plotly.graph_objects as go
from   plotly.subplots import make_subplots

from scipy.stats import norm
from scipy.stats import spearmanr

from markdown_pdf import Section
from markdown_pdf import MarkdownPdf


from libs.Basic import pdwritecsv, pdreadcsv, create_dir, all_equal_list, echo_print, read_txt, write_txt, dumpdic, loaddic
from libs.gene_lib import *
from libs.config_lib import *
from libs.stat_lib import *
from libs.reactome_lib import *
from libs.biomart_lib import *

from libs.graphic_lib import plotly_colors_proteins

from libs.calc_degs_lib import CALC_DEGS
from libs.tcga_gdc_lib import GDC

from project_context_GDC import load_project_context

ctx = load_project_context()

SUBTYPE_GENES=ctx.SUBTYPE_GENES
HISTOLOGY_GENES=ctx.HISTOLOGY_GENES
TUMOR_CLASS=ctx.TUMOR_CLASS
GLOBAL_SUBTYPE=ctx.GLOBAL_SUBTYPE
HISTOLOGY=ctx.HISTOLOGY
SITE_MAP=ctx.SITE_MAP
colors=ctx.colors


# print('recursionlimit', sys.getrecursionlimit())
# sys.setrecursionlimit(20000)
# print('recursionlimit', sys.getrecursionlimit())

class MTD(object):
	def __init__(self, disease:str, gene_protein:str, s_omics:str, project:str, s_project:str, 
			     root0:Path, root0_data:Path,
				 case_list:List, dic_case_list:dict, has_age:bool=True, has_gender:bool=True, exp_normalization:bool=False, 
				 std_filename:str='', std_filename_list:list=[],
				 geneset_num:int=0, ptw_min_num_of_degs_cut:int=3,
				 tolerance_pPMI:float=.15, s_pathw_enrichm_method:str='enricher',
				 LFC_cut_inf:float=0.40, fdr_ptw_cutoff_list:List=[],
				 num_of_genes_list:List=[3], lfc_list = [], fdr_list = [],
				 min_lfc_modulation:float=0.40, type_sat_ptw_index:str='linear_sat', 
				 saturation_lfc_param:float=5., enr_db_list:List=[], pPMI_normalized:bool=False):

		self.root0 = Path(root0)
		self.root_colab = create_dir(root0, 'colab')
		self.root_src   = create_dir(root0, 'src')

		self.project = project
		self.s_project = s_project		

		self.root0_data = Path(root0_data)
		self.root_project = create_dir(root0_data, s_project)
		self.disease = title_replace(disease)
		self.root_disease = create_dir(self.root_project, self.disease)

		self.root_result   = create_dir(self.root_disease, 'results')
		self.root_lfc      = create_dir(self.root_disease, 'lfc')
		self.root_data	   = create_dir(self.root_disease, 'data')
		self.root_enrich   = create_dir(self.root_disease, 'enrichment_analyses')
		self.root_ressum   = create_dir(self.root_disease, 'res_summ')
		self.root_figure   = create_dir(self.root_disease, 'figures')
		self.root_pathway  = create_dir(self.root_disease, 'pathway_summaries')
		self.root_config   = create_dir(self.root_disease, 'config')
		self.root_llm	   = create_dir(self.root_disease, 'llm')
		self.root_nc	   = create_dir(self.root_disease, 'non_coding')
		self.root_pdf_summ = create_dir(self.root_disease, 'pdf_summary')
		self.root_enrich_sampling = create_dir(self.root_disease, 'enrich_sampling')
		self.root_enrich_random  = create_dir(self.root_disease, 'enrichment_random')
		self.root_ptw_modulation = create_dir(self.root_disease, 'pathway_modulation')
		
		self.cfg  = Config(root0=root0, root_disease=self.root_disease, disease=self.disease, case_list=case_list)


		self.gene = Gene(root0=root0)
		self.med_max_ptw = 'median'

		self.s_pathw_enrichm_method = s_pathw_enrichm_method
		self.selected_pivot_pathway_list = []
		self.selected_pivot_symb_list = []

		lfc_list = list(lfc_list)
		if lfc_list == []:
			lfc_list = np.round(np.arange(1.0, -0.01, -.025), 3)
			lfc_list[-1] = 0

		self.lfc_list = lfc_list

		fdr_list = list(fdr_list)
		if fdr_list == []:
			fdr_list = np.arange(0.05, 0.76, .01)
		self.fdr_list = fdr_list

		self.dfsim = pd.DataFrame()

		self.type_sat_ptw_index   = type_sat_ptw_index
		self.saturation_lfc_param = saturation_lfc_param
		self.min_lfc_modulation   = min_lfc_modulation
		self.pPMI_normalized  = pPMI_normalized
		self.dff_pPMI = pd.DataFrame()

		self.max_pathways = 10
		self.abs_mod_diff_cutoff = 1.
		self.min_mod_diff_cutoff = 0.4
		self.highly_mod_cutoff = 3.0
		self.min_corr_sig = 0.7

		self.dfpiv_high_inter, self.dfpiv_high = pd.DataFrame(), pd.DataFrame()
		self.dfpiv_high_female, self.dfpiv_high_male = pd.DataFrame(), pd.DataFrame()

		self.gene_protein = gene_protein
		self.s_omics = s_omics
		self.s_gene_protein = 'Protein' if gene_protein == 'protein' else 'Gene'
		self.s_deg_dap = 'DAP' if gene_protein == 'protein' else 'DEG'

		self.has_gender = has_gender
		self.unique_gender_list = ['Female', 'Male'] if self.has_gender else ['unique']

		self.has_age = has_age

		self.case, self.group, self.age, self.gender = '','','',''
		self.icase = 0

		''' renaming bad symbols in dflfc_ori '''
		self.symbols2_list = []
		self.locs_list = []
		self.dfnot = None

		# fname_given_lfc: is the original LFC table given by the study
		# fname_final_lfc_ori: is the final LFC table, corrected with no duplications
		# open(fname_given_lfc) --> dflfc_all
		self.fname_given_lfc_table0 = f"{self.disease}_ALL_LFC_%s_x_CTRL_%s.tsv"
		# open(fname_final_lfc_ori) - corrected --> dflfc_ori
		self.fname_final_lfc_table0 = f"{self.disease}_final_LFC_%s_x_CTRL_%s.tsv"

		# lfc columns for given lfc table
		self.lfc_cols_default = ['ensembl_id', 'symbol', 'biotype', 'description', 'lfc', 'abs_lfc', 'pval', 'fdr']

		self.fname_lfc_nodup_table0 = f"{self.disease}_NO_DUP_LFC_%s_x_CTRL_%s.tsv"
		self.fname_final_lfc_ori = ''
		self.fname_all_dfenr = 'pathway_all_list_for_%s.tsv'

		self.fname_lfc_mod_summary = "lfc_modulation_summary_lfc_threshold_%.2f_diff_cutoff_%.2f.txt"

		''' enricher_GO_Biological_Process_2021_taubate_covid19_proteomics_for_g3_female_elder_x_ctrl_not_normalized_cutoff_lfc_0.200_fdr_0.700_pathway_pval_0.050_fdr_0.500_num_genes_3.tsv '''
		fname_enrich_table0 = f"{s_pathw_enrichm_method}_%s_{self.disease}_{s_omics}_for_%s_x_ctrl_%s"
		self.fname_enrich_table0 = fname_enrich_table0

		self.fname_enr_simulation = f"enr_simulation_{self.disease}_%s_%s.tsv"

		self.fname_stringdb  = ''

		self.case_list = case_list
		self.dic_case_list = dic_case_list

		self.my_colors=['navy', 'red', 'darkcyan', 'darkgreen', 'orange', 'brown', 'darksalmon',
		'magenta', 'darkturquoise', 'orange', 'darkred', 'indigo', 'magenta', 'maroon', 'black',
		'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgrey', 'olivedrab', 'navy'] + plotly_colors_proteins

		if has_gender:
			self.group_female_list = [x for x in case_list if '_female' in x]
			self.group_male_list   = [x for x in case_list if '_male' in x]
			self.list_order = self.group_female_list + self.group_male_list
			self.group_list = [x.replace('_female', '') for x in self.group_female_list]
		else:
			self.list_order = case_list
			self.group_list = case_list		
			self.group_female_list, self.group_male_list = [], []	

		self.group_colors = ['gray', 'blue', 'orange', 'red'] + plotly_colors_proteins
		self.group_colors = self.group_colors[:len(self.group_list)]

		''' expression normalization: if None -> not normalized '''
		self.exp_normalization = exp_normalization
		self.normalization = 'quantile_norm' if exp_normalization == True else 'not_normalized'

		self.std_filename = std_filename
		self.std_filename_list = std_filename_list

		self.set_db(geneset_num)
		self.geneset_num = geneset_num

		self.quantile_list  = np.round(np.arange(0, 1, 0.05), 2)
		self.LFC_cut_inf = LFC_cut_inf

		fdr_ptw_cutoff_list =  [np.round(x,2) for x in fdr_ptw_cutoff_list]
		if fdr_ptw_cutoff_list == []:
			fdr_ptw_cutoff_list = list(np.arange(0.05, 0.80, 0.05))
		self.fdr_ptw_cutoff_list = fdr_ptw_cutoff_list

		if num_of_genes_list is None or not isinstance(num_of_genes_list, list):
			num_of_genes_list = [3]

		self.num_of_genes_list = num_of_genes_list

		self.reset_degs_and_df_enr()

		self.df_summ = pd.DataFrame()
		self.dfsimb_perc =  pd.DataFrame()
		self.dfpiv, self.dfpiv_symbs_cases =  pd.DataFrame(),  pd.DataFrame()

		self.group_index = -1

		self.figure_svg, self.figure_stringdb = '', ''

		self.toggle_viewer = True
		self.hidden_urls = True

		self.fname_all_fdr_lfc_correlation = 'all_fdr_lfc_correlation.tsv_LFC_cut_inf_%.3f.json'
		self.fname_dic_fdr_lfc_correlation = 'dic_fdr_lfc_correlation_case_%s_LFC_cut_inf_%.3f.json'

		self.fname_degs_and_pathways_summary = "degs_and_pathways_summary.tsv"
		self.fname_degs_simulation = 'degs_simulation_for_%s_lfc.tsv'

		self.fname_txt_updegs = "upreg_degs_for_case_%s_lfc_%.2f_fdr_%.2f_%s.txt"
		self.fname_txt_dwdegs = "downreg_degs_for_case_%s_lfc_%.2f_fdr_%.2f_%s.txt"

		self.fname_enr_gene_stat = 'enrichment_pathways_gene_statistics_for_case_%s_and_database_%s_%s.tsv'

		self.fname_summ_comparing_2cases = 'comparing_2cases_pathw_%s_comparing_%s_x_%s.tsv'
		self.fname_pPMI = 'pPMI_%s_type_%s_saturation_%s_index_%s.tsv'
		self.fname_pathway_case_index = 'pathway_index_%s_%s_%s.tsv'
		self.fname_one_pathway_symb_LFC = 'pathway_%s_gene_modulations_%s_%s.tsv'

		self.fname_big_summary_txt = "%s_summary_all_cases_degs_geneset_%d_%s.txt"
		self.fname_degs_summary = 'summary_genes_per_case.tsv'
		self.fname_pathway_summary  = 'summary_genes_in_pathways_per_case_%s.tsv'
		self.fname_pathways_per_case = 'summary_pathways_per_case.tsv'
		self.fname_cutoff_table = 'cutoff_table_all_cases_quantiles_for_col_%s_for_geneset_%d_%s_%s.tsv'

		self.fname_summ_deg_ptw = 'summary_bca_x_default_cutoffs_result.tsv'

		self.max_iloop_plot_lines = 20

		self.param_defaults = 1, 0.05, -1, -1, -1, -1

		self.root_assets = create_dir(root0, 'assets')

		root_fig_md = self.root_disease / 'figures'
		if root_fig_md.exists():
			try:
				shutil.rmtree(root_fig_md)
			except FileNotFoundError:
				print(f"Folder {root_fig_md} not found.")
			except OSError as e:
				print(f"Error deleting folder {root_fig_md}: {e}")

		self.root_fig_md   = create_dir(self.root_disease, 'figures')
		self.curr_figname  = 'unknown_md.png'
		


		'''----- self.root_colab = where data is -------'''
		self.root_bioplanet = create_dir(self.root_colab, 'bioplanet')
		self.root_kegg		= create_dir(self.root_colab, 'kegg')
		self.root_refseq	= create_dir(self.root_colab, 'refseq')
		self.root_hgnc		= create_dir(self.root_colab, 'hgnc')

		''' ---- Affymetrix ---'''
		self.fname_affy = 'Human_Agilent_WholeGenome_4x44k_v2_MSigDB_v71.tsv'
		self.root_affymetrix = create_dir(self.root_colab, 'affymetrix')
		self.root_affy	   = create_dir(self.root_disease, 'affy')

		''' affymetrix experiment table - probe x symbols '''
		self.df_gpl = pd.DataFrame()

		''' ---- KEGG ---'''
		self.kegg_fname	= 'kegg_pathways.tsv'
		self.fname_kegg_pathways  = 'kegg_pathways.tsv'
		self.fname_kegg_gene_comp = 'kegg_gene_compound.tsv'

		self.root_owl= create_dir(self.root0_data, 'owl')
		self.root_reactome = create_dir(self.root0_data, 'reactome')

		self.pathway_id, self.pathway = 'DUMMY', 'DUMMY'

		self.df_lfc, self.dflfc_ori  = pd.DataFrame(), pd.DataFrame()
		self.dic_lfc, self.fig_lfc = {}, {}

		self.label, self.symbol = '', ''

		self.mim_gene_accession, self.mim_gene_description  = '', ''

		self.enr_db_list = enr_db_list

		self.ptw_min_num_of_degs_cut = ptw_min_num_of_degs_cut

		self.tolerance_pPMI = tolerance_pPMI


		self.reactome = Reactome(root_owl=self.root_owl, root_reactome=self.root_reactome)

		self.dfr = None
		''' below, df reactome with pathway_original '''
		self.df_reactome_gmt, self.df_enr_reactome = pd.DataFrame(), pd.DataFrame()

		# https://www.w3schools.com/css/tryit.asp?filename=trycss_table_border-spacing
		self.s_css  = """
h1 {text-align:center;font-family: Arial; color: DarkRed;} 
h2 {text-align:center;font-family: Arial; color:navy;} 
h3, h4 {text-align:left; font-family: Arial; color: black;} 
body {text-align:justify; font-family: Arial; color: black;} 
table {font-family: Arial; border-collapse: collapse; margin-left: auto; margin-right: auto; margin-top: 20px; margin-bottom: 20px;} 
th, td {border: 1px solid #ccc; padding: 8px; text-align: right; } 
th {background-color: #f2f2f2; font-weight: bold;} 
"""

	def translate_case(self, case:str) -> str:
		try:
			return self.dic_case_list[case]
		except:
			return case


	def set_which_model(self, which_model:str):
		self.which_model = which_model

	def open_reactome_gmt_for_pathway_analysis(self, verbose:bool=False) -> bool:
		if self.df_reactome_gmt is not None and not self.df_reactome_gmt.empty:
			return True

		df_gmt = self.reactome.open_reactome_gmt(verbose=verbose)
		if df_gmt is None or df_gmt.empty:
			self.df_reactome_gmt = pd.DataFrame()
			return False

		cols = ['pathway_id', 'pathway', 'genes', 'n']
		df_gmt = df_gmt[cols]
		df_gmt.columns =  ['pathway_id', 'pathway_original', 'genes_pathway', 'ngenes_pathway']
		self.df_reactome_gmt = df_gmt

		return True


	def merge_reactome(self, df_enr:pd.DataFrame) -> pd.DataFrame:

		if isinstance(self.df_enr_reactome, pd.DataFrame) and not self.df_enr_reactome.empty:
			return self.df_enr_reactome

		ret = self.open_reactome_gmt_for_pathway_analysis()
		if ret:
			dfn = pd.merge(df_enr, self.df_reactome_gmt, how="inner", on='pathway_id')
			self.df_enr_reactome = dfn
		else:
			self.df_enr_reactome = pd.DataFrame()
			dfn = pd.DataFrame()

		return dfn


	def reactome_find_genes_in_pathway(self, pathway:str, _type:str='pathway') -> Tuple[str, List]:
		'''
			_types: pathway, pathway_id
			return pathway -> pathway or pathway_id -> pathway and all genes
			if no results or > 1: return '', []
		'''

		ret = self.open_reactome_gmt_for_pathway_analysis()
		if not ret:
			return '', []

		if _type == 'pathway':
			dfa = self.df_reactome_gmt[self.df_reactome_gmt.pathway_original == pathway]

			if dfa.empty:
				dfa = self.df_reactome_gmt[self.df_reactome_gmt.pathway_original.str.contains(pathway)]
		else:
			dfa = self.df_reactome_gmt[self.df_reactome_gmt.pathway_id == pathway]

		if dfa.empty:
			print(f"Nothing found for {_type}=='{pathway}'")
			return '', []

		if len(dfa) > 1:
			print(f"There are {len(dfa)} results for {_type}=='{pathway}'")
			print(">>>", "; ".join(dfa.pathway_original))
			return '', []

		genes_in_pathway = dfa.iloc[0].genes_pathway
		if _type == 'pathway':
			s_pathway = dfa.iloc[0].pathway_id
		else:
			s_pathway = dfa.iloc[0].pathway_original

		if isinstance(genes_in_pathway, str):
			genes_in_pathway = eval(genes_in_pathway)

		genes_in_pathway = list(genes_in_pathway)
		return s_pathway, genes_in_pathway


	def rgb_to_hex(self, r, g, b):
		return '#{:02x}{:02x}{:02x}'.format(int(r), int(g), int(b))

	def calc_color(self, lfc, max_red = 4, min_blue = -4):

		perc = (lfc+max_red)/(max_red-min_blue)
		if perc > 1:
			perc = 1
		elif perc < 0:
			perc = 0

		pct_diff = 1.0 - perc
		blue_color = min(255, pct_diff*2 * 255)
		red_color = min(255, perc*2 * 255)

		pcg_green = perc if perc <= 0.5 else 1-perc
		green_color = 255 * 2*pcg_green

		# print(perc, red_color, green_color, blue_color)

		return self.rgb_to_hex(red_color, green_color, blue_color)



	def check_lfc_names(self, verbose:bool=False):
		
		for case in self.case_list:
			_, _, _, _ = self.open_case(case, verbose=False)
			if verbose:
				print("\nEcho Parameters:")
				self.echo_parameters()
				print("")
			
			fname_final, _, _ = self.set_lfc_names()
			filename = osjoin(self.root_lfc, fname_final)

			if exists(filename):
				print(">>> lfc file exists:", filename, '\n')
				dfq = pdreadcsv(fname_final, self.root_lfc)
				cols = list(dfq.columns)

				print("Checking columns:")
				for col in self.lfc_cols_default:
					if col in cols:
						print('\t', col, 'ok')
					else:
						print('\t', col, '???')
				print("\n")
			else:
				print("Error: lfc file does not exists",  filename, '\n')


	def set_lfc_names(self):

		if self.std_filename:
			fname_given_lfc		= self.fname_given_lfc_table0%(self.case, self.normalization)
			fname_final_lfc_ori = self.fname_final_lfc_table0%(self.case, self.normalization)
		else:
			fname_given_lfc = self.std_filename_list[self.icase]
			fname_final_lfc_ori = fname_given_lfc.replace('_ALL_', '_final_')

		self.fname_given_lfc = fname_given_lfc
		self.fname_final_lfc_ori = fname_final_lfc_ori

		if (self.age == '') and (self.gender == ''):
			title = self.group

		elif self.age == '':
			title = f'{self.group} {self.gender}'

		elif self.gender == '':
			title = f'{self.group} {self.age}'
		else:
			title = f'{self.group} {self.gender} {self.age}'

		title += f" ({self.normalization})"

		return fname_final_lfc_ori, fname_given_lfc, title

	def review_proteomics_table(self, fname_final_lfc_ori:str, fname_given_lfc:str, verbose:bool=False) -> pd.DataFrame:
		'''
			fname_given_lfc: is the original LFC table given by the study
			fname_final_lfc_ori: is the final LFC table, corrected with no duplications
		'''
		filename = self.root_lfc / fname_given_lfc

		if not filename.exists():
			print(f"Error: could not find table '{filename}'")
			raise Exception('stop: fix dup_dflfc_ori()')

		dflfc_ori = pdreadcsv(fname_given_lfc, self.root_lfc, verbose=verbose)

		if dflfc_ori is None or dflfc_ori.empty:
			print(f"Error: could not find data for table '{fname_given_lfc}'")
			raise Exception('stop: fix dup_dflfc_ori()')

		cols = list(dflfc_ori.columns)
		if "abs_lfc" not in cols:
			dflfc_ori["abs_lfc"] = np.abs(dflfc_ori.lfc)

		if "symbol_pipe" not in cols:
			dflfc_ori["symbol_pipe"] = None

		dflfc_ori = self.remove_repeated_according_max_LFC(dflfc_ori)

		if 'ensembl_id' in cols and 'biotype' in cols:
			pass
		else:
			dflfc_ori = self.biomart_fix_dflfc_ori(dflfc_ori)

		pdwritecsv(dflfc_ori, fname_final_lfc_ori, self.root_lfc, verbose=True)

		return dflfc_ori

	def remove_repeated_according_max_LFC(self, dflfc_ori:pd.DataFrame) -> pd.DataFrame:
		'''
		remove_repeated_according_max_LFC()
			add: symbol_prev
			fix symbol with replace_symbol_to_synonym()
			remove repeated symbols, conserve the highest abs(LFC)
		'''

		if 'symbol_prev' not in dflfc_ori.columns:
			dflfc_ori['symbol_prev'] = dflfc_ori.symbol
			dflfc_ori.loc[:, 'symbol'] = [self.gene.replace_symbol_to_synonym(x) for x in dflfc_ori.symbol]

		dflfc_ori = dflfc_ori.sort_values(['symbol', 'abs_lfc'], ascending=[True, False]).copy()
		dflfc_ori.reset_index(inplace=True, drop=True)

		''' remove duplicates - stay with the biggest LFC'''
		previous = ''; goods = []
		for i in range(len(dflfc_ori)):

			if not isinstance(dflfc_ori.iloc[i].symbol, str):
				goods.append(False)
			elif dflfc_ori.iloc[i].symbol != previous:
				previous = dflfc_ori.iloc[i].symbol
				goods.append(True)
			else:
				goods.append(False)

		dflfc_ori = dflfc_ori[goods].copy()
		dflfc_ori.reset_index(inplace=True, drop=True)

		return dflfc_ori


	def biomart_fix_dflfc_ori(self, dflfc_ori:pd.DataFrame) -> pd.DataFrame:

		for symbol_new, symbol_old in self.symbols2_list:
			dflfc_ori.loc[ dflfc_ori.symbol == symbol_old, 'symbol'] = symbol_new

		print(">>> running biomart_fix_dflfc_ori() ...")
		bm = Biomart()

		dfbm = bm.open_biomart_hsapiens()

		if dfbm is None or dfbm.empty:
			print("Please run biomart_wget.ipynb -> download_biomart_hsapiens()")
			raise Exception("stop: biomart_fix_dflfc_ori()")

		cols_dfbm = ['symbol', 'biotype', 'description', 'ensembl_transcript_id', 'ensembl_gene_id',
					 'chromosome', 'start_position', 'end_position']

		# remove description
		cols_dflfc = ['entry_id', 'symbol', 'uniprot_name', 'lfc', 'abs_lfc', 'pval', 'fdr',
					  'mean_exp', 't', 'B', 'symbol_pipe', 'symbol_prev']

		dfn = pd.merge(dflfc_ori[cols_dflfc], dfbm[cols_dfbm], how='inner', on='symbol')

		''' could not map '''
		dfnot = dflfc_ori[~dflfc_ori.symbol.isin(dfn.symbol)].copy()

		if len(dfnot) > 0:
			print(f"There are {len(dfnot)} genes not mapped in dfnot.")

			''' add new columns to merge '''
			cols = ['biotype', 'ensembl_transcript_id', 'ensembl_gene_id', 'chromosome', 'start_position', 'end_position']
			for col in cols:
				dfnot[col] = None

			dflfc_new = pd.concat([dfn, dfnot])
		else:
			print(f"All symbols were mapped in df biomart.")
			dflfc_new = dfn

		dflfc_new.reset_index(inplace=True, drop=True)

		# renaming to ensembl_id
		cols_new = ['entry_id', 'symbol', 'uniprot_name', 'lfc', 'abs_lfc', 'pval', 'fdr', 'mean_exp', 't', 'B',
					'symbol_pipe', 'symbol_prev', 'biotype', 'description', 'ensembl_transcript_id',
					'ensembl_id', 'chromosome', 'start_position', 'end_position']

		dflfc_new.columns = cols_new

		cols_order = ['entry_id', 'symbol', 'ensembl_id', 'uniprot_name', 'biotype', 'description', 'lfc', 'abs_lfc', 'pval', 'fdr',
					  'mean_exp', 't', 'B', 'symbol_pipe', 'symbol_prev', 'ensembl_transcript_id',
					  'chromosome', 'start_position', 'end_position']
		dflfc_new = dflfc_new[cols_order]

		''' fixing locs_list -= mannually mapped LOCs '''
		for symbol, ensembl_id, biotype in self.locs_list:
			print(">>> changing", symbol, ensembl_id, biotype)

			n = len( dflfc_new[dflfc_new.symbol == symbol] )
			dflfc_new.loc[ dflfc_new.symbol == symbol, ('ensembl_id', 'biotype')] = [ensembl_id, biotype]*n

		dfa = dflfc_new[ pd.isnull(dflfc_new.ensembl_id)]
		if dfa.empty:
			print("All genes are well mapped according to ensembl_id")

			dflfc_new = self.remove_repeated_according_max_LFC(dflfc_new)
		else:
			print(f"It remaines {len(dfa)} genes not well mapped according to ensembl_id")

		self.dfnot = dfa

		return dflfc_new

	def open_dflfc_ori(self, verbose:bool=False):
		'''
			read final dflfc_ori table --> the corrected one

			if not found:
				if self.s_omics == 'microarray':
					_ = self.review_LFC_table_with_affy_annot_wo_case(force=False, calc_interm_tables=False, verbose=False)

				elif self.s_omics == 'proteomics':
					_ = self.review_proteomics_table(fname_final_lfc_ori, fname_given_lfc, verbose=verbose)
			
			in the end:
				dflfc_ori = pdreadcsv(fname_final_lfc_ori, self.root_lfc, verbose=verbose)
		'''

		fname_final_lfc_ori, fname_given_lfc, _ = self.set_lfc_names()
		filename = self.root_lfc / fname_final_lfc_ori

		if not filename.exists():
			if self.s_omics == 'microarray':
				_ = self.review_LFC_table_with_affy_annot_wo_case(force=False, calc_interm_tables=False, verbose=False)

			elif self.s_omics == 'proteomics':
				_ = self.review_proteomics_table(fname_final_lfc_ori, fname_given_lfc, verbose=verbose)

			elif self.s_omics == 'RNA-Seq':
				pass

			else:
				print(f"There is no fix duplication method to {self.s_omics}")

		if not filename.exists():
			_ = self.import_from_GDC(prog_id="TCGA", force=False, verbose=verbose)

			if not filename.exists():
				print(f"Error: could not find {filename}")
				self.dflfc_ori = pd.DataFrame()
				return False

		# read final dflfc_ori table --> the corrected one
		dflfc_ori = pdreadcsv(fname_final_lfc_ori, self.root_lfc, verbose=verbose)

		changed = False
		if 'abs_lfc' not in dflfc_ori.columns:
			dflfc_ori['abs_lfc'] = [np.abs(x) if not pd.isnull(x) else x for x in dflfc_ori.lfc]
			changed = True

		if 'gene_id' in dflfc_ori.columns:
			dflfc_ori = dflfc_ori.rename(columns={"gene_id": "ensembl_id"})
			changed = True

		if 'geneid' in dflfc_ori.columns:
			dflfc_ori = dflfc_ori.rename(columns={"geneid": "ensembl_id"})
			changed = True

		if 'gene_type' in dflfc_ori.columns:
			dflfc_ori = dflfc_ori.rename(columns={"gene_type": "biotype"})
			changed = True

		if changed:
			_ = pdwritecsv(dflfc_ori, fname_final_lfc_ori, self.root_lfc, verbose=True)

		if verbose: print(f">>> dflfc_ori: contains {len(dflfc_ori)} probes.")

		goods = [True if isinstance(x, str) else False for x in dflfc_ori.symbol]
		dflfc_ori = dflfc_ori[goods].copy()
		dflfc_ori.reset_index(inplace=True, drop=True)

		self.dflfc_ori = dflfc_ori
		self.valid_genes = list(dflfc_ori.symbol)
		self.n_total_genes = len(dflfc_ori)

		if verbose: print(f">>> dflfc_ori: has {self.n_total_genes} valid symbols.")

		ret = self.check_dflfc_ori_duplicates()
		if not ret:
			print(f"Error: problems with table: '{fname_final_lfc_ori}'")


		return ret

	def check_dflfc_ori_duplicates(self):
		dfg = self.dflfc_ori.groupby('symbol').count().reset_index().iloc[:, :2]
		dfg.columns = ['symbol', 'n']
		dfn = dfg[dfg.n > 1]

		self.dfn = dfn
		if not dfn.empty:
			print("There are repeated symbols in dflfc_ori, see self.dfn")
			print("Rerun: review_LFC_table_with_affy_annot_wo_case(force=True, calc_interm_tables=True, verbose=False)")
			return False

		return True


	def open_df_lfc_complete(self, verbose:bool=False):

		self.dflfc_ori = pd.DataFrame()

		fname_final_lfc_ori, _, title = self.set_lfc_names()
		self.title = title

		filename = osjoin(self.root_lfc, fname_final_lfc_ori)

		if not exists(filename):
			df = pd.DataFrame({})
			print('LFC table does not exists: %s'%(filename))
			return False

		if verbose: print("reading df_lfc ori", self.group, self.gender, self.age)
		df = pdreadcsv(fname_final_lfc_ori, self.root_data)

		''' microarray columns from limma'''
		cols1 = ['probe', 'symbol', 'geneid', 'description', 'logFC', 'meanExpr',
				 't.stat', 'p-value', 'fdr', 'B', 'chr.range', 'org.chromosome',
				 'forward.reverse', 'nuc.sequence', 'gemmaid', 'go.term']

		cols3 = ['entry_id', 'symbol', 'uniprot_name', 'description', 'lfc', 'pval', 'fdr', 'mean_exp', 't', 'B']

		if all_equal_list(df.columns, cols1):

			cols = ['probe', 'symbol', 'gene_id', 'description', 'lfc', 'mean_expr',
					'tstat', 'pval', 'fdr', 'B', 'chr_ange', 'org_hromosome',
					'forward_everse', 'nuc_equence', 'gemma_id', 'go_term']
			df.columns = cols

			df['abs_lfc'] = [np.abs(x) if not pd.isnull(x) else x for x in df.lfc]

			df['multi_symbol'] = [x for x in df.symbol]
			df['symbol'] = [x.split('|')[0] if not pd.isnull(x) and isinstance(x, str) else x for x in df.symbol]

			cols = ['probe', 'symbol', 'gene_id', 'description', 'fdr', 'abs_lfc', 'lfc', 'pval', 'multi_symbol', 'mean_expr',
				   'tstat', 'B', 'chr_ange', 'org_hromosome', 'forward_everse', 'nuc_equence', 'gemma_id', 'go_term']
			df = df[cols]

			df = df.sort_values(['fdr', 'abs_lfc'], ascending=[True, False])
			df.reset_index(inplace=True, drop=True)

			pdwritecsv(df, fname_final_lfc_ori, self.root_data, verbose=verbose)

		elif all_equal_list(df.columns, cols3):

			df['abs_lfc'] = [np.abs(x) if not pd.isnull(x) else x for x in df.lfc]
			df['symbol_pipe'] = df.symbol
			df['symbol'] = [x.split('|')[0] if isinstance(x, str) else None for x in df.symbol]

			cols = ['entry_id', 'symbol', 'uniprot_name', 'description', 'lfc', 'abs_lfc', 'pval', 'fdr', 'mean_exp', 't', 'B', 'symbol_pipe']
			df = df[cols]

			df = df.sort_values(['fdr', 'abs_lfc'], ascending=[True, False])
			df.reset_index(inplace=True, drop=True)

			pdwritecsv(df, fname_final_lfc_ori, self.root_data, verbose=verbose)
		else:
			if 'abs_lfc' not in df.columns:
				print(f"Review the LFC table, please: '{fname_final_lfc_ori}'")
				raise Exception('stop: LFC table')

		if df is None or df.empty:
			df  = None; ret = False
			self.dic_lfc, self.fig_lfc = {}, {}
			self.dflfc_ori = pd.DataFrame()
			print('lfc table is empty')
		else:
			self.dic_lfc  = df.to_dict('records')
			self.fig_lfc  = px.bar(df, x='symbol', y='lfc')

			df = df[ ~pd.isnull(df.symbol)].copy()

			''' duplicated probes '''
			df = df.sort_values(['symbol', 'abs_lfc'], ascending=[True, False])
			df = df.drop_duplicates('symbol')

			df = df.sort_values(['fdr', 'abs_lfc'], ascending=[True, False])
			df.reset_index(inplace=True, drop=True)
			if verbose: print(f'lfc table length = {len(df)}')

			self.dflfc_ori = df

		return True



	def set_which_db(self, geneset_lib):
		self.geneset_lib = geneset_lib
		mat = [i for i in range(len(self.enr_db_list)) if self.enr_db_list[i] == geneset_lib ]

		if mat == []:
			mat = [i for i in range(len(self.enr_db_list)) if geneset_lib.lower() in self.enr_db_list[i].lower() ]

		try:
			self.geneset_num = mat[0]
		except:
			self.geneset_num = -1

		self.set_db(self.geneset_num)

		return self.geneset_num

	def set_db(self, geneset_num: int, verbose:bool=False):
		self.geneset_num = geneset_num

		if geneset_num == -1:
			self.geneset_lib = 'Undefined'
		else:
			try:
				self.geneset_lib = self.enr_db_list[geneset_num]
			except:
				self.geneset_lib = 'Undefined'

		if verbose: print(f">>> {self.geneset_lib}")

		return

	def set_enrichment_name(self) -> Tuple[str, str]:
		if self.geneset_lib is None:
			self.set_db(self.geneset_num)

		fname = self.fname_enrich_table0%(self.geneset_lib, self.case, self.normalization)
		fname = title_replace(fname)

		s_LFC_cut = f"{self.LFC_cut:.3f}"
		s_lfc_FDR_cut  = f"{self.lfc_FDR_cut:.3f}"
		s_pval_pathway = f"{self.ptw_pval_cut:.3f}"
		s_fdr_pathway  = f"{self.ptw_FDR_cut:.3f}"

		fname += f"_cutoff_lfc_{s_LFC_cut}_fdr_{s_lfc_FDR_cut}.tsv"
		self.fname_enrich_table = fname

		fname_cutoff = fname.replace('.tsv', f"_pathway_pval_{s_pval_pathway}_fdr_{s_fdr_pathway}_num_genes_{self.ptw_min_num_of_degs_cut}.tsv")

		return fname, fname_cutoff

	def open_enrichment_analysis(self, save_EP_xls:bool=False, 
								 force:bool=False, verbose:bool=False) -> bool:
		self.df_enr, self.pathway_in = pd.DataFrame(), False

		fname, _ = self.set_enrichment_name()
		filename = osjoin(self.root_enrich, fname)
		if verbose: print(">>>> open enrichment_analysis(): ", filename)

		if not exists(filename):
			if verbose: print(f"Warning: EA table for {self.case} does not exist: '{filename}'")
			return False

		df_enr0 = pdreadcsv(fname, self.root_enrich, verbose=verbose)

		if df_enr0 is None or df_enr0.empty:
			print(f"EA table is empty: '{filename}'")
			return False

		cols_ori = list(df_enr0.columns)

		if not 'num_of_genes' in df_enr0.columns or force:
			if 'num_genes' in df_enr0.columns:
				cols = ['pathway', 'pathway_id', 'pval', 'fdr', 'odds_ratio', 'combined_score', 'genes', 'num_genes']

				if all_equal_list(cols, cols_ori[:len(cols)]):
					df_enr0 = df_enr0[cols]
					cols = ['pathway', 'pathway_id', 'pval', 'fdr', 'odds_ratio', 'combined_score', 'genes', 'num_of_genes']
					df_enr0.columns = cols
					ret = pdwritecsv(df_enr0, fname, self.root_enrich, verbose=verbose)
			else:
				cols = ['pathway', 'overlap', 'pval', 'fdr', 'old_pval', 'old_fdr', 'odds_ratio', 'combined_score', 'genes']

				if len(df_enr0.columns) == len(cols):
					df_enr0.columns = cols
				elif len(df_enr0.columns) == len(cols)-2:
					cols = ['pathway', 'overlap', 'pval', 'fdr', 'odds_ratio', 'combined_score', 'genes']
					df_enr0.columns = cols

				try:
					if ";" in df_enr0.iloc[0].genes:
						mat = [0 if not isinstance(x, str) or x == '' else len(x.split(';')) for x in df_enr0.genes]
					else:
						mat = [0 if not isinstance(x, str) or x == '[]' else len(eval(x)) for x in df_enr0.genes]
				except:
					print("Error: please review df_enr0 structure and columns genes")
					raise Exception('Stop: open enrichment_analysis, column genes')

				df_enr0['num_of_genes'] = mat
				ret = pdwritecsv(df_enr0, fname, self.root_enrich, verbose=verbose)

		else:
			ret = 1
			for col in ['pathway', 'pval', 'fdr', 'genes', 'num_of_genes' ]:
				ret *= col in df_enr0.columns

			if ret != 1:
				cols = list(df_enr0.columns)
				print(f"Rename manually the enriched table, wrong columns: {', '.join(cols)}")
				raise Exception('stop')

		self.df_enr0 = df_enr0
		''' calculates degs in pathways '''
		self.df_enr = self.calc_sig_enriched_pathway(save_EP_xls, verbose)

		return True


	def calc_sig_enriched_pathway(self, save_EP_xls:bool=False, verbose:bool=False) -> pd.DataFrame:
		''' calculates degs in pathways '''
		df_enr = self.df_enr0.copy()

		df_enr = df_enr[ (df_enr.pval < self.ptw_pval_cut) &
						 (df_enr.fdr  < self.ptw_FDR_cut) &
						 (df_enr.num_of_genes >= self.ptw_min_num_of_degs_cut) ]
	
		df_enr.reset_index(inplace=True, drop=True)

		if df_enr.empty:
			return df_enr

		self.calc_enrichment_parameters(df_enr)

		if save_EP_xls:
			_, fname = self.set_enrichment_name()
			pdwritecsv(df_enr, fname, self.root_result, verbose=verbose)

			'''
			pip install openpyxl
			from openpyxl import Workbook
			'''
			fname = self.fname_enrich_table0%(self.geneset_lib, self.case, self.normalization) + '.xlsx'
			filename = osjoin(self.root_result, fname)
			df_enr.to_excel(filename, sheet_name=self.case, index=False)

		return df_enr

	def calc_enrichment_parameters(self, df_enr:pd.DataFrame):
		self.df_enr = df_enr

		if df_enr is None or df_enr.empty:
			self.clear_degs_vars()
			return

		degs_in_pathways = []
		pathway_list, pathway_id_list, pathway_fdr_list = [], [], []
		num_of_genes_in_pathway_list = []
		for i in range(len(df_enr)):
			row = df_enr.iloc[i]
			genes = row.genes
			if isinstance(genes, str):
				genes = eval(genes)
			degs_in_pathways += list(genes)

			pathway_list.append(row.pathway)
			num_of_genes_in_pathway_list.append(row.num_of_genes)

			try:
				# now, only Ractome has pathway_id
				pathway_id_list.append(row.pathway_id)
			except:
				pass
			pathway_fdr_list.append(row.fdr)

		degs_in_pathways = list(np.unique(degs_in_pathways))

		self.degs_in_pathways		  = degs_in_pathways
		self.degs_ensembl_in_pathways = [x for x in degs_in_pathways if x in self.degs_ensembl]
		self.degs_not_in_pathways	  = [x for x in self.degs		 if x not in degs_in_pathways]
		self.degs_not_in_pathways.sort()
		self.degs_ensembl_not_in_pathways = [x for x in self.degs_ensembl if x not in degs_in_pathways]
		self.degs_ensembl_not_in_pathways.sort()

		# without ensembl
		self.degs_up_not_ensembl = [x for x in self.degs_up if x not in self.degs_up_ensembl]
		self.degs_up_not_ensembl.sort()
		self.degs_dw_not_ensembl = [x for x in self.degs_dw if x not in self.degs_dw_ensembl]
		self.degs_dw_not_ensembl.sort()


		self.pathway_list	  = pathway_list
		self.pathway_id_list  = pathway_id_list
		self.pathway_fdr_list = pathway_fdr_list
		self.num_of_genes_in_pathway_list = num_of_genes_in_pathway_list

		self.degs_up_ensembl_in_pathways = [x for x in self.degs_up_ensembl  if x	 in degs_in_pathways]
		self.degs_up_ensembl_in_pathways.sort()
		self.degs_dw_ensembl_in_pathways = [x for x in self.degs_dw_ensembl  if x	 in degs_in_pathways]
		self.degs_dw_ensembl_in_pathways.sort()

		self.degs_up_ensembl_not_in_pathways = [x for x in self.degs_up_ensembl  if x not in degs_in_pathways]
		self.degs_up_ensembl_not_in_pathways.sort()
		self.degs_dw_ensembl_not_in_pathways = [x for x in self.degs_dw_ensembl  if x not in degs_in_pathways]
		self.degs_dw_ensembl_not_in_pathways.sort()

		self.n_pathways					 = len(df_enr)
		self.n_degs_in_pathways			 = len(degs_in_pathways)
		self.n_degs_ensembl_in_pathways	 = len(self.degs_ensembl_in_pathways )
		self.n_degs_not_in_pathways		 = len(self.degs_not_in_pathways)
		self.n_degs_ensembl_not_in_pathways = len(self.degs_ensembl_not_in_pathways)

		self.n_degs_up_not_ensembl = len(self.degs_up_not_ensembl)
		self.n_degs_dw_not_ensembl = len(self.degs_dw_not_ensembl)
		self.n_degs_not_ensembl = self.n_degs_up_not_ensembl + self.n_degs_dw_not_ensembl

		self.n_degs_up_ensembl_in_pathways	 = len(self.degs_up_ensembl_in_pathways)
		self.n_degs_dw_ensembl_in_pathways	 = len(self.degs_dw_ensembl_in_pathways)
		self.n_degs_up_ensembl_not_in_pathways = len(self.degs_up_ensembl_not_in_pathways)
		self.n_degs_dw_ensembl_not_in_pathways = len(self.degs_dw_ensembl_not_in_pathways)


	def set_pathway_cutoff_params(self, ptw_FDR_cut:float, 
								  ptw_pval_cut:float, ptw_min_num_of_degs_cut:int):

		self.ptw_FDR_cut = ptw_FDR_cut
		self.ptw_pval_cut = ptw_pval_cut
		self.ptw_min_num_of_degs_cut = ptw_min_num_of_degs_cut


	def split_case(self, case):
		if case != self.case_list[self.group_index]:
			try:
				self.group_index = [i for i in range(len(self.case_list)) if case == self.case_list[i]][0]
			except:
				self.group_index = 0
				case = self.case_list[0]

		self.case = case

		if not self.has_age and not self.has_gender:
			group = case
			gender, age = '', ''
		else:
			mat = case.split('_')
			if len(mat) == 1:
				group = case
				gender, age = '', ''

			elif len(mat) == 2:
				'''
					g2b_male
					g3_male_elder
				'''
				if self.has_gender and self.has_age:
					group = mat[0]
					gender = mat[1]
					age = ''
				else:
					if self.has_gender:
						group  = mat[0]
						gender = mat[1]
						age	= ''
					else:
						group  = mat[0]
						age	= mat[1]
						gender = ''
			elif len(mat) == 3:
				group  = mat[0]
				gender = mat[1]
				age	= mat[2]
			else:
				print("Error: Houstou we have problems.")
				print("Problems spliting case?", case)
				raise Exception('Stop: case')

		self.group  = group
		self.gender = gender
		self.age	= age

		return

	def get_best_ptw_cutoff(self, med_max_ptw:str='median', verbose:bool=False):

		aux_geneset_num = self.geneset_num

		'''
			row['quantile'], row.LFC_cut, row.lfc_FDR_cut, \
			row.ptw_pval_cut, row.ptw_FDR_cut, row.ptw_min_num_of_degs_cut, \
			row.n_pathways, row.n_degs_in_pathways, \
			row.n_degs_in_pathways_mean, row.n_degs_in_pathways_median, row.n_degs_in_pathways_std, \
			row.toi1_median, row.toi2_median, row.toi3_median, row.toi4_median
		'''
		self.quantile, self.LFC_cut, self.lfc_FDR_cut, \
		self.ptw_pval_cut, self.ptw_FDR_cut, self.ptw_min_num_of_degs_cut, \
		self.n_pathways_best, self.n_degs_in_pathways_best, \
		self.n_degs_in_pathways_mean, self.n_degs_in_pathways_median, self.n_degs_in_pathways_std, \
		self.toi1_median, self.toi2_median, self.toi3_median, self.toi4_median  = \
		self.cfg.get_cfg_best_ptw_cutoff(self.case, self.normalization, self.geneset_num, med_max_ptw=med_max_ptw, verbose=verbose)

		if self.geneset_num == -1:
			self.geneset_num = aux_geneset_num

		self.set_db(self.geneset_num, verbose=verbose)

	def echo_parameters(self, want_echo_default:bool=False, jump_line:bool=True, echo:bool=False):
		stri = ''
		if want_echo_default:
			stri += self.echo_default(echo=echo)
			if jump_line: stri += '\n'

		stri += self.echo_degs_all(echo=echo)
		if jump_line: stri += '\n'
		stri += self.echo_enriched_pathways(echo=echo)
		return stri


	def echo_default(self, echo:bool=False) -> str:
		stri  = f"geneset lib '{self.geneset_lib}' num={self.geneset_num}\n"
		stri += f"Normalization={self.normalization}; has age={self.has_age} and has gender={self.has_gender}\n"
		if echo: print(stri)
		return stri


	def echo_degs(self, echo:bool=False) -> str:
		stri  = f"For case {self.icase} '{self.case}' ('{self.translate_case(self.case)}'), there are {self.n_degs}/{self.n_degs_ensembl} {self.s_deg_dap}s/{self.s_deg_dap}s with ensembl_id\n"
		stri += f"{self.s_deg_dap}'s cutoffs: abs(LFC)={self.LFC_cut:.3f}; FDR={self.lfc_FDR_cut:.3f}\n"
		if echo: print(stri)
		return stri

	def echo_degs_all(self, echo:bool=False) -> str:
		self.echo_degs()

		stri  = f"\t{self.n_degs}/{self.n_degs_ensembl} {self.s_deg_dap}s/ensembl.\n"
		stri += f"\t\tUp {self.n_degs_up}/{self.n_degs_up_ensembl} {self.s_deg_dap}s/ensembl.\n"
		stri += f"\t\tDw {self.n_degs_dw}/{self.n_degs_dw_ensembl} {self.s_deg_dap}s/ensembl.\n"
		if echo: print(stri)
		return stri

	def echo_enriched_pathways(self, echo:bool=False) -> str:

		stri  = f"Found {self.n_pathways} (best={self.n_pathways_best}) pathways for geneset num={self.geneset_num} '{self.geneset_lib}'\n"
		stri += f"Pathway cutoffs p-value={self.ptw_pval_cut:.3f} fdr={self.ptw_FDR_cut:.3f} min genes={self.ptw_min_num_of_degs_cut}"

		if self.df_enr is not None and not self.df_enr.empty:
			stri += f"{self.s_deg_dap}s found in enriched pathways:\n"
			stri += f"\tThere are {self.n_degs_ensembl} {self.s_deg_dap}s found in pathways\n"
			stri += f"\t{self.n_degs_in_pathways} (best={self.n_degs_in_pathways_best}) {self.s_deg_dap}s in pathways and {self.n_degs_not_in_pathways}/{self.n_degs_ensembl_not_in_pathways} {self.s_deg_dap}s/ensembl not in pathways\n"
			stri += "\n"

			stri += f"\t{self.n_degs_up_ensembl_in_pathways} {self.s_deg_dap}s ensembl Up in pathways\n"
			stri += f"\t{self.n_degs_up_ensembl_not_in_pathways} {self.s_deg_dap}s Up ensembl not in pathways\n"

			stri += "\n"

			stri += f"\t{self.n_degs_dw_ensembl_in_pathways} {self.s_deg_dap}s ensembl Dw in pathways\n"
			stri += f"\t{self.n_degs_dw_ensembl_not_in_pathways} {self.s_deg_dap}s Dw ensembl not in pathways"

		else:
			stri += "No enrichment analysis was calculated."

		if echo: print(stri)
		return stri


	def summary_degs_and_pathways(self, check:bool=False, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		filename = osjoin(self.root_ressum, self.fname_summ_deg_ptw)

		if exists(filename) and not force:
			dfsum = pdreadcsv(self.fname_summ_deg_ptw, self.root_ressum, verbose=verbose)
			col0 = list(dfsum.columns)[0]
			dfsum = dfsum.set_index(col0)
			dfsum.index.names = ['index']
			dfsum = dfsum.infer_objects(copy=False).fillna(0)

			cols5 = list(dfsum.columns)
			ncols = np.arange(0, len(cols5), 2)
			for ncol in ncols:
				dfsum[cols5[ncol]] = dfsum[cols5[ncol]].astype(int)

			return dfsum

		dic = {}; icount = -1
		for case in self.case_list:
			ret, _, _, _ = self.open_case(case, verbose=False)
			if not ret: continue

			''' Best Cutoff Algorithm params '''
			n_degs_bca	= self.n_degs
			n_degs_up_bca = self.n_degs_up
			n_degs_dw_bca = self.n_degs_dw

			n_degs_ensembl_bca	= self.n_degs_ensembl
			n_degs_up_ensembl_bca = self.n_degs_up_ensembl
			n_degs_dw_ensembl_bca = self.n_degs_dw_ensembl

			n_pathways_bca					 = self.n_pathways
			n_degs_in_pathways_bca			 = self.n_degs_in_pathways
			n_degs_not_in_pathways_bca		 = self.n_degs_not_in_pathways
			n_degs_ensembl_in_pathways_bca	 = self.n_degs_ensembl_in_pathways
			n_degs_ensembl_not_in_pathways_bca = self.n_degs_ensembl_not_in_pathways

			n_degs_up_ensembl_in_pathways_bca	 = self.n_degs_up_ensembl_in_pathways
			n_degs_dw_ensembl_in_pathways_bca	 = self.n_degs_dw_ensembl_in_pathways
			n_degs_up_ensembl_not_in_pathways_bca = self.n_degs_up_ensembl_not_in_pathways
			n_degs_dw_ensembl_not_in_pathways_bca = self.n_degs_dw_ensembl_not_in_pathways

			''' Default params '''
			ret, _, _, _ = self.open_case_params(case, LFC_cut=1, lfc_FDR_cut=0.05, ptw_FDR_cut=0.05)
			if not ret: continue

			n_degs_default = self.n_degs
			n_degs_div	 = n_degs_default if n_degs_default > 0 else 1
			n_degs_up_default = self.n_degs_up
			n_degs_dw_default = self.n_degs_dw

			n_degs_ensembl_default	= self.n_degs_ensembl
			n_degs_ensembl_div		= n_degs_ensembl_default if n_degs_ensembl_default > 0 else 1
			n_degs_up_ensembl_default = self.n_degs_up_ensembl
			n_degs_dw_ensembl_default = self.n_degs_dw_ensembl

			n_pathways_default		 = self.n_pathways
			n_degs_in_pathways_default = self.n_degs_in_pathways
			n_degs_ensembl_in_pathways_default	 = self.n_degs_ensembl_in_pathways
			n_degs_ensembl_not_in_pathways_default = self.n_degs_ensembl_not_in_pathways
			n_degs_not_in_pathways_default		 = self.n_degs_not_in_pathways

			n_degs_up_ensembl_in_pathways_default	 = self.n_degs_up_ensembl_in_pathways
			n_degs_dw_ensembl_in_pathways_default	 = self.n_degs_dw_ensembl_in_pathways
			n_degs_up_ensembl_not_in_pathways_default = self.n_degs_up_ensembl_not_in_pathways
			n_degs_dw_ensembl_not_in_pathways_default = self.n_degs_dw_ensembl_not_in_pathways

			'''----------- BCA frequencies & probabilities  --------------------'''
			for perc in [False, True]:
				icount += 1
				dic[icount] = {}
				dic2 = dic[icount]

				if not perc:
					dic2['case'] = case

					'''------------- BCA - best cutoff algorithm -----------------------'''
					dic2['n_degs_bca']	 = n_degs_bca
					if check:
						soma = (n_degs_up_bca + n_degs_dw_bca)
						dic2['n_degs_bca_sum'] = soma
					dic2['n_degs_up_bca']  = n_degs_up_bca
					dic2['n_degs_dw_bca']  = n_degs_dw_bca

					dic2['n_degs_ensembl_bca'] = n_degs_ensembl_bca
					if check:
						soma = (n_degs_up_ensembl_bca + n_degs_dw_ensembl_bca)
						dic2['n_degs_ensembl_bca_sum'] = soma
					dic2['n_degs_up_ensembl_bca'] = n_degs_up_ensembl_bca
					dic2['n_degs_dw_ensembl_bca'] = n_degs_dw_ensembl_bca

					dic2['n_pathways_bca'] = n_pathways_bca

					'''----------- in and out pathway --------------'''
					n_all_degs_bca = n_degs_in_pathways_bca + n_degs_not_in_pathways_bca
					dic2['n_all_degs_bca']		= n_all_degs_bca
					dic2['n_degs_in_pathways_bca']	 = n_degs_in_pathways_bca
					dic2['n_degs_not_in_pathways_bca'] = n_degs_not_in_pathways_bca

					'''----------- Ensembl in and out pathway --------------'''
					n_all_degs_ensembl_bca = n_degs_ensembl_in_pathways_bca + n_degs_ensembl_not_in_pathways_bca
					dic2['n_all_degs_ensembl_bca'] = n_all_degs_ensembl_bca
					if check: dic2['n_degs_ensembl_bca_again'] = n_degs_ensembl_bca
					dic2['n_degs_ensembl_in_pathways_bca']	 = n_degs_ensembl_in_pathways_bca
					if check: dic2['n_degs_ensembl_in_pathways_bca_rep'] = None

					'''----------- Up/Dw in and out pathway --------------'''
					if check:
						soma = n_degs_up_ensembl_in_pathways_bca + n_degs_dw_ensembl_in_pathways_bca
						dic2['n_degs_in_pathways_bca_sum'] = soma
					dic2['n_degs_up_ensembl_in_pathways_bca'] = n_degs_up_ensembl_in_pathways_bca
					dic2['n_degs_dw_ensembl_in_pathways_bca'] = n_degs_dw_ensembl_in_pathways_bca

					'''----------- Up/Dw NOT in and out pathway --------------'''
					dic2['n_degs_ensembl_not_in_pathways_bca'] = n_degs_ensembl_not_in_pathways_bca
					dic2['n_degs_not_in_pathways_bca']	= n_degs_not_in_pathways_bca
					if check:
						dic2['n_degs_not_in_pathways_bca_again'] = n_degs_not_in_pathways_bca
						soma = n_degs_up_ensembl_not_in_pathways_bca + n_degs_dw_ensembl_not_in_pathways_bca
						dic2['n_degs_not_in_pathways_bca_sum'] = soma
					dic2['n_degs_up_ensembl_not_in_pathways_bca'] = n_degs_up_ensembl_not_in_pathways_bca
					dic2['n_degs_dw_ensembl_not_in_pathways_bca'] = n_degs_dw_ensembl_not_in_pathways_bca

					'''---------- default frequencies & probabilities  -----------------'''
					dic2['n_degs_default']	 = n_degs_default
					if check:
						soma = (n_degs_up_default + n_degs_dw_default)
						dic2['n_degs_default_sum'] = soma
					dic2['n_degs_up_default']  = n_degs_up_default
					dic2['n_degs_dw_default']  = n_degs_dw_default

					dic2['n_degs_ensembl_default'] = n_degs_ensembl_default
					if check:
						soma = (n_degs_up_ensembl_default + n_degs_dw_ensembl_default)
						dic2['n_degs_ensembl_default_sum'] = soma
					dic2['n_degs_up_ensembl_default'] = n_degs_up_ensembl_default
					dic2['n_degs_dw_ensembl_default'] = n_degs_dw_ensembl_default

					dic2['n_pathways_default'] = n_pathways_default

					'''----------- in and out pathway --------------'''
					if n_pathways_default == 0:
						dic2['n_all_degs_default'] = None
						dic2['n_degs_in_pathways_default']	 = None
						dic2['n_degs_not_in_pathways_default'] = None

						'''----------- Ensembl in and out pathway --------------'''
						dic2['n_all_degs_ensembl_default'] = None
						if check: dic2['n_degs_ensembl_default_again'] = None
						dic2['n_degs_ensembl_in_pathways_default']	 = None
						if check: dic2['n_degs_ensembl_in_pathways_default_rep'] = None

						'''----------- Up/Dw in and out pathway --------------'''
						if check:
							dic2['n_degs_in_pathways_default_sum'] = None
						dic2['n_degs_up_ensembl_in_pathways_default'] = None
						dic2['n_degs_dw_ensembl_in_pathways_default'] = None

						'''----------- Up/Dw NOT in and out pathway --------------'''
						dic2['n_degs_ensembl_not_in_pathways_default'] = None
						dic2['n_degs_not_in_pathways_default'] = None
						if check:
							dic2['n_degs_not_in_pathways_default_again'] = None
							dic2['n_degs_not_in_pathways_default_sum'] = None
						dic2['n_degs_up_ensembl_not_in_pathways_default'] = None
						dic2['n_degs_dw_ensembl_not_in_pathways_default'] = None

					else:
						n_all_degs_default = n_degs_in_pathways_default + n_degs_not_in_pathways_default
						dic2['n_all_degs_default'] = n_all_degs_default
						dic2['n_degs_in_pathways_default']	 = n_degs_in_pathways_default
						dic2['n_degs_not_in_pathways_default'] = n_degs_not_in_pathways_default

						'''----------- Ensembl in and out pathway --------------'''
						n_all_degs_ensembl_default = n_degs_ensembl_in_pathways_default + n_degs_ensembl_not_in_pathways_default
						dic2['n_all_degs_ensembl_default'] = n_all_degs_ensembl_default
						if check: dic2['n_degs_ensembl_default_again'] = n_degs_ensembl_default
						dic2['n_degs_ensembl_in_pathways_default']	 = n_degs_ensembl_in_pathways_default
						if check: dic2['n_degs_ensembl_in_pathways_default_rep'] = None

						'''----------- Up/Dw in and out pathway --------------'''
						if check:
							soma = n_degs_up_ensembl_in_pathways_default + n_degs_dw_ensembl_in_pathways_default
							dic2['n_degs_in_pathways_default_sum'] = soma
						dic2['n_degs_up_ensembl_in_pathways_default'] = n_degs_up_ensembl_in_pathways_default
						dic2['n_degs_dw_ensembl_in_pathways_default'] = n_degs_dw_ensembl_in_pathways_default

						'''----------- Up/Dw NOT in and out pathway --------------'''
						dic2['n_degs_ensembl_not_in_pathways_default'] = n_degs_ensembl_not_in_pathways_default
						dic2['n_degs_not_in_pathways_default']	= n_degs_not_in_pathways_default
						if check:
							dic2['n_degs_not_in_pathways_default_again'] = n_degs_not_in_pathways_default
							soma = n_degs_up_ensembl_not_in_pathways_default + n_degs_dw_ensembl_not_in_pathways_default
							dic2['n_degs_not_in_pathways_default_sum'] = soma
						dic2['n_degs_up_ensembl_not_in_pathways_default'] = n_degs_up_ensembl_not_in_pathways_default
						dic2['n_degs_dw_ensembl_not_in_pathways_default'] = n_degs_dw_ensembl_not_in_pathways_default

				else:
					dic2['case'] = case + '_perc'

					'''------------- BCA - best cutoff algorithm -----------------------'''
					dic2['n_degs_bca']	 = 1.00
					if check:
						soma = (n_degs_up_bca + n_degs_dw_bca)
						dic2['n_degs_bca_sum'] = soma / n_degs_bca
					dic2['n_degs_up_bca']  = n_degs_up_bca / n_degs_bca
					dic2['n_degs_dw_bca']  = n_degs_dw_bca / n_degs_bca

					dic2['n_degs_ensembl_bca']	= n_degs_ensembl_bca / n_degs_bca
					if check:
						soma = (n_degs_up_ensembl_bca + n_degs_dw_ensembl_bca)
						dic2['n_degs_ensembl_bca_sum'] = soma / n_degs_bca
					dic2['n_degs_up_ensembl_bca'] = n_degs_up_ensembl_bca / n_degs_ensembl_bca
					dic2['n_degs_dw_ensembl_bca'] = n_degs_dw_ensembl_bca / n_degs_ensembl_bca

					dic2['n_pathways_bca'] = None

					'''----------- in and out pathway --------------'''
					dic2['n_all_degs_bca']	= 1.00
					n_all_degs_bca = n_degs_in_pathways_bca + n_degs_not_in_pathways_bca
					dic2['n_degs_in_pathways_bca']	 = n_degs_in_pathways_bca / n_all_degs_bca
					dic2['n_degs_not_in_pathways_bca'] = n_degs_not_in_pathways_bca / n_all_degs_bca

					'''----------- Ensembl in and out pathway --------------'''
					n_all_degs_ensembl_bca = n_degs_ensembl_in_pathways_bca + n_degs_ensembl_not_in_pathways_bca
					if check: dic2['n_degs_ensembl_bca_again'] = 1
					dic2['n_all_degs_ensembl_bca']		= n_all_degs_ensembl_bca / n_degs_ensembl_bca
					dic2['n_degs_ensembl_in_pathways_bca']	 = n_degs_ensembl_in_pathways_bca / n_degs_in_pathways_bca
					if check: dic2['n_degs_ensembl_in_pathways_bca_rep'] = n_degs_ensembl_in_pathways_bca / n_all_degs_ensembl_bca

					'''----------- Up/Dw in and out pathway --------------'''
					if check:
						soma = (n_degs_up_ensembl_in_pathways_bca + n_degs_dw_ensembl_in_pathways_bca)
						dic2['n_degs_in_pathways_bca_sum'] = soma / n_degs_in_pathways_bca
					dic2['n_degs_up_ensembl_in_pathways_bca']  = n_degs_up_ensembl_in_pathways_bca / n_degs_in_pathways_bca
					dic2['n_degs_dw_ensembl_in_pathways_bca']  = n_degs_dw_ensembl_in_pathways_bca / n_degs_in_pathways_bca

					'''----------- Up/Dw NOT in and out pathway --------------'''
					dic2['n_degs_ensembl_not_in_pathways_bca'] = n_degs_ensembl_not_in_pathways_bca / n_all_degs_ensembl_bca
					if check:
						dic2['n_degs_not_in_pathways_bca_again'] = 1
						soma = n_degs_up_ensembl_not_in_pathways_bca + n_degs_dw_ensembl_not_in_pathways_bca
						dic2['n_degs_not_in_pathways_bca_sum'] = soma / n_degs_not_in_pathways_bca
					dic2['n_degs_up_ensembl_not_in_pathways_bca']  = n_degs_up_ensembl_not_in_pathways_bca / n_degs_ensembl_not_in_pathways_bca
					dic2['n_degs_dw_ensembl_not_in_pathways_bca']  = n_degs_dw_ensembl_not_in_pathways_bca / n_degs_ensembl_not_in_pathways_bca

					'''---------- default frequencies & probabilities  -----------------'''

					dic2['n_degs_default'] = 1.00
					if check:
						soma = (n_degs_up_default + n_degs_dw_default)
						dic2['n_degs_default_sum'] = soma / n_degs_div
					dic2['n_degs_up_default']  = n_degs_up_default / n_degs_div
					dic2['n_degs_dw_default']  = n_degs_dw_default / n_degs_div

					dic2['n_degs_ensembl_default']	= n_degs_ensembl_default / n_degs_div
					if check:
						soma = (n_degs_up_ensembl_default + n_degs_dw_ensembl_default)
						dic2['n_degs_ensembl_default_sum'] = soma / n_degs_div
					dic2['n_degs_up_ensembl_default'] = n_degs_up_ensembl_default / n_degs_ensembl_div
					dic2['n_degs_dw_ensembl_default'] = n_degs_dw_ensembl_default / n_degs_ensembl_div


					'''----------- in and out pathway --------------'''
					if n_pathways_default == 0:
						dic2['n_pathways_default'] = None
						dic2['n_all_degs_default'] = None
						dic2['n_degs_in_pathways_default']	 = None
						dic2['n_degs_not_in_pathways_default'] = None

						'''----------- Ensembl in and out pathway --------------'''
						if check: dic2['n_degs_ensembl_default_again'] = None
						dic2['n_all_degs_ensembl_default']		= None
						dic2['n_degs_ensembl_in_pathways_default']	 = None
						if check: dic2['n_degs_ensembl_in_pathways_default_rep'] = None

						'''----------- Up/Dw in and out pathway --------------'''
						if check:
							dic2['n_degs_in_pathways_default_sum'] = None
						dic2['n_degs_up_ensembl_in_pathways_default'] = None
						dic2['n_degs_dw_ensembl_in_pathways_default'] = None

						'''----------- Up/Dw NOT in and out pathway --------------'''
						dic2['n_degs_ensembl_not_in_pathways_default'] = None
						if check:
							dic2['n_degs_not_in_pathways_default_again'] = None
							dic2['n_degs_not_in_pathways_default_sum'] = None
						dic2['n_degs_up_ensembl_not_in_pathways_default'] = None
						dic2['n_degs_dw_ensembl_not_in_pathways_default'] = None
					else:
						dic2['n_pathways_default'] = None
						dic2['n_all_degs_default'] = 1

						n_all_degs_default = n_degs_in_pathways_default + n_degs_not_in_pathways_default
						n_all_genes_in_pathways_div = n_all_degs_default if n_all_degs_default > 0 else 1
						dic2['n_degs_in_pathways_default']	 = n_degs_in_pathways_default / n_all_genes_in_pathways_div
						dic2['n_degs_not_in_pathways_default'] = n_degs_not_in_pathways_default / n_all_genes_in_pathways_div

						'''----------- Ensembl in and out pathway --------------'''
						n_all_degs_ensembl_default = n_degs_ensembl_in_pathways_default + n_degs_ensembl_not_in_pathways_default
						if check: dic2['n_degs_ensembl_default_again'] = 1
						dic2['n_all_degs_ensembl_default']		= n_all_degs_ensembl_default / n_all_genes_in_pathways_div
						dic2['n_degs_ensembl_in_pathways_default']	 = n_degs_ensembl_in_pathways_default / n_all_degs_ensembl_default
						if check: dic2['n_degs_ensembl_in_pathways_default_rep'] = n_degs_ensembl_in_pathways_default / n_degs_in_pathways_default

						'''----------- Up/Dw in and out pathway --------------'''
						soma = (n_degs_up_ensembl_in_pathways_default + n_degs_dw_ensembl_in_pathways_default)
						if check:
							dic2['n_degs_in_pathways_default_sum'] = soma / n_degs_in_pathways_default
						dic2['n_degs_up_ensembl_in_pathways_default']  = n_degs_up_ensembl_in_pathways_default / soma
						dic2['n_degs_dw_ensembl_in_pathways_default']  = n_degs_dw_ensembl_in_pathways_default / soma

						'''----------- Up/Dw NOT in and out pathway --------------'''
						dic2['n_degs_ensembl_not_in_pathways_default'] = n_degs_ensembl_not_in_pathways_default / n_all_degs_ensembl_default
						if check:
							dic2['n_degs_not_in_pathways_default_again'] = 1
							soma = n_degs_up_ensembl_not_in_pathways_default + n_degs_dw_ensembl_not_in_pathways_default
							dic2['n_degs_not_in_pathways_default_sum'] = soma / n_degs_not_in_pathways_default
						dic2['n_degs_up_ensembl_not_in_pathways_default']  = n_degs_up_ensembl_not_in_pathways_default / n_degs_ensembl_not_in_pathways_default
						dic2['n_degs_dw_ensembl_not_in_pathways_default']  = n_degs_dw_ensembl_not_in_pathways_default / n_degs_ensembl_not_in_pathways_default

		dfsum = pd.DataFrame(dic)
		dfsum = dfsum.T
		dfsum = dfsum.set_index('case')
		dfsum = dfsum.T
		pdwritecsv(dfsum, self.fname_summ_deg_ptw, self.root_ressum,  index=True, verbose=verbose)

		return dfsum


	def open_case(self, case:str, save_file:bool=False, force:bool=False, save_EP_xls:bool=False,
				  prompt_verbose:bool=False, verbose:bool=False) -> Tuple[bool, List, List, pd.DataFrame]:
		self.reset_degs_and_df_enr()
		self.case = case
		self.set_icase(case)
		self.split_case(case)
		
		if prompt_verbose: print(">>> case", case)
		self.get_best_ptw_cutoff()

		ret_lfc = self.open_dflfc_ori(verbose=verbose)

		degs, degs_ensembl, dflfc = self.list_of_degs(save_file=save_file, force=force, prompt_verbose=prompt_verbose, verbose=verbose)

		ret_enr = self.open_enrichment_analysis(force=force, save_EP_xls=save_EP_xls, verbose=verbose)

		if verbose:
			if ret_lfc * ret_enr != 1:
				print("Warning: Enrichment Analysis was not performed.")
				print(f"For case {case}, pathways cutoffs are: pval={self.ptw_pval_cut:.3f} and fdr={self.ptw_FDR_cut:.3f}")

		return ret_lfc, degs, degs_ensembl, dflfc


	def open_case_params(self, case:str, LFC_cut:float=1.0, lfc_FDR_cut:float=0.05,
						 ptw_pval_cut:float=0.05, ptw_FDR_cut:float=0.05, ptw_min_num_of_degs_cut:int=3,
						 force:bool=False, prompt_verbose:bool=False, verbose:bool=False) -> Tuple[bool, List, List, pd.DataFrame]:

		self.reset_degs_and_df_enr()
		self.case = case
		self.set_icase(case)
		self.split_case(case)

		if prompt_verbose: print(">>> case", case)

		self.LFC_cut	  = LFC_cut
		self.lfc_FDR_cut  = lfc_FDR_cut
		self.ptw_pval_cut = ptw_pval_cut
		self.ptw_FDR_cut  = ptw_FDR_cut
		self.ptw_min_num_of_degs_cut = ptw_min_num_of_degs_cut

		ret_lfc = self.open_dflfc_ori(verbose=verbose)
		degs, degs_ensembl, dflfc = self.list_of_degs(force=force, prompt_verbose=prompt_verbose, verbose=verbose)

		ret_enr = self.open_enrichment_analysis(force=force, verbose=verbose)

		if ret_lfc * ret_enr != 1:
			if verbose:
				print("Warning: Enrichment Analysis were not performed.")
				print(f"For case {case}, pathways cutoffs are: pval={self.ptw_pval_cut:.3f} and fdr={self.ptw_FDR_cut:.3f}")

		return ret_lfc, degs, degs_ensembl, dflfc


	def open_case_params_wo_pathways(self, case:str, LFC_cut:float, lfc_FDR_cut:float,
									 prompt_verbose:bool=False, verbose:bool=False) -> Tuple[bool, List, List, pd.DataFrame]:
		self.reset_degs_and_df_enr()
		self.case = case
		self.set_icase(case)
		self.split_case(case)

		self.LFC_cut = LFC_cut
		self.lfc_FDR_cut = lfc_FDR_cut

		ret = self.open_dflfc_ori(verbose=verbose)
		degs, degs_ensembl, dflfc = self.list_of_degs(force=False, prompt_verbose=prompt_verbose, verbose=verbose)

		return ret, degs, degs_ensembl, dflfc

	def open_case_simple(self, case:str, verbose:bool=False) -> bool:
		'''
			open_case_simple:
				set case
				open dflfc_ori
		'''
		self.reset_degs_and_df_enr()
		self.case = case
		self.set_icase(case)

		self.split_case(case)
		ret = self.open_dflfc_ori(verbose=verbose)

		return ret

	def set_icase(self, case:str):
		mat = [i for i in range(len(self.case_list)) if self.case_list[i] == case]

		if len(mat) != 1:
			print(f"Error: {self.case} not in {self.case_list}")
			raise Exception(f"\n\n---------- Stop: mat {mat} ------------")
		
		self.icase = mat[0]

	def cols_omics(self) -> List:

		if self.s_omics == 'microarray':
			cols = ['probe', 'symbol', 'symbol_prev', 'symb_or_syn', 'biotype', '_type',
					'lfc', 'abs_lfc', 'fdr',  'description',
					'desc_gff', 'description_prev', 'accession', 'ensembl_id', ]
		elif self.s_omics == 'proteomics':
			cols = ['symbol', 'biotype', 'lfc', 'abs_lfc', 'fdr',  'description', 'ensembl_id', ]
		else:
			cols = ['symbol', 'biotype', 'lfc', 'abs_lfc', 'fdr',  'description', 'ensembl_id', ]

		return cols


	def calc_all_genes_in_pubmed_per_case(self, force:bool=False, prompt_verbose:bool=False,
										  verbose:bool=False) -> pd.DataFrame:

		filename = osjoin(self.root_ressum, self.fname_degs_summary)

		if exists(filename) and not force:
			return pdreadcsv(self.fname_degs_summary, self.root_ressum, verbose=verbose)

		dic={}; icount=-1

		for case in self.case_list:
			if prompt_verbose: print(f">>> case {case}")

			ret, _, _, _ = self.open_case(case, verbose)
			if not ret or self.degs is None:
				if verbose: print(f"There are no {self.s_deg_dap} for {case}")
				continue

			if self.dfsimb_perc is None or self.dfsimb_perc.empty:
				degs_in_pubmed	 = []
				degs_not_in_pubmed = []
			else:
				degs_in_pubmed	   = [x for x in self.dfsimb_perc.symbol if x	 in self.degs]
				degs_not_in_pubmed = [x for x in self.dfsimb_perc.symbol if x not in degs_in_pubmed]

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			dic2['case'] = case
			dic2['normalization']  = self.normalization
			dic2['geneset_num']	= self.geneset_num
			dic2['quantile']	   = self.quantile
			dic2['LFC_cut'] = self.LFC_cut
			dic2['lfc_FDR_cut'] = self.lfc_FDR_cut

			dic2['n_degs_in_pubmed']	 = len(degs_in_pubmed)
			dic2['n_degs_not_in_pubmed'] = len(degs_not_in_pubmed)

			dic2['n_degs']	= self.n_degs
			dic2['n_degs_up'] = self.n_degs_up
			dic2['n_degs_dw'] = self.n_degs_dw
			dic2['n_degs_up_ensembl'] = self.n_degs_up_ensembl
			dic2['n_degs_dw_ensembl'] = self.n_degs_dw_ensembl

			dic2['n_pathways']				 = self.n_pathways
			dic2['n_all_genes_in_pathways']	= self.n_degs_in_pathways + self.n_degs_not_in_pathways
			dic2['n_degs_in_pathways']		 = self.n_degs_in_pathways
			dic2['n_degs_in_pathways_ensembl'] = self.n_degs_ensembl_in_pathways
			dic2['n_degs_not_in_pathways']	 = self.n_degs_not_in_pathways

			dic2['n_degs_up_in_pathways']	 = self.n_degs_up_ensembl_in_pathways
			dic2['n_degs_dw_in_pathways']	 = self.n_degs_dw_ensembl_in_pathways
			dic2['n_degs_up_not_in_pathways'] = self.n_degs_up_ensembl_not_in_pathways
			dic2['n_degs_dw_not_in_pathways'] = self.n_degs_dw_ensembl_not_in_pathways

			dic2['degs_in_pubmed']	 = degs_in_pubmed
			dic2['degs_not_in_pubmed'] = degs_not_in_pubmed

			dic2['degs']	= self.degs
			dic2['degs_up'] = self.degs_up
			dic2['degs_dw'] = self.degs_dw
			dic2['degs_up_ensembl'] = self.degs_up_ensembl
			dic2['degs_dw_ensembl'] = self.degs_dw_ensembl

			dic2['degs_in_pathways']			= self.degs_in_pathways
			dic2['degs_ensembl_in_pathways']	= self.degs_ensembl_in_pathways
			dic2['degs_up_ensembl_in_pathways'] = self.degs_up_ensembl_in_pathways
			dic2['degs_dw_ensembl_in_pathways'] = self.degs_dw_ensembl_in_pathways

			dic2['degs_not_in_pathways']		 = self.degs_not_in_pathways
			dic2['degs_ensembl_not_in_pathways'] = self.degs_ensembl_not_in_pathways
			dic2['degs_up_ensembl_not_in_pathways']	  = self.degs_up_ensembl_not_in_pathways
			dic2['degs_dw_ensembl_not_in_pathways']	  = self.degs_dw_ensembl_not_in_pathways


			all_text=''
			df_enr = self.df_enr

			if df_enr is None or df_enr.empty:
				n_enr = 0
			else:
				n_enr = len(df_enr)
				for i_enr in range(len(df_enr)):
					row = df_enr.iloc[i_enr]
					text = f"{i_enr+1}) {row.pathway} ({row.pathway_id}) fdr: {row.fdr:.2e} {self.s_deg_dap}s ({row.num_of_genes}) {row.genes}\n"
					all_text += text

					all_text = all_text[:-1]


			dic2['pathway_list'] = all_text
			dic2['pathway_id_list']  = self.pathway_id_list
			dic2['pathway_fdr_list'] = self.pathway_fdr_list
			dic2['num_of_genes_in_pathway_list'] = self.num_of_genes_in_pathway_list
			dic2['n_enr'] = n_enr


		dfa = pd.DataFrame(dic).T
		# dfa = dfa.set_index('case')

		ret = pdwritecsv(dfa, self.fname_degs_summary, self.root_ressum, verbose=verbose)
		return dfa

	def calc_only_genes_in_pathway_per_case_DEPRECATED(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		if self.df_summ is None or self.df_summ.empty:
			self.open_enriched_pathways_summary()

		fname = self.fname_pathway_summary%(self.normalization)
		filename = osjoin(self.root_ressum, fname)

		if exists(filename) and not force:
			dfa = pdreadcsv(fname, self.root_ressum, verbose=verbose)
			return dfa

		dic = {}
		for i in range(len(self.df_summ)):
			case  = self.df_summ.iloc[i].case
			genes = self.df_summ.iloc[i].all_degs
			genes = eval(genes)

			if not isinstance(genes, list):
				self.genes = genes
				print("genes must be a list")
				raise Exception('stop: genes as list')

			if case not in dic:
				dic[case] = []

			dic[case] += genes

		dicq = {}
		i = -1
		for case in dic.keys():
			genes = np.unique(dic[case])
			i+= 1
			dicq[i] = {}
			dic2 = dicq[i]

			dic2['case'] = case
			dic2['genes'] = genes
			dic2['n'] = len(genes)

		dfa = pd.DataFrame(dicq).T
		dfa = dfa.set_index('case')

		dfa = dfa.loc[self.list_order]
		dfa = dfa.reset_index()

		ret = pdwritecsv(dfa, fname, self.root_ressum, verbose=verbose)

		return dfa

	def save_up_and_down_degs(self, verbose:bool=False):

		''' list_of_degs calculates:
				self.degs, self.degs_ensembl, self.dflfc, self.n_degs,
				self.dfdegs_up, self.dfdegs_dw
				self.degs_up, self.degs_dw
		'''
		'''------------- degs_up ------------------------------------------------'''
		fname = self.fname_txt_updegs%(case, self.LFC_cut, self.lfc_FDR_cut, self.normalization)

		text = "\n".join(self.degs_up)
		write_txt(text, fname, self.root_ressum, verbose=verbose)

		fname = fname.replace('.txt', '_geneid.txt')
		degs_geneid_up = [self.gene.find_mygene_geneid(x) for x in self.degs_up]
		degs_geneid_up = [str(x) for x in degs_geneid_up if x is not None]
		text = "\n".join(degs_geneid_up)
		write_txt(text, fname, self.root_ressum, verbose=verbose)

		if len(self.dfdegs_up) > 150:
			degs_up150 = [x for x in self.dfdegs_up.iloc[:150].symbol if isinstance(x, str)]
			text = "\n".join(degs_up150)
			fname2 = fname.replace('.txt', '_150limit.txt')
			write_txt(text, fname2, self.root_ressum, verbose=verbose)

			degs_geneid_up150 = [self.gene.find_mygene_geneid(x) for x in degs_up150]
			text = "\n".join(degs_geneid_up150)
			fname2 = fname.replace('.txt', '_150limit.txt')
			write_txt(text, fname2, self.root_ressum, verbose=verbose)

		'''------------- degs_up ------- end ----------------------------------'''

		'''------------- degs down ---------------------------------------------'''
		fname = self.fname_txt_dwdegs%(case, self.LFC_cut, self.lfc_FDR_cut, self.normalization)

		text = "\n".join(self.degs_dw)
		write_txt(text, fname, self.root_ressum, verbose=verbose)

		fname = fname.replace('.txt', '_geneid.txt')
		degs_geneid_dw = [self.gene.find_mygene_geneid(x) for x in degs_dw]
		degs_geneid_dw = [str(x) for x in degs_geneid_dw if x is not None]
		text = "\n".join(degs_geneid_dw)
		write_txt(text, fname, self.root_ressum, verbose=verbose)

		if len(self.dfdegs_dw) > 150:
			degs_dw150 = [x for x in self.dfdegs_dw.iloc[:150].symbol if isinstance(x, str)]
			text = "\n".join(degs_dw150)
			fname2 = fname.replace('.txt', '_150limit.txt')
			write_txt(text, fname2, self.root_ressum, verbose=verbose)

			degs_geneid_dw150 = [self.gene.find_mygene_geneid(x) for x in degs_dw150]
			text = "\n".join(degs_geneid_dw150)
			fname2 = fname.replace('.txt', '_150limit.txt')
			write_txt(text, fname2, self.root_ressum, verbose=verbose)
		'''------------- degs down ------- end --------------------------------'''

		return

	def test_DEGs(self, prompt_verbose:bool=True, verbose:bool=False):
		all_degs, all_degs_ensembl, degs_in_pathways = [], [], []

		for case in self.case_list:
			if prompt_verbose: print(f"<<< case: {case}")
			ret, _, _, _ = self.open_case(case, verbose)

			if not ret:
				continue

			all_degs += self.degs
			all_degs_ensembl += self.degs_ensembl
			degs_in_pathways += self.degs_in_pathways

		all_degs = list(np.unique(all_degs))
		all_degs_ensembl = list(np.unique(all_degs_ensembl))
		degs_in_pathways = list(np.unique(degs_in_pathways))

		if prompt_verbose:
			print(f">> {self.s_deg_dap}s calc {len(all_degs)}: '{', '.join(all_degs)}'")
			print(f">> {self.s_deg_dap}s with ensembl_id calc {len(all_degs_ensembl)}: '{', '.join(all_degs_ensembl)}'")
			print(f">> {self.s_deg_dap}s in pathways {len(degs_in_pathways)}: '{', '.join(degs_in_pathways)}'")


	def calc_degs_and_pathways_summary(self, force:bool=False, prompt_verbose:bool=False,
									   verbose:bool=False) -> bool:

		filename = osjoin(self.root_ressum, self.fname_degs_and_pathways_summary)

		if exists(filename) and not force:
			return self.open_enriched_pathways_summary(verbose=verbose)

		dic_summ = {}; icount = -1

		for case in self.case_list:
			if prompt_verbose: print(f">>> case: {case}")

			ret, _, _, _ = self.open_case(case, verbose=verbose)

			if not ret or self.degs is None or self.degs == []:
				continue

			icount += 1
			dic_summ[icount] = {}
			dic2 = dic_summ[icount]

			dic2['case'] = case
			dic2['normalization'] = self.normalization
			dic2['geneset_num'] = self.geneset_num
			dic2['quantile'] = self.quantile

			dic2['LFC_cut']	  = self.LFC_cut
			dic2['lfc_FDR_cut']	  = self.lfc_FDR_cut
			dic2['ptw_pval_cut'] = self.ptw_pval_cut
			dic2['ptw_FDR_cut']  = self.ptw_FDR_cut
			dic2['ptw_min_num_of_degs_cut'] = self.ptw_min_num_of_degs_cut

			dic2['med_max_ptw'] = self.med_max_ptw
			dic2['toi1_median'] = self.toi1_median
			dic2['toi2_median'] = self.toi2_median
			dic2['toi3_median'] = self.toi3_median
			dic2['toi4_median'] = self.toi4_median

			'''----------- number of DEGs -------------------'''
			dic2['n_degs'] = self.n_degs
			dic2['n_degs_up'] = self.n_degs_up
			dic2['n_degs_dw'] = self.n_degs_dw
			dic2['n_degs_up_ensembl'] = self.n_degs_up_ensembl
			dic2['n_degs_dw_ensembl'] = self.n_degs_dw_ensembl

			dic2['n_pathways']			 = self.n_pathways
			dic2['n_degs_in_pathways']	 = self.n_degs_in_pathways
			dic2['n_degs_not_in_pathways'] = self.n_degs_not_in_pathways

			dic2['n_degs_up_ensembl_in_pathways'] = self.n_degs_up_ensembl_in_pathways
			dic2['n_degs_dw_ensembl_in_pathways'] = self.n_degs_dw_ensembl_in_pathways

			dic2['n_degs_up_ensembl_not_in_pathways'] = self.n_degs_up_ensembl_not_in_pathways
			dic2['n_degs_dw_ensembl_not_in_pathways'] = self.n_degs_dw_ensembl_not_in_pathways

			'''------------- Pathways ------------------'''
			dic2['pathway_list']	 = self.pathway_list
			dic2['pathway_id_list']  = self.pathway_id_list
			dic2['pathway_fdr_list'] = self.pathway_fdr_list
			dic2['num_of_genes_in_pathway_list'] = self.num_of_genes_in_pathway_list


			'''--------------- DEGS ------------------'''
			dic2['degs'] = self.degs
			dic2['degs_up'] = self.degs_up
			dic2['degs_dw'] = self.degs_dw
			dic2['degs_up_ensembl'] = self.degs_up_ensembl
			dic2['degs_dw_ensembl'] = self.degs_dw_ensembl

			dic2['degs_in_pathways'] = self.degs_in_pathways
			dic2['degs_ensembl_in_pathways'] = self.degs_ensembl_in_pathways
			dic2['degs_not_in_pathways'] = self.degs_not_in_pathways
			dic2['degs_ensembl_not_in_pathways'] = self.degs_ensembl_not_in_pathways

			dic2['degs_up_ensembl_in_pathways'] = self.degs_up_ensembl_in_pathways
			dic2['degs_up_ensembl_not_in_pathways'] = self.degs_up_ensembl_not_in_pathways
			dic2['degs_dw_ensembl_in_pathways'] = self.degs_dw_ensembl_in_pathways
			dic2['degs_dw_ensembl_not_in_pathways'] = self.degs_dw_ensembl_not_in_pathways

		self.df_summ = pd.DataFrame(dic_summ).T

		ret2 = pdwritecsv(self.df_summ, self.fname_degs_and_pathways_summary, self.root_ressum, verbose=verbose)

		return ret2

	def open_enriched_pathways_summary(self, verbose:bool=False) -> bool:
		self.df_summ = pdreadcsv(self.fname_degs_and_pathways_summary, self.root_ressum, verbose=verbose)
		return self.df_summ is not None and not self.df_summ.empty


	def list_of_degs_set_params(self, LFC_cut:float, lfc_FDR_cut:float, force:bool=False,
								save_file:bool=False, prompt_verbose:bool=False, 
								verbose:bool=False) -> Tuple[List, List, pd.DataFrame]:

		self.reset_degs_and_df_enr()

		self.LFC_cut = LFC_cut
		self.lfc_FDR_cut = lfc_FDR_cut

		# degs, degs_ensembl, dflfc
		return self.list_of_degs(force=force, save_file=save_file, prompt_verbose=prompt_verbose, verbose=verbose)


	def reset_degs_and_df_enr(self):
		self.dflfc, self.dfdegs_up, self.dfdegs_dw = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
		self.dflfc_ensembl, self.dfdegs_up_ensembl, self.dfdegs_dw_ensembl = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

		self.degs, self.n_degs  = [], 0
		self.degs_up, self.n_degs_up  = [], 0
		self.degs_dw, self.n_degs_dw  = [], 0

		self.degs_ensembl, self.n_degs_ensembl  = [], 0
		self.degs_up_ensembl, self.n_degs_up_ensembl  = [], 0
		self.degs_dw_ensembl, self.n_degs_dw_ensembl  = [], 0

		self.biotype_list = []
		''' group by biotype '''
		self.dfg_dflc_ori , self.dfg_dflc = None, None

		self.valid_genes, self.annotated_genes = [], []
		self.n_total_genes, self.n_annotated_genes = -1, -1

		self.n_pathways_best, self.n_degs_in_pathways_best = 0, 0
		self.all_genes_annotatted_in_pathway = []

		self.df_enr, self.pathway_in, self.df_enr_reactome = pd.DataFrame(), False, pd.DataFrame()
		
		self.clear_degs_vars()


	def clear_degs_vars(self):
		self.degs_in_pathways		  = []
		self.degs_ensembl_in_pathways = []
		self.degs_not_in_pathways	 = []
		self.degs_up_not_in_pathways  = []
		self.degs_dw_not_in_pathways  = []
		self.degs_ensembl_not_in_pathways = []

		self.pathway_list	  = []
		self.pathway_id_list  = []
		self.pathway_fdr_list = []
		self.num_of_genes_in_pathway_list = []

		self.degs_up_ensembl = []
		self.degs_dw_ensembl = []

		self.degs_up_ensembl_in_pathways = []
		self.degs_dw_ensembl_in_pathways = []
		self.degs_up_ensembl_not_in_pathways = []
		self.degs_dw_ensembl_not_in_pathways = []
		self.degs_up_not_ensembl = []
		self.degs_dw_not_ensembl = []

		self.n_pathways			= 0
		self.n_degs_in_pathways	= 0
		self.n_degs_ensembl_in_pathways	= 0
		self.n_degs_not_in_pathways		= 0
		self.n_degs_ensembl_not_in_pathways = 0

		self.n_degs_not_ensembl = 0
		self.n_degs_up_not_ensembl = 0
		self.n_degs_dw_not_ensembl = 0

		self.n_degs_up_ensembl_in_pathways = 0
		self.n_degs_dw_ensembl_in_pathways = 0
		self.n_degs_up_ensembl_not_in_pathways = 0
		self.n_degs_dw_ensembl_not_in_pathways = 0		


	def list_of_degs(self, force:bool=False, save_file:bool=False,
					 prompt_verbose:bool=False, verbose:bool=False) -> Tuple[List, List, pd.DataFrame]:
		''' calculates:
				self.degs, self.dflfc, self.n_degs,
				self.dfdegs_up, self.dfdegs_dw
				self.degs_up, self.degs_dw

				self.degs_ensembl, self.n_degs_ensembl
				self.degs_up_ensembl, self.n_degs_up_ensembl
				self.degs_dw_ensembl, self.n_degs_dw_ensembl
		'''
		stri = '_DAP_' if self.gene_protein == 'protein' else '_DEG_'

		if self.dflfc_ori is None or self.dflfc_ori.empty:
			print(f"No dflfc table was calculated for this case {self.case}")
			return [], [], pd.DataFrame()

		dflfc = self.dflfc_ori.copy()
		dflfc = dflfc[ (dflfc.abs_lfc >= self.LFC_cut) & (dflfc.fdr < self.lfc_FDR_cut)].copy()

		if dflfc.empty:
			self.dflfc = dflfc
			if verbose: print(f"There are no {self.s_deg_dap}s for dflfc fdr < {self.lfc_FDR_cut:.3f} and lfc >= {self.LFC_cut:.3f}")
			return [], [], pd.DataFrame()

		dflfc = dflfc.sort_values('fdr', ascending=True)
		dflfc.reset_index(inplace=True, drop=True)

		if save_file:
			fname = f'bca_{self.s_deg_dap}s_for_{self.disease}_case_{self.case}_best_cutoff_abs_lfc_{self.LFC_cut:.3f}_fdr_{self.lfc_FDR_cut:.3f}.tsv'
			ret = pdwritecsv(dflfc, fname, self.root_result, verbose=verbose)

			fname = f'bca_{self.s_deg_dap}s_for_{self.disease}_case_{self.case}_best_cutoff_abs_lfc_{self.LFC_cut:.3f}_fdr_{self.lfc_FDR_cut:.3f}.txt'
			text = "\n".join(dflfc.symbol)
			write_txt(text, fname, self.root_result, verbose=verbose)

		# if verbose: print(f"Found %d {self.s_deg_dap}"%(len(dflfc)))
		''' official HUGO names '''
		degs = list(dflfc.symbol)
		self.dflfc = dflfc

		self.degs = degs
		self.n_degs = len(degs)

		'''
			TEC - transcritps to be validate by experiment
			link: https://www.ensembl.org/info/genome/genebuild/biotypes.html
		'''
		dfdegs_ensembl = dflfc[ ~pd.isnull(dflfc.ensembl_id) ].copy()
		dfdegs_ensembl.reset_index(inplace=True, drop=True)
		degs_ensembl = list(dfdegs_ensembl.symbol)
		self.dfdegs_ensembl = dfdegs_ensembl
		self.degs_ensembl = degs_ensembl
		self.n_degs_ensembl = len(degs_ensembl)

		self.biotype_list = np.unique(dfdegs_ensembl.biotype)

		self.dfg_dflc_ori = self.get_dflfc_ori_biotypes()
		self.dfg_dflc	 = self.get_dflfc_biotypes()

		''' ---------- Up ------------------------'''
		self.dfdegs_up = dflfc[dflfc.lfc > 0].copy()
		if self.dfdegs_up.empty:
			self.degs_up = []
			self.n_degs_up = 0
		else:
			self.dfdegs_up.reset_index(inplace=True, drop=True)
			self.degs_up = list(self.dfdegs_up.symbol)
			self.degs_up.sort()
			self.n_degs_up = len(self.degs_up)

		self.dfdegs_up_ensembl = dfdegs_ensembl[dfdegs_ensembl.lfc > 0].copy()
		if self.dfdegs_up_ensembl.empty:
			self.degs_up_ensembl = []
			self.n_degs_up_ensembl = 0
		else:
			self.dfdegs_up_ensembl.reset_index(inplace=True, drop=True)
			self.degs_up_ensembl = list(self.dfdegs_up_ensembl.symbol)
			self.degs_up_ensembl.sort()
			self.n_degs_up_ensembl = len(self.degs_up_ensembl)

		''' ---------- Down ------------------------'''
		self.dfdegs_dw = dflfc[dflfc.lfc < 0].copy()
		if self.dfdegs_dw.empty:
			self.degs_dw = []
			self.n_degs_dw = 0
		else:
			self.dfdegs_dw.reset_index(inplace=True, drop=True)
			self.degs_dw = list(self.dfdegs_dw.symbol)
			self.degs_dw.sort()
			self.n_degs_dw = len(self.degs_dw)

		self.dfdegs_dw_ensembl = dfdegs_ensembl[dfdegs_ensembl.lfc < 0].copy()
		if self.dfdegs_dw_ensembl.empty:
			self.degs_dw_ensembl = []
			self.n_degs_dw_ensembl = 0
		else:
			self.dfdegs_dw_ensembl.reset_index(inplace=True, drop=True)
			self.degs_dw_ensembl = list(self.dfdegs_dw_ensembl.symbol)
			self.degs_dw_ensembl.sort()
			self.n_degs_dw_ensembl = len(self.degs_dw_ensembl)

		fname_final_lfc_ori, fname_given_lfc, title_dummy = self.set_lfc_names()
		fname0 = fname_final_lfc_ori

		biotype_list = list(self.dfg_dflc.biotype)
		biotype_list.sort()

		self.dfg_up_biotype, self.dfg_dw_biotype = None, None

		if self.n_degs_up > 0:
			for biotype in biotype_list:
				dfa = self.dfdegs_up[self.dfdegs_up.biotype == biotype]

				if save_file and not dfa.empty:
					fname_up = fname0.replace('_final_', f"{stri}UP_{biotype}_").replace('.tsv', '.txt')
					filename = osjoin(self.root_result, fname_up)

					if not exists(filename) or force:
						text = "\n".join(dfa.symbol)
						write_txt(text, fname_up, self.root_result, verbose=verbose)

		if self.n_degs_dw > 0:
			for biotype in biotype_list:
				dfa = self.dfdegs_dw[self.dfdegs_dw.biotype == biotype]

				if save_file and not dfa.empty:
					fname_dw = fname0.replace('_final_', f"{stri}DW_{biotype}_").replace('.tsv', '.txt')
					filename = osjoin(self.root_result, fname_dw)

					if not exists(filename) or force:
						text = "\n".join(dfa.symbol)
						write_txt(text, fname_dw, self.root_result, verbose=verbose)


		if prompt_verbose:
			print(f"\t{self.s_deg_dap}s {len(self.degs)}")
			print(f"\t\tUp (#{len(self.degs_up)})")
			print(f"\t\tDw (#{len(self.degs_dw)})")

			print("\nUp-regulated per biotype")
			if self.n_degs_up > 0:
				dfg = self.get_df_biotypes(self.dfdegs_up)
				self.dfg_up_biotype = dfg
				print(dfg)

			print("\nDown-regulated per biotype")
			if self.n_degs_dw > 0:
				dfg = self.get_df_biotypes(self.dfdegs_dw)
				self.dfg_dw_biotype = dfg
				print(dfg)


		return degs, degs_ensembl, dflfc

	def set_dfsim(self):
		if self.dfsim is None or self.dfsim.empty:
			self.fdr_list = []
			self.lfc_list = []
		else:
			self.fdr_list = np.unique(self.dfsim.lfc_FDR_cut)
			self.lfc_list = np.unique(self.dfsim.LFC_cut)
	
	def open_simulation_table(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		if self.dfsim is not None and not self.dfsim.empty and not force:
			self.set_dfsim()
			return self.dfsim

		filename = osjoin(self.root_config, self.cfg.fname_lfc_cutoff)

		if exists(filename):
			dfsim = pdreadcsv(self.cfg.fname_lfc_cutoff, self.root_config, verbose=verbose)
		else:
			print(f"Could not found simulation table: {filename}")
			dfsim = pd.DataFrame()

		self.dfsim = dfsim
		self.set_dfsim()

		return dfsim


	'''
	cutoff_list = [(1, 0.05), ...]
	'''
	def calc_degs_cutoff_simulation(self, cutoff_list:List, force:bool=False, save_file:bool=False, 
									n_echo:int=-1, verbose:bool=False) -> pd.DataFrame:

		filename = self.root_config / self.cfg.fname_lfc_cutoff

		if filename.exists() and not force:
			dfsim = pdreadcsv(self.cfg.fname_lfc_cutoff, self.root_config)
			self.dfsim = dfsim
			return dfsim

		if len(cutoff_list) == 0:
			cutoff_list = [(x, y) for x in self.lfc_list for y in self.fdr_list]

		icount=-1
		dic= {}
		for case in self.case_list:

			print(">>>", case)
			if not self.open_case_simple(case):
				continue

			for LFC_cut, lfc_FDR_cut in cutoff_list:
				_, _, _ = self.list_of_degs_set_params(LFC_cut, lfc_FDR_cut, force=force, save_file=save_file, verbose=verbose)

				icount += 1
				dic[icount] = {}
				dic2 = dic[icount]

				dic2['case'] = case
				dic2['normalization'] = self.normalization

				dic2['cutoff'] = f"{LFC_cut:.3f} - {lfc_FDR_cut:.3f}"

				dic2['LFC_cut'] = LFC_cut
				dic2['lfc_FDR_cut'] = lfc_FDR_cut

				dic2['degs'] = self.degs
				dic2['n_degs'] = self.n_degs
				dic2['degs_ensembl'] = self.degs_ensembl
				dic2['n_degs_ensembl'] = self.n_degs_ensembl

				dic2['degs_up'] = self.degs_up
				dic2['n_degs_up'] = self.n_degs_up
				dic2['degs_up_ensembl'] = self.degs_up_ensembl
				dic2['n_degs_up_ensembl'] = self.n_degs_up_ensembl

				dic2['degs_dw'] = self.degs_dw
				dic2['n_degs_dw'] =self.n_degs_dw
				dic2['degs_dw_ensembl'] = self.degs_dw_ensembl
				dic2['n_degs_dw_ensembl'] = self.n_degs_dw_ensembl

				if n_echo > 9 and icount%n_echo==0:

					print(f"{self.s_deg_dap} cutoff: lfc={self.LFC_cut:.3f}; lfc_fdr={self.lfc_FDR_cut:.3f}")
					print(f"\tthere are {len(self.degs)} {self.s_deg_dap}s")
					print(f"\t\t{self.n_degs_up} Up: ({', '.join(self.degs_up)})")
					print(f"\t\t{self.n_degs_dw} Dw: ({', '.join(self.degs_dw)})")
					print("")

		dfsim = pd.DataFrame(dic).T
		self.dfsim = dfsim
		self.cfg.save_best_lfc_cutoff(dfsim, verbose=True)

		return dfsim


	def prepare_df_two_conditions(self, case, dfn=None):

		self.split_case(case)
		ret_lfc = self.open_dflfc_ori()

		df_lfc = self.dflfc_ori

		cols = ['symbol', 'lfc', 'pval', 'fdr']
		if dfn is None:
			dfn = pd.merge(self.dfsun, df_lfc[cols], how='outer', on='symbol')
		else:
			dfn = pd.merge(dfn, df_lfc[cols], how='outer', on='symbol')

		dfn = dfn[~pd.isnull(dfn.biopax_num)]
		dfn.biopax_num = dfn.biopax_num.astype(int)
		dfn.num_parent = dfn.num_parent.astype(int)

		cols = list(dfn.columns)
		dfn.lfc = [0 if pd.isnull(x) else x for x in dfn.lfc]
		cols[-1] = 'fdr_%s'%(case)
		cols[-2] = 'pval_%s'%(case)
		cols[-3] = 'lfc_%s'%(case)
		dfn.columns = cols

		return dfn


	def compare_two_conditions(self, verbose=False):
		if verbose: print("Comparing 2 conditions ...")
		case1, case2 = self.group1, self.group2

		df_comp = self.prepare_df_two_conditions(case1, None)
		df_comp = self.prepare_df_two_conditions(case2, df_comp)
		df_comp.reset_index(inplace=True, drop=True)

		exc_list=[]; diff_list=[]
		for i in range(len(df_comp)):
			lfc1 = df_comp.iloc[i]['lfc_%s'%(case1)]
			lfc2 = df_comp.iloc[i]['lfc_%s'%(case2)]
			diff = lfc2 - lfc1

			exc = True if lfc1*lfc2 < 0 else False

			exc_list.append(exc)
			diff_list.append(diff)

		df_comp['exchange_sign'] = exc_list
		df_comp['lfc_diff'] = diff_list
		df_comp['abs_lfc_diff'] = np.abs(diff_list)

		df_comp = df_comp.sort_values('abs_lfc_diff', ascending=False)
		df_comp.reset_index(inplace=True, drop=True)

		self.df_comp = df_comp
		return df_comp


	def calc_up_donw_reg(self, val2, val1, cutoff=0.4):

		diff = val2 - val1

		if val2 > cutoff:
			stri2 = 'up'
		elif val2 < -cutoff:
			stri2 = 'dw'
		else:
			stri2 = '--'

		if val1 > cutoff:
			stri1 = 'up'
		elif val1 < -cutoff:
			stri1 = 'dw'
		else:
			stri1 = '--'

		upup_dwdw = stri2 + '/' + stri1

		if stri2 == 'up' and stri1 == 'dw' or stri1 == 'up' and stri2 == 'dw':
			signal = ' and inverted'
		else:
			signal = '			 '

		if diff >= cutoff:
			modulated = 'upmodulated	'
		elif diff <= -cutoff:
			modulated = 'downmodulated  '
		else:
			modulated = '~ not modulated'

		return upup_dwdw, modulated+signal

	def treat_three_values(self, val):
		if val == 0:
			stri = ' %.2f'%(val)
		elif val > 0:
			stri = '+%.2f'%(val)
		else:
			stri = '%.2f'%(val)

		return stri

	def calc_venn_df_summ_between_two_cases(self, case1:str, case2:str, col:str, first150:bool=False,
											use_ensembl_degs:bool=True, in_pathway:bool=False,
											print_inverted:bool=True, figsize:tuple=(8,6),
											title_font_size=15, dpi:int=300, format:str='png', 
											facecolor:str='white', save_fig:bool=True, 
											verbose:bool=False) -> Tuple[object, str, str, str, str, str, List, List, List, List]:


		if self.df_summ is None or self.df_summ.empty:
			print("Impossible to continue: self.df_summ is None or empty")
			return None, '', '', '', '', '', [], [], [], ''

		if in_pathway:
			if col.startswith('degs_up'):
				col = 'degs_up_ensembl_in_pathways'
			else:
				col = 'degs_dw_ensembl_in_pathways'
		else:
			if use_ensembl_degs:
				if col.startswith('degs_up'):
					col = 'degs_up_ensembl'
				else:
					col = 'degs_dw_ensembl'

		# print(f">>>> {col} in_pathway {in_pathway} use_ensembl_degs {use_ensembl_degs}")
		print(f">>>> Selected column: {col}")

		if '_dw' in col:
			col_inv = col.replace('_dw', '_up')
		else:
			col_inv = col.replace('_up', '_dw')

		vals1 = self.df_summ.loc[self.df_summ.case == case1].iloc[0][col]
		if isinstance(vals1, str): vals1 = eval(vals1)
		if first150 and len(vals1) > 150:
			vals1 = vals1[:150]
		set1 = set(vals1)

		vals1_inv = self.df_summ.loc[self.df_summ.case == case1].iloc[0][col_inv]
		if isinstance(vals1_inv, str): vals1_inv = eval(vals1_inv)

		# all_vals1 = vals1 + vals1_inv

		vals2 = self.df_summ.loc[self.df_summ.case == case2].iloc[0][col]
		if isinstance(vals2, str): vals2 = eval(vals2)
		if first150 and len(vals2) > 150:
			vals2 = vals2[:150]
		set2 = set(vals2)

		vals2_inv = self.df_summ.loc[self.df_summ.case == case2].iloc[0][col_inv]
		if isinstance(vals2_inv, str): vals2_inv = eval(vals2_inv)

		# all_vals2 = vals2 + vals2_inv

		vals1_inverted = [x for x in vals1 if x in vals2_inv]
		vals2_inverted = [x for x in vals2 if x in vals1_inv]

		inverted_degs = list(np.unique(vals1_inverted + vals2_inverted))

		''' --- without contraries -------------'''
		all_vals = list(np.unique(vals1 + vals2))

		fig = plt.figure(figsize=figsize)

		venn2([set2, set1], (case2, case1))

		col2 = "_".join(col.split('_')[:2])
		if in_pathway:
			title = f"{case2} x {case1} for {col2} {self.s_deg_dap}s in pathways\ntotal #{len(all_vals)}"
		else:
			if use_ensembl_degs:
				title = f"{case2} x {case1} for {col2} - all {self.s_deg_dap}s with ensembl\ntotal #{len(all_vals)}"
			else:
				title = f"{case2} x {case1} for {col2} - all {self.s_deg_dap}s\ntotal #{len(all_vals)}"

		if first150:
			title += ' (first 150 genes)'
		plt.title(title, size=title_font_size)

		commons = list(set1.intersection(set2))
		commons.sort()

		only1 = [x for x in vals1 if x not in commons]
		only1.sort()

		only2 = [x for x in vals2 if x not in commons]
		only2.sort()

		s_all	 = f"all {self.s_gene_protein}s (#{len(all_vals)}): {', '.join(all_vals)}"
		s_commons = f"commons (#{len(commons)}): {', '.join(commons)}"
		s_case1   = f"only {case1} (#{len(only1)}): {', '.join(only1)}"
		s_case2   = f"only {case2} (#{len(only2)}): {', '.join(only2)}"


		if len(inverted_degs) == 0:
			s_inverted = 'No inversions'
		else:
			s_inverted = f"inverted (#{len(inverted_degs)}): {', '.join(inverted_degs)}"


		if verbose:
			print(f'\n{s_commons}\n')
			print(s_case1, '\n')
			print(s_case2, '\n')

			if print_inverted: print(s_inverted,'\n')
		else:
			if print_inverted: print(s_inverted)

		if save_fig:
			fnamefig = title_replace(title)
			filefull = osjoin(self.root_figure, fnamefig + '.' + format)
			plt.savefig(filefull, dpi=dpi, format=format, facecolor=facecolor)

		plt.show();

		return fig, s_all, s_commons, s_case1, s_case2, s_inverted, commons, only1, only2, inverted_degs

	# calc_venn_dfenrich_between_two_cases
	def calc_venn_pathways_between_two_cases(self, case1:str, case2:str, figsize:tuple=(8,6),
											title_font_size=15, dpi:int=300, format:str='png', 
											facecolor:str='white', 
											save_fig:bool=True, verbose:bool=False) -> Tuple[object, str, str, str, str, list, list, list]:

		if self.df_summ is None or self.df_summ.empty:
			self.open_enriched_pathways_summary()

		dfqq = self.df_summ.loc[self.df_summ.case == case1]
		if dfqq.empty:
			print("No summary for {case1}")
			return None, '', '', '', '', [], [], []
		
		vals1 = eval(dfqq.iloc[0].pathway_list)
		set1 = set(vals1)

		dfqq = self.df_summ.loc[self.df_summ.case == case2]
		if dfqq.empty:
			print("No summary for {case2}")
			return None, '', '', '', '', [], [], []
		
		vals2 = eval(dfqq.iloc[0].pathway_list)
		set2 = set(vals2)

		all_vals = list(np.unique(vals1 + vals2))

		fig = plt.figure(figsize=figsize)

		venn2([set1, set2], (case1, case2))

		title = f"{case1} x {case2}: total enriched pathways = #{len(all_vals)}"
		plt.title(title, size=title_font_size)

		commons = list(set1.intersection(set2))
		commons.sort()

		only1 = [x for x in vals1 if x not in commons]
		only1.sort()

		only2 = [x for x in vals2 if x not in commons]
		only2.sort()

		s_all	 = f"all pathways (#{len(all_vals)}): {', '.join(all_vals)}"
		s_commons = f"commons (#{len(commons)}): {', '.join(commons)}"
		s_case1   = f"only {case1} (#{len(only1)}): {', '.join(only1)}"
		s_case2   = f"only {case2} (#{len(only2)}): {', '.join(only2)}"

		if verbose:
			print(s_commons, '\n')
			print(s_case1, '\n')
			print(s_case2)

		if save_fig:
			fnamefig = title_replace(title)
			filefull = osjoin(self.root_figure, fnamefig + '.' + format)
			plt.savefig(filefull, dpi=dpi, format=format, facecolor=facecolor)

		plt.show()

		return fig, s_all, s_commons, s_case1, s_case2, commons, only1, only2


	def summary_degs_up_down(self, per_biotype:bool=False, ensembl:bool=False,
							 save_file:bool=False, before_best_cutoff:bool=False,
							 verbose:bool=False) -> pd.DataFrame:

		dic = {}; i=-1
		for case in self.case_list:
			if before_best_cutoff:
				ret, _, _, _ = self.open_case_params(case, LFC_cut=1, lfc_FDR_cut=0.05, ptw_FDR_cut=0.05, verbose=False)
			else:
				ret, _, _, _ = self.open_case(case, verbose=False)
			
			if not ret:
				continue

			i+=1
			dic[i] = {}
			dic2 = dic[i]
			dic2['case'] = case
			dic2['tot_measured'] = len(self.dflfc_ori)

			if per_biotype:
				''' upregulated all or ensembl '''
				if ensembl:
					df = self.dfdegs_up_ensembl
					sufix = 'up_ens'
				else:
					df = self.dfdegs_up
					sufix = 'up'

				if df is not None and not df.empty:
					dfg = self.get_df_biotypes(df)

					for j in range(len(dfg)):
						row = dfg.iloc[j]
						dic2[f"{row.biotype}_{sufix}"] = row['n']

				''' dwregulated all or ensembl '''
				if ensembl:
					df = self.dfdegs_dw_ensembl
					sufix = 'dw_ens'
				else:
					df = self.dfdegs_dw
					sufix = 'dw'

				if df is not None and not df.empty:
					dfg = self.get_df_biotypes(df)

					for j in range(len(dfg)):
						row = dfg.iloc[j]
						dic2[f"{row.biotype}_{sufix}"] = row['n']

			else:
				dic2['n_degs']	= self.n_degs
				dic2['n_degs_up'] = self.n_degs_up
				dic2['n_degs_up_ensembl'] = self.n_degs_up_ensembl
				dic2['n_degs_dw'] = self.n_degs_dw
				dic2['n_degs_dw_ensembl'] = self.n_degs_dw_ensembl
				dic2['n_ptw'] = len(self.df_enr)

		dfa = pd.DataFrame(dic).T
		dfa = dfa.infer_objects().fillna(0)

		cols = list(dfa.columns)[1:]
		dfa[cols] = dfa[cols].astype(int)

		dfa = dfa.set_index('case')
		dfa = dfa.T

		if save_file:
			if per_biotype:
				fname = f'resume_of_degs_per_biotype_{per_biotype}_ensembl_{ensembl}_BCA_{not before_best_cutoff}.tsv'
			else:
				fname = f'resume_of_degs_BCA_{not before_best_cutoff}.tsv'

			pdwritecsv(dfa.T, fname, self.root_ressum, index=True, verbose=verbose)

		return dfa

	def barplot_up_down_genes_per_case(self, per_biotype:bool=False, ensembl:bool=False,
									   before_best_cutoff:bool=False, title:str='',
									   log10:bool=False, show_val:bool=False,
									   yaxis_title:str='', width:int=800, height:int=600,
									   plot_bgcolor:str='lightgray', nround:int=2, verbose:bool=False):

		

		dfa = self.summary_degs_up_down(per_biotype=per_biotype, ensembl=ensembl, 
										before_best_cutoff=before_best_cutoff, verbose=verbose)
		dfa = dfa.reset_index()
		cols = list(dfa.columns)
		cols[0] = 'deg_col'
		dfa.columns = cols
		dfa.index.name = 'index'

		dfa = dfa[dfa.deg_col != 'tot_measured']
		dfa.reset_index(inplace=True, drop=True)

		before = 'before best cutoff' if before_best_cutoff else 'with best cutoff'

		if title == '':
			if not per_biotype:
				title = f'Up and Down {self.s_deg_dap}s {before}'
			else:
				if not ensembl:
					title = f'Up and Down {self.s_deg_dap}s per biotype {before}'
				else:
					title = f'Up and Down {self.s_deg_dap}s per biotype having ensembl_id {before}'


		if yaxis_title == '':
			yaxis_title = f'# {self.s_deg_dap}s'

		if log10:
			yaxis_title = 'log10 ' + yaxis_title

		fig = go.Figure()

		cols = list(dfa.columns)
		
		pos_plot = -1
		x_tick_num_list, x_tick_label_list = [], []
		for case in self.list_order:
			pos_plot += 1

			general_max = dfa[case].max()
			if log10:
				general_max = np.log10(general_max)

			icase=0
			for i in range(len(dfa)):
				icase += 1
				if icase==3:
					x_tick_label_list.append(case)
				else:
					x_tick_label_list.append('')

				row = dfa.iloc[i]

				color = self.my_colors[i]
				col   = row.deg_col
				val0  = float(row[case])

				if log10:
					val = np.log10(val0) if isfloat(val0) and val0 > 0 else 0
				else:
					val = val0

				pos_plot += 1
				x_tick_num_list.append(pos_plot)

				fig.add_trace(go.Bar(x=[pos_plot], y=[val], marker_color=color, name=f"{case}: {col}"))
				# fig.add_trace(go.Bar(x=[case], y=[val],  marker_color=color, name=f"{case}: {col}"))

				if show_val:
					if pd.notnull(val) and val != 0:
						# in graphic_lib
						delta_y = define_delta_y(val, general_max)

						fig.add_annotation(text=str(np.round(val0,nround)), x=pos_plot, y=val+delta_y, showarrow=False)


		fig.update_layout(
					autosize=True,
					title=title,
					width=width,
					height=height,
					plot_bgcolor=plot_bgcolor,
					xaxis_title="cases",
					# xaxis_showticklabels=False,
					xaxis=dict(
						tickmode='array',
						tickvals=x_tick_num_list,
						ticktext=x_tick_label_list
					),
					yaxis_title=yaxis_title,
					showlegend=True,
					font=dict(
						family="Arial",
						size=14,
						color="Black"
					)

		)


		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		return fig, dfa

	def count_sampling(self, geneset_num_list:List=[0, 1, 2, 4, 5, 7], prompt_verbose:bool=False):

		dic={}; i=-1
		for geneset_num in geneset_num_list:
			self.set_db(geneset_num, verbose=True)

			s_start = f"enricher_{self.geneset_lib}"

			for case in self.case_list:
				files = [x for x in os.listdir(self.root_enrich_sampling) if x.startswith(s_start) and case in x]
				if prompt_verbose: print("\tcase", case, len(files))

				i+=1
				dic[i] = {}
				dic2 = dic[i]

				dic2['geneset_num'] = self.geneset_num
				dic2['geneset_lib'] = self.geneset_lib
				dic2['case'] = case
				dic2['n'] = len(files)

			if prompt_verbose: print('')

		dfa = pd.DataFrame(dic).T
		return dfa

	def plot_degs_in_pathways_vs_toi_per_case(self, selected_toi_col:str='toi4_median', title:str='',
								 width:int=1100, height:int=600, plot_all_dfi:bool=True, sel_colors:List=[],
								 plot_bgcolor:str='lightgray', verbose:bool=False):
		if title is None or title == '':
			title=f"The {selected_toi_col} space"

		if sel_colors is None or sel_colors == []:
			sel_colors = self.my_colors

		yaxis_title1 = f"# {self.s_deg_dap}s in pathways"
		yaxis_title2 = "# of pathways"

		dfcut = self.build_all_cutoffs_table(selected_toi_col, force=False, verbose=verbose)

		dfcut = dfcut.sort_values(['case', selected_toi_col], ascending=[True, True])
		dfcut.reset_index(inplace=True, drop=True)

		subplot_titles=[f"# {self.s_deg_dap}s in pathways", "# pathways", f"# pathways per {self.s_deg_dap}s in pathways"]
		fig = make_subplots(rows=3, cols=1, subplot_titles=subplot_titles)

		df_all_fdr = self.open_all_fdr_lfc_correlation()

		'''------------ calc max -------------'''
		maxi_x = 0
		for icase, case in enumerate(self.case_list):

			df_fdr = df_all_fdr[df_all_fdr.case == case]
			if df_fdr.empty:
				print(f"Error: no correlation was calculated for case '{case}'")
				raise Exception('stop: plot_degs_in_pathways_vs_toi_per_case()')

			for fdr in df_fdr.fdr:
				df2 = dfcut[ (dfcut.case == case) & (dfcut.lfc_FDR_cut == fdr)  & (dfcut.med_max_ptw == 'median')]
				if df2.empty:
					continue

				maxi_x2 = np.round(df2[selected_toi_col].max(), 3) + 0.005
				if maxi_x2 > maxi_x:
					maxi_x = maxi_x2


		'''------------ plot loop ------------'''
		dic_visible = {}

		if plot_all_dfi:
			for icase, case in enumerate(self.case_list):
				 
				dfi = self.calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(case)

				df_fdr = df_all_fdr[df_all_fdr.case == case]
				if df_fdr.empty:
					print(f"Error: no correlation was calculated for case '{case}'")
					raise Exception('stop: plot_degs_in_pathways_vs_toi_per_case()')

				dic_visible[case] = 0

				is_visible = True if icase == 0 else False
				i = -1;
				for fdr in df_fdr.fdr:
					i += 1

					fdr = np.round(fdr, 2)

					dfi2 = dfi[ (dfi.case == case) & (dfi.lfc_FDR_cut == fdr) ].copy()
					if dfi2.empty:
						continue

					dfi2['pathways_per_degs'] = dfi2.n_pathways / dfi2.n_degs_in_pathways
					dfi2 = dfi2.sort_values(selected_toi_col, ascending=True)
					dfi2.reset_index(inplace=True, drop=True)

					text_ini = f'case {case}<br>lfc_FDR_cut={fdr:.2f}'

					hovertext_list = []
					for _, row in dfi2.iterrows():
						text =  f'LFC_cut={row.LFC_cut:.2f}<br>ptw_FDR_cut={row.ptw_FDR_cut:.2f}<br>toi4={row.toi4_median:.3f}'
						text += f'<br>---------<br># {self.s_deg_dap}s in ptws={row.n_degs_in_pathways}<br># pathways={row.n_pathways}'
						hovertext_list.append(text_ini + '<br>' + text)

					'''--- mask of trues and falses '''
					dic_visible[case] += 3

					name1 = f"{case} fdr={fdr:.3f} for {self.s_deg_dap}s"
					name2 = f"{case} fdr={fdr:.3f} for pathways"
					name3 = f"{case} fdr={fdr:.3f} for pathways/{self.s_deg_dap}s"

					color = sel_colors[i]

					fig.add_trace(go.Scatter(x=dfi2[selected_toi_col], y=dfi2.n_degs_in_pathways,  
											hovertext=hovertext_list, line=dict(dash='dash'), 
											marker_color=color, name=name1), row=1, col=1 )
					
					fig.add_trace(go.Scatter(x=dfi2[selected_toi_col], y=dfi2.n_pathways,  
											hovertext=hovertext_list, line=dict(dash='dash'), 
											marker_color=color, name=name2), row=2, col=1 )
					
					fig.add_trace(go.Scatter(x=dfi2[selected_toi_col], y=dfi2.pathways_per_degs, 
											hovertext=hovertext_list, line=dict(dash='dash'),
											marker_color=color, name=name3), row=3, col=1 )
					fig.update_layout(title=title,
									width=width,
									height=height,
									xaxis_title=selected_toi_col,
									yaxis_title=yaxis_title1,
									yaxis2_title=yaxis_title2,
									xaxis2_title=selected_toi_col,
									xaxis=dict(range=[0, maxi_x]),
									plot_bgcolor=plot_bgcolor,
									legend_title="cases",
									showlegend=True)
		else:
			for icase, case in enumerate(self.case_list):

				df_fdr = df_all_fdr[df_all_fdr.case == case]
				if df_fdr.empty:
					print(f"Error: no correlation was calculated for case '{case}'")
					raise Exception('stop: plot_degs_in_pathways_vs_toi_per_case()')

				dic_visible[case] = 0

				is_visible = True if icase == 0 else False
				i = -1;
				for fdr in df_fdr.fdr:
					i += 1

					fdr = np.round(fdr, 2)

					df2 = dfcut[(dfcut.case == case) & (dfcut.lfc_FDR_cut == fdr)  & 
								(dfcut.med_max_ptw == 'median') &
								(dfcut.n_pathways > 0) & (dfcut.n_degs_in_pathways > 0) ].copy()
					df2.reset_index(inplace=True, drop=True)
					if df2.empty:
						continue

					df2['pathways_per_degs'] = df2.n_pathways / df2.n_degs_in_pathways

					text_ini = f'case {case}<br>lfc_FDR_cut={fdr:.2f}'

					hovertext_list = []
					for _, row in df2.iterrows():
						text =  f'LFC_cut={row.LFC_cut:.2f}<br>ptw_FDR_cut={row.ptw_FDR_cut:.2f}<br>toi4={row.toi4_median:.3f}'
						text += f'<br>---------<br># {self.s_deg_dap}s in ptws={row.n_degs_in_pathways}<br># pathways={row.n_pathways}'
						hovertext_list.append(text_ini + '<br>' + text)

					'''--- mask of trues and falses '''
					dic_visible[case] += 3

					name1 = f"{case} fdr={fdr:.3f} for {self.s_deg_dap}s"
					name2 = f"{case} fdr={fdr:.3f} for pathways"
					name3 = f"{case} fdr={fdr:.3f} for pathways/{self.s_deg_dap}s"

					color = sel_colors[i]

					fig.add_trace(go.Scatter(x=df2[selected_toi_col], y=df2.n_degs_in_pathways,  
											hovertext=hovertext_list, hoverinfo="text",
											marker_color=color, visible=is_visible, name=name1), row=1, col=1 )
					
					fig.add_trace(go.Scatter(x=df2[selected_toi_col], y=df2.n_pathways, 
											hovertext=hovertext_list, hoverinfo="text",
											marker_color=color, visible=is_visible, name=name2), row=2, col=1 )
					
					fig.add_trace(go.Scatter(x=df2[selected_toi_col], y=df2.pathways_per_degs, 
											hovertext=hovertext_list, hoverinfo="text",
											marker_color=color, visible=is_visible, name=name3), row=3, col=1 )

					fig.update_layout(title=title,
									width=width,
									height=height,
									xaxis_title=selected_toi_col,
									yaxis_title=yaxis_title1,
									yaxis2_title=yaxis_title2,
									xaxis2_title=selected_toi_col,
									xaxis=dict(range=[0, maxi_x]),
									plot_bgcolor=plot_bgcolor,
									legend_title="cases",
									showlegend=True)
					
		# add dropdown menus to the figure
		buttons=[]
		for case in self.case_list:
			buttons.append(dict(method='update',
								label=case,
								visible=True,
								args=[ {'visible': list(sum( [tuple([True]  * dic_visible[case2]) if case == case2 else \
															  tuple([False] * dic_visible[case2]) for case2 in self.case_list], () ))} ]
								)
						  )

		# some adjustments to the updatemenus
		updatemenu = []
		your_menu = dict()
		updatemenu.append(your_menu)

		updatemenu[0]['buttons'] = buttons
		updatemenu[0]['direction'] = 'down'
		updatemenu[0]['showactive'] = True
		updatemenu[0]['showactive'] = True
		updatemenu[0]['x'] = 1
		updatemenu[0]['y'] = 1.2

		fig.update_layout(title=title,
						  width=width,
						  height=height,
						  xaxis_title=selected_toi_col,
						  yaxis_title=yaxis_title1,
						  yaxis2_title=yaxis_title2,
						  xaxis2_title=selected_toi_col,
						  xaxis=dict( range=[0, maxi_x]),
						  xaxis2=dict( range=[0, maxi_x]),
						  legend_title="cases",
						  showlegend=True,
						  plot_bgcolor=plot_bgcolor,
						  updatemenus=updatemenu)

		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		return fig

	def plot_degs_vs_lfc_per_fdr_per_case(self, title:str='',
								 width:int=1100, height:int=600, sel_colors:List=[],
								 plot_bgcolor:str='lightgray', verbose:bool=False):
		'''
			removed: selected_toi_col:str='toi4_median', plot_all_dfi:bool=True, 
		'''

		if title is None or title == '':
			title = f'scatter plot - {self.s_deg_dap}s versus LFC_cut per FDR'

		if sel_colors is None or sel_colors == []:
			sel_colors = self.my_colors

		xaxis_title = f"abs LFC"
		yaxis_title = f"# {self.s_deg_dap}s"

		dfsim = self.open_simulation_table()

		dfsim = dfsim.sort_values(['case', 'lfc_FDR_cut', 'LFC_cut'], ascending=[True, False, False])

		fig = go.Figure()

		dic_visible = {}
		for icase in range(len(self.case_list)):
			case = self.case_list[icase]

			dic_visible[case] = 0
			is_visible = True if icase == 0 else False
			i = -1;
			for i in range(len(self.fdr_list)):
				lfc_FDR_cut = self.fdr_list[i]
				color = sel_colors[i]
				name = f"{lfc_FDR_cut:.3f}"

				dfsim2 = dfsim[ (dfsim.case == case) & (dfsim.lfc_FDR_cut == lfc_FDR_cut)]
				if dfsim2.empty:
					continue

				dic_visible[case] += 1

				text_ini = f'case {case}<br>lfc_FDR cutoff={lfc_FDR_cut:.3f}'

				hovertext_list = []
				for j in range(len(dfsim2)):
					row = dfsim2.iloc[j]
					text =  f'LFC_cut={row.LFC_cut:.3f}'
					text += f'# {self.s_deg_dap}s {row.n_degs}<br># Up={row.n_degs_up} Down={row.n_degs_dw}'
					hovertext_list.append(text_ini + '<br>' + text)

				fig.add_trace(go.Scatter(x=dfsim2.LFC_cut, y=dfsim2.n_degs, hovertext=hovertext_list, hoverinfo="text",
										 mode='markers', marker={'color':color}, visible=is_visible, name=name))

			fig.update_layout(
						autosize=True,
						title=title,
						width=width,
						height=height,
						xaxis_title=xaxis_title,
						yaxis_title=yaxis_title,
						showlegend=True,
						legend_title='lfc_FDR_cutoff',
						plot_bgcolor=plot_bgcolor,
						font=dict(
							family="Arial",
							size=14,
							color="Black"
						)
			)

		# add dropdown menus to the figure
		buttons=[]
		for case in self.case_list:
			buttons.append(dict(method='update',
								label=case,
								visible=True,
								args=[ {'visible': list(sum( [tuple([True]  * dic_visible[case2]) if case == case2 else \
															  tuple([False] * dic_visible[case2]) for case2 in self.case_list], () ))} ]
								)
						  )

		# some adjustments to the updatemenus
		updatemenu = []
		your_menu = dict()
		updatemenu.append(your_menu)

		updatemenu[0]['buttons'] = buttons
		updatemenu[0]['direction'] = 'down'
		updatemenu[0]['showactive'] = True
		updatemenu[0]['showactive'] = True
		updatemenu[0]['x'] = 1
		updatemenu[0]['y'] = 1.2

		fig.update_layout(
			autosize=True,
			title=title,
			width=width,
			height=height,
			xaxis_title=xaxis_title,
			yaxis_title=yaxis_title,
			showlegend=True,
			legend_title='lfc_FDR_cutoff',
			font=dict(
				family="Arial",
				size=14,
				color="Black"
			),
			plot_bgcolor=plot_bgcolor,
			updatemenus=updatemenu
		)

		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		return fig


	def barplot_sampling_cutoffs(self, title:str='Sampling cutoffs per geneset', yaxis_title:str='number of samples',
								 width:int=1100, height:int=600, geneset_num_list:List=[0, 1, 2, 4, 5, 7],
								 colors:List=['green', 'red', 'blue', 'brown', 'yellow', 'cyan', 'lightgreen', 'pink', 'gray', 'lightblue'],
								 plot_bgcolor:str='lightgray', prompt_verbose:bool=False, verbose:bool=False):

		dfa = self.count_sampling(geneset_num_list=geneset_num_list, prompt_verbose=prompt_verbose)

		geneset_lista = dfa.geneset_lib.unique()
		colors_geneset = colors[:len(geneset_lista)]

		_n_rows = int(np.ceil(len(self.case_list)/2))
		fig = make_subplots(rows=_n_rows, cols=2, subplot_titles=self.case_list)

		nrow=1; ncol=0
		for case in self.case_list:

			dfa2 = dfa[dfa.case == case].copy()
			dfa2 = dfa2.sort_values('geneset_lib')
			dfa2.reset_index(inplace=True, drop=True)

			ncol += 1
			if ncol == 3:
				ncol = 1
				nrow += 1

			fig.add_trace(go.Bar(x=dfa2.geneset_lib, y=dfa2.n, marker_color=colors_geneset, name=case), row=nrow, col=ncol)

		fig.update_layout(
					autosize=True,
					title=title,
					width=width,
					height=height*_n_rows,
					plot_bgcolor=plot_bgcolor,
					xaxis_title="",
					xaxis2_title="",
					yaxis_title=yaxis_title,
					showlegend=False,
					font=dict(
						family="Arial",
						size=14,
						color="Black"
					)
		)

		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		return fig, dfa


	def barplot_degs_summary(self, title:str='', yaxis_title:str='',
							 colors:List=[], width:int=900, height:int=600,
							 verbose:bool=False) -> object:


		if self.df_summ is None or self.df_summ.empty:
			_ = self.calc_degs_and_pathways_summary(force=False, prompt_verbose=False, verbose=verbose)

		if self.df_summ is None or self.df_summ.empty:
			print("Error: no df_summ")
			return None

		if title == '':
			title = f"Number of {self.s_deg_dap}s and pathways per case"
		if yaxis_title == '':
			yaxis_title = f"# {self.s_deg_dap}s and pathways"

		if colors == []:
			colors = ['darkgreen', 'darkred', 'navy', 'olivedrab', 'red', 'blue', 'darkcyan']
			'''
			'orange', 'brown', 'darksalmon', 'marron',
			'magenta', 'darkturquoise', 'orange', 'darkred', 'indigo', 'magenta',  'black',
			'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgrey', 'olivedrab', 'navy']
			'''

		fig = go.Figure()

		cols = ['n_degs', 'n_degs_up', 'n_degs_dw', 'n_degs_in_pathways',
				'n_degs_up_ensembl_in_pathways', 'n_degs_dw_ensembl_in_pathways', 'n_pathways']
		for i in range(len(cols)):
			col = cols[i]
			color = colors[i]

			fig.add_trace(go.Bar(x=self.df_summ.case, y=self.df_summ[col], name=col, marker_color=color))

		fig.update_layout(
					autosize=True,
					title=title,
					width=width,
					height=height,
					xaxis_title="cases",
					yaxis_title=yaxis_title,
					showlegend=True,
					font=dict(
						family="Arial",
						size=14,
						color="Black"
					)
		)

		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		return fig

	def calc_pPMI_summary_and_pivot_tables(self, pathway_list:List=[], title_head:str='', with_pPMI_obs:bool=False, 
								do_case_translate:bool=True, echo_pathway:bool=True, 
								force:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, str, pd.DataFrame]:
		'''
			pPMI - pseudo-Pathway Modulation Index
		'''

		dff, text_all, df_sum_ptw = \
			self.calc_pPMI_summary(pathway_list=pathway_list, title_head=title_head, with_pPMI_obs=with_pPMI_obs, 
						  		   do_case_translate=do_case_translate, echo_pathway=echo_pathway, 
								   force=force, verbose=verbose)

		if dff is None or dff.empty or df_sum_ptw is None or df_sum_ptw.empty:
			self.dfpiv = pd.DataFrame()
			self.selected_pivot_pathway_list = []
			return pd.DataFrame(), dff, text_all, df_sum_ptw

		col_pPMI = 'pPMI_norm' if self.pPMI_normalized else 'pPMI_total'

		dfpiv = pd.pivot_table(dff, values=col_pPMI, index='pathway', columns='case', fill_value=None)
		self.dfpiv = dfpiv

		# save in memory list of pathways
		selected_pivot_pathway_list = []

		for pathway in dfpiv.index:
			selected_pivot_pathway_list.append(pathway)

		self.selected_pivot_pathway_list = selected_pivot_pathway_list

		return dfpiv, dff, text_all, df_sum_ptw


	def save_all_summaries_in_markdown(self, pathway_list:List, title_head:str,
				lfc_threshold:float=1, corr_threshold:float=0.6,
				modulations:List=[ ('all', 'up'),('all', 'dw'),('all', 'none'), ('corr', 'up'), ('corr', 'dw'), ('corr', 'none') ],
				width:int=1000, height:int=700, line_width:int=2, 
				plot_bgcolor:str='lightgray', vertical_spacing:float=0.25, verbose:bool=False):

		if verbose: print(">>> saving pathway summary\n")

		if len(pathway_list) > self.max_pathways:
			print(f"Warning: it is not recomended more than {self.max_pathways}, you sent {len(pathway_list)}")
			return

		for pathway in pathway_list:
			if verbose: print("\t", pathway)
			_, _, _, _ =  self.save_pathway_summary_in_markdown(pathway=pathway, title_head=title_head, 
									with_pPMI_obs=True,  do_case_translate=True, echo_pathway=True, 
									lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
									modulations=modulations, width=width, height=height, line_width=line_width,
									fontsize=18, leg_fontsize=18, 
									tick_fontsize=18, title_fontsize=20,
									plot_bgcolor=plot_bgcolor, vertical_spacing=vertical_spacing, verbose=verbose)

		if verbose: 
			print(f"\nAll summaries saved in {self.root_pathway}\n")


	def save_pathway_summary_in_markdown(self, pathway:str, title_head:str,
								with_pPMI_obs:bool=False, do_case_translate:bool=True, echo_pathway:bool=True, 
								lfc_threshold:float=1, corr_threshold:float=0.6,
								modulations:List=[ ('all', 'up'),('all', 'dw'),('all', 'none'), ('corr', 'up'), ('corr', 'dw'), ('corr', 'none') ],
								width:int=1000, height:int=700, line_width:int=2, 
								fontsize:int=16, leg_fontsize:int=16, 
								tick_fontsize:int=16, title_fontsize:int=18,
								plot_bgcolor:str='lightgray', vertical_spacing:float=0.20,
								force_calc_pPMI_summary:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, str, pd.DataFrame]:
		'''
			pPMI - pseudo-Pathway Modulation Index
		'''
		pathway_list=[pathway]

		dff, text_all, df_sum_ptw = \
			self.calc_pPMI_summary(pathway_list=pathway_list, title_head=title_head, with_pPMI_obs=with_pPMI_obs, 
						  		   do_case_translate=do_case_translate, echo_pathway=echo_pathway, 
								   force=force_calc_pPMI_summary, verbose=verbose)

		if dff is None or dff.empty or df_sum_ptw is None or df_sum_ptw.empty:
			self.dfpiv = pd.DataFrame()
			self.selected_pivot_pathway_list = []
			return pd.DataFrame(), dff, text_all, df_sum_ptw

		col_pPMI = 'pPMI_norm' if self.pPMI_normalized else 'pPMI_total'

		dfpiv = pd.pivot_table(dff, values=col_pPMI, index='pathway', columns='case', fill_value=None)

		figname_list = self.plot_all_pathways_LFC_modulations(pathway_list=pathway_list, 
										lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
										modulations=modulations, width=width, height=height, line_width=line_width,
										fontsize=fontsize, leg_fontsize=leg_fontsize, 
										tick_fontsize=tick_fontsize, title_fontsize=title_fontsize,
										plot_bgcolor=plot_bgcolor, vertical_spacing=vertical_spacing, verbose=verbose)

		root_fig= self.root_figure.replace(self.root_project, '')
		root_fig = root_fig[1:] if root_fig.startswith('/') else root_fig
		root_fig = osjoin('..', root_fig)

		for figname, pathway in zip(figname_list, pathway_list):
			# filename = osjoin(root_fig, figname) + '.png'
			filename = osjoin('figures/', figname) + '.png'
			# print(">>> ", filename)
			# s_picture = f'\n\n<iframe width="100%" height="100%" src="{filename}"" title="{pathway}" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>\n'
			s_picture = f'\n\n![{pathway}]({filename})\n'
			text_all += s_picture
	
		fname = title_replace(pathway) + '.txt'
		write_txt(text_all, fname, self.root_pathway)

		return dfpiv, dff, text_all, df_sum_ptw
	

	# calc_one_pathway_case_gene_modulations
	def calc_one_pathway_case_symb_LFC(self, pathway_and_id:List, case:str, verbose:bool=False) -> pd.DataFrame:

		if not self.open_case_simple(case, verbose=verbose):
			return pd.DataFrame()

		dflfc = self.dflfc_ori

		pathway2, pathway_id = pathway_and_id
		pathway, genes_in_pathway = self.reactome_find_genes_in_pathway(pathway_id, _type='pathway_id')

		assert pathway.lower().strip() == pathway2.lower().strip(), f"Diferent pathways '{pathway}' and '{pathway2}'"

		if len(genes_in_pathway) == 0:
			print(f"No genes were found for '{pathway_id} - {pathway}'")
			return pd.DataFrame()

		dflfc_ptw = dflfc[dflfc.symbol.isin(genes_in_pathway)]
		if len(dflfc_ptw) == 0:
			print(f"No genes were found for '{pathway_id} - {pathway}' in LFC table")
			return pd.DataFrame()

		cols = ['symbol', 'lfc', 'abs_lfc', 'fdr']
		dflfc_ptw = dflfc_ptw[cols].copy()
		dflfc_ptw.reset_index(inplace=True, drop=True)

		dflfc_ptw['pathway'] = pathway
		dflfc_ptw['pathway_id'] = pathway_id
		dflfc_ptw['case'] = case

		cols = ['case', 'pathway_id', 'pathway', 'symbol', 'lfc', 'abs_lfc', 'fdr']
		dflfc_ptw = dflfc_ptw[cols]

		return dflfc_ptw
	
	# calc_one_pathway_all_cases_gene_modulations
	def calc_one_pathway_all_cases_symb_LFC(self, pathway_and_id:list, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		# pathway, pathway_id = pathway_and_id
		pathway, _ = pathway_and_id

		fname = self.fname_one_pathway_symb_LFC%(pathway, self.geneset_lib, self.normalization)
		fname = title_replace(fname)
		filename = osjoin(self.root_ressum, fname)

		if exists(filename) and not force:
			dff = pdreadcsv(fname, self.root_ressum, verbose=verbose)
			return dff

		df_list = []
		for case in self.case_list:
			dfa =  self.calc_one_pathway_case_symb_LFC(pathway_and_id, case, verbose=verbose)
			if dfa is None or dfa.empty: continue
			df_list.append(dfa)

		if df_list == []:
			print(f"No gene modulations was found for {pathway}.")
			return pd.DataFrame()

		dff = pd.concat(df_list)
		dff = dff.sort_values(['case', 'pathway_id', 'symbol'])
		dff.reset_index(inplace=True, drop=True)

		_ = pdwritecsv(dff, fname, self.root_ressum, verbose=verbose)

		return dff

	def calc_pivot_one_pathway_LFC(self, pathway_and_id:List, force:bool=False, 
								   verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame]:
		'''
			dfpiv is for one pathway
			pivot tabel: symbs x cases
			cells: LFC

			returns: dfpiv, dff
		'''

		dff = self.calc_one_pathway_all_cases_symb_LFC(pathway_and_id, force=force, verbose=verbose)

		if dff is None or dff.empty:
			dfpiv = pd.DataFrame()

			self.dfpiv_symbs_cases = dfpiv
			self.selected_pivot_symb_list = []

			return dfpiv, dff

		# dfpiv symbs x cases
		dfpiv = pd.pivot_table(dff, values='lfc', index='symbol', columns='case', fill_value=None)

		self.dfpiv_symbs_cases = dfpiv
		self.selected_pivot_symb_list = dfpiv.index.to_list

		return dfpiv, dff


	def df_to_plotly(self, df):
		return {'z': df.values.tolist(),
				'x': df.columns.tolist(),
				'y': df.index.tolist()}

	def plot_pathway_heatmap_simple(self, pathway_list:List, sel_list:List, title_head:str='', all_cases:bool=True,
									maxLen:int=50,  with_pPMI_obs:bool=False, do_case_translate:bool=True,
									width:int=1100, height0:int=350, del_height:int=75,
									line_width:int=2, plot_bgcolor:str='lightgray',  
									lfc_threshold:float=1, corr_threshold:float=0.6,
									modulations:List=[ ('all', 'up'),('all', 'dw'),('all', 'none'), ('corr', 'up'), ('corr', 'dw'), ('corr', 'none') ],   
									plot_bars:bool=True, plot_lines:bool=False,
									fontsize_circos:int=14, space_factor:float=3.,
									verbose:bool=False) -> Tuple[object, pd.DataFrame, str, str, pd.DataFrame, list]:

		if not isinstance(pathway_list, list) or pathway_list == []:
			print("Pathway list is empty or not a list.")
			return None, pd.DataFrame(), '', '', pd.DataFrame(), []
		
		if len(pathway_list) > self.max_pathways:
			print(f"Warning: it is not recomended more than {self.max_pathways}, you sent {len(pathway_list)}")
			return None, pd.DataFrame(), '', '', pd.DataFrame(), []

		dfpiv, dff, text_all, _ = \
			self.calc_pPMI_summary_and_pivot_tables(pathway_list=pathway_list, title_head=title_head, with_pPMI_obs=with_pPMI_obs, 
										   			do_case_translate=do_case_translate, verbose=verbose)

		self.dff = dff

		maxi = int(np.max(np.ceil( dff.pPMI_total.to_list())))
		mini = int(np.min(np.floor(dff.pPMI_total.to_list())))
		if mini < 0:
			mini = -mini
		zlim = maxi if maxi > mini else mini

		zmin = -zlim
		zmax = +zlim

		not_in = [x for x in pathway_list if x not in self.selected_pivot_pathway_list]

		if len(not_in) > 0:
			print("Some pathways look wrong, please fix the pathway list.")
			print(">>>", "; ".join(not_in))
			return None, pd.DataFrame(), '', '', pd.DataFrame(), []

		fig = go.Figure()

		if all_cases:
			dfpiv2 =  pd.DataFrame(dfpiv.loc[pathway_list, :].copy())
			title = f"pPMI heatmap for '{title_head}'"
		else:
			cols = self.case
			dfpiv2 =  pd.DataFrame(dfpiv.loc[pathway_list, cols].copy())
			title = f"Heatmap for '{title_head}' for {self.case}"

		lista = [str(x) for x in dfpiv2.index.to_list()]
		new_indexes = [break_line_per_length(x, sep='<br>', maxLen=maxLen) for x in lista]
		dfpiv2 = dfpiv2.reset_index(drop=True)
		dfpiv2 = dfpiv2.set_index(pd.Series(new_indexes))

		stri = 'normalized' if self.pPMI_normalized else 'not normalized'
		title += " - " + stri

		df_vectors = self.df_to_plotly(dfpiv2)

		text_list=[]; lista_ratio=[]
		for i in range(len(dfpiv2)):
			row  = dfpiv2.iloc[i]
			pathway = row.name

			lista =[]; lista2=[]
			for j in range(len(row.index)):
				col = row.index[j]
				case = col
				# pPMI = row[col]

				row_dff = dff[ (dff.case == case) & (dff.pathway == pathway)]
				if row_dff.empty:
					# print(i, 'empty', pathway, case, row.index)
					lista.append(None)
					# lista.append(f'<br>x:{case}<br>y:{pathway}<br>pPMI:{pPMI:.3f}')
				else:
					# print(i, 'Ok', pathway, case, row.index)
					row_dff = row_dff.iloc[0]
					mod_up_list = row_dff['mod_up_in_pathway']
					lfc_up_list = row_dff['lfc_up']

					mod_up_list = mod_up_list if isinstance(mod_up_list, list) else eval(mod_up_list)
					lfc_up_list = lfc_up_list if isinstance(lfc_up_list, list) else eval(lfc_up_list)

					# up_mod_list = self.build_modulation_list(mod_up_list, lfc_up_list)

					mod_dw_list = row_dff['mod_dw_in_pathway']
					lfc_dw_list = row_dff['lfc_dw']

					mod_dw_list = mod_dw_list if isinstance(mod_dw_list, list) else eval(mod_dw_list)
					lfc_dw_list = lfc_dw_list if isinstance(lfc_dw_list, list) else eval(lfc_dw_list)

					n_mod_in_pathway	= row_dff['n_mod_in_pathway']
					n_mod_up_in_pathway = row_dff['n_mod_up_in_pathway']
					n_mod_dw_in_pathway = row_dff['n_mod_dw_in_pathway']

					if mod_up_list == []:
						text_up_in_pathway = '---'
					else:
						mat_modulation = self.build_modulation_list(mod_up_list, lfc_up_list)
						text_up_in_pathway = break_list(mat_modulation, n_elems=5, CR='<br>')


					mod_dw_in_pathway = row_dff['mod_dw_in_pathway']
					if mod_dw_in_pathway == '[]':
						text_dw_in_pathway = '---'
					else:
						mat_modulation = self.build_modulation_list(mod_dw_list, lfc_dw_list)
						text_dw_in_pathway = break_list(mat_modulation, n_elems=5, CR='<br>')

					s_ratio = f"<br>#Up={n_mod_up_in_pathway} #Dw={n_mod_dw_in_pathway} #Modulated={n_mod_in_pathway}"
					stri_up_dw = s_ratio+f'<br><br>Up({n_mod_up_in_pathway}): {text_up_in_pathway}<br><br>Dw({n_mod_dw_in_pathway}): {text_dw_in_pathway}'
					lista.append(stri_up_dw)
					lista2.append(stri_up_dw)

			text_list.append(lista)
			lista_ratio.append(lista2)

		fig.add_trace(go.Heatmap(df_vectors, zmin=zmin, zmax=zmax, text=text_list, colorscale='RdBu_r') )

		n_lines = len(dfpiv2)
		n_lines -= 3
		if n_lines < 0: n_lines=0

		height = height0 + del_height*n_lines

		fig.update_layout(
					autosize=True,
					title=title,
					width=width,
					height=height,
					plot_bgcolor=plot_bgcolor,
					xaxis_title="cases",
					yaxis_title='pathways',
					showlegend=False,
					font=dict(
						family="Arial",
						size=14,
						color="Black"
					)
		)

		figname0 = title_replace(title)
		figname = osjoin(self.root_figure, figname0+'.html')

		fig.write_html(figname)
		figname = figname.replace('.html', '.png')
		fig.write_image(figname)
		if verbose: print(">>> HTML and png saved:", figname)		

		figname = osjoin(self.root_fig_md, figname0+'.png')
		fig.write_image(figname)
		self.curr_figname = figname

		self.dic_curr_circo = {}
		if sel_list != []:
			dfpiv_cases, figname_circos_list = \
				self.calc_mod_and_plot_circle_sel_list(sel_list, title_head=title_head,
							plot_bars=plot_bars, plot_lines=plot_lines,
							fontsize_circos=fontsize_circos, lfc_threshold=lfc_threshold, 
							space_factor=space_factor, with_pPMI_obs=with_pPMI_obs, 
							do_case_translate=do_case_translate, verbose=verbose)

			self.save_all_summaries_in_markdown(pathway_list=pathway_list, title_head=title_head,
							lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
							modulations=modulations, width=width, height=height, line_width=line_width,
							verbose=verbose)

			text_table_modulation = self.to_text_df_mod_inv_high()
			self.build_dic_curr_circo(figname_circos_list, sel_list)

		else:
			dfpiv_cases, figname_circos_list = pd.DataFrame(), []
			self.dfpiv_mod, self.dfpiv_inv = pd.DataFrame(), pd.DataFrame()
			self.dfpiv_high, self.dfpiv_high_inter = pd.DataFrame(), pd.DataFrame()
			self.dfpiv_high_female, self.dfpiv_high_male = pd.DataFrame(), pd.DataFrame()
			text_table_modulation = ''

		self.dfpiv_cases = dfpiv_cases

		return fig, dfpiv2, text_all, text_table_modulation, dfpiv_cases, figname_circos_list

	# self.build_modulation_list(mod_dw_list, lfc_dw_list)
	def build_modulation_list(self, mod_symb_list:List, lfc_list:List, max_regs:int=30) -> List:
		if len(mod_symb_list) != len(lfc_list):
			print("Error: len(mod_symb_list) !=  len(lfc_list)")
			return [", ".join(mod_symb_list) + ' LFC: ' + ",".join(lfc_list)]
		

		if len(mod_symb_list) > max_regs:
			dfa = pd.DataFrame({'symb': mod_symb_list, 'lfc': lfc_list})

			dfa['abs_lfc'] = np.abs(dfa['lfc'])
			dfa = dfa.sort_values('abs_lfc', ascending=False)
			dfa = dfa.head(max_regs)
			dfa = dfa.sort_values('symb', ascending=True)

			mod_symb_list = dfa.symb.to_list()
			lfc_list =  dfa.lfc.to_list()
		
		lista=[]
		for i in range(len(mod_symb_list)):
			symb = mod_symb_list[i]
			lfc  = lfc_list[i]
			lista.append(f"{symb} ({lfc:.2f})")

		return lista
	

	def study_name(self, title_head:str) -> Tuple[str, str]:
		fname = f"study_{title_head}.pdf"
		fname = title_replace(fname)

		filename_pdf = osjoin(self.root_pdf_summ, fname)

		return fname, filename_pdf


	def modulation_name(self, title:str) -> Tuple[str, str]:
		fname = title + ".pdf"
		fname = title_replace(fname)
		filename_pdf = osjoin(self.root_pdf_summ, fname)

		return fname, filename_pdf

	
	def write_pdf_simple(self, title:str, s_text:str, toc_level:int=5, verbose:bool=False):

		pdf = MarkdownPdf(toc_level=toc_level, optimize=True)
		
		_, filename_pdf = self.modulation_name(title) 


		try:
			pdf.add_section( Section(s_text), user_css=self.s_css)
			pdf.save(filename_pdf)
			if verbose: print(f"File saved: '{filename_pdf}'")
			ret = True

		except:
			print(f"Error: could not save: '{filename_pdf}'")
			ret = False

		return ret



	def save_degs_modulation_table(self, gem, title_head:str, chosen_model:int, context:str,
						temperature:float=0.1, topK:int=50, topP:float=0.10, maxOutputTokens:int=-1,
						stopSequences:List=[], start_dummy_word:str='Of course. ',
						force_ia:bool=False, force:bool=False,
						toc_level:int=5, prompt:bool=False, verbose:bool=False) -> bool:
		
		title_analysis = f"{self.gene_protein} modulations table for {title_head}"
		fname = title_analysis + '.tsv'
		fname = title_replace(fname)

		filefull = osjoin(self.root_pdf_summ, fname)
		if exists(filefull) and not force:
			if verbose: print(f"\t\tAlready calculated {fname}.")
			return True

		ret = self.calc_degs_modulation_table_using_AI(title_analysis=title_analysis, fname=fname,
					gem=gem, chosen_model=chosen_model, context=context,
					temperature=temperature, topK=topK, topP=topP, maxOutputTokens=maxOutputTokens,
					stopSequences=stopSequences, start_dummy_word=start_dummy_word,
					force_ia=force_ia, toc_level=toc_level, prompt=prompt, verbose=verbose)

		return ret	



	def run_ia_analysis(self, gem, title_head:str, chosen_model:int, 
						prefix_question:str, text_all:str, text_table_modulation:str, context:str,
						temperature:float=0.1, topK:int=50, topP:float=0.10, maxOutputTokens:int=-1,
						stopSequences:List=[], 
						force:bool=False, verbose:bool=False) -> Tuple[str, str]:
		
		
		ia_analysis, ia_model = \
		gem.large_text_pathways_analysis(head=title_head, root_pathway=self.root_pathway, 
										chosen_model=chosen_model, prefix_question=prefix_question, 
										text_all=text_all, text_table_modulation=text_table_modulation, context=context, 
										temperature=temperature, topK=topK, topP=topP, maxOutputTokens=maxOutputTokens,
										stopSequences=stopSequences,  force=force, verbose=verbose)

		self.ia_analysis = ia_analysis
		self.ia_model = ia_model
		self.chosen_model = chosen_model

		if ia_analysis is None or ia_analysis == '' or ia_analysis.startswith('Error') or ia_analysis == 'None':
			print(f"Error: ia_analysis '{ia_analysis}'")
			return 'Error', ia_model

		return ia_analysis, ia_model	



	def calc_degs_modulation_table_using_AI(self, title_analysis:str, fname:str, 
						gem, chosen_model:int, context:str,
						temperature:float=0.1, topK:int=50, topP:float=0.10, maxOutputTokens:int=-1,
						stopSequences:List=[], start_dummy_word:str='Of course. ',
						force_ia:bool=False, toc_level:int=5, prompt:bool=False, verbose:bool=False) -> bool:
		
		if self.dfpiv_mod is None or self.dfpiv_mod.empty:
			print("Error: did not calculate self.dfpiv_mod")
			return False
	
		if prompt: print(f"\t\tCalculating {self.s_deg_dap} modulations ...")

		dfpiv_mod2 = self.dfpiv_mod.reset_index().copy()

		lista = self.dfpiv_high_female.index.to_list()
		dfpiv_mod2 = dfpiv_mod2[~dfpiv_mod2.symbol.isin(lista)]

		lista = self.dfpiv_high_male.index.to_list()
		dfpiv_mod2 = dfpiv_mod2[~dfpiv_mod2.symbol.isin(lista)]

		lista_mod = dfpiv_mod2.symbol.to_list()
		lista_inv = self.dfpiv_inv.index.to_list()

		if self.has_gender:
			lista_high_inter = self.dfpiv_high_inter.index.to_list()
			lista_high_inter = [x for x in lista_high_inter if x not in lista_inv]
			lista_high_female = self.dfpiv_high_female.index.to_list()
			lista_high_male = self.dfpiv_high_male.index.to_list()
		else:
			lista_high_inter, lista_high_female, lista_high_male = [], [], []

		text_mod  = f"There are {len(lista_mod)} slightly modulated {self.s_deg_dap}s\n\n"
		text_mod += f"There are {len(lista_inv)} {self.s_deg_dap}s having LFC inversion between cases.\n\n"

		if self.has_gender:
			text_mod += f"There are {len(lista_high_inter)} {self.s_deg_dap}s having highly LFC differences between genders.\n\n"
			text_mod += f"There are {len(lista_high_female)} {self.s_deg_dap}s highly modulated LFC for females.\n\n"
			text_mod += f"There are {len(lista_high_male)} {self.s_deg_dap}s highly modulated LFC for males.\n\n"


		prefix_question = f"""
Analyze all the following list of {self.gene_protein}s. 
Search in the scientific literature, only for {self.disease}, for each sybmol, if:
  1) it was defined a biomarker.
  2) it was patented as a biomarker.
  3) it was seleted as target for therapeutics.
  4) it was patented as target for therapeutics.
"""

		lista_all = lista_mod + lista_inv

		if self.has_gender:
			lista_all += lista_high_inter + lista_high_female + lista_high_male
	
		lista_all = np.unique(lista_all)

		text_degs = f"""List of {self.gene_protein}s (#{len(lista_all)}): {"; ".join(lista_all)}."""

		text_mod += text_degs + '\n\n'
		if verbose: print(text_mod)

		prefix_question2 = prefix_question + """
Build a table markdown table, having as columns: each symbol and all 4 searchers.
Return comments and tables in markdown format. Do not embed a tsv file text.
"""

		if prompt: print(f"\t\tStarting IA analysis: model {chosen_model} ...")

		ia_analysis, ia_model = gem.large_text_pathways_analysis(head=title_analysis, root_pathway=self.root_pathway, 
											chosen_model=chosen_model, prefix_question=prefix_question2, 
											text_all=text_degs, text_table_modulation='', context=context, 
											temperature=temperature, topK=topK, topP=topP, 
											maxOutputTokens=maxOutputTokens, stopSequences=stopSequences,
											force=force_ia, verbose=verbose)
		
		
		if prompt: print(f"\t\tEnd IA analysis: model '{ia_model}'.")

		self.ia_analysis = ia_analysis
		new_text = self.replace_headers2(ia_analysis, start_dummy_word=start_dummy_word)
	
		title_analysis += f"_model {ia_model}"
		self.title_analysis = title_analysis

		header1 = f'# Modulations and {ia_model} analysis\n\n'
		s_text = header1 + text_mod + '\n' + new_text
		self.s_text = s_text

		ret = self.write_pdf_simple(title=title_analysis, s_text=s_text, toc_level=toc_level, verbose=verbose)

		if not ret:
			print(f"Error writing in {title_analysis} pdf.")
			return False
		
		if prompt: print("\t\tIA building table tsv ....")

		prefix_question2 = prefix_question + """\n
Build a table having as rows the symbols and as columns the 4 questions.
Answer the 4 questions with Yes or No.
The first column is the 'symbol' of the {self.gene_protein}.
The second column is 'selected as biomarker' as 'sel_biomarker'
The third column is 'patented as biomarker' as 'pat_biomarker'
The fourth column is 'selected for therapeutic' as 'sel_therapeutic'
The fifths (last) column is 'patentend for therapeutic' as 'pat_therapeutic'
Return a tsv file with respective header, separate char as '\t', and nothing more."""

		title_analysis_tsv = title_analysis + '_tsv_table'

		self.prefix_question2 = prefix_question2
		self.title_analysis_tsv = title_analysis_tsv
		self.text_degs = text_degs

		text_tsv, ia_model = gem.large_text_pathways_analysis(head=title_analysis_tsv, root_pathway=self.root_pathway, 
											chosen_model=chosen_model, prefix_question=prefix_question2, 
											text_all=text_degs, text_table_modulation='', context=context, 
											temperature=temperature, topK=topK, topP=topP, 
											maxOutputTokens=maxOutputTokens, stopSequences=stopSequences,
											force=force_ia, verbose=verbose)


		ret = write_txt(text_tsv, fname, self.root_pdf_summ)

		if not ret:
			print(f"Error writing in {fname}.")
			return False

		dfg = pdreadcsv(fname, self.root_pdf_summ)
		if dfg is None or dfg.empty:
			print(f"Error: could not read {fname}")
			return False
	
		cols = [x.lower().replace(' ', '_').replace('(','').replace(')','') for x in dfg.columns]
		cols[0] = 'symbol'
		dfg.columns = cols

		# no selection, patent neither therapeutics
		dfg['not_proposed'] = [True if dfg.iloc[i,1] == 'No' and dfg.iloc[i,2] == 'No' and
									   dfg.iloc[i,3] == 'No' and dfg.iloc[i,4] == 'No' else False for i in range(len(dfg))]

		# may be selected, but not  patent neither therapeutics
		dfg['weak_proposed'] = [True if dfg.iloc[i,1] != 'No' and dfg.iloc[i,2] == 'No' and
										dfg.iloc[i,3] == 'No' and dfg.iloc[i,4] == 'No' else False for i in range(len(dfg))]

		dfg['slightly_mod']	  = [True if x in lista_mod  else False for x in dfg.symbol]

		if self.has_gender:
			dfg['highly_inter_mod']  = [True if x in self.dfpiv_high_inter.index.to_list()  else False for x in dfg.symbol]
			dfg['highly_female_mod'] = [True if x in self.dfpiv_high_female.index.to_list() else False for x in dfg.symbol]
			dfg['highly_male_mod']   = [True if x in self.dfpiv_high_male.index.to_list()   else False for x in dfg.symbol]

		dfa = self.dfpiv_cases.reset_index()
		dfn = pd.merge(dfg, dfa, how="inner", on='symbol')

		if self.has_gender:
			seq = np.arange(0, len(self.group_female_list))

			dfa_female = dfn[self.group_female_list]
			dfa_male   = dfn[self.group_male_list]

			dfn['spearman_corr_female'] = [spearmanr(list(dfa_female.iloc[i]), seq)[0] for i in range(len(dfa_female))]
			dfn['spearman_corr_male']   = [spearmanr(list(  dfa_male.iloc[i]), seq)[0] for i in range(len(dfa_male))]

			dfn['corr_female_sig'] = dfn.spearman_corr_female.abs() >= self.min_corr_sig
			dfn['corr_male_sig']   = dfn.spearman_corr_male.abs()   >= self.min_corr_sig
			dfn['corr_both_sig']   = [ (dfn.iloc[i].corr_female_sig == True) and (dfn.iloc[i].corr_male_sig == True) for i in range(len(dfn))]

			lista_ini = [cols[0]] + ['not_proposed', 'weak_proposed'] + cols[1:9]
							
			lista_cor = ['corr_both_sig',  'corr_female_sig', 'corr_male_sig',  'spearman_corr_female', 'spearman_corr_male']
			cols2 = lista_ini + lista_cor + self.group_female_list + self.group_male_list
			dfn = dfn[cols2]
		else:
			seq = np.arange(0, len(self.case_list))

			dfa = dfn[self.case_list]

			dfn['spearman_corr'] = [spearmanr(list(dfa.iloc[i]), seq)[0] for i in range(len(dfa))]
			dfn['corr_sig']	  = dfn.spearman_corr.abs() >= self.min_corr_sig

			lista_ini = [cols[0]] + ['not_proposed', 'weak_proposed'] + cols[1:5]
							
			lista_cor = ['corr_sig', 'spearman_corr']
			cols2 = lista_ini + lista_cor + self.case_list
			dfn = dfn[cols2]

		dfn = dfn.sort_values(['not_proposed', 'weak_proposed', 'symbol'], ascending=[False, False, True])

		ret = pdwritecsv(dfn, fname, self.root_pdf_summ)
		if ret:
			if prompt: print(f"\t\tSaved table '{fname}' ....")
		else:
			print(f"\t\tError: saving table '{fname}' ....")

		return ret


	def plot_pathway_heatmap_by_gender(self, pathway_list:List, sel_list:List, title_head:str, maxLen:int=50,
				with_pPMI_obs:bool=False, do_case_translate:bool=True,
				width:int=1100, height0:int=350, del_height:int=75, line_width:int=2, 
				plot_bgcolor:str='lightgray', vertical_spacing=0.15, 
				lfc_threshold:float=1, corr_threshold:float=0.6,
				modulations:List=[ ('all', 'up'),('all', 'dw'),('all', 'none'), ('corr', 'up'), ('corr', 'dw'), ('corr', 'none') ],
				plot_bars:bool=True, plot_lines:bool=False,
				fontsize_circos:int=14, space_factor:float=3.,
				verbose:bool=False) -> Tuple[object, pd.DataFrame, pd.DataFrame, str, str, pd.DataFrame, list]:
		
		if not isinstance(pathway_list, list) or pathway_list == []:
			print("Pathway list is empty or not a list.")
			return None, pd.DataFrame(), pd.DataFrame(), '', '', pd.DataFrame(), []
	
		if len(pathway_list) > self.max_pathways:
			print(f"Warning: it is not recomended more than {self.max_pathways}, you sent {len(pathway_list)}")
			return None, pd.DataFrame(), pd.DataFrame(), '', '', pd.DataFrame(), []

		dfpiv, dff, text_all, _ = \
			self.calc_pPMI_summary_and_pivot_tables(pathway_list=pathway_list, title_head=title_head, with_pPMI_obs=with_pPMI_obs, 
													do_case_translate=do_case_translate, verbose=verbose)

		self.dff = dff

		maxi = int(np.max(np.ceil( dff.pPMI_total.to_list())))
		mini = int(np.min(np.floor(dff.pPMI_total.to_list())))
		if mini < 0:
			mini = -mini
		zlim = maxi if maxi > mini else mini

		zmin = -zlim
		zmax = +zlim

		not_in = [x for x in pathway_list if x not in self.selected_pivot_pathway_list]

		if len(not_in) > 0:
			print("Some pathways look wrong, please fix the pathway list.")
			print(">>>", "; ".join(not_in))
			return None, pd.DataFrame(), pd.DataFrame(), '', '', pd.DataFrame(), []

		gender_list = ['female', 'male', 'male-female']
		fig = make_subplots(rows=3, cols=1, subplot_titles=gender_list, vertical_spacing=vertical_spacing)

		dfpiv_female = dfpiv.loc[pathway_list, self.group_female_list].copy()
		lista = [str(x) for x in dfpiv_female.index.to_list()]
		new_indexes = [break_line_per_length(x, sep='<br>', maxLen=maxLen) for x in lista]
		dfpiv_female = dfpiv_female.reset_index(drop=True)
		dfpiv_female = dfpiv_female.set_index(pd.Series(new_indexes))
		cols_female = list(dfpiv_female.columns)
		cols = [x.replace('_female','') for x in dfpiv_female.columns.tolist()]
		dfpiv_female.columns = cols

		dfpiv_male = dfpiv.loc[pathway_list, self.group_male_list].copy()
		lista = [str(x) for x in dfpiv_male.index.tolist()]
		new_indexes = [break_line_per_length(x, sep='<br>', maxLen=maxLen) for x in lista]
		dfpiv_male = dfpiv_male.reset_index(drop=True)
		dfpiv_male = dfpiv_male.set_index(pd.Series(new_indexes))
		cols_male = list(dfpiv_male.columns)
		dfpiv_male.columns = cols

		dfpiv3 = pd.DataFrame()
		nrows = 0

		for gender in gender_list:
			if gender == 'female':
				# cols = self.group_female_list
				dfpiv3 = dfpiv_female
				nrow = 1

			elif gender == 'male':
				# cols = self.group_male_list
				dfpiv3 = dfpiv_male
				nrow = 2
			else:
				dfpiv3 = dfpiv_male.infer_objects().fillna(0) - dfpiv_female.infer_objects().fillna(0)
				nrow = 3

			if gender == 'male-female':
				text_list=None
			else:
				text_list=[]; lista_ratio=[]
				for i in range(len(dfpiv3)):
					row = dfpiv3.iloc[i]
					pathway = row.name

					lista =[]; lista2=[]
					for j in range(len(row.index)):
						# col = row.index[j]
						# pPMI = row[col]
						if gender == 'male':
							case = cols_male[j]
						else:
							case = cols_female[j]

						row_dff = dff[ (dff.case == case) & (dff.pathway == pathway)]
						if row_dff.empty:
							# print(i, 'empty', pathway, case, row.index)
							lista.append(None)
							# lista.append(f'<br>x:{case}<br>y:{pathway}<br>pPMI:{pPMI:.3f}')
						else:
							# print(i, 'Ok', pathway, case, row.index)
							row_dff = row_dff.iloc[0]
							mod_up_list = row_dff['mod_up_in_pathway']
							lfc_up_list = row_dff['lfc_up']

							mod_up_list = mod_up_list if isinstance(mod_up_list, list) else eval(mod_up_list)
							lfc_up_list = lfc_up_list if isinstance(lfc_up_list, list) else eval(lfc_up_list)

							# up_mod_list = self.build_modulation_list(mod_up_list, lfc_up_list)

							mod_dw_list = row_dff['mod_dw_in_pathway']
							lfc_dw_list = row_dff['lfc_dw']

							mod_dw_list = mod_dw_list if isinstance(mod_dw_list, list) else eval(mod_dw_list)
							lfc_dw_list = lfc_dw_list if isinstance(lfc_dw_list, list) else eval(lfc_dw_list)

							n_mod_in_pathway	= row_dff['n_mod_in_pathway']
							n_mod_up_in_pathway = row_dff['n_mod_up_in_pathway']
							n_mod_dw_in_pathway = row_dff['n_mod_dw_in_pathway']

							if mod_up_list == []:
								text_up_in_pathway = '---'
							else:
								mat_modulation = self.build_modulation_list(mod_up_list, lfc_up_list)
								text_up_in_pathway = break_list(mat_modulation, n_elems=5, CR='<br>')


							mod_dw_in_pathway = row_dff['mod_dw_in_pathway']
							if mod_dw_in_pathway == '[]':
								text_dw_in_pathway = '---'
							else:
								mat_modulation = self.build_modulation_list(mod_dw_list, lfc_dw_list)
								text_dw_in_pathway = break_list(mat_modulation, n_elems=5, CR='<br>')

							s_ratio = f"<br>#Up={n_mod_up_in_pathway} #Dw={n_mod_dw_in_pathway} #Modulated={n_mod_in_pathway}"
							stri_up_dw = s_ratio+f'<br><br>Up({n_mod_up_in_pathway}): {text_up_in_pathway}<br><br>Dw({n_mod_dw_in_pathway}): {text_dw_in_pathway}'
							lista.append(stri_up_dw)
							lista2.append(stri_up_dw)

					text_list.append(lista)
					lista_ratio.append(lista2)

			df_vectors = self.df_to_plotly(dfpiv3)
			nrows += len(dfpiv3)

			fig.add_trace(go.Heatmap(df_vectors, zmin=zmin, zmax=zmax, text=text_list, 
									 colorscale='RdBu_r', name=gender), row=nrow, col=1)

		stri = 'normalized' if self.pPMI_normalized else 'not normalized'
		title = f"pPMI heatmap for '{title_head}' - " + stri

		n_lines = nrows
		n_lines -= 9
		if n_lines < 0: n_lines=0

		height = height0 + del_height*n_lines

		fig.update_layout(
					autosize=True,
					title=title,
					width=width,
					height=height,
					plot_bgcolor=plot_bgcolor,
					xaxis3_title="cases",
					yaxis_title='pathways',
					showlegend=True,
					font=dict(
						family="Arial",
						size=14,
						color="Black"
					)
		)

		figname0 = title_replace(title)
		figname = osjoin(self.root_figure, figname0+'.html')

		fig.write_html(figname)
		figname = figname.replace('.html', '.png')
		fig.write_image(figname)
		if verbose: print(">>> HTML and png saved:", figname)

		figname = osjoin(self.root_fig_md, figname0+'.png')
		fig.write_image(figname)
		self.curr_figname = figname
		
		if sel_list != []:
			dfpiv_cases, figname_circos_list = \
				self.calc_mod_and_plot_circle_sel_list(sel_list, title_head=title_head, 
							  		plot_bars=plot_bars, plot_lines=plot_lines,
									fontsize_circos=fontsize_circos, lfc_threshold=lfc_threshold, 
									space_factor=space_factor, with_pPMI_obs=with_pPMI_obs, 
									do_case_translate=do_case_translate, verbose=verbose)
			
			self.save_all_summaries_in_markdown(pathway_list=pathway_list, title_head=title_head, 
									lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
									modulations=modulations, width=width, height=height, line_width=line_width,
									verbose=verbose)
			
			text_table_modulation = self.to_text_df_mod_inv_high()
			self.build_dic_curr_circo(figname_circos_list, sel_list)
		else:
			dfpiv_cases, figname_circos_list = pd.DataFrame(), []
			self.dfpiv_mod, self.dfpiv_inv = pd.DataFrame(), pd.DataFrame()
			self.dfpiv_high, self.dfpiv_high_inter = pd.DataFrame(), pd.DataFrame()
			self.dfpiv_high_female, self.dfpiv_high_male = pd.DataFrame(), pd.DataFrame()
			text_table_modulation = ''
			self.dic_curr_circo = {}

		self.dfpiv_cases = dfpiv_cases

		return fig, dfpiv_female, dfpiv_male, text_all, text_table_modulation, dfpiv_cases, figname_circos_list
	
	def build_dic_curr_circo(self, figname_circos_list:List, sel_list:List):
		
		self.dic_curr_circo = {}

		if self.has_gender:
			case_list = self.group_female_list + self.group_male_list
		else:
			case_list = self.case_list
		
		for ifig, filename in enumerate(figname_circos_list):
			case = case_list[ifig]
			self.dic_curr_circo[case] = {}
			self.dic_curr_circo[case]['filename'] = filename
			self.dic_curr_circo[case]['sel_list']  = sel_list
			

	def to_text_df_mod_inv_high(self) -> str:

		text = '## Modulations and Inversions\n'

		if self.has_gender:
			lista = [ ('female', self.group_female_list), ('male', self.group_male_list) ]
		else:
			lista = [ ('', self.case_list) ]

		if self.dfpiv_mod is None or self.dfpiv_mod.empty:
			text += f'There are no {self.s_deg_dap}s with significative modulations between cases.\n'
		else:

			for term_gender, case_list in lista:

				if term_gender == '':
					text += f'### {self.s_deg_dap}s modulations between cases.\n'
				else:
					text += f'### {self.s_deg_dap}s modulations between {term_gender} cases.\n'

				cols = case_list
				rows = self.dfpiv_mod.index.to_list()
				row_name = self.s_deg_dap + 's'
				text += df_to_md_table(self.dfpiv_mod, cols, rows, row_name, nround=2)
				text += '\n'

		text += '\n'

		if self.dfpiv_inv is None or self.dfpiv_inv.empty:
			text += f'There are no {self.s_deg_dap}s with inverted modulations between cases.\n'
		else:
			for term_gender, case_list in lista:

				if term_gender == '':
					text += f'### {self.s_deg_dap}s inverted modulations between cases.\n'
				else:
					text += f'### {self.s_deg_dap}s inverted modulations between {term_gender} cases.\n'

				cols = case_list
				rows = self.dfpiv_inv.index.to_list()
				row_name = self.s_deg_dap + 's'
				text += df_to_md_table(self.dfpiv_inv, cols, rows, row_name, nround=2)
				text += '\n'
		text += '\n'

		#--------------------- highly modulated --------------------
		if self.has_gender:
			#-------------- between genders -------------------------------------------
			if self.dfpiv_high_inter is None or self.dfpiv_high_inter.empty:
				text += f'There are no {self.s_deg_dap}s highly modulated between genders.\n'
			else:
				text += f'### {self.s_deg_dap}s highly modulated for between genders.\n'

				cols = self.group_list
				rows = self.dfpiv_high_inter.index.to_list()
				row_name = self.s_deg_dap + 's'
				text += df_to_md_table(self.dfpiv_high_inter, cols, rows, row_name, nround=2)
				text += '\n'
			text += '\n'

			#-------------- females -----------------------------------------------------
			if self.dfpiv_high_female is None or self.dfpiv_high_female.empty:
				text += f'There are no {self.s_deg_dap}s highly modulated for females.\n'
			else:
				text += f'### {self.s_deg_dap}s highly modulated for females.\n'

				cols = self.group_female_list
				rows = self.dfpiv_high_female.index.to_list()
				row_name = self.s_deg_dap + 's'
				text += df_to_md_table(self.dfpiv_high_female, cols, rows, row_name, nround=2)
				text += '\n'
			text += '\n'

			#-------------- males -----------------------------------------------------
			if self.dfpiv_high_male is None or self.dfpiv_high_male.empty:
				text += f'There are no {self.s_deg_dap}s highly modulated for males.\n'
			else:
				text += f'### {self.s_deg_dap}s highly modulated for males.\n'

				cols = self.group_male_list
				rows = self.dfpiv_high_male.index.to_list()
				row_name = self.s_deg_dap + 's'
				text += df_to_md_table(self.dfpiv_high_male, cols, rows, row_name, nround=2)
				text += '\n'
			text += '\n'

		else:
			if self.dfpiv_high is None or self.dfpiv_high.empty:
				text += f'There are no {self.s_deg_dap}s highly modulated between cases.\n'
			else:
				text += f'### {self.s_deg_dap}s highly modulated between cases.\n'

				cols = self.case_list
				rows = self.dfpiv_high.index.to_list()
				row_name = self.s_deg_dap + 's'
				text += df_to_md_table(self.dfpiv_high, cols, rows, row_name, nround=2)
				text += '\n'
			text += '\n'

		return text


	def study_exists(self, title_head:str, verbose:bool=False) -> bool:
		_, filename_pdf = self.study_name(title_head)

		if exists(filename_pdf):
			if verbose: print(f"File already exists: {filename_pdf}")
			ret = True
		else:
			ret = False

		return ret

	def calc_LFC_df_all_pathways(self, pathway_list:List, verbose:bool=False) -> pd.DataFrame:

		df_list=[]
		for pathway in pathway_list:
			pathway_id, _ = self.reactome_find_genes_in_pathway(pathway, _type='pathway')
			pathway_and_id = [pathway, pathway_id]

			dfpiv, _ = self.calc_pivot_one_pathway_LFC(pathway_and_id, verbose=verbose)

			cols = list(dfpiv.columns)
			dfpiv = dfpiv.reset_index()
			cols = ['symbol'] + cols
			dfpiv.columns = cols
			
			df_list.append(dfpiv)

		dfpiv = pd.concat(df_list)
		dfpiv = dfpiv.drop_duplicates('symbol')

		return dfpiv

	def is_modulated(self, dfpiv:pd.DataFrame, lfc_threshold=1.0) -> pd.DataFrame:
		goods=[]
		for i in range(len(dfpiv)):
			goods.append( np.sum([True if abs(val) >= lfc_threshold else False for val in dfpiv.iloc[i] ] ) > 0 )

		return dfpiv[goods]

	def get_availabel_memory_MB(self):
		memory_info = psutil.virtual_memory()

		return round(memory_info.available / (1024**2), 2)

	def plot_all_lfc_path(self, df_anal, diff:float=1., corr_cutoff:float=0.6, 
					      margin:dict=dict(l=80, r=80, t=100, b=80), 
						  show_figure:bool=True) -> pd.DataFrame:

		all_dics = {}; icount=-1
		for (want_female, want_male) in [(True, False), (False, True), (True, True), (False, False)]:
			
			if want_female and want_male:
				similar_modulation_for_both_list = [True, False]
			else:
				similar_modulation_for_both_list = [False]

			if not want_female and not want_male:
				# do not see correlations
				high_correlation_list = [False]
			else:
				high_correlation_list = [True, False]
			
			for high_correlation in high_correlation_list:
				for similar_modulation_for_both in similar_modulation_for_both_list:

					dic_fig = self.plot_genes_lfc_path(df_anal, want_female=want_female, want_male=want_male, 
									 	diff=diff, corr_cutoff=corr_cutoff, 
										high_correlation=high_correlation, similar_modulation_for_both=similar_modulation_for_both,
										width=800, height=600, margin=margin, verbose=False)
					
					if dic_fig == {}:
						continue
	
					for _, dic2 in dic_fig.items():
						icount += 1
						all_dics[icount] = dic2

					if show_figure:
						show_title = True
						for _, dic in dic_fig.items():
							symbol = dic['symbol']
							fig = dic['fig']
							title = dic['title']
							if show_title:
								print(">>>", title)
								show_title = False
							print("##", symbol)
							fig.show()
							print("\n")

		if all_dics == {}:
			return pd.DataFrame()
	
		dffig = pd.DataFrame(all_dics).T
		dffig.reset_index(inplace=True, drop=True)
		return dffig

	def write_pdf_lfc_path(self, dffig:pd.DataFrame, title_head:str, toc_level:int=5, verbose:bool=False):

		fname = F'lfc_path_plots_for_{title_head}.pdf'
		fname = title_replace(fname)
		filename_pdf = osjoin(self.root_pdf_summ, fname)

		pdf = MarkdownPdf(toc_level=toc_level, optimize=True)

		s_md_plots = f"# LFC-path plot for {title_head}\n"

		cols = ['want_female', 'want_male', 'high_correlation','similar_modulation_for_both']
		dfa = dffig[cols].drop_duplicates()

		for i in range(len(dfa)):
			row = dfa.iloc[i]

			want_female = row['want_female']
			want_male = row['want_male']
			high_correlation = row['high_correlation']
			similar_modulation_for_both = row['similar_modulation_for_both']

			df2 = dffig[(dffig.want_female == want_female) & (dffig.want_male == want_male) &
			   			(dffig.high_correlation == high_correlation) & (dffig.similar_modulation_for_both == similar_modulation_for_both) ]
			
			if df2.empty:
				continue

			first = True
			for j in range(len(df2)):
				row2 = df2.iloc[j]
				symbol = row2.symbol
				title = row2.title

				if first:
					first = False

					if want_female and want_male:
						want_gender = 'female and male'
					elif not want_female and not want_male:
						want_gender = 'not female and not male'
					else:
						if want_female:
							want_gender = 'only female'
						else:
							want_gender = 'only male'

					if not want_female and not want_male:
						high_correlation2 = 'independ of correlations'
					else:
						high_correlation2 = 'high correlation' if high_correlation else 'low correlation'

					similar_modulation_for_both2 = ' and similar modulation' if similar_modulation_for_both else ''

					s_md_plots += f"## {want_gender} and {high_correlation2}{similar_modulation_for_both2}\n"

				figname0 = title_replace(title)
				if len(figname0) > 250:
					figname0 = figname0[:250]
				figname = figname0 + '.png'

				filename = self.root_fig_md / figname

				if filename.exists():
					s_md_plots += f"### {symbol}\n"
					s_md_plots += f"""![{symbol})]({filename})""" + '\n'
				else:
					print(f"Error: could not find {filename}")

				s_md_plots += '\n' + title + '\n'
				pdf.add_section( Section(s_md_plots), user_css=self.s_css)
				s_md_plots = ''

		try:
			pdf.save(filename_pdf)
			if verbose: print(f"File saved: '{filename_pdf}'")
			ret = True

		except:
			print(f"Error: could not save: '{filename_pdf}'")
			ret = False

		return ret


	def write_pdf_study(self, title_head:str, pathway_list:List, text_all:str, 
						text_ia:str, ia_model:str,_type_fig:str='png', toc_level:int=4,
						want_circle_plots:bool=True, want_line_plots:bool=False, 
						want_heatmap_plots:bool=True, want_expression_tables:bool=True,
						want_ia:bool=True, want_modulation_summary:bool=True,
						start_dummy_word:str='Of course. ',
						type_modulation:str='all', up_down_all:str='all', 
						lfc_threshold:float=1.0, diff_cutoff:float=1.0,
						force_summary:bool=False, prompt:bool=True, verbose:bool=False) -> bool:

		if len(pathway_list) > self.max_pathways:
			print(f"Warning: it is not recomended more than {self.max_pathways}, you sent {len(pathway_list)}")
			return False

		fname, filename_pdf = self.study_name(title_head)

		if exists(filename_pdf) and not force_summary:
			if verbose: print(f"File already exists: {filename_pdf}")
			return False
		
		dfpiv_mod = self.dfpiv_mod
		dfpiv_inv = self.dfpiv_inv
		dfpiv_high = self.dfpiv_high
		dfpiv_high_inter = self.dfpiv_high_inter
		dfpiv_high_female = self.dfpiv_high_female
		dfpiv_high_male = self.dfpiv_high_male

		pdf = MarkdownPdf(toc_level=toc_level, optimize=True)

		if prompt: print(f"Starting: memory {self.get_availabel_memory_MB()} MBytes.")

		s_md_head = f"# Study: {title_head}\n"
		s_md_heatmap = f"## Heatmap\n![Heatmap: {title_head}]({self.curr_figname})\n"
		s_md = s_md_head + '\n' + s_md_heatmap + '\n'
		pdf.add_section( Section(s_md), user_css=self.s_css)

		if want_circle_plots:
			s_md_circos = '## Circle-plot\n\n'
			if len(self.dic_curr_circo) > 0:
				if self.has_gender:
					case_list = self.group_female_list + self.group_male_list
				else:
					case_list = self.case_list
				
				for case in case_list:

					dic = self.dic_curr_circo[case]
					filename = dic['filename']
					sel_list = dic['sel_list']

					s_md_head = f"### {case}\n"
					s_md_in = f"![{case}]({filename})\n"

					for i_sector, pathway in enumerate(sel_list):
						s_md_in += f"\nSector {i_sector+1}: {pathway}\n"
			
					s_md_circos += s_md_head + s_md_in + '\n'
					pdf.add_section( Section(s_md_circos), user_css=self.s_css)
					s_md_circos = ''

			del(s_md_circos)
			if prompt: print(f"Finished Circos: memory {self.get_availabel_memory_MB()} MBytes.")


		DEG_DAPs = self.s_deg_dap+'s'

		if want_line_plots:
			s_md_plots = f"## Lineplots\n"

			term_corrs = ['_positive_correlation', '_negative_correlation', '_no_correlation', '_at_least_one_upregulated','_at_least_one_downregulated','_all_non-regulated']
			n_loops = len(pathway_list) * len(term_corrs)

			#------------------- n_loops: max is 5 pathways * 6 -----------------------------------------
			if n_loops > self.max_iloop_plot_lines:
				term_corrs = ['_positive_correlation', '_negative_correlation']

			i_loop = 0
			for pathway in pathway_list:
				pathway_und = title_replace(pathway)

				s_md_plots += f"### {pathway}\n"

				for term in term_corrs:
					term2 = term.replace('_',' ').strip()

					if prompt: print(">>>", pathway, ">>>", term)

					mat = [x for x in os.listdir(self.root_figure) if pathway_und in x and DEG_DAPs in x and _type_fig in x and term in x]

					fname = ''
					if len(mat) != 1:

						if len(mat) == 0:
							if verbose: print(f"Warning: searching for os.listdir() in {self.root_figure} pathway {pathway_und} having {DEG_DAPs} and {term} and {_type_fig}")
							continue
						
						#--- get the one with the shortest name ------
						L = 9999
						for fname2 in mat:
							if len(fname2) < L:
								fname = fname2
								L = len(fname2)
					else:
						fname = mat[0]

					file_ori  = osjoin(self.root_figure, fname)

					if not exists(file_ori):
						continue

					file_dest = osjoin(self.root_fig_md, fname)

					if exists(file_dest):
						s_md_plots += f"#### {term2}\n"
						s_md_plots += f"""![{pathway}]({file_dest})""" + '\n'
						pdf.add_section( Section(s_md_plots), user_css=self.s_css)

						i_loop += 1
						if i_loop >= self.max_iloop_plot_lines:
							print(f"Plot_lines overflow: max {self.max_iloop_plot_lines} ")
							break
					
					s_md_plots = ''

				if i_loop >= self.max_iloop_plot_lines:
					break

			if prompt: print(f"Finished Lineplots - i_loop {i_loop}: memory {self.get_availabel_memory_MB()} MBytes.")


		#------------------ heatmaps -------------------
		if want_heatmap_plots:
			s_md_plots = f"## {self.s_deg_dap}s heatmap\n"

			'''
			it saves image in the ./figures - auxiliar dir for markdown

			_ = self.plot_set_pathways_genes_heatmap(set_pathway=set_pathway, 
										type_modulation='all', up_down_all='all',
										lfc_threshold=lfc_threshold, diff_cutoff=diff_cutoff,
										width=1100, row_height=30, header_height=200,
										plot_bgcolor='lightgray', verbose=False)
			'''

			print_header3 = True
			s_filtered = self.set_s_filtered(type_modulation, up_down_all, diff_cutoff, lfc_threshold)

			for pathway in pathway_list:
				pathway_id, _ = self.reactome_find_genes_in_pathway(pathway, _type='pathway')

				for unique_gender in self.unique_gender_list:

					title = self.get_title_heatmap_set_pathway(unique_gender, pathway_id, pathway,
															   title_head,  s_filtered)
					figname0 = title_replace(title)
					if len(figname0) > 250:
						figname0 = figname0[:250]
					figname = figname0 + '.png'

					filename = osjoin(self.root_fig_md, figname)

					if exists(filename):
						if print_header3:
							s_md_plots += f"### {title_head}\n"
							print_header3 = False
						
						if unique_gender == 'unique':
							s_md_plots += f"#### {pathway} ({pathway_id})\n"
						else:
							s_md_plots += f"#### {unique_gender}: {pathway} ({pathway_id})\n"
							
						s_md_plots += f"""![{pathway})]({filename})""" + '\n'
					else:
						print(f"Error: could not find {filename}")
				
				pdf.add_section( Section(s_md_plots), user_css=self.s_css)
				s_md_plots = ''

			if prompt: print(f"Finished Heatmap: memory {self.get_availabel_memory_MB()} MBytes.")

		#-------------------------------- Summary -----------------------------------------
		s_md_text = f"## Summary (text)\n\n{text_all}\n"
		pdf.add_section( Section(s_md_text), user_css=self.s_css)


		#-------------------------------- datatable -----------------------------------------
		if want_expression_tables:
			s_md_mod = '## Modulations and Inversions\n\n'

			#---------- Modulations --------------
			row_name = 'symbol'

			if not dfpiv_mod.empty:
				if self.has_gender:
					#------------- female --------------
					self.dfpiv_mod = dfpiv_mod
					dfpiv_female = dfpiv_mod.loc[:, self.group_female_list]
					self.dfpiv_female = dfpiv_female
					cols = self.group_female_list
					# cols = [x.replace('_female', '') for x in cols]

					rows = dfpiv_female.index.to_list()
					s_markdown_female = df_to_md_table(dfpiv_female, cols, rows, row_name, nround=2)

					#------------- male --------------
					dfpiv_male = dfpiv_mod.loc[:, self.group_male_list]
					cols = self.group_male_list
					# cols = [x.replace('_male', '') for x in cols]

					rows = dfpiv_male.index.to_list()
					s_markdown_male = df_to_md_table(dfpiv_male, cols, rows, row_name, nround=2)

					s_md_mod += f"### Female modulated {self.s_deg_dap}s LFC table\n\n{s_markdown_female}\n"

					s_md_mod += f"### Male modulated {self.s_deg_dap}s LFC table\n\n{s_markdown_male}\n"
				else:
					cols = self.case_list

					rows = dfpiv_mod.index.to_list()
					s_markdown = df_to_md_table(dfpiv_mod, cols, rows, row_name, nround=2)

					s_md_mod += f"### Modulated {self.s_deg_dap}s LFC table\n\n{s_markdown}\n"
			else:
					s_md_mod += f"There are no significative modulated {self.s_deg_dap}s\n"


			s_md_mod += f'\nModulation: abs(lfc) >= {self.abs_mod_diff_cutoff}.\n'
			pdf.add_section( Section(s_md_mod), user_css=self.s_css)
			s_md_mod = ''

			#---------- Inversions --------------
			if not dfpiv_inv.empty:
				if self.has_gender:
					#------------- female --------------
					dfpiv_female = dfpiv_inv.loc[:, self.group_female_list]
					cols = self.group_female_list
					# cols = [x.replace('_female', '') for x in cols]

					rows = dfpiv_female.index.to_list()
					s_markdown_female = df_to_md_table(dfpiv_female, cols, rows, row_name, nround=2)

					#------------- male --------------
					dfpiv_male = dfpiv_inv.loc[:, self.group_male_list]
					cols = self.group_male_list
					# cols = [x.replace('_male', '') for x in cols]

					rows = dfpiv_male.index.to_list()
					s_markdown_male = df_to_md_table(dfpiv_male, cols, rows, row_name, nround=2)

					s_md_mod += f"### Female inversions {self.s_deg_dap}s LFC table\n\n{s_markdown_female}\n"

					s_md_mod += f"### Male inversions {self.s_deg_dap}s LFC table\n\n{s_markdown_male}\n"

				else:
					cols = self.case_list

					rows = dfpiv_inv.index.to_list()
					s_markdown = df_to_md_table(dfpiv_inv, cols, rows, row_name, nround=2)

					s_md_mod += f"### Inversions {self.s_deg_dap}s LFC table\n\n{s_markdown}\n"
			else:
					s_md_mod += f"There are no significative inversions {self.s_deg_dap}s\n"

			if self.has_gender:
				s_md_mod += '\nInversion: at least one inversion in one of the genders.\n'

			pdf.add_section( Section(s_md_mod), user_css=self.s_css)
			s_md_mod = ''

			#---------------------- highly modulated -------------------

			if self.has_gender:
				#----------- dfpiv_high_inter, dfpiv_high_female, dfpiv_high_male --------------
				s_md_mod += '\n'
				if not dfpiv_high_inter.empty:
					#------------- female --------------
					cols = self.group_list
					rows = dfpiv_high_inter.index.to_list()
					s_markdown_high_inter = df_to_md_table(dfpiv_high_inter, cols, rows, row_name, nround=2)

					s_md_mod += f"### Highly modulated between genders: female-male {self.s_deg_dap}s LFC table\n\n{s_markdown_high_inter}\n"
				else:
					s_md_mod += f"There are no highly modulated {self.s_deg_dap}s between genders.\n"

				s_md_mod += '\n'
				if not dfpiv_high_female.empty:
					#------------- female --------------
					cols = self.group_female_list
					rows = dfpiv_high_female.index.to_list()
					s_markdown_high_female = df_to_md_table(dfpiv_high_female, cols, rows, row_name, nround=2)

					s_md_mod += f"### Female highly modulated: {self.s_deg_dap}s LFC table\n\n{s_markdown_high_female}\n"
				else:
					s_md_mod += f"There are no highly modulated {self.s_deg_dap}s for female.\n"


				s_md_mod += '\n'
				if not dfpiv_high_male.empty:
					#------------- female --------------
					cols = self.group_male_list
					rows = dfpiv_high_male.index.to_list()
					s_markdown_high_male = df_to_md_table(dfpiv_high_male, cols, rows, row_name, nround=2)

					s_md_mod += f"### Male highly modulated: {self.s_deg_dap}s LFC table\n\n{s_markdown_high_male}\n"
				else:
					s_md_mod += f"There are no highly modulated {self.s_deg_dap}s for male.\n"


			else:
				s_md_mod += '\n'
				if not dfpiv_high.empty:
					#------------- female --------------
					cols = self.case_list
					rows = dfpiv_high.index.to_list()
					s_markdown_high = df_to_md_table(dfpiv_high, cols, rows, row_name, nround=2)

					s_md_mod += f"### Highly modulated {self.s_deg_dap}s LFC table\n\n{s_markdown_high}\n"
				else:
					s_md_mod += f"There are no highly modulated {self.s_deg_dap}s.\n"

			pdf.add_section( Section(s_md_mod), user_css=self.s_css)
			del(s_md_mod)

			
		if want_ia:
			new_text = self.replace_headers2(text_ia, start_dummy_word=start_dummy_word)
			s_md_text = f"# {ia_model} analysis\n\n{new_text}\n"

			pdf.add_section( Section(s_md_text), user_css=self.s_css)

			if prompt: print(f"Finished IA1: memory {self.get_availabel_memory_MB()} MBytes.")
			del(s_md_text)
			if prompt: print(f"Finished IA2: memory {self.get_availabel_memory_MB()} MBytes.")

		if want_modulation_summary:
			set_pathway = [title_head, pathway_list, []]

			s_md_text = self.one_pathway_LFC_modulation(set_pathway, lfc_threshold=lfc_threshold, diff_cutoff=diff_cutoff)
			pdf.add_section( Section(s_md_text), user_css=self.s_css)

			# del(s_md_text)
		
		try:
			pdf.save(filename_pdf)
			if verbose or prompt: print(f"PDF saved at '{filename_pdf}'")
			ret = True
		except Exception as e:
			ret = False
			print(f"An unexpected error occurred saving '{filename_pdf}': '{e}'")
		finally:
			del(pdf)

		return ret

	def replace_headers2(self, text_ia:str, start_dummy_word:str='Of course. ') -> str:
		'''
			remove dummy comments
			if theres is no header 1 and 2
				replace all headers to ##, ###, ####
		'''
		ret1 = text_starts_with_word(text_ia, word='# ')
		ret2 = text_starts_with_word(text_ia, word='## ')

		if start_dummy_word is not None and start_dummy_word != '':
			new_text = re.sub(start_dummy_word, '', text_ia)
		else:
			new_text = text_ia
		
		if not ret1 and not ret2:
			new_text = re.sub(r'###', r'##', new_text)
			new_text = re.sub(r'####', r'###', new_text)
			new_text = re.sub(r'#####', r'####', new_text)
			new_text = re.sub(r'######', r'#####', new_text)

		return new_text

	def which_reacotme_id(self, pathway):
		if self.dfr is None:
			dfr =  self.reactome.open_reactome(verbose=False)
			if dfr is None or dfr.empty:
				return None
			self.dfr = dfr

		dfq = self.dfr[self.dfr.pathway == pathway]
		if dfq.empty:
			return None
		return dfq.iloc[0].pathway_id

	def rebuild_column(self, group_age, gender):
		mat = group_age.split('_')

		if len(mat) == 1:
			return group_age + '_' + gender

		return mat[0] + '_' + gender + '_' + mat[1]


	def plot_all_pathways_LFC_modulations(self, pathway_list:List, lfc_threshold:float=1, corr_threshold:float=0.6,
						   modulations:List=[ ('all', 'up'),('all', 'dw'),('all', 'none'), ('corr', 'up'), ('corr', 'dw'), ('corr', 'none') ],
						   width:int=1000, height:int=700, line_width:int=2, 
						   fontsize:int=16, leg_fontsize:int=16, 
						   tick_fontsize:int=16, title_fontsize:int=18,
						   plot_bgcolor:str='lightgray', vertical_spacing:float=0.20, verbose:bool=False) -> List:

		figname_list = []
		
		for pathway in pathway_list:
			pathway_id, _ = self.reactome_find_genes_in_pathway(pathway, _type='pathway')
			
			pathway_and_id = [pathway, pathway_id]
			
			for (type_modulation, up_down_all) in modulations:

				_, _, _, figname0 = \
					self.lineplot_one_pathway_LFCs(pathway_and_id, type_modulation=type_modulation, up_down_all=up_down_all,
												lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
												width=width, height=height, line_width=line_width,
												fontsize=fontsize, leg_fontsize=leg_fontsize, 
												tick_fontsize=tick_fontsize, title_fontsize=title_fontsize,
												plot_bgcolor=plot_bgcolor, vertical_spacing=vertical_spacing, verbose=verbose)
				
				if figname0 != '': figname_list.append(figname0)

		return figname_list
				

	def plot_one_pathway_LFC_modulations(self, pathway:str, lfc_threshold:float=1, corr_threshold:float=0.6,
						   modulations:List=[ ('all', 'up'),('all', 'dw'),('all', 'none'), ('corr', 'up'), ('corr', 'dw'), ('corr', 'none') ],
						   width:int=1000, height:int=700, line_width:int=2, 
						   fontsize:int=16, leg_fontsize:int=16, 
						   tick_fontsize:int=16, title_fontsize:int=18,
						   plot_bgcolor:str='lightgray', vertical_spacing:float=0.20, verbose:bool=False) -> List:

		pathway_id, _ = self.reactome_find_genes_in_pathway(pathway, _type='pathway')
		
		pathway_and_id = [pathway, pathway_id]
		if verbose: print(pathway, pathway_id)
		
		figname_list = []
		for (type_modulation, up_down_all) in modulations:

			_, _, _, figname = \
				self.lineplot_one_pathway_LFCs(pathway_and_id, type_modulation=type_modulation, up_down_all=up_down_all,
											   lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
											   width=width, height=height, line_width=line_width,
											   plot_bgcolor=plot_bgcolor, vertical_spacing=vertical_spacing, verbose=verbose)
			if figname != '': figname_list.append(figname)

		return figname_list


	def lineplot_one_pathway_LFCs(self, pathway_and_id:List, 
								  type_modulation:str='all', up_down_all:str='up', 
								  lfc_threshold:float=1., corr_threshold:float=0.4,
								  width:int=1100, height:int=700, line_width:int=3,
								  fontsize:int=16, leg_fontsize:int=16, 
								  tick_fontsize:int=16, title_fontsize:int=18,
								  plot_bgcolor:str='lightgray',  vertical_spacing:float=0.20,
								  verbose:bool=False) -> Tuple[object, pd.DataFrame, pd.DataFrame, str]:
		
		if self.has_gender:
			fig, dff, dfpiv, figname0 = \
			self.lineplot_one_pathway_LFCs_by_gender(pathway_and_id, type_modulation=type_modulation, 
								up_down_all=up_down_all, lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
								width=width, height=height, line_width=line_width,
								fontsize=16, leg_fontsize=16, tick_fontsize=16, title_fontsize=18,
								plot_bgcolor=plot_bgcolor, vertical_spacing=vertical_spacing, verbose=verbose)
		else:
			fig, dff, dfpiv, figname0 = \
			self.lineplot_one_pathway_LFCs_wo_gender(pathway_and_id, type_modulation=type_modulation, 
								up_down_all=up_down_all, lfc_threshold=lfc_threshold, corr_threshold=corr_threshold,
								width=width, height=height, line_width=line_width,
								fontsize=16, leg_fontsize=16, tick_fontsize=16, title_fontsize=18,
								plot_bgcolor=plot_bgcolor, vertical_spacing=vertical_spacing, verbose=verbose)
			
		return fig, dff, dfpiv, figname0


	def lineplot_one_pathway_LFCs_by_gender(self, pathway_and_id:List, 
							type_modulation:str='all', up_down_all:str='up', 
							lfc_threshold:float=1., corr_threshold:float=0.4,
							width:int=1100, height:int=700, line_width:int=3,
							fontsize:int=16, leg_fontsize:int=16, 
							tick_fontsize:int=16, title_fontsize:int=18,
							plot_bgcolor:str='lightgray',  vertical_spacing:float=0.20,
							verbose:bool=False) -> Tuple[object, pd.DataFrame, pd.DataFrame, str]:
		'''
			LFC modulation, line plot
		'''

		if not isinstance(pathway_and_id, list) or pathway_and_id == []:
			print("Pathway and ID is an empty list.")
			return None, pd.DataFrame(), pd.DataFrame(), ''

		pathway, _ = pathway_and_id

		dfpiv, dff = self.calc_pivot_one_pathway_LFC(pathway_and_id, verbose=verbose)
		found = False

		symbols = list(dfpiv.index)
		
		dfp = pd.DataFrame()
		gender_list = []; gender=''; correlations=[]
		fig = go.Figure()

		def calc_all_modulation(dfp:pd.DataFrame, up_down_all:str, lfc_threshold:float):
			if up_down_all == 'up':
				dfp = dfp[ [ np.sum( dfp.iloc[i].apply(lambda x: 1 if x >= lfc_threshold else 0) ) > 0  for i in range(len(dfp))] ]
				s_filtered = f"at least one upregulated: LFC >= {lfc_threshold}"
	
			elif up_down_all == 'dw':
				dfp = dfp[ [ np.sum( dfp.iloc[i].apply(lambda x: 1 if x <= -lfc_threshold else 0) ) > 0  for i in range(len(dfp))] ]
				s_filtered = f"at least one downregulated: LFC <= {-lfc_threshold}"
	
			else:
				dfp = dfp[ [ (np.sum( dfp.iloc[i].apply(lambda x: 1 if x >= lfc_threshold else 0) ) +
							np.sum( dfp.iloc[i].apply(lambda x: 1 if x <= -lfc_threshold else 0) ) ) == 0  for i in range(len(dfp))] ]
	
				s_filtered = f"all non-regulated: abs any LFC < {lfc_threshold}"

			return dfp, s_filtered			


		def calc_correlations(dfp:pd.DataFrame, up_down_all:str, seqx:List, corr_threshold:float):

			corr_list = [  spearmanr(seqx, dfp.iloc[i,:])[0] for i in range(len(dfp))]

			if up_down_all == 'up':
				correlations = [ corr_list[i] >=  corr_threshold for i in range(len(dfp))]
				dfp = dfp.loc[ correlations ]
				correlations = [x for x in corr_list if x >= corr_threshold]
				s_filtered = f'positive correlation: Spearman corr >= {corr_threshold}'
	
			elif up_down_all == 'dw':
				correlations = [ corr_list[i] <= -corr_threshold for i in range(len(dfp))]
				dfp = dfp.loc[ correlations ]
				correlations = [x for x in corr_list if x <= -corr_threshold]
				s_filtered = f'negative correlation: Spearman corr <= {-corr_threshold}'
	
			else:
				correlations = [ abs(corr_list[i]) < corr_threshold for i in range(len(dfp))]
				dfp = dfp.loc[ correlations ]
				correlations = [x for x in corr_list if abs(x) < corr_threshold]
				s_filtered = f'no correlation: Spearman abs(corr) < {corr_threshold}'

			return dfp, s_filtered, correlations		

		dfp_female, dfp_male = pd.DataFrame(), pd.DataFrame()
		corr_female, corr_male = [], []
		s_filtered = '???'

		for (gender, gender_list) in [ ('female', self.group_female_list), ('male', self.group_male_list)]:
			dfp = dfpiv[gender_list].copy()
			dfp = dfp.fillna(0)

			dfp2 = pd.DataFrame()

			if type_modulation == 'all':
				dfp2, s_filtered = calc_all_modulation(dfp, up_down_all, lfc_threshold)
			else:
				seqx = list(np.arange(len(gender_list)))
				dfp2, s_filtered, correlations = calc_correlations(dfp, up_down_all, seqx, corr_threshold)

			if gender == 'female':
				dfp_female = dfp2.copy()
				corr_female = correlations
			else:
				dfp_male = dfp2.copy()
				corr_male = correlations

		if dfp_female.empty and dfp_male.empty:
			if verbose: print(f"No data found for {type_modulation} {up_down_all}")
			return None, pd.DataFrame(), pd.DataFrame(), ''

		leg_y_offset = -0.25

		if not dfp_female.empty and not dfp_male.empty:
			fig = make_subplots(rows=2, cols=1, subplot_titles=['female', 'male'], vertical_spacing=vertical_spacing)
		elif not dfp_female.empty:
			fig = make_subplots(rows=1, cols=1, subplot_titles=['female'], vertical_spacing=vertical_spacing)
			height = int(height * 0.7)
			leg_y_offset = -0.35
		else:
			fig = make_subplots(rows=1, cols=1, subplot_titles=['male'], vertical_spacing=vertical_spacing)
			height = int(height * 0.7)
			leg_y_offset = -0.35

		nrow = 0
		num_legends = 0
		for (gender, dfp, correlations) in [('female', dfp_female, corr_female), ('male', dfp_male, corr_male)]:
			if dfp.empty: continue

			nrow += 1
			cols = list(dfp.columns)
			df2 = dfp.reset_index()
			symbols = list(df2.symbol)

			df_list=[]
			for col in cols:
				dfa = df2.loc[:, ['symbol', col] ].copy()
				dfa['case'] = col
				dfa.columns = ['symbol', 'lfc', 'case']
				df_list.append(dfa)
			
			df = pd.concat(df_list)
			df.reset_index(inplace=True, drop=True)
			# self.df = df
			found=True

			dic_symb={}
			for icolor, symbol in enumerate(symbols):
				dic_symb[symbol] = self.my_colors[icolor]

			i=-1
			for symbol in symbols:
				df2 = df[df.symbol == symbol].copy()
				df2.reset_index(inplace=True, drop=True)
				# self.df2 = df2
				color = dic_symb[symbol]

				if type_modulation == 'all':
					name = f"{gender:6} - {symbol}"
				else:
					i+=1
					corr = correlations[i]
					name = f"{gender:6} - {symbol:8} corr: {corr:.3f}"

				num_legends += len(df2)
				fig.add_trace(go.Scatter(x=df2.case, y=df2.lfc,
										mode='lines', line=dict(color=color, width=line_width), name=name), row=nrow, col=1)

		if not found:
			if verbose: print(f"Warning: no genes are modulated using the filter {type_modulation} - {up_down_all}")
			return None, pd.DataFrame(), pd.DataFrame(), ''

		title = f"{self.s_deg_dap}s' LFC line-plot for '{pathway}'<br>{s_filtered}"
		title = title.strip()

		num_legends -= 15
		if num_legends > 0:
			inc_heigh = np.ceil(  num_legends / 5 )
			height += inc_heigh * 30
			leg_y_offset -= inc_heigh * 0.025


		fig.update_layout(
					autosize=True,
					title=dict( text=title, font=dict(size=title_fontsize, family="Arial", color="black")),
					width=width,
					height=height,
					plot_bgcolor=plot_bgcolor,
					xaxis_title="cases",
					yaxis_title=f"{self.s_deg_dap}s LFC",
					showlegend=True,
					font   = dict(family="Arial", size=fontsize, color="black" ),
					legend = dict( orientation="h", yanchor="bottom", y=leg_y_offset, xanchor="left", x=0, 
				   				   # entrywidth=100, entrywidthmode="pixels"
				   				   font=dict(family="Arial", size=leg_fontsize, color="black" ) ),
		)

		# Update font size for all subplot titles
		fig.update_annotations(font_size=title_fontsize) 

		fig.update_xaxes(tickfont=dict(size=tick_fontsize, family='Arial', color='black'))
		fig.update_yaxes(tickfont=dict(size=tick_fontsize, family='Arial', color='black'))

		figname0 = title_replace(title)
		figname0 = figname0.replace("'","").replace('"','')

		figname = osjoin(self.root_figure, figname0+'.html')
		if verbose: print(">>> HTML and png saved:", figname)

		try:
			fig.write_html(figname)
			fig.write_image(figname.replace('.html', '.png'))
		except:
			print(f"Error writing {figname}")

		figname = osjoin(self.root_fig_md, figname0+'.png')
		try:
			fig.write_image(figname)
		except:
			print(f"Error writing {figname}")

		return fig, dff, dfpiv, figname0
	

	def lineplot_one_pathway_LFCs_wo_gender(self, pathway_and_id:List, 
						type_modulation:str='all', up_down_all:str='up', 
						lfc_threshold:float=1., corr_threshold:float=0.4,
						width:int=1100, height:int=700, line_width:int=3,
						fontsize:int=16, leg_fontsize:int=16, 
						tick_fontsize:int=16, title_fontsize:int=18,
						plot_bgcolor:str='lightgray',  vertical_spacing:float=0.20,
						verbose:bool=False) -> Tuple[object, pd.DataFrame, pd.DataFrame, str]:
		'''
			LFC modulation, line plot only by cases: wo gender and age
		'''

		if not isinstance(pathway_and_id, list) or pathway_and_id == []:
			print("Pathway and ID is an empty list.")
			return None, pd.DataFrame(), pd.DataFrame(), ''

		pathway, _ = pathway_and_id

		dfpiv, dff = self.calc_pivot_one_pathway_LFC(pathway_and_id, verbose=verbose)

		def calc_all_modulation(dfp:pd.DataFrame, up_down_all:str, lfc_threshold:float):
			if up_down_all == 'up':
				dfp = dfp[ [ np.sum( dfp.iloc[i].apply(lambda x: 1 if x >= lfc_threshold else 0) ) > 0  for i in range(len(dfp))] ]
				s_filtered = f"at least one upregulated: LFC >= {lfc_threshold}"
	
			elif up_down_all == 'dw':
				dfp = dfp[ [ np.sum( dfp.iloc[i].apply(lambda x: 1 if x <= -lfc_threshold else 0) ) > 0  for i in range(len(dfp))] ]
				s_filtered = f"at least one downregulated: LFC <= {-lfc_threshold}"
	
			else:
				dfp = dfp[ [ (np.sum( dfp.iloc[i].apply(lambda x: 1 if x >= lfc_threshold else 0) ) +
							np.sum( dfp.iloc[i].apply(lambda x: 1 if x <= -lfc_threshold else 0) ) ) == 0  for i in range(len(dfp))] ]
	
				s_filtered = f"all non-regulated: abs any LFC < {lfc_threshold}"

			return dfp, s_filtered			


		def calc_correlations(dfp:pd.DataFrame, up_down_all:str, seqx:List, corr_threshold:float):

			corr_list = [  spearmanr(seqx, dfp.iloc[i,:])[0] for i in range(len(dfp))]

			if up_down_all == 'up':
				correlations = [ corr_list[i] >=  corr_threshold for i in range(len(dfp))]
				dfp = dfp.loc[ correlations ]
				correlations = [x for x in corr_list if x >= corr_threshold]
				s_filtered = f'positive correlation: Spearman corr >= {corr_threshold}'
	
			elif up_down_all == 'dw':
				correlations = [ corr_list[i] <= -corr_threshold for i in range(len(dfp))]
				dfp = dfp.loc[ correlations ]
				correlations = [x for x in corr_list if x <= -corr_threshold]
				s_filtered = f'negative correlation: Spearman corr <= {-corr_threshold}'
	
			else:
				correlations = [ abs(corr_list[i]) < corr_threshold for i in range(len(dfp))]
				dfp = dfp.loc[ correlations ]
				correlations = [x for x in corr_list if abs(x) < corr_threshold]
				s_filtered = f'no correlation: Spearman abs(corr) < {corr_threshold}'

			return dfp, s_filtered, correlations		

		correlations=[]

		dfp = dfpiv.copy()
		dfp = dfp.fillna(0)

		if type_modulation == 'all':
			dfp2, s_filtered = calc_all_modulation(dfp, up_down_all, lfc_threshold)
		else:
			seqx = list(np.arange(len(self.case_list)))
			dfp2, s_filtered, correlations = calc_correlations(dfp, up_down_all, seqx, corr_threshold)

		if dfp2.empty:
			if verbose: print(f"No data found for {type_modulation} {up_down_all}")
			return None, pd.DataFrame(), pd.DataFrame(), ''

		leg_y_offset = -0.25

		fig = make_subplots(rows=1, cols=1, subplot_titles=['female', 'male'], vertical_spacing=vertical_spacing)

		num_legends = 0

		cols = list(dfp2.columns)

		if len(dfp2) > 100:
			dfp2 = self.is_modulated(dfp2, lfc_threshold=1)

		df2 = dfp2.reset_index()

		n_colors = len(self.my_colors)
		if len(df2) > n_colors:
			print("To many symbols", len(df2))
			df2 = df2.head(n_colors)
			
		symbols = list(df2.symbol)

		# self.df2 = df2
		# self.symbols = symbols

		df_list=[]
		for col in cols:
			dfa = df2.loc[:, ['symbol', col] ].copy()
			dfa['case'] = col
			dfa.columns = ['symbol', 'lfc', 'case']
			df_list.append(dfa)
		
		df = pd.concat(df_list)
		df.reset_index(inplace=True, drop=True)

		dic_symb={}
		for icolor, symbol in enumerate(symbols):
			dic_symb[symbol] = self.my_colors[icolor]

		i=-1
		for symbol in symbols:
			df2 = df[df.symbol == symbol].copy()
			df2.reset_index(inplace=True, drop=True)
			color = dic_symb[symbol]

			if type_modulation == 'all':
				name = symbol
			else:
				i+=1
				corr = correlations[i]
				name = f"{symbol:8} corr: {corr:.3f}"

			num_legends += len(df2)
			fig.add_trace(go.Scatter(x=df2.case, y=df2.lfc,
									mode='lines', line=dict(color=color, width=line_width), name=name), row=1, col=1)

		title = f"{self.s_deg_dap}s' LFC line-plot for '{pathway}'<br>{s_filtered}"
		title = title.strip()

		num_legends -= 15
		if num_legends > 0:

			mini = df2.lfc.min()
			maxi= df2.lfc.max()

			ampli = maxi-mini
			delta = ampli/3.
			
			inc_heigh = np.ceil(  num_legends / 5 )
			height	   += inc_heigh * 30
			leg_y_offset -= inc_heigh * delta


		fig.update_layout(
					autosize=True,
					title=dict( text=title, font=dict(size=title_fontsize, family="Arial", color="black")),
					width=width,
					height=height,
					plot_bgcolor=plot_bgcolor,
					xaxis_title="cases",
					yaxis_title=f"{self.s_deg_dap}s LFC",
					showlegend=True,
					font   = dict(family="Arial", size=fontsize, color="black" ),
					legend = dict( orientation="h", yanchor="bottom", y=leg_y_offset, xanchor="left", x=0, 
				   				   # entrywidth=100, entrywidthmode="pixels"
				   				   font=dict(family="Arial", size=leg_fontsize, color="black" ) ),
		)

		# Update font size for all subplot titles
		fig.update_annotations(font_size=title_fontsize) 

		fig.update_xaxes(tickfont=dict(size=tick_fontsize, family='Arial', color='black'))
		fig.update_yaxes(tickfont=dict(size=tick_fontsize, family='Arial', color='black'))

		figname0 = title_replace(title)
		figname0 = figname0.replace("'","").replace('"','')

		figname = osjoin(self.root_figure, figname0+'.html')
		if verbose: print(">>> HTML and png saved:", figname)

		try:
			fig.write_html(figname)
			fig.write_image(figname.replace('.html', '.png'))
		except:
			print(f"Error writing {figname}")

		figname = osjoin(self.root_fig_md, figname0+'.png')
		try:
			fig.write_image(figname)
		except:
			print(f"Error writing {figname}")

		return fig, dff, dfpiv, figname0

	def calc_pathway_modulation_up_down(self, df3:pd.DataFrame, regulation:str, type_ptw_mod:str, case:str,
										with_pPMI_obs:bool, do_case_translate:bool,
										echo_pathway:bool=True) -> Tuple[str, dict]:
		dic = {}

		if df3.empty:
			return '', dic

		text = '\n'

		for i in range(len(df3)):
			row = df3.iloc[i]
			_ident = '  * '

			if echo_pathway:
				if pd.isnull(row.fdr):
					text += echo_print(f"{_ident}**{row.pathway}** - not enriched - pPMI is {regulation}")
				else:
					text += echo_print(f"{_ident}**{row.pathway}** (FDR={row.fdr:.2e}) - pPMI is {regulation}")

			if with_pPMI_obs:
				text += '\n'
				text += echo_print(f"pPMI_total = {row.pPMI_total:.3f} ratio_up_dw = {row.ratio_up_dw:.3f}")
				text += echo_print(f"All {self.s_deg_dap}s #{row.n_mod_in_pathway} ({row.perc_mod_in_pathway*100:.2f}%): Up={row.n_mod_up_in_pathway} and Down={row.n_mod_dw_in_pathway}")
				
			text += '\n\n'

			_ident = '	* '

			_up = 'Upregulated' if do_case_translate else 'Up'
			if row.n_mod_up_in_pathway > 0:
				mod_up_list = row.mod_up_in_pathway
				lfc_up_list = row.lfc_up

				mod_up_list = mod_up_list if isinstance(mod_up_list, list) else eval(mod_up_list)
				lfc_up_list = lfc_up_list if isinstance(lfc_up_list, list) else eval(lfc_up_list)

				text_modulation = self.build_modulation_list(mod_up_list, lfc_up_list)
				text += echo_print(f"{_ident}{_up} #{row.n_mod_up_in_pathway} = {', '.join(text_modulation)}\n")


			_dw = 'Downregulated' if do_case_translate else 'Dw'
			if row.n_mod_dw_in_pathway > 0:
				mod_dw_list = row.mod_dw_in_pathway
				lfc_dw_list = row.lfc_dw

				mod_dw_list = mod_dw_list if isinstance(mod_dw_list, list) else eval(mod_dw_list)
				lfc_dw_list = lfc_dw_list if isinstance(lfc_dw_list, list) else eval(lfc_dw_list)

				text_modulation = self.build_modulation_list(mod_dw_list, lfc_dw_list)
				text += echo_print(f"{_ident}{_dw} #{row.n_mod_dw_in_pathway} = {', '.join(text_modulation)}\n")

			dic[i] = {}
			dic2 = dic[i]

			dic2['case']		 = case
			dic2['type_ptw_mod'] = type_ptw_mod
			dic2['pathway']	  = row.pathway
			dic2['fdr']		  = row.fdr
			dic2['pPMI_total'] = row.pPMI_total
			dic2['ratio_up_dw']	   = row.ratio_up_dw

			dic2['n_mod_in_pathway']	= row.n_mod_in_pathway
			dic2['perc_mod_in_pathway'] = row.perc_mod_in_pathway
			dic2['n_mod_up_in_pathway'] = row.n_mod_up_in_pathway
			dic2['n_mod_dw_in_pathway'] = row.n_mod_dw_in_pathway

			dic2['mod_up_in_pathway'] = row.mod_up_in_pathway
			dic2['mod_dw_in_pathway'] = row.mod_dw_in_pathway

			dic2['lfc_up'] = row.lfc_up
			dic2['lfc_dw'] = row.lfc_dw

			text += echo_print("")

		return text, dic

	def calc_pPMI_summary(self, pathway_list:List=[], title_head:str='', with_pPMI_obs:bool=False, 
						  do_case_translate:bool=True, 
						  echo_pathway:bool=True, force:bool=False, 
						  verbose:bool=False) -> Tuple[pd.DataFrame, str, pd.DataFrame]:

		dff = self.calc_all_pPMI(force=force, verbose=verbose)

		if dff is None or dff.empty:
			self.dff = dff
			return pd.DataFrame(), '', pd.DataFrame()
		
		if pathway_list != []:
			dff = dff[dff.pathway.isin(pathway_list)].copy()
			
		dff = dff.sort_values('pathway', ascending=True)
		dff.reset_index(inplace=True, drop=True)
		self.dff = dff

		if len(pathway_list) == 1:
			pathway = pathway_list[0]
			text_all = f"## Study of the pathway '{pathway}'\n"
		else:
			text_all = ''

		df_list = []
		for case in self.case_list:
			if df_list == []:
				title_head = f"### Study: **{title_head}**"
			else:
				title_head = ''
				
			text, dfa = self.to_text_pathway_modulation(title_head, dff, case=case, with_pPMI_obs=with_pPMI_obs, 
											   			echo_pathway=echo_pathway, do_case_translate=do_case_translate)
			
			if dfa is None or dfa.empty:
				if verbose: print(f"Could not transform to_text_pathway_modulation() for case {case}")
				continue
	
			df_list.append(dfa)
			if text_all == '':
				text_all = text
			else:
				if len(pathway_list) == 1:
					text_all += text
				else:
					text_all += '\n' + text

		if pathway_list == []:
			fname = f"summary_pathway_modulation_degs_{self.normalization}_geneset_{self.geneset_num}_{self.geneset_lib}.txt"
			write_txt(text_all, fname, self.root_ressum, verbose=verbose)

		if df_list == []:
			print("Nothing found in calc_pPMI_summary()")
			df_sum_ptw = pd.DataFrame()
		else:
			df_sum_ptw = pd.concat(df_list)
			df_sum_ptw.reset_index(inplace=True, drop=True)

		return dff, text_all, df_sum_ptw

	def to_text_pathway_modulation(self, title_head:str, dff:pd.DataFrame, case:str, with_pPMI_obs:bool=False,
								   do_case_translate:bool=True, echo_pathway:bool=True) -> Tuple[str, pd.DataFrame]:

		text = title_head

		if do_case_translate:
			text += f"\n\n#### case {case} ({self.translate_case(case)})\n"
		else:
			text += f"\n\n#### case {case}\n"


		regulation = "upregulated"
		type_ptw_mod = 'up'
		df2 = dff[dff.case == case]
		df3 = df2[df2.pPMI_total >= self.tolerance_pPMI]
		text2, dic_up = self.calc_pathway_modulation_up_down(df3, regulation, type_ptw_mod, case, with_pPMI_obs, do_case_translate, echo_pathway)
		if text2 != '': text += text2

		regulation = "balanced"
		type_ptw_mod = 'not'
		df3 = df2[ (df2.pPMI_total > -self.tolerance_pPMI) & (df2.pPMI_total < self.tolerance_pPMI) ]
		text2, dic_not = self.calc_pathway_modulation_up_down(df3, regulation, type_ptw_mod, case, with_pPMI_obs, do_case_translate, echo_pathway)
		if text2 != '': text += text2

		regulation = "downregulated"
		type_ptw_mod = 'down'
		df3 = df2[df2.pPMI_total <= -self.tolerance_pPMI]
		text2, dic_dw = self.calc_pathway_modulation_up_down(df3, regulation, type_ptw_mod, case, with_pPMI_obs, do_case_translate, echo_pathway)
		if text2 != '': text += text2

		df_list = []
		if dic_up != {}:
			df_up  = pd.DataFrame(dic_up).T
			df_list.append(df_up)

		if dic_not != {}:
			df_not = pd.DataFrame(dic_not).T
			df_list.append(df_not)

		if dic_dw != {}:
			df_dw  = pd.DataFrame(dic_dw).T
			df_list.append(df_dw)

		if df_list == []:
			return text, pd.DataFrame()
		
		dfa = pd.concat(df_list)
		dfa.reset_index(inplace=True, drop=True)

		return text, dfa


	def calc_all_pPMI(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		if self.type_sat_ptw_index == 'discrete':
			s_saturation = 'none'
		elif self.type_sat_ptw_index == 'linear':
			s_saturation = 'none'
		elif self.type_sat_ptw_index == 'linear_sat':
			s_saturation = f'{self.saturation_lfc_param:.3f}'
		else:
			print("Please choose: discrete, linear, or linear_sat as type_sat_ptw_index")
			raise Exception("stop: calc_all_pPMI()")

		fname = self.fname_pPMI%(self.project, self.type_sat_ptw_index, s_saturation, self.normalization)
		fname = title_replace(fname)
		filename = osjoin(self.root_ressum, fname)

		if exists(filename) and not force:
			dff = pdreadcsv(fname, self.root_ressum, verbose=verbose)
			return dff

		self.df_enr_all = self.list_all_df_erichments()

		df_list = []
		for case in self.case_list:
			print(">>> case", case)
			dfa =  self.calc_pPMI_per_case(case, verbose=verbose)

			if dfa is None or dfa.empty:
				continue

			fname2 = self.fname_pathway_case_index%(self.project, case, self.normalization)
			_ = pdwritecsv(dfa, fname2, self.root_ressum, verbose=verbose)

			df_list.append(dfa)

		if df_list == []:
			print("Nothing was found.")
			return pd.DataFrame()

		dff = pd.concat(df_list)
		dff.reset_index(inplace=True, drop=True)

		_ = pdwritecsv(dff, fname, self.root_ressum, verbose=verbose)

		return dff

	def list_all_df_erichments(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		fname_all_dfenr = self.fname_all_dfenr%(self.disease)
		filename = osjoin(self.root_result, fname_all_dfenr)

		if exists(filename) and not force:
			dfall = pdreadcsv(fname_all_dfenr, self.root_result, verbose=verbose)
			return dfall

		df_list = []

		cols = ['pathway', 'pathway_id']
		for case in self.case_list:
			ret, _, _, _ = self.open_case(case, verbose=False)
			if not ret:
				continue

			try:
				dfa = self.df_enr[ cols]
			except:
				dfa = self.df_enr[ ['pathway'] ]

			df_list.append(dfa)

		dfall = pd.concat(df_list)
		dfall = dfall.drop_duplicates()
		dfall = dfall.sort_values('pathway')
		dfall.reset_index(inplace=True, drop=True)

		pdwritecsv(dfall, fname_all_dfenr, self.root_result, verbose=verbose)
		return dfall


	def calc_pPMI_per_case(self, case:str, verbose:bool=False) -> pd.DataFrame:

		ret, _, _, _ = self.open_case(case, verbose=verbose)
		if not ret:
			return pd.DataFrame()

		icount = -1; dic = {}
		for i in range(len(self.df_enr_all)):
			row = self.df_enr_all.iloc[i]

			pathway_id = row.pathway_id
			pathway, genes_in_pathway = self.reactome_find_genes_in_pathway(pathway_id, _type='pathway_id')

			dfa = self.df_enr0[self.df_enr0.pathway_id == pathway_id]
			if dfa.empty:
				fdr, pval = None, None
			else:
				row = self.df_enr0.iloc[0]
				fdr  = row.fdr
				pval = row.pval


			if len(genes_in_pathway) == 0:
				# there are pathways related to other cases!
				# print(f"No genes were found for '{pathway_id} - {pathway}'")
				continue

			# print(f"###Found '{pathway_id} - {pathway}' genes {num_of_genes} / {genes_in_pathway}")

			#---------------- Look for all genes in dflfc ------------------------
			dflfc_ptw = self.dflfc_ori[self.dflfc_ori.symbol.isin(genes_in_pathway)].copy()
			if dflfc_ptw.empty:
				if verbose: print(f"No genes were found for '{pathway_id} - {pathway}' in LFC table for case {case}")
				continue
		
			# minimum modulation like 0.4 else is noise - remove them
			dflfc_ptw = dflfc_ptw[dflfc_ptw.abs_lfc >= self.min_lfc_modulation]

			if dflfc_ptw.empty:
				if verbose: print(f"No genes were found for '{pathway_id} - {pathway}' in LFC table for case {case}")
				continue

			dflfc_ptw.reset_index(inplace=True, drop=True)

			#------------- modulated genes for determined pathway ------------------------------
			dflfc_ptw_up = dflfc_ptw[dflfc_ptw.lfc > 0].copy()
			dflfc_ptw_up.reset_index(inplace=True, drop=True)

			dflfc_ptw_dw = dflfc_ptw[dflfc_ptw.lfc < 0].copy()
			dflfc_ptw_dw.reset_index(inplace=True, drop=True)

			mod_up_in_pathway = list(dflfc_ptw_up.symbol)
			mod_dw_in_pathway = list(dflfc_ptw_dw.symbol)

			n_mod_up_in_pathway = len(mod_up_in_pathway)
			n_mod_dw_in_pathway = len(mod_dw_in_pathway)

			# here aaa

			#-------------- sum all modulated genes to calculate the pPMI positive/negative index ---------------
			lfc_pPMI_up, lfc_pPMI_dw = 0, 0

			if self.type_sat_ptw_index == 'discrete':
				lfc_pPMI_up = n_mod_up_in_pathway
				lfc_pPMI_dw = n_mod_dw_in_pathway

			elif self.type_sat_ptw_index == 'linear':
				lfc_pPMI_up = dflfc_ptw_up.lfc.sum()
				lfc_pPMI_dw = dflfc_ptw_dw.lfc.sum()

			elif self.type_sat_ptw_index == 'linear_sat':
				if dflfc_ptw_up.empty:
					lfc_up = []
				else:
					lfc_up = [x if x < self.saturation_lfc_param else self.saturation_lfc_param for x in dflfc_ptw_up.lfc]

				if dflfc_ptw_dw.empty:
					lfc_dw = []
				else:
					lfc_dw = [-self.saturation_lfc_param if x < -self.saturation_lfc_param else x for x in dflfc_ptw_dw.lfc]

				lfc_pPMI_up = np.sum(lfc_up)
				lfc_pPMI_dw = np.sum(lfc_dw)

			# --------------- total pPMI ------------------------------------------
			# calc ratio: if down == 0 like down = 0.01; if up == 0; like up = 0.01
			
			if lfc_pPMI_dw == 0:
				ratio = lfc_pPMI_up * self.saturation_lfc_param
			else:
				if lfc_pPMI_up == 0:
					ratio = 1 / (self.saturation_lfc_param * lfc_pPMI_dw)
				else:
					ratio = lfc_pPMI_up / lfc_pPMI_dw

			ratio = abs(ratio)
			if ratio > 0:
				pPMI_total = np.log2(ratio)
			else:
				pPMI_total = None

			#--------------------- normalized pPMI --------------------------------------
			lfc_pPMI_up_norm = lfc_pPMI_up / n_mod_up_in_pathway if n_mod_up_in_pathway > 0 else 0
			lfc_pPMI_dw_norm = lfc_pPMI_dw / n_mod_dw_in_pathway if n_mod_dw_in_pathway > 0 else 0

			# min modulario = 0.40 ... now 0.20 as a threshold if no modulatio
			min2_lfc_modulation = self.min_lfc_modulation/2

			if lfc_pPMI_dw_norm == 0:
				ratio_norm = lfc_pPMI_up_norm / min2_lfc_modulation
			else:
				if lfc_pPMI_up_norm == 0:
					ratio_norm = min2_lfc_modulation / lfc_pPMI_dw_norm
				else:
					ratio_norm = lfc_pPMI_up_norm / lfc_pPMI_dw_norm

			ratio_norm = abs(ratio_norm)
			if ratio_norm > 0:
				pPMI_norm = np.log2(ratio_norm)
			else:
				pPMI_norm = None
			
			'''
			if pathway == 'Platelet Aggregation (Plug Formation)':
													print(">>>", case, lfc_pPMI_up, lfc_pPMI_dw, ratio, pPMI_total)
													print("up", dflfc_ptw_up.symbol, dflfc_ptw_up.lfc)
													print("dw", dflfc_ptw_dw.symbol, dflfc_ptw_dw.lfc, '\n\n')
			'''

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			dic2['case'] = case
			dic2['pathway_id'] = pathway_id
			dic2['pathway']	= pathway

			dic2['pPMI_total'] = pPMI_total
			dic2['ratio_up_dw'] = ratio

			dic2['pPMI_norm'] = pPMI_norm
			dic2['ratio_norm_up_dw'] = ratio_norm

			dic2['lfc_pPMI_up'] = lfc_pPMI_up
			dic2['lfc_pPMI_dw'] = lfc_pPMI_dw

			n_mod_in_pathway = len(mod_up_in_pathway) + len(mod_dw_in_pathway)

			N = len(genes_in_pathway)
			dic2['n_genes_annot_in_pathway'] = N

			dic2['n_mod_in_pathway']	= n_mod_in_pathway
			dic2['perc_mod_in_pathway'] = n_mod_in_pathway / N
			dic2['n_mod_up_in_pathway'] = len(mod_up_in_pathway)
			dic2['n_mod_dw_in_pathway'] = len(mod_dw_in_pathway)

			dic2['fdr'] = fdr
			dic2['pval'] = pval

			dic2['mod_up_in_pathway'] = mod_up_in_pathway
			dic2['mod_dw_in_pathway'] = mod_dw_in_pathway

			dic2['lfc_up'] = list( np.round(dflfc_ptw_up.lfc, 3) )
			dic2['lfc_dw'] = list( np.round(dflfc_ptw_dw.lfc, 3) )

			dic2['genes_annot_in_pathway'] = genes_in_pathway

		if len(dic) == 0:
			return pd.DataFrame()

		df = pd.DataFrame(dic).T
		df = df.sort_values("pPMI_total", ascending=False)
		df.reset_index(inplace=True, drop=True)

		return df

	def set_s_filtered(self, type_modulation:str, up_down_all:str, diff_cutoff:float, lfc_threshold:float) -> str:

		if self.has_gender:
			s_filtered = 'comparison between genders:<br>'
		else:
			s_filtered = 'comparing last-first case:<br>'

		s_filtered += f"type of LFC modulation = '{type_modulation}'"
		if type_modulation == 'all':
			pass
		elif type_modulation == 'similar':
			s_filtered += f' (abs_diff < {diff_cutoff})'
		else:
			#--------- different: not all, not similar ------------
			s_filtered += f' (abs diff >= {diff_cutoff})'
		
		s_filtered += ' and LFC diff '

		if up_down_all == 'is up':
			s_filtered += f'up (diff >= {lfc_threshold})'
		elif up_down_all == 'dw':
			s_filtered += f'is dw (diff <= {-lfc_threshold})'
		else:
			s_filtered += f'is up or dw (abs_diff >= {lfc_threshold})'

		return s_filtered



	def calc_one_pathway_genes_modulation(self, pathway_and_id:List, 
									type_modulation:str='all', up_down_all:str='all',
									lfc_threshold:float=1.,  diff_cutoff:float=1,
									verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame]:
		'''
			type_modulation = all, similar, different - uses diff_cutoff
				all: no filter
				similar: abs_diff < diff_cutoff]
				different: abs_diff >= diff_cutoff]

				up: means female have diff >= than male diff
				dw: means female have diff < =than male diff

			up_down_all uses lfc_threshold:
				if has_gender:
					up: means female have diff >= than male diff  
					dw: means female have diff < =than male diff
					all: both cases accepted
				else:
					up: means last case >= first case
					dw: means female have diff < than male diff
					all: both cases accepted
		'''

		if not isinstance(pathway_and_id, list) or pathway_and_id == []:
			print("Pathway and ID is an empty list.")
			return pd.DataFrame(), pd.DataFrame()

		# pathway, pathway_id = pathway_and_id
		pathway, _ = pathway_and_id

		dfpiv, dff = self.calc_pivot_one_pathway_LFC(pathway_and_id, verbose=verbose)

		dfpiv = dfpiv.fillna(0)

		if self.has_gender:
			# female
			cols = self.group_female_list
			df2 = dfpiv[cols]
			dfpiv['diff_female'] = df2[cols[-1]] - df2[cols[0]]
			dfpiv['abs_diff_female'] = np.abs(dfpiv['diff_female'].to_list())  

			# male
			cols = self.group_male_list
			df2 = dfpiv[cols]
			dfpiv['diff_male'] = df2[cols[-1]] - df2[cols[0]]
			dfpiv['abs_diff_male'] = np.abs(dfpiv['diff_male'].to_list())  

			dfpiv['diff'] = dfpiv['diff_male'] - dfpiv['diff_female']
		else:
			cols = self.case_list
			dfpiv['diff'] = dfpiv[cols[-1]] - dfpiv[cols[0]]

		dfpiv['abs_diff'] = np.abs(dfpiv['diff'].to_list())

		if type_modulation == 'all':
			pass
		elif type_modulation == 'similar':
			dfpiv = dfpiv[dfpiv.abs_diff < diff_cutoff]
		else:
			#--------- different: not all, not similar ------------
			dfpiv = dfpiv[dfpiv.abs_diff >= diff_cutoff]


		if up_down_all == 'up':
			# dfpiv = dfpiv[ [is_greater_equal(dfpiv.iloc[i], lfc_threshold) for i in range(len(dfpiv))] ]
			dfpiv = dfpiv[ dfpiv['diff'] >= lfc_threshold ]
		elif up_down_all == 'dw':
			# dfpiv = dfpiv[ [is_less_equal(dfpiv.iloc[i], lfc_threshold) for i in range(len(dfpiv))] ]
			dfpiv = dfpiv[ dfpiv['diff'] <= -lfc_threshold ]
		else:
			#--- all - o not filter --------
			pass

		if dfpiv.empty:
			if verbose: print(f"No genes are modulated using the filter {type_modulation} {up_down_all} for '{pathway}'")
			return pd.DataFrame(), pd.DataFrame()

		dfpiv = dfpiv.sort_values('diff', ascending=True)

		return dff, dfpiv


	def get_fname_pathway_modulation(self, title_head:str, type_modulation:str, up_down_all:str, 
									 lfc_threshold:float, diff_cutoff:float) -> Tuple[str, str, str, str]:
		fname_summ = f'pPMI_table_{title_head}_summary_filtered_{type_modulation}_up_down_all_{up_down_all}_at_least_{lfc_threshold:.1f}_diff_cutoff_{diff_cutoff:.1f}.tsv'
		fname_summ = title_replace(fname_summ)
		fullname_summ = osjoin(self.root_ptw_modulation, fname_summ)

		fname_anal = f'pPMI_table_{title_head}_analytics_filtered_{type_modulation}_up_down_all_{up_down_all}_at_least_{lfc_threshold:.1f}_diff_cutoff_{diff_cutoff:.1f}.tsv'
		fname_anal = title_replace(fname_anal)
		fullname_anal = osjoin(self.root_ptw_modulation, fname_anal)

		return fname_summ, fullname_summ,fname_anal, fullname_anal

	def get_pathway_modulation(self, title_head:str, type_modulation:str, up_down_all:str, 
									 lfc_threshold:float, diff_cutoff:float,
									 verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame]:

		fname_summ, fullname_summ,fname_anal, fullname_anal = \
		self.get_fname_pathway_modulation(title_head, type_modulation, up_down_all, lfc_threshold, diff_cutoff)


		if not exists(fullname_summ):
			df_summ = pd.DataFrame()
		else:
			df_summ = pdreadcsv(fname_summ, self.root_ptw_modulation, verbose=verbose)

		if not exists(fullname_anal):
			df_anal = pd.DataFrame()
		else:
			df_anal = pdreadcsv(fname_anal, self.root_ptw_modulation, verbose=verbose)

		return df_summ, df_anal


	# plot_and_calc_pathway_modulation
	def calc_pathway_modulation(self,set_pathway:List, 
							type_modulation:str='all', up_down_all:str='all',
							lfc_threshold:float=1, diff_cutoff:float=2., 
							force:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, str]:

		self.type_modulation = type_modulation
		self.lfc_threshold	= lfc_threshold
		self.diff_cutoff = diff_cutoff

		s_filtered = self.set_s_filtered(type_modulation, up_down_all, diff_cutoff, lfc_threshold)

		try:
			title_head, pathway_list, _ = set_pathway
		except:
			print(f"Error: wrong set_pathway: {set_pathway}")
			return pd.DataFrame(), pd.DataFrame(), s_filtered

		pathway_and_id_list = [[x, self.which_reacotme_id(x)] for x in pathway_list]

		fname_summ, fullname_summ,fname_anal, fullname_anal = \
			self.get_fname_pathway_modulation(title_head, type_modulation, up_down_all, lfc_threshold, diff_cutoff)

		if exists(fullname_summ) and  exists(fullname_anal) and not force:
			df_summ, df_anal = \
			self.get_pathway_modulation(title_head, type_modulation, up_down_all, lfc_threshold, diff_cutoff, verbose=verbose)
			return df_summ, df_anal, s_filtered
		
		dic={}; icount=-1
		# dfpiv_female_list, dfpiv_male_list,
		df_anal_list=[]

		for pathway_and_id in pathway_and_id_list:
			pathway, pathway_id = pathway_and_id
			if verbose: print(">>>", pathway, pathway_id, '\n')

			dff, dfpiv = self.calc_one_pathway_genes_modulation(pathway_and_id,
								type_modulation=type_modulation, up_down_all=up_down_all,
								lfc_threshold=lfc_threshold, diff_cutoff=diff_cutoff, verbose=verbose)

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			dic2['title_head'] = title_head
			dic2['pathway_id'] = pathway_id
			dic2['pathway'] = pathway
			dic2['type_modulation'] = type_modulation
			dic2['up_down_all'] = up_down_all
			dic2['lfc_threshold'] = lfc_threshold
			dic2['diff_cutoff'] = diff_cutoff
			dic2['filter'] = s_filtered
			dic2['dff'] = dff
			dic2['dfpiv'] = dfpiv

			df2 = dfpiv.reset_index().copy()
			df2['pathway_id'] = pathway_id
			df2['pathway'] = pathway

			df_anal_list.append(df2)

		if df_anal_list == []:
			print("Nothing found.")
			return pd.DataFrame(), pd.DataFrame(), s_filtered

		df_anal = pd.concat(df_anal_list)
		df_anal.reset_index(inplace=True, drop=True)

		if self.has_gender:
			cols2 = self.group_female_list
			seqx = np.arange(len(cols2))

			corr_list = [ spearmanr(seqx, df_anal.loc[i, cols2].to_list())[0] for i in range(len(df_anal))]
			df_anal['spearman_cor_female'] = corr_list

			cols2 = self.group_male_list
			seqx = np.arange(len(cols2))

			corr_list = [ spearmanr(seqx, df_anal.loc[i, cols2].to_list())[0] for i in range(len(df_anal))]
			df_anal['spearman_cor_male'] = corr_list
			
			if not df_anal.empty:
				df_anal = df_anal.sort_values(['spearman_cor_female', 'spearman_cor_male', 'diff'], ascending=[False, False, False])
		else:
			cols2 = self.case_list
			seqx = np.arange(len(cols2))

			corr_list = [ spearmanr(seqx, df_anal.loc[i, cols2].to_list())[0] for i in range(len(df_anal))]
			df_anal['spearman_cor'] = corr_list
			if not df_anal.empty:
				df_anal = df_anal.sort_values(['spearman_cor', 'diff'], ascending=[False, False])

		df_anal.reset_index(inplace=True, drop=True)
		df_summ = pd.DataFrame(dic).T

		_ = pdwritecsv(df_summ, fname_summ, self.root_ptw_modulation, verbose=verbose)
		_ = pdwritecsv(df_anal, fname_anal, self.root_ptw_modulation, verbose=verbose)

		return df_summ, df_anal, s_filtered

	def get_pPMI(self, case:str, pathway_id:str, verbose:bool=False) -> object:
		if self.dff_pPMI.empty:
			self.dff_pPMI = self.calc_all_pPMI()
		
		dfa = self.dff_pPMI[ (self.dff_pPMI.case == case) & (self.dff_pPMI.pathway_id == pathway_id)]
		
		if dfa.empty:
			pPMI = None
			if verbose: print(f"Error: get_pPMI {case} {pathway_id}")
		else:
			pPMI = dfa.iloc[0].pPMI_total
			pPMI = float(pPMI)
			# print(f"Ok: get_pPMI {case} {pathway_id} {pPMI}")

		return pPMI

	def get_title_heatmap_set_pathway(self, unique_gender:str, pathway_id:str, pathway:str,
									  title_head:str, s_filtered:str) -> str:
		if unique_gender == 'unique':
			title = f"{self.s_deg_dap} heatmap for {pathway} ({pathway_id})"
		else:
			gender = unique_gender
			title = f"{self.s_deg_dap} heatmap for {gender} and {pathway} ({pathway_id})"
			
		title += f"<br>in {title_head} " + s_filtered
		return title


	def plot_set_pathways_genes_heatmap(self, set_pathway:List, maxLen:int=75, limit_length:int=42,
									type_modulation:str='all', up_down_all:str='all',
									lfc_threshold:float=1,  diff_cutoff:float=2.,
									width:int=1100, row_height:int=70, header_height:int=200,
									plot_bgcolor:str='lightgray',  margin:dict=dict(l=100, r=80, t=120, b=80), 
									verbose:bool=False) -> dict:
		'''
			plot_set_pathways_genes_heatmap?
				. plot all modulated genes/proteins in a heatmap
				. type_modulation: all, similar, different
				. up_down_all: up, dw or all

			it saves image in the ./figures - auxiliar dir for markdown
		'''

		try:
			title_head, _, _ = set_pathway
		except:
			print(f"Error: wrong set_pathway: {set_pathway}")
			return {}

		_, df_anal,s_filtered = self.calc_pathway_modulation(set_pathway, 
								type_modulation=type_modulation, up_down_all=up_down_all,
								lfc_threshold=lfc_threshold, diff_cutoff=diff_cutoff, 
								verbose=verbose)

		if df_anal is None or df_anal.empty:
			print(f"Pathway list {set_pathway[0]} has no data.")
			return {}
		
		pathway_id_list = np.unique(df_anal.pathway_id)

		dic_fig={}; icount=-1
		for pathway_id in pathway_id_list:
			df2 = df_anal[df_anal.pathway_id == pathway_id].copy()
			if df2.empty:
				print(f"Error: ??? no pathway {pathway_id}???")
				continue

			df2.reset_index(inplace=True, drop=True)
			pathway = df2.iloc[0].pathway
			self.df2 = df2

			maxi = df2[self.case_list].max().max()
			mini = df2[self.case_list].min().min()

			zlim = int(maxi-mini)+1
			zmin = -zlim
			zmax = +zlim

			for unique_gender in self.unique_gender_list:
				if unique_gender == 'unique':
					cols = self.case_list
					col_diff = "diff"
					col_abs_diff = "abs_diff"
					col_corr = 'spearman_cor'
				else:
					if unique_gender == 'Female':
						cols = self.group_female_list
						col_diff = "diff_female"
						col_abs_diff = "abs_diff_female"
						col_corr = 'spearman_cor_female'
					else:
						cols = self.group_male_list
						col_diff = "diff_male"
						col_abs_diff = "abs_diff_male"
						col_corr = 'spearman_cor_male'

				sel_cols = ['symbol'] + cols + [col_diff, col_abs_diff, col_corr]
				df3 = df2[sel_cols].copy()

				new_cols = ['symbol'] + cols + ['diff', 'abs_diff', 'spearman_cor']
				df3.columns = new_cols

				lista = [f"diff={df3.loc[i,'diff']:.2f} cor={df3.loc[i,'spearman_cor']:.2f} - {df3.loc[i,'symbol'] + '  '}" for i in range(len(df3))]
				lista = [break_line_per_length(x, sep='<br>', maxLen=maxLen) for x in lista]
				df3 = df3.set_index(pd.Series(lista))
				df3 = df3.sort_values(['spearman_cor', 'abs_diff'], ascending=[True, False])

				plot_cols = cols + ['diff', 'spearman_cor']
				df4 = df3[plot_cols].copy()

				# artificial scale, for color propositions
				df4.spearman_cor *= zlim

				cols_label = []
				for case in cols:
					pPMI = self.get_pPMI(case, pathway_id)

					if pPMI is None:
						cols_label.append(f"{case}<br>pPMI=---")
					else:
						cols_label.append(f"{case}<br>pPMI={pPMI:.2f}")

				df4.columns = cols_label + ['diff', 'spearman_cor']

				if len(df4) > limit_length:
					center = int(len(df4)/2)
					third = int(limit_length/3)
					sixth = int(limit_length/6)

					df4_up  = df4.head(third)
					df4_dw  = df4.tail(third)
					df4_mid = df4.iloc[(center-sixth):(center+sixth+1)]

					df_empty = create_empty_df(df4, vals=[])
					
					df4 = pd.concat([df4_up, df_empty, df4_mid, df_empty, df4_dw])

					lista_up  = df4_up.index.to_list()
					lista_dw  = df4_dw.index.to_list()
					lista_mid = df4_mid.index.to_list()
					lista = lista_up + ['..'] + lista_mid + ['...'] + lista_dw
					df4 = df4.set_index(pd.Series(lista))

				#----------- heatmap ----------------------
				fig = go.Figure(data=go.Heatmap(
					z=df4.values,
					x=list(df4.columns),
					y=df4.index, zmin=zmin, zmax=zmax, colorscale='RdBu_r'))

				title = self.get_title_heatmap_set_pathway(unique_gender, pathway_id, pathway,
														   title_head, s_filtered)
				
				nrows = len(df4)
				height = row_height * nrows + header_height

				fig.update_layout(
							autosize=True,
							title=title,
							width=width,
							height=height,
							plot_bgcolor=plot_bgcolor,
							xaxis_title="cases",
							yaxis_title=f"{self.s_deg_dap}s",
							showlegend=False,
							margin=margin,
							font=dict(
								family="Arial",
								size=14,
								color="Black"
							)
				)

				icount += 1
				dic_fig[icount] = {}
				dic2 = dic_fig[icount]

				dic2['pathway_id'] = pathway_id
				dic2['pathway'] = pathway
				dic2['zlim'] = zlim
				dic2['fig'] = fig
			
	
				figname0 = title_replace(title)
				if len(figname0) > 250:
					figname0 = figname0[:250]
				figname = osjoin(self.root_figure, figname0+'.html')
				if verbose: print(">>> HTML and png saved:", figname)

				try:
					fig.write_html(figname)
					fig.write_image(figname.replace('.html', '.png'))

					figname = osjoin(self.root_fig_md, figname0+'.png')
					fig.write_image(figname)

				except:
					print(f"Error writing {figname0}")

		return dic_fig

	def degs_to_text_all_cases_summary(self, force:bool=False, verbose:bool=False) -> str:

		fname = self.fname_big_summary_txt%(self.disease, self.geneset_num, self.geneset_lib)
		filename = osjoin(self.root_ressum, fname)

		if exists(filename) and not force:
			text = read_txt(fname, self.root_ressum, verbose=verbose)
			return text

		now = datetime.now()

		s_date_time = now.strftime("%d/%m/%Y %H:%M:%S")

		text  = echo_print(f"# Project {self.project} for {self.disease}: {self.s_omics}")
		text += echo_print(f"# Calculated in {s_date_time}")
		text += echo_print(f"# geneset {self.geneset_num} - {self.geneset_lib}\n")

		for case in self.case_list:
			ret, _, _, _ = self.open_case(case, verbose=False)
			if not ret:
				continue

			stri = f"{self.icase}) {case} ('{self.translate_case(self.case)}') \n"
			text += echo_print(stri)

			text += echo_print(f"BCA's LFC {self.s_deg_dap} cutoffs: abs(lfc)={self.LFC_cut:.3f}; fdr={self.lfc_FDR_cut:.3f}")
			text += echo_print("")

			text += echo_print(f"\tThere are a total of {self.n_degs} {self.s_deg_dap}s, and {self.n_degs_ensembl} have ensembl_id")
			text += echo_print("")
			text += echo_print(f"\t\tAll {self.s_deg_dap}s Up ({self.n_degs_up}): {', '.join(self.degs_up)}")
			text += echo_print("")
			text += echo_print(f"\t\tAll {self.s_deg_dap}s Dw ({self.n_degs_dw}): {', '.join(self.degs_dw)}")

			text += echo_print("\n")
			text += echo_print(f"\tWith ensembl_id ({self.n_degs_ensembl})")
			text += echo_print("")
			text += echo_print(f"\t\t{self.s_deg_dap}s Up ({self.n_degs_up_ensembl}): {', '.join(self.degs_up_ensembl)}")
			text += echo_print("")
			text += echo_print(f"\t\t{self.s_deg_dap}s Dw ({self.n_degs_dw_ensembl}): {', '.join(self.degs_dw_ensembl)}")


			text += echo_print("\n")
			text += echo_print(f"\tWithout ensembl_id ({self.n_degs_not_ensembl})")
			text += echo_print("")
			text += echo_print(f"\t\t{self.s_deg_dap}s Up ({self.n_degs_up_not_ensembl}): {', '.join(self.degs_up_not_ensembl)}")
			text += echo_print("")
			text += echo_print(f"\t\t{self.s_deg_dap}s Dw ({self.n_degs_dw_not_ensembl}): {', '.join(self.degs_dw_not_ensembl)}")
			text += echo_print("\n")

			if self.df_enr is not None and not self.df_enr.empty:
				text += echo_print("")

				text += echo_print(f"\tBCA's pathway cutoffs: pval={self.ptw_pval_cut:.3f} fdr={self.ptw_FDR_cut:.3f} minimum number of genes={self.ptw_min_num_of_degs_cut}")
				text += echo_print(f"\tFound {self.n_pathways} (best={self.n_pathways_best}) pathways according to '{self.geneset_lib}'\n")

				for iptw in range(len(self.pathway_list)):
					text += echo_print(f"{iptw+1}) {self.pathway_id_list[iptw]:18} FDR={self.pathway_fdr_list[iptw]:.2e} - {self.pathway_list[iptw]} (# {self.num_of_genes_in_pathway_list[iptw]})")

				text += echo_print("")
				text += echo_print(f"\t{self.s_deg_dap}s found in enriched pathways:")
				text += echo_print(f"\t\t{self.n_degs_in_pathways} (best={self.n_degs_in_pathways_best}) in pathways {self.s_deg_dap}s; out of {self.n_degs} {self.s_deg_dap}s")
				text += echo_print("")

				text += echo_print(f"\t\t{self.s_deg_dap}s in pathways - possibly only {self.s_deg_dap}s with ensembl_id are annotated in pathways")
				text += echo_print(f"\t\tUp and Dw {self.s_deg_dap}s in pathways:")
				text += echo_print("")
				text += echo_print(f"\t\t\tUp {self.n_degs_up_ensembl_in_pathways}: {', '.join(self.degs_up_ensembl_in_pathways)}")
				text += echo_print("")
				text += echo_print(f"\t\t\tDw {self.n_degs_dw_ensembl_in_pathways}: {', '.join(self.degs_dw_ensembl_in_pathways)}")
				text += echo_print("\n")

				text += echo_print(f"\t\tUp and Dw {self.s_deg_dap}s not in pathways:")
				text += echo_print(f"\t\t\tAll {self.s_deg_dap}s not in pathways {self.n_degs_not_in_pathways}")
				text += echo_print(f"\t\t\t{self.s_deg_dap}s with ensembl_id not in pathways {self.n_degs_ensembl_not_in_pathways}")
				text += echo_print("")
				text += echo_print(f"\t\t\tUp ensembl {self.n_degs_up_ensembl_not_in_pathways}: {', '.join(self.degs_up_ensembl_not_in_pathways)}")
				text += echo_print("")
				text += echo_print(f"\t\t\tDw ensembl {self.n_degs_dw_ensembl_not_in_pathways}: {', '.join(self.degs_dw_ensembl_not_in_pathways)}")
				text += echo_print("")

				if self.n_degs == self.n_degs_in_pathways + self.n_degs_not_in_pathways:
					# stri = "All right:"
					# text += echo_print(f"\t\t{stri} {self.n_degs} all {self.s_deg_dap}s == {self.n_degs_in_pathways} in + {self.n_degs_not_in_pathways:} not in pathways")
					pass
				else:
					stri = "There are problems:"
					text += echo_print(f"\t\t{stri} {self.n_degs} all {self.s_deg_dap}s != {self.n_degs_in_pathways} in + {self.n_degs_not_in_pathways:} not in pathways")

			else:
				text += echo_print("\nNo enrichment analysis was found.")

			text += echo_print("")

		write_txt(text, fname, self.root_ressum, verbose=verbose)

		return text


	def plot_cutoff_simulation_histograms(self, col:str='toi1_median', width:int=1100, height:int=700):

		lista=[]
		local_case_list = []
		for case in self.case_list:
			dfi = self.calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(case)
			if dfi is None or dfi.empty:
				continue
			median = np.median(dfi[col])
			mean = np.mean(dfi[col])
			std = np.std(dfi[col])

			stri = f'{median:.1f}/{mean:.1f} ({std:.1f})'
			lista.append(f"{case}<br>{stri}")
			local_case_list.append(case)

		_n_rows = int(np.ceil(len(local_case_list)/4))
		fig = make_subplots(rows=_n_rows, cols=4, subplot_titles=local_case_list)

		nrow=1
		ncol=0

		for case in local_case_list:
			dfi = self.calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(case)
			if dfi is None or dfi.empty:
				continue

			ncol += 1
			if ncol == 5:
				ncol = 1
				nrow += 1

			fig.add_trace(go.Histogram(x=dfi[col], name=case), row=nrow, col=ncol)

		fig.update_layout(title=f"{self.project} - frequency of {col} in pathways'",
						  width=width,
						  height=height*_n_rows,
						  xaxis_title=col,
						  yaxis_title="frequency",
						  showlegend=False)

		return fig


	def plot_all_degs_up_down_per_cutoffs(self, width:int=1100, height:int=400,
										  title:str='', color_up:str='darkred', color_dw:str='navy',
										  plot_bgcolor:str='whitesmoke', y_anchor:float=1.01, verbose:bool=False):

		dfsim = self.open_simulation_table()

		if title is None or title == '':
			title=f"{self.project} - Number of Up/Downregulated {self.s_deg_dap}s per LFC-FDR cutoffs"

		cases = self.case_list

		fig = make_subplots(rows=len(cases), cols=1, subplot_titles=cases)

		row = 0
		col = 1
		for case in cases:
			row+=1

			showlegend = row == 1

			df2 = dfsim[dfsim.case == case].copy()
			df2.reset_index(inplace=True, drop=True)

			fig.add_trace(go.Bar(x=df2.cutoff, y=df2.n_degs_up, marker={'color':color_up}, name='Up',   showlegend=showlegend), row=row, col=col)
			fig.add_trace(go.Bar(x=df2.cutoff, y=df2.n_degs_dw, marker={'color':color_dw}, name='Down', showlegend=showlegend), row=row, col=col)


		yaxis_title = f"number of Up/Down {self.s_deg_dap}s"

		fig.update_layout(title=title,
						  width=width,
						  height=height*len(cases),
						  xaxis_title="",
						  xaxis2_title="LFC-FDR cutoffs",
						  yaxis_title=yaxis_title,
						  yaxis2_title=yaxis_title,
						  legend_title=f"{self.s_deg_dap}s regulation",
						  legend=dict( orientation="h",  yanchor="top", y=y_anchor, xanchor="right", x=1),
						  plot_bgcolor=plot_bgcolor,
						  showlegend=True)

		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		return fig


	def calc_all_annoted_genes(self):
		if self.geneset_num != 0:
			print("Only defined for Reactome - geneset_num == 0")
			return False

		df_enr = self.merge_reactome(self.df_enr)

		annotated_genes = []
		for i in range(len(df_enr)):
			symbs = df_enr.iloc[i].genes_pathway
			symbs = eval(symbs) if isinstance(symbs, str) else list(symbs)
			annotated_genes += symbs

		self.annotated_genes = list(np.unique(annotated_genes))
		self.n_annotated_genes = len(self.annotated_genes)

		return self.n_annotated_genes > 0


	def calc_all_enrichment_cutoff_params_and_ndxs_and_geneset_lib(self, force:bool=False,
																   verbose:bool=False) -> bool:
		
		print(">>> geneset_lib", self.geneset_lib)
		ret = True
		for case in self.case_list:
			print(">>>", case)
			dfa = self.calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(case, force=force, verbose=verbose)
			if dfa.empty:
				print(f"Error: no results for {case}")
				ret = False

		return ret
					
		
	''' old: calc_best_cutoff_parameters_by_case_geneset '''
	def calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(self, case:str, force:bool=False,
																	   verbose:bool=False) -> pd.DataFrame:
		icount=-1

		fname = self.fname_enr_gene_stat%(case, self.geneset_lib, self.normalization)
		filename = osjoin(self.root_ressum, fname)

		if exists(filename) and not force:
			dfi = pdreadcsv(fname, self.root_ressum)
			self.dfi = dfi
			return dfi

		# ret_lfc, degs, dfdegs = self.open_case(case, verbose=verbose)
		_, _, _, _ = self.open_case(case, verbose=verbose)
		if verbose: print(self.geneset_lib, case, self.normalization)

		files = [x for x in os.listdir(self.root_enrich) if '_pathway_pval_' in x \
				 and case in x and self.geneset_lib in x and self.normalization in x and not '~lock' in x]
		if files == []:
			if verbose: print("No files found in calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib()")
			return pd.DataFrame()

		if verbose: print(f"For case {case} there are {len(files)} files")

		dic = {}
		icount = -1
		for fname_enr in files:
			df_enr = pdreadcsv(fname_enr, self.root_enrich)

			if df_enr is None:
				print(f"Could not read: '{fname_enr}'")
				continue

			if len(df_enr) < 3:
				if verbose: print(f"Warning: there is no sufficient data (#df_enr) in: '{fname_enr}'")
				continue

			'''
			name = 'enricher_Reactome_2022_medulloblastoma_microarray_for_WNT_x_ctrl_not_normalized_cutoff_lfc_0.950_fdr_0.200_pathway_pval_0.050_fdr_0.450_num_genes_3.tsv'
			np.array(name.split('_'))
			'''
			mat = fname_enr.split('_')
			try:
				LFC_cut = float(mat[-11])
				lfc_FDR_cut = float(mat[-9])
				ptw_pval_cut = float(mat[-6])
				ptw_FDR_cut = float(mat[-4])
				ptw_min_num_of_degs_cut = int(mat[-1][:-4])
			except:
				print(f"Error: in {case} split={str(mat)}")
				return pd.DataFrame()

			self.df_enr_reactome = None
			df_enr = self.merge_reactome(df_enr)

			toi1_list, toi2_list, toi3_list, toi4_list = [], [], [], []
			all_genes, n_degs_in_pathways_list, all_genes_annotatted_in_pathway = [], [], []
			for i in range(len(df_enr)):
				row		 = df_enr.iloc[i]
				ptw_FDR_cut = row.fdr
				# pathway	 = row.pathway
				# pathway_id  = row.pathway_id

				genes = row.genes if isinstance(row.genes, list) else eval(row.genes)
				n_degs_in_pathways = len(genes)

				genes_annotatted_in_pathway = row.genes_pathway if isinstance(row.genes_pathway, list) else eval(row.genes_pathway)

				n_degs_in_pathways_list.append(n_degs_in_pathways)
				all_genes += list(genes)
				all_genes_annotatted_in_pathway += [str(x) for x in genes_annotatted_in_pathway]

				n_genes_annotatted_in_pathway = row.ngenes_pathway
				ratio = n_degs_in_pathways / n_genes_annotatted_in_pathway

				toi1 = np.sqrt(-np.log10(ptw_FDR_cut) * ratio)
				toi2 = np.sqrt(-np.log10(lfc_FDR_cut) * -np.log10(ptw_FDR_cut))
				toi3 = (-np.log10(lfc_FDR_cut) * -np.log10(ptw_FDR_cut) * ratio) ** 1/3.
				toi4 = (LFC_cut * -np.log10(lfc_FDR_cut) * -np.log10(ptw_FDR_cut) * ratio) ** 1/4.

				toi1_list.append(toi1)
				toi2_list.append(toi2)
				toi3_list.append(toi3)
				toi4_list.append(toi4)

			n_pathways = len(df_enr)
			toi1_mean = np.mean(toi1_list)
			toi1_median = np.median(toi1_list)
			toi1_std = np.std(toi1_list)

			toi2_mean = np.mean(toi2_list)
			toi2_median = np.median(toi2_list)
			toi2_std = np.std(toi2_list)

			toi3_mean = np.mean(toi3_list)
			toi3_median = np.median(toi3_list)
			toi3_std = np.std(toi3_list)

			toi4_mean = np.mean(toi4_list)
			toi4_median = np.median(toi4_list)
			toi4_std = np.std(toi4_list)

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			all_genes = list(np.unique(all_genes))
			n_degs_in_pathways = len(all_genes)
			all_genes_annotatted_in_pathway = list(np.unique(all_genes_annotatted_in_pathway))

			dic2['case'] = case
			dic2['cutoff'] = f'{LFC_cut:.3f} - {lfc_FDR_cut:.3f}'
			dic2['LFC_cut'] = LFC_cut
			dic2['lfc_FDR_cut'] = lfc_FDR_cut
			dic2['ptw_pval_cut'] = ptw_pval_cut
			dic2['ptw_FDR_cut'] = ptw_FDR_cut
			dic2['ptw_min_num_of_degs_cut'] = ptw_min_num_of_degs_cut

			dic2['n_pathways'] = n_pathways
			dic2['all_genes_annotatted_in_pathway'] = len(all_genes_annotatted_in_pathway)
			dic2['n_degs_in_pathways'] = n_degs_in_pathways

			dic2['n_degs_in_pathways_mean']   = np.mean(n_degs_in_pathways_list)
			dic2['n_degs_in_pathways_median'] = np.median(n_degs_in_pathways_list)
			dic2['n_degs_in_pathways_std']	= np.std(n_degs_in_pathways_list)

			dic2['toi1_mean'] = toi1_mean
			dic2['toi1_median'] = toi1_median
			dic2['toi1_std'] = toi1_std

			dic2['toi2_mean'] = toi2_mean
			dic2['toi2_median'] = toi2_median
			dic2['toi2_std'] = toi2_std

			dic2['toi3_mean'] = toi3_mean
			dic2['toi3_median'] = toi3_median
			dic2['toi3_std'] = toi3_std

			dic2['toi4_mean'] = toi4_mean
			dic2['toi4_median'] = toi4_median
			dic2['toi4_std'] = toi4_std

			dic2['all_genes'] = all_genes
			dic2['all_genes_annotatted_in_pathway'] = all_genes_annotatted_in_pathway

		dfi = pd.DataFrame(dic).T
		dfi = dfi.sort_values(['case', 'ptw_FDR_cut', 'lfc_FDR_cut', 'LFC_cut',], ascending=[True, True, True, False])

		_ = pdwritecsv(dfi, fname, self.root_ressum, verbose=verbose)
		self.dfi = dfi

		return dfi


	def best_cutoff_quantiles(self, case:str, selected_col:str='toi1_median', med_max_ptw:str='median',
							  verbose:bool=False) -> dict:
		''' best_cutoff_quantiles calculates the quantiles related to each toi = selected_toi_col '''

		self.med_max_ptw = med_max_ptw
		dfi = self.calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(case, verbose=verbose)

		''' impossible to calculate quantiles with less than 3 values '''
		if dfi is None or len(dfi) < 3:
			return {}

		quantile_vals = np.quantile(dfi[selected_col], self.quantile_list)

		dic = {}
		for i in range(len(self.quantile_list)):
			quantile = self.quantile_list[i]

			if quantile < 0.95:
				lim_inf = quantile_vals[i]
				lim_sup = quantile_vals[i+1]

				df2 = dfi[ (dfi[selected_col] > lim_inf) & (dfi[selected_col] <= lim_sup) ].copy()
			else:
				lim_inf = quantile_vals[i]
				lim_sup = None
				df2 = dfi[ (dfi[selected_col] > lim_inf) ].copy()

			if df2.empty:
				if verbose: print(f"Warning: empty quantile: {quantile} and case {self.case} - {self.normalization}")

				mat = ['lim_inf', 'lim_sup', 'LFC_cut', 'lfc_FDR_cut', 'ptw_pval_cut',
					   'ptw_FDR_cut', 'ptw_min_num_of_degs_cut', 'n_pathways', 'n_degs_in_pathways',
					   'n_degs_in_pathways_mean', 'n_degs_in_pathways_median ', 'n_degs_in_pathways_std',
					   'toi1_mean', 'toi1_median', 'toi1_std',
					   'toi2_mean', 'toi2_median', 'toi2_std',
					   'toi3_mean', 'toi3_median', 'toi3_std',
					   'toi4_mean', 'toi4_median', 'toi4_std'  ]

				dic[quantile] = [None] * len(mat)
			else:

				if med_max_ptw == 'pathway':
					df2 = df2.sort_values(['n_pathways', 'n_degs_in_pathways'], ascending=[False, False])
					df2.reset_index(inplace=True, drop=True)
					row = df2.iloc[0]
				else:
					df2 = df2.sort_values(selected_col, ascending=False)
					if med_max_ptw == 'median':
						med_max = df2[selected_col].median()
					else:
						''' maximum '''
						med_max = df2[selected_col].max()

					'''------ filter med_max, if there are more than 1 ro2 -----'''
					df3 = df2[ (df2[selected_col] == med_max)].copy()
					df3 = df3.sort_values(['n_pathways', 'n_degs_in_pathways'], ascending=[False, False])
					df3.reset_index(inplace=True, drop=True)
					row = df2.iloc[0]

				dic[quantile] = [lim_inf, lim_sup, row.LFC_cut, row.lfc_FDR_cut, row.ptw_pval_cut,
								 row.ptw_FDR_cut, row.ptw_min_num_of_degs_cut, row.n_pathways, row.n_degs_in_pathways,
								 row.n_degs_in_pathways_mean, row.n_degs_in_pathways_median, row.n_degs_in_pathways_std,
								 row.toi1_mean, row.toi1_median, row.toi1_std,
								 row.toi2_mean, row.toi2_median, row.toi2_std,
								 row.toi3_mean, row.toi3_median, row.toi3_std,
								 row.toi4_mean, row.toi4_median, row.toi4_std  ]

				if verbose:
					print(f"For case {case} {row.n_degs_in_pathways_median:.1f} median genes for len={len(df2)}")
					print(f"\tthe best cutoff LFC = {row.LFC_cut:.3f} and FDR = {row.lfc_FDR_cut:.3f} and Pathway FDR = {row.ptw_FDR_cut:.3f}")

		return dic

	def test_build_all_cutoffs_table(self, case:str='WNT',  selected_toi_col:str='toi4_median',
									 med_max_ptw:str='median', lista_fdr:List=[0.3, 0.35, 0.4],
									 verbose:bool=False):

		dfcut = self.build_all_cutoffs_table(selected_toi_col, verbose=verbose)

		df2 = dfcut[(dfcut.case == case) & (dfcut.med_max_ptw == med_max_ptw) & (dfcut.lfc_FDR_cut.isin(lista_fdr))]
		df2 = df2.sort_values(selected_toi_col, ascending=False)
		return df2

	def build_all_cutoffs_table(self, selected_toi_col:str='toi4_median',
								force:bool=False, verbose:bool=False) -> pd.DataFrame:

		# self.fname_cutoff_table = 'cutoff_table_all_cases_quantiles_for_col_%s_for_geneset_%d_%s_%s.tsv'
		fname = self.fname_cutoff_table%(selected_toi_col, self.geneset_num, self.geneset_lib, self.normalization)
		filename = osjoin(self.root_ressum, fname)

		if exists(filename) and not force:
				return pdreadcsv(fname, self.root_ressum, verbose=verbose)

		dic={}; icount=-1
		for case in self.case_list:
			'''
			 dictionary with following values:

			 quantile, LFC_cut, lfc_FDR_cut, ptw_pval_cut, ptw_FDR_cut,
			 ptw_min_num_of_degs_cut, n_pathways, n_degs_in_pathways, n_degs_in_pathways_mean, n_degs_in_pathways_median,
			 n_degs_in_pathways_std, toi1_mean, toi1_median, toi1_std, 2, 3, 4
			'''

			for med_max_ptw in ['median', 'maximum', 'pathway']:
				print(f">>> case {case} column {selected_toi_col} - {med_max_ptw}")
				self.med_max_ptw = med_max_ptw

				''' best_cutoff_quantiles calculates the quantiles related to each toi = selected_toi_col '''
				dic_quant = self.best_cutoff_quantiles(case, selected_toi_col, med_max_ptw, verbose=verbose)

				if dic_quant is None or len(dic_quant) == 0:
					print(f"Error: nohting found for {case}, {med_max_ptw}. and {selected_toi_col}")
					continue

				i=-1
				for quantile, (lim_inf, lim_sup, LFC_cut, lfc_FDR_cut,
							   ptw_pval_cut, ptw_FDR_cut,
							   ptw_min_num_of_degs_cut, n_pathways, n_degs_in_pathways,
							   n_degs_in_pathways_mean, n_degs_in_pathways_median, n_degs_in_pathways_std,
							   toi1_mean, toi1_median, toi1_std,
							   toi2_mean, toi2_median, toi2_std,
							   toi3_mean, toi3_median, toi3_std,
							   toi4_mean, toi4_median, toi4_std) in dic_quant.items():

					icount += 1
					i += 1
					dic[icount] = {}
					dic2 = dic[icount]

					dic2['case'] = case
					dic2['geneset_num'] = self.geneset_num
					dic2['normalization'] = self.normalization
					dic2['med_max_ptw'] = med_max_ptw
					dic2['parameter'] = selected_toi_col
					dic2['quantile'] = quantile
					dic2['quantile_val_inf'] = lim_inf
					dic2['quantile_val_sup'] = lim_sup

					dic2['LFC_cut'] = LFC_cut
					dic2['lfc_FDR_cut'] = lfc_FDR_cut
					dic2['ptw_pval_cut'] = ptw_pval_cut
					dic2['ptw_FDR_cut']  = ptw_FDR_cut
					dic2['ptw_min_num_of_degs_cut'] = ptw_min_num_of_degs_cut
					dic2['n_pathways'] = n_pathways
					dic2['n_degs_in_pathways'] = n_degs_in_pathways
					dic2['n_degs_in_pathways_mean'] = n_degs_in_pathways_mean
					dic2['n_degs_in_pathways_median'] = n_degs_in_pathways_median
					dic2['n_degs_in_pathways_std'] = n_degs_in_pathways_std

					dic2['toi1_mean'] = toi1_mean
					dic2['toi1_median'] = toi1_median
					dic2['toi1_std'] = toi1_std

					dic2['toi2_mean'] = toi2_mean
					dic2['toi2_median'] = toi2_median
					dic2['toi2_std'] = toi2_std

					dic2['toi3_mean'] = toi3_mean
					dic2['toi3_median'] = toi3_median
					dic2['toi3_std'] = toi3_std

					dic2['toi4_mean'] = toi4_mean
					dic2['toi4_median'] = toi4_median
					dic2['toi4_std'] = toi4_std

		dfcut = pd.DataFrame(dic).T
		_ = pdwritecsv(dfcut, fname, self.root_ressum, verbose=verbose)

		return dfcut

	def calc_best_cutoffs_params(self, selected_toi_col:str='toi4_median', n_best_sample:int=3,
								 save_config:bool=False, verbose:bool=False) -> pd.DataFrame:

		fname = self.fname_cutoff_table%(selected_toi_col, self.geneset_num, self.geneset_lib, self.normalization)
		filename = osjoin(self.root_ressum, fname)

		if not exists(filename):
			print(f"Could not find: '{filename}'")
			return pd.DataFrame()

		dfcut = pdreadcsv(fname, self.root_ressum, verbose=verbose)

		if dfcut is None or dfcut.empty:
			print(f"Error: could not read: '{filename}'")
			return pd.DataFrame()

		df_all_fdr = self.open_all_fdr_lfc_correlation()

		df_list = []

		for case in self.case_list:
			df_fdr = df_all_fdr[df_all_fdr.case == case]
			if df_fdr.empty:
				print(f"Error: no correlation was calculated for case '{case}'")
				raise Exception('stop: plot_degs_in_pathways_vs_toi_per_case()')

			for med_max_ptw in ['median', 'maximum', 'pathway']:

				df2 = dfcut[(dfcut.case == case) & (dfcut.med_max_ptw == med_max_ptw) &
							(dfcut.lfc_FDR_cut.isin(df_fdr.fdr))].copy()

				if df2.empty: continue

				if med_max_ptw == 'pathway':
					df2 = df2.sort_values(['n_pathways', 'n_degs_in_pathways'], ascending=[False, False])
				else:
					df2 = df2.sort_values(selected_toi_col, ascending=False)

				df2.reset_index(inplace=True, drop=True)

				if len(df2) >= n_best_sample:
					dfa = pd.DataFrame(df2.iloc[n_best_sample-1]).T
				else:
					dfa = pd.DataFrame(df2.iloc[len(df2)-1]).T

				df_list.append(dfa)

		if df_list == []:
			dfconfig = pd.DataFrame()
		else:
			dfconfig = pd.concat(df_list)
			dfconfig.reset_index(inplace=True, drop=True)

			if save_config:
				self.cfg.save_best_ptw_cutoff(dfconfig, verbose=True)
				print(f"Chosen n_best_sample: {n_best_sample}")

		return dfconfig

	def display_best_cutoff_params(self, npoints:int=10, selected_toi_col:str='toi4_median', 
								   med_max_ptw:str='median') -> pd.DataFrame:
		
		dic = {}; icount=-1
		for ipoint in range(1,npoints+1):
			dfconfig = self.calc_best_cutoffs_params(selected_toi_col=selected_toi_col, n_best_sample=ipoint, save_config=False, verbose=False)
			if dfconfig is None or dfconfig.empty:
				continue
			dfconfig = dfconfig[dfconfig.med_max_ptw == med_max_ptw]

			for j in range(len(dfconfig)):
				row = dfconfig.iloc[j]
				
				icount += 1
				dic[icount] = {}
				dic2 = dic[icount]
				dic2['case'] = row.case
				dic2['ipoint'] = ipoint
					
				dic2[f'n_ptws'] = row.n_pathways
				dic2[f'n_degs'] = row.n_degs_in_pathways

				dic2[f'lfc_cut'] = row.LFC_cut
				dic2[f'lfc_FDR_cut'] = row.lfc_FDR_cut
				dic2[f'ptw_FDR_cut'] = row.ptw_FDR_cut

				dic2[f'toi4_median'] = row.toi4_median
				dic2[f'toi3_median'] = row.toi3_median
				dic2[f'toi2_median'] = row.toi2_median
				dic2[f'toi1_median'] = row.toi1_median

		dfpoints = pd.DataFrame(dic).T
		dfpoints = dfpoints.sort_values(['case', 'ipoint'])
		dfpoints.reset_index(inplace=True, drop=True)

		return dfpoints


	def show_multiple_best_cutoffs_params(self, selected_toi_col:str='toi4_median', n_best_sample_list:List=[],
								 save_config:bool=False, verbose:bool=False) -> pd.DataFrame:

		if len(n_best_sample_list) != len(self.case_list):
			print("Error: n_best_sample_list is a list of same dimension of case_list.")
			return pd.DataFrame()

		fname = self.fname_cutoff_table%(selected_toi_col, self.geneset_num, self.geneset_lib, self.normalization)
		filename = osjoin(self.root_ressum, fname)

		if not exists(filename):
			print(f"Could not find: '{filename}'")
			return pd.DataFrame()

		dfcut = pdreadcsv(fname, self.root_ressum, verbose=verbose)

		if dfcut is None or dfcut.empty:
			print(f"Error: could not read: '{filename}'")
			return pd.DataFrame()

		df_all_fdr = self.open_all_fdr_lfc_correlation()

		df_list = []

		for i in range(len(self.case_list)):
			case = self.case_list[i]
			n_best_sample = n_best_sample_list[i]

			df_fdr = df_all_fdr[df_all_fdr.case == case]
			if df_fdr.empty:
				print(f"Error: no correlation was calculated for case '{case}'")
				raise Exception('stop: plot_degs_in_pathways_vs_toi_per_case()')

			for med_max_ptw in ['median', 'maximum', 'pathway']:

				df2 = dfcut[(dfcut.case == case) & (dfcut.med_max_ptw == med_max_ptw) &
							(dfcut.lfc_FDR_cut.isin(df_fdr.fdr))]

				if df2.empty: continue

				if med_max_ptw == 'pathway':
					df2 = df2.sort_values(['n_pathways', 'n_degs_in_pathways'], ascending=[False, False])
				else:
					df2 = df2.sort_values(selected_toi_col, ascending=False)
					''' get the first 3 (n_best_sample) rows
						there are cases that changing a little LFC_cut or lfc_FDR_cut
						the number of n_pathways increases!
					'''
					df2 = df2.iloc[:n_best_sample]
					df2 = df2.sort_values(['n_pathways', 'n_degs_in_pathways'], ascending=[False, False])

				dfa = pd.DataFrame(df2.iloc[0]).T
				df_list.append(dfa)

		if df_list == []:
			dfconfig = pd.DataFrame()
		else:
			dfconfig = pd.concat(df_list)
			dfconfig.reset_index(inplace=True, drop=True)

			if save_config:
				self.cfg.save_best_ptw_cutoff(dfconfig, verbose=True)
				print(f"Chosen n_best_sample: {n_best_sample_list}")

		return dfconfig


	def plot_genes_and_pathways_frequecies_per_cases(self, selected_toi_col:str='toi4_median',
													 med_max_ptw:str='median', width:int=800,
													 height:int=400, verbose:bool=False) -> List:

		dfcut = self.build_all_cutoffs_table(selected_toi_col, force=False, verbose=verbose)

		dfcut = dfcut.sort_values(['case','quantile'])
		dfcut.reset_index(inplace=True, drop=True)
		self.dfcut = dfcut

		if dfcut is None or dfcut.empty:
			return []

		fig_list = []
		for _plot in ['genes', 'pathways']:

			fig = go.Figure()

			for i in range(len(self.quantile_list)):
				quantile = self.quantile_list[i]
				name=f'{quantile}'
				color = self.my_colors[i]

				df2 = dfcut[ (dfcut['quantile'] == quantile) & (dfcut.med_max_ptw == med_max_ptw)].copy()
				df2.reset_index(inplace=True, drop=True)

				if _plot == 'genes':
					fig.add_trace(go.Bar(x=df2.case, y=df2.n_degs_in_pathways, marker_color=color, name=name)) # marker_color=color,
					plot_name = f"'Best' number of {self.s_deg_dap}s in pathways"
					yaxis_title = f"# {self.s_deg_dap}s in pathways"
				else:
					fig.add_trace(go.Bar(x=df2.case, y=df2.n_pathways, marker_color=color, name=name))
					plot_name = f"'Best' number of enriched pathways"
					yaxis_title = "# of pathways"

			title=f"{plot_name} for {selected_toi_col}  {med_max_ptw} values in quantiles"
			fig.update_layout(title=title,
							  width=width,
							  height=height,
							  xaxis_title='cases x quantiles',
							  yaxis_title=yaxis_title,
							  legend_title="Quantiles",
							  showlegend=True)

			figname = title_replace(title)
			figname = osjoin(self.root_figure, figname+'.html')

			fig.write_html(figname)
			if verbose: print(">>> HTML and png saved:", figname)
			fig.write_image(figname.replace('.html', '.png'))

			fig_list.append(fig)

		return fig_list

	def open_all_fdr_lfc_correlation(self, verbose:bool=False) -> pd.DataFrame:
		fname = self.fname_all_fdr_lfc_correlation%(self.LFC_cut_inf)
		filename = osjoin(self.root_config, fname)

		if exists(filename):
			df_all_fdr = pdreadcsv(fname, self.root_config, verbose=verbose)
		else:
			print(f"File all_fdr_lfc_correlation not found: {filename}")
			print("Please run: 'pubmed_taubate_new03_up_down_simulation'")
			df_all_fdr = pd.DataFrame()

		return df_all_fdr


	def open_fdr_lfc_correlation(self, case:str, LFC_cut_inf:float=2, verbose:bool=False) -> pd.DataFrame:

		if LFC_cut_inf is None:
			LFC_cut_inf = self.LFC_cut_inf

		fname = self.fname_dic_fdr_lfc_correlation%(case, LFC_cut_inf)
		filename = osjoin(self.root_config, fname)

		if exists(filename):
			dic = loaddic(fname, self.root_config, verbose=verbose)
		else:
			dic = None

		if dic is None or len(dic) == 0:
			print(f"Could not find {filename}, thus df_fdr is empty")
			df_fdr = pd.DataFrame()
		else:
			df_fdr = dic['df_fdr']

		return df_fdr

	def calc_all_LFC_FDR_cutoffs(self, cols2:List=['n_degs', 'LFC_cut'], corr_cutoff:float=-0.90, 
								 nregs_fdr:int=20, want_improve:bool=False,
								 method:str='spearman', force:bool=False, verbose:bool=False):

		self.df_all_fdr = None

		fname = self.fname_all_fdr_lfc_correlation%(self.LFC_cut_inf)
		filename = osjoin(self.root_config, fname)

		if exists(filename) and not force:
			self.df_all_fdr = pdreadcsv(fname, self.root_config, verbose=verbose)
			return self.df_all_fdr


		df_list = []
		for case in self.case_list:
			if verbose: print(">>>", case)

			ret, dic_return = self.calc_nDEG_curve_per_LFC_FDR(case=case, cols2=cols2, corr_cutoff=corr_cutoff, 
															   nregs_fdr=nregs_fdr, want_improve=want_improve,
															   method=method, force=force, verbose=verbose)

			if not ret:
				print(">>>", case, "could not calculate correlation table.")
				continue

			df_fdr = dic_return['df_fdr']
			df_list.append(df_fdr)


		if df_list == []:
			df_all_fdr = None
		else:
			df_all_fdr = pd.concat(df_list)
			self.df_all_fdr = df_all_fdr
			df_all_fdr.reset_index(inplace=True, drop=True)
			ret = pdwritecsv(df_all_fdr, fname, self.root_config, verbose=verbose)

		return df_all_fdr


	def plot_all_LFC_FDR_cutoffs(self, width:int=1100, height:int=700, title:str=None,
								 cols2:List=['n_degs', 'LFC_cut'],
								 corr_cutoff:float=-0.90, nregs_fdr:int=5, method:str='spearman',
								 verbose:bool=False) -> dict:

		dic_fig = {}
		for case in self.case_list:

			ret, dic_fig_return, _ = self.plot_nDEG_curve_per_LFC_FDR(case,
										  width=width, height=height, title=title, cols2=cols2,
										  corr_cutoff=corr_cutoff, nregs_fdr=nregs_fdr, method=method, verbose=verbose)

			if not ret:
				continue

			dic_fig[case] = {}
			dic2 = dic_fig[case]

			for key, fig in dic_fig_return.items():
				dic2[key] = fig

		return dic_fig


	'''
		calc_and_plot_nDEG_curve_per_LFC_FDR was splited in
			1) calc_nDEG_curve_per_LFC_FDR
				   return ret, dic

			2) plot_nDEG_curve_per_LFC_FDR
				   return ret, dic_fig, df_fdr
	'''

	def calc_nDEG_curve_per_LFC_FDR(self, case:str,
									cols2:List=['n_degs', 'LFC_cut'],
									corr_cutoff:float=-0.90, nregs_fdr:int=20, want_improve:bool=False,
									method:str='spearman', force:bool=False, verbose:bool=False) -> Tuple[bool, dict]:
		'''
			calc_nDEG_curve_per_LFC_FDR
				calculates lfc_FDR cut when Spearman's correlation (e.g.) <= -.90
				correlation must always decrease
				no more than ireg_fdr, exception repeated correlations

				the DAP x LFC x fdr lanscape looks very irregular
				at least 5 regs are necessary aftare found a corr <= corr_cutoff

				somentime the curves jumps! other, the look the same (superposition)
				this ocurrs in ou first studies os PBMC proteomics COVID-19 and microarray Medulloblastoma

				saves nregs_fdr values
				dic contains df_fdr, name_list, fdrs --> necessary do draw the plot

				return ret, dic

			calc per case, fdr, lfc, correlations
		'''
		dfsim = self.open_simulation_table()

		fname = self.fname_dic_fdr_lfc_correlation%(case, self.LFC_cut_inf)
		filename = osjoin(self.root_config, fname)

		if exists(filename) and not force:
			dic = loaddic(fname, self.root_config, verbose=verbose)
			return True, dic

		dfsim = dfsim[dfsim.case == case].copy()
		dfsim = dfsim.sort_values(['lfc_FDR_cut', 'LFC_cut'], ascending=[True, False])
		dfsim.reset_index(inplace=True, drop=True)

		fdrs = dfsim.lfc_FDR_cut.unique()

		name_list, fdr_list = [], []

		dic = {}; icount=-1; found=False; ireg_fdr=0; corr_previous=0
		for fdr in fdrs:
			dfsim2 = dfsim[ (dfsim.lfc_FDR_cut == fdr) & (dfsim.LFC_cut >= self.LFC_cut_inf) ]
			if len(dfsim2) < 3:
				continue

			corr = dfsim2[cols2].corr(method=method).iloc[0,1]

			n_degs_min = dfsim2.n_degs.min()
			dfsim2_min = dfsim2[dfsim2.n_degs == n_degs_min]

			n_degs_max = dfsim2.n_degs.max()
			dfsim2_max = dfsim2[dfsim2.n_degs == n_degs_max]

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			dic2['case'] = case
			dic2['fdr'] = fdr
			dic2['corr'] = corr

			''' correlation must be negative '''
			if pd.isnull(corr) or corr > corr_cutoff:
				name = f"fdr={fdr:.2f} not found corr."
				name_list.append(name)

				# print(">>> NOT", icount, case, fdr, n_degs_min, n_degs_max)
				dic2['first'] = False
				dic2['chosen'] = False

			else:
				''' correlation is negative, usually cor <= -0.90 '''
				if not found:
					found = True

					ireg_fdr = 1

					dic2['first'] = True
					dic2['chosen'] = True
					corr_previous = corr

					name = f"fdr={fdr:.2f} corr={corr:.3f} ***"

					if verbose: print(">>> ***", icount, case, fdr, n_degs_min, n_degs_max)
				else:

					''' correlation must improve '''
					if want_improve:
						if corr < corr_previous or corr <= -0.99:
							corr_previous = corr
							ireg_fdr += 1
					else:
						ireg_fdr += 1

					dic2['first'] = False
					dic2['chosen'] = True
					corr_previous = corr

					name = f"fdr={fdr:.2f} corr={corr:.3f}"


			dic2['label'] = name
			dic2['method'] = method

			dic2['n_degs_min'] = n_degs_min
			dic2['n_degs_max'] = n_degs_max

			dic2['LFC_cut_inf'] = self.LFC_cut_inf

			dic2['degs_min'] = dfsim2_min.iloc[0].degs
			dic2['degs_max'] = dfsim2_max.iloc[0].degs


			name_list.append(name)
			fdr_list.append(fdr)

			''' no more than ireg_fdr, exception repeated correlations'''
			if ireg_fdr >= nregs_fdr:
				break

		if len(dic) == 0:
			print("Nothing found in calc_nDEG_curve_per_LFC_FDR()")
			return False, {}

		if verbose:
			print(">>> ", case, icount, case, ireg_fdr, nregs_fdr)

		df_fdr = pd.DataFrame(dic).T

		dic_return = {}
		dic_return['df_fdr']	= df_fdr
		dic_return['name_list'] = name_list
		dic_return['fdrs']	  = fdr_list

		_ = dumpdic(dic_return, fname, self.root_config, verbose=verbose)

		return True, dic_return

	def plot_all_dfsim(self, dfsim:pd.DataFrame, case:str, fdr_list:List, which_list:list=['deg', 'up', 'down'],
					   width:int=1100, height:int=700, title:str='', verbose:bool=False):

		dfsim2 = dfsim[dfsim.case == case].copy()
		dfsim2.reset_index(inplace=True, drop=True)

		if title is None or title == '':
			title = f"{self.s_deg_dap}s curve per fdr x lfc for {case} - all FDR's"

		title2 = title

		yaxis_title = f"num of {self.s_deg_dap}s"
		yaxis_title2 = yaxis_title
		xaxis_title = "LFC_cut"
		legend_title = 'lfc_FDR cut'

		dic_fig = {}
		for which in which_list:
			''' many figures: one case, deg, up and down '''
			fig = go.Figure()

			for fdr in fdr_list:

				dfsim3 = dfsim2[dfsim2.lfc_FDR_cut == fdr]
				if dfsim3.empty:
					continue
			
				name = f'{fdr:.2f}'

				if which == 'deg':
					col = 'n_degs'
					yaxis_title2 = yaxis_title
					title2 = title
				else:
					col = 'n_degs_up' if which == 'up' else 'n_degs_dw'
					stri = 'Up' if which == 'up' else 'Down'
					yaxis_title2 = yaxis_title + ' ' + stri
					title2 = title.replace('s curve', f's {stri} curve')

				legend_title = 'lfc_FDR cut'

				fig.add_trace(go.Scatter(x=dfsim3.LFC_cut, y=dfsim3[col], name=name))


			fig.update_layout(
						autosize=True,
						title=title2,
						width=width,
						height=height,
						xaxis_title=xaxis_title,
						yaxis_title=yaxis_title2,
						legend_title=legend_title,
						showlegend=True,
						font=dict(
							family="Arial",
							size=14,
							color="Black"
						)
			)
			dic_fig[which] = fig

			figname = title_replace(title2)
			figname = osjoin(self.root_figure, figname+'.html')

			fig.write_html(figname)
			if verbose: print(">>> HTML and png saved:", figname)
			fig.write_image(figname.replace('.html', '.png'))

		return dic_fig


	def plot_nDEG_curve_per_LFC_FDR(self, case:str, which_list:list=['deg', 'up', 'down'],
									width:int=1100, height:int=700, title:str='',
									cols2:List=['n_degs', 'LFC_cut'], corr_cutoff:float=-0.90, 
									nregs_fdr:int=5, method:str='spearman', 
									verbose:bool=False) -> Tuple[bool, dict, pd.DataFrame]:

		'''
			plot_nDEG_curve_per_LFC_FDR

			call calc_nDEG_curve_per_LFC_FDR
					   return ret, dic

			return ret, dic_fig, df_fdr
		'''

		ret, dic_return = self.calc_nDEG_curve_per_LFC_FDR(case=case, cols2=cols2, 
														   corr_cutoff=corr_cutoff, nregs_fdr=nregs_fdr,
														   method=method, verbose=verbose)

		if not ret:
			print("Could not get data from calc_nDEG_curve_per_LFC_FDR()")
			return ret, {}, pd.DataFrame()


		df_fdr = dic_return['df_fdr']
		df_fdr = df_fdr[df_fdr.chosen == True]
		# name_list = dic_return['name_list']
		# fdrs	  = df_fdr.fdr.to_list()

		if title is None or title == '':
			title = f'{self.s_deg_dap}s curve per fdr x lfc for {case}<br>for corr_cutoff={corr_cutoff:.3f} and correlation dependes on LFC_cut_inf={self.LFC_cut_inf:.3f}'

		title2 = title

		yaxis_title = f"num of {self.s_deg_dap}s"
		yaxis_title2 = yaxis_title
		xaxis_title = "LFC_cut"
		legend_title = 'lfc_FDR cut'

		dic_fig = {}
		for which in which_list:
			''' many figures: one case, deg, up and down '''
			fig = go.Figure()

			for i, row in df_fdr.iterrows():
				fdr  = row['fdr']
				corr = row['corr']

				if which == 'deg':
					try:
						name = row.label
					except:
						name = f'{fdr:.2f} \t {corr:.3f}'
				else:
					name = f'{fdr:.2f}'


				dfsim = self.dfsim[self.dfsim.case == case]
				dfsim = dfsim[ (dfsim.lfc_FDR_cut == fdr) & (dfsim.LFC_cut >= self.LFC_cut_inf) ]

				if which == 'deg':
					col = 'n_degs'
					yaxis_title2 = yaxis_title
					title2 = title
					legend_title = 'lfc_FDR cut - Spearman corr.'
				else:
					col = 'n_degs_up' if which == 'up' else 'n_degs_dw'
					stri = 'Up' if which == 'up' else 'Down'
					yaxis_title2 = yaxis_title + ' ' + stri
					title2 = title.replace('s curve', f's {stri} curve')
					legend_title = 'lfc_FDR cut'

				fig.add_trace(go.Scatter(x=dfsim.LFC_cut, y=dfsim[col], name=name))


			fig.update_layout(
						autosize=True,
						title=title2,
						width=width,
						height=height,
						xaxis_title=xaxis_title,
						yaxis_title=yaxis_title2,
						legend_title=legend_title,
						showlegend=True,
						font=dict(
							family="Arial",
							size=14,
							color="Black"
						)
			)
			dic_fig[which] = fig

			figname = title_replace(title2)
			figname = osjoin(self.root_figure, figname+'.html')

			fig.write_html(figname)
			if verbose: print(">>> HTML and png saved:", figname)
			fig.write_image(figname.replace('.html', '.png'))

		return True, dic_fig, df_fdr


	def plot_toi_versus_genes_and_pathways(self, case:str, selected_toi_col:str, _plot:str='genes',
				width:int=1100, height:int=500, plot_all_dfi:bool=True,
				colors:List=['olivedrab', 'navy', 'red', 'darkcyan', 'darkgreen', 'orange', 'brown', 'darksalmon',
							 'magenta', 'darkturquoise', 'orange', 'darkred', 'indigo', 'magenta', 'maroon', 'black',
							 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgrey', 'darkgreen', 'navy'],
				verbose:bool=False) -> List:
		'''
			_plot = ['genes', 'pathways']
		'''

		dfcut = self.build_all_cutoffs_table(selected_toi_col, force=False, verbose=verbose)
		# fig = make_subplots(rows=2, cols=1, subplot_titles=['genes', 'pathways'])

		if _plot == 'genes':
			title=f"Best number of {self.s_deg_dap}s x {selected_toi_col} for {case}"
			yaxis_title = f"# {self.s_deg_dap}s in pathways"
		else:
			title = f"Best number of enriched pathways x {selected_toi_col} for {case}"
			yaxis_title = "# of pathways"

		fig = go.Figure()
		i = -1

		df2  = dfcut[ (dfcut.case == case) ]
		fdr_list = np.unique(df2.lfc_FDR_cut)

		if plot_all_dfi: 
			dfi = self.calc_enrichment_cutoff_params_and_ndxs_per_case_and_geneset_lib(case)
		else:
			dfi = pd.DataFrame()

		for fdr in fdr_list:
			i += 1

			df2  = dfcut[ (dfcut.case == case) & (dfcut.lfc_FDR_cut == fdr) ]
			if plot_all_dfi: 
				dfi2 = dfi[   (dfi.case == case)   & (dfi.lfc_FDR_cut == fdr) ]
			else:
				dfi2 = pd.DataFrame()

			if df2.empty:
				continue

			name1 = f"{case} fdr={fdr:.3f} for {_plot}"
			name2 = name1 + '_all'
			df2  = df2.sort_values( selected_toi_col, ascending=True)
			if plot_all_dfi: dfi2 = dfi2.sort_values(selected_toi_col, ascending=True)
			
			color = colors[i]

			if _plot == 'genes':
				fig.add_trace(go.Scatter(x=df2[selected_toi_col],   y=df2.n_degs_in_pathways,  marker_color=color, name=name1) ) # , row=1, col=1)
				if plot_all_dfi: fig.add_trace(go.Scatter(x=dfi2[selected_toi_col], y=dfi2.n_degs_in_pathways, line=dict(dash='dash'), marker_color=color, name=name2)) # , row=1, col=1)
			else:
				fig.add_trace(go.Scatter(x=df2[selected_toi_col],   y=df2.n_pathways, marker_color=color, name=name1)) #, row=2, col=1)
				if plot_all_dfi: fig.add_trace(go.Scatter(x=dfi2[selected_toi_col], y=dfi2.n_pathways, line=dict(dash='dash'), marker_color=color, name=name2)) #, row=2, col=1)

		fig.update_layout(title=title,
						  width=width,
						  height=height,
						  xaxis_title=selected_toi_col,
						  yaxis_title=yaxis_title,
						  # xaxis2_title=selected_toi_col,
						  legend_title="cases",
						  showlegend=True)

		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))

		#fig['layout']['yaxis']['title']=yaxis_title1
		# fig['layout']['yaxis2']['title']=yaxis_title2
		return fig


	def calc_min_max(self, vals):
		mini = np.min(vals)
		maxi = np.max(vals)

		if mini * maxi < 0:
			all_max0 = maxi if np.abs(maxi) >= np.abs(mini) else mini
			all_max0_abs = np.abs(all_max0)
			all_max = int(np.round(all_max0_abs, 0))

			if all_max < all_max0_abs:
				all_max += 1

			ret = [-all_max, +all_max]
		else:
			mini_abs = np.abs(mini)
			maxi_abs = np.abs(maxi)

			all_max0 = maxi_abs if maxi_abs >= mini_abs else mini_abs
			all_max = int(np.round(all_max0, 0))

			if all_max < all_max0:
				all_max += 1

			if mini < 0:
				ret = [-all_max, 0]
			else:
				ret = [0, all_max]

		return ret

	''' old: plot_all_lfc_path'''
	def plot_all_dpiv_lfc_path(self, pathway_and_id:List, width:int=500, height:int=300, 
							   margin:dict=dict(l=80, r=80, t=100, b=80), verbose:bool=False) -> dict:

		# pathway, pathway_id = pathway_and_id
		dfpiv, _  = self.calc_pivot_one_pathway_LFC(pathway_and_id, verbose=verbose)

		dic = {}
		for i in range(len(dfpiv)):
			row  = dfpiv.iloc[i]
			symbol = str(row.name)
			vals_female = list(row[self.group_female_list])
			vals_male   = list(row[self.group_male_list])
			fig = self.plot_lfc_path(vals_female, vals_male, title=symbol, 
							width=width, height=height, margin=margin, verbose=verbose)

			dic[symbol] = fig

		return dic


	def plot_genes_lfc_path(self, df:pd.DataFrame, want_female:bool=True, want_male:bool=True, 
						 	diff:float=1, corr_cutoff:float=0.6, 
							high_correlation:bool=True, similar_modulation_for_both:bool=False,
							width:int=500, height:int=300, 
							margin:dict=dict(l=80, r=80, t=100, b=80), verbose:bool=False) -> dict:
		
		title0 = f"LFC path for {self.gene_protein} '%s'<br>"

		''' removing multiple pathways - 'pathway_id', 'pathway','''
		cols = ['symbol', 'g2a_female', 'g2a_male', 'g2b_female', 'g2b_male', 'g3_female_adult',
				'g3_female_elder', 'g3_male_adult', 'g3_male_elder', 'diff_female', 'abs_diff_female',
				'diff_male', 'abs_diff_male', 'diff', 'abs_diff', 'spearman_cor_female', 'spearman_cor_male']
		
		df = df[cols]
		df = df.drop_duplicates('symbol')
	
		if want_female and not want_male:
			df2 = df[ (df.abs_diff_female >= diff) & (df.abs_diff_male < diff) ]
			title0 += f"female modulation diff >= {diff:.1f} and male modulation diff < {diff:.1f}"

			if high_correlation:
				title0 += f"<br>and female correlation >= {corr_cutoff:.1f}"
				df2 = df2[ (df2.spearman_cor_female >= corr_cutoff) | (df2.spearman_cor_female <= -corr_cutoff) ]
			else:
				title0 += f"<br>and female correlation < {corr_cutoff:.1f}"
				df2 = df2[ (df2.spearman_cor_female > -corr_cutoff) &  (df2.spearman_cor_female < corr_cutoff) ]

		elif want_male and not want_female:
			df2 = df[ (df.abs_diff_male >= diff) & (df.abs_diff_female < diff) ]
			title0 += f"male modulation diff >= {diff:.1f} and female modulation diff < {diff:.1f}"

			if high_correlation:
				title0 += f"<br>and male correlation >= {corr_cutoff:.1f}"
				df2 = df2[ (df2.spearman_cor_male >= corr_cutoff) | (df2.spearman_cor_male <= -corr_cutoff) ]
			else:
				title0 += f"<br>and male correlation < {corr_cutoff:.1f}"
				df2 = df2[ (df2.spearman_cor_male > -corr_cutoff) &  (df2.spearman_cor_male < corr_cutoff) ]
		
		elif not want_male and not want_female:
			df2 = df[ (df.abs_diff_male < diff) & (df.abs_diff_female < diff) ]
			title0 += f"male modulation diff < {diff:.1f} and female modulation diff < {diff:.1f}"

			title0 += f"<br>any correlations"
		
		else:
			if similar_modulation_for_both:
				df2 = df[ (df.abs_diff_male < diff) & (df.abs_diff_female < diff) ]
				title0 += f"similar modulation for female and male - diffs < {diff:.1f}"
			else:
				df2 = df[ (df.abs_diff_male >= diff) & (df.abs_diff_female >= diff) ]
				title0 += f"female and male modulation diffs >= {diff:.1f}"

			if high_correlation:
				title0 += f"<br>and both high correlations >= {corr_cutoff:.1f}"
				df2 = df2[ ( (df2.spearman_cor_female >= corr_cutoff) | (df2.spearman_cor_female <= -corr_cutoff)) &
			               ( (df2.spearman_cor_male   >= corr_cutoff) | (df2.spearman_cor_male <= -corr_cutoff)) ]
			else:
				title0 += f"<br>and both low correlations < {corr_cutoff:.1f}"
				df2 = df2[ (df2.spearman_cor_female < corr_cutoff) &  (df2.spearman_cor_male < corr_cutoff) ]
				df2 = df2[ ( (df2.spearman_cor_female > -corr_cutoff) | (df2.spearman_cor_female < corr_cutoff)) &
			               ( (df2.spearman_cor_male   > -corr_cutoff) | (df2.spearman_cor_male   < corr_cutoff)) ]
		if df2.empty:
			if verbose: print(f"Warning: nothing found for '{title0}'")
			return {}
	
		symbol_list = df2.symbol.to_list()

		dic = {}; icount=-1
		for symbol in symbol_list:
			dfsymb  = df2[df2.symbol == symbol]
			if dfsymb.empty: continue

			row = dfsymb.iloc[0]

			vals_female = list(row[self.group_female_list])
			vals_male   = list(row[self.group_male_list])

			title = title0%(symbol)
			fig = self.plot_lfc_path(vals_female, vals_male, title, width=width, height=height, margin=margin, verbose=verbose)

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			dic2['want_female'] = want_female
			dic2['want_male'] = want_male
			dic2['high_correlation'] = high_correlation
			dic2['similar_modulation_for_both'] = similar_modulation_for_both
			dic2['symbol'] = symbol
			dic2['fig'] = fig
			dic2['title'] = title
			dic2['df'] = df2

		return dic

	def plot_lfc_path(self, vals_female:List, vals_male:List, title:str, 
				      width:int=500, height:int=300, 
					  margin:dict=dict(l=80, r=80, t=100, b=80), verbose:bool=False):
		
		groups = self.group_list
		colors = self.group_colors

		line_traces = []
		for i in range(len(groups)-1):
			name  = f"{groups[i]}-{groups[i+1]}"
			line_traces.append(go.Scatter(x=vals_female[i:(i+1)], y=vals_male[i:(i+1)],
										  marker_color=colors[i+1], name=name, mode="lines",  ))

		arrow_traces = []
		for i in range(len(vals_female) - 1):
			arrow_trace = go.layout.Annotation(
				dict( x=vals_female[i+1], y=vals_male[i+1],
					  xref="x", yref="y", text="",
					  showarrow=True, axref="x", ayref="y",
					  ax=vals_female[i],  ay=vals_male[i],
					  arrowhead=3,  arrowwidth=1, arrowcolor=colors[i+1],
				)
			)
			arrow_traces.append(arrow_trace)

		fig = go.Figure(data=line_traces)

		x_range = self.calc_min_max(vals_female)
		y_range = self.calc_min_max(vals_male)

		fig.update_layout(title=title,
					  width=width,
					  height=height,
					  annotations=arrow_traces,
					  xaxis_title="female",
					  yaxis_title="male",
					  showlegend=True,
					  xaxis_range = x_range,
					  yaxis_range = y_range,
					  margin=margin)

		figname0 = title_replace(title)
		if len(figname0) > 250:
			figname0 = figname0[:250]

		figname = osjoin(self.root_figure, figname0+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)

		figname = osjoin(self.root_figure, figname0+'.png')
		fig.write_image(figname)

		figname = osjoin(self.root_fig_md, figname0+'.png')
		fig.write_image(figname)

		return fig


	def open_affymetrix_human_table(self, fname_affy:str='Human_Agilent_WholeGenome_4x44k_v2_MSigDB_v71.tsv',
									verbose:bool=False) -> pd.DataFrame:
		self.fname_affy = fname_affy
		df_affy = pdreadcsv(fname_affy, self.root_affymetrix, verbose=verbose)

		try:
			df_affy.columns = ['probe', 'symbol', 'description']
		except:
			print("Review df_affy columns, please.")

		self.df_affy = df_affy

		return df_affy

	def test_affymetrix_probe(self, probe:str) -> bool:
		return not self.df_gpl[self.df_gpl.probe == probe].empty


	def calc_unique_LFC_replace_symbol_to_synonym(self, df):

		df = df.sort_values(['symbol', 'abs_lfc'], ascending=[True, False])
		df.reset_index(inplace=True, drop=True)

		previous = ''; goods = []
		for i in range(len(df)):

			if not isinstance(df.iloc[i].symbol, str) or df.iloc[i].symbol == '':
				goods.append(True)
			elif df.iloc[i].symbol != previous:
				previous = df.iloc[i].symbol
				goods.append(True)
			else:
				goods.append(False)

		df = df[goods]
		df['symb_or_syn'] = df.symbol
		df.loc[:, 'symbol'] = [self.gene.replace_symbol_to_synonym(x) for x in df.symbol]

		return df


	def review_LFC_table_with_affymetrix_annotation(self, case, force:bool=False, calc_interm_tables:bool=False, verbose:bool=False) -> bool:
		'''
			set case
			call review_LFC_table_with_affy_annot_wo_case
		'''

		ret, _, _, _ = self.open_case(case, verbose=verbose)
		if not ret:
			return False

		ret = self.review_LFC_table_with_affy_annot_wo_case(force=force, calc_interm_tables=calc_interm_tables, verbose=verbose)

	def review_LFC_table_with_affy_annot_wo_case(self, force:bool=False, calc_interm_tables:bool=False,
												 verbose:bool=False) -> bool:
		'''
		There are 3 final tables:
			self.df_good, self.df_new_symbol, self.df_empty_gpl_lnc_new, self.df_empty

			self.df_good: all annotated RNAs
			self.df_empty_gpl_lnc_new:  new not annotated LNC - df_empty_gpl_lnc_new['_type'] = 'LNC'
			self.df_empty: all empty symbols and not LNC

			all data is stored in dflfc_new

		'''
		case = self.case
		print(">>>", case)

		self.dflfc_all, self.df_good, self.df_empty = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
		self.df_good_gpl, self.df_empty_gpl = pd.DataFrame(), pd.DataFrame()
		self.df_empty_gpl_lnc_new, self.df_empty_gpl_lnc, self.df_empty_gpl_not_lnc = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

		fname_final_lfc_ori = self.fname_final_lfc_table0%(case, self.normalization)
		filename = osjoin(self.root_data, fname_final_lfc_ori)

		if not calc_interm_tables and exists(filename) and not force:
			dflfc_new = pdreadcsv(fname_final_lfc_ori, self.root_data, verbose=verbose)
			self.dflfc_new = dflfc_new
			return True

		fname = self.fname_given_lfc_table0%(self.case, self.normalization)
		dflfc_all = pdreadcsv(fname, self.root_data)
		if verbose: print(">>> Before clean duplicates: len(dflfc_all)", len(dflfc_all))
		dflfc_all = self.calc_unique_LFC_replace_symbol_to_synonym(dflfc_all)
		if verbose: print(">>> After clean duplicates: len(dflfc_all)", len(dflfc_all))
		self.dflfc_all = dflfc_all

		if dflfc_all is None or dflfc_all.empty:
			return False

		df_good = dflfc_all[ ~pd.isnull(dflfc_all.symbol) ].copy()
		df_good.reset_index(inplace=True, drop=True)
		df_good['_type'] = 'hasSymbol'
		self.df_good = df_good
		cols_good_rename = ['probe', 'symbol_prev', 'lfc', 'abs_lfc', 'pval', 'fdr', 'description_prev', 'symbol_pipe',
							'mean_exp', 't', 'B', 'symb_or_syn', '_type']
		df_good.columns = cols_good_rename

		df_empty = dflfc_all[ pd.isnull(dflfc_all.symbol) ].copy()
		self.df_empty = df_empty
		df_empty['_type'] = 'UNK'
		df_empty.reset_index(inplace=True, drop=True)
		df_empty.columns = cols_good_rename

		if len(dflfc_all) != len(df_good)+len(df_empty):
			print(f"Problems with lens: dflfc_all:{len(dflfc_all)}, df_good={len(df_good)}, and df_empty={len(df_empty)}")
			raise Exception("Stop: review_LFC_table_with_affymetrix_annotation()")

		if self.df_gpl is None or self.df_gpl.empty:
			_ = self.open_affymetrix_final_table(verbose=True)
			if self.df_gpl is None or self.df_gpl.empty:
				return False

		df_good_gpl = pd.merge(df_good, self.df_gpl, how='inner', on='probe')
		df_good_gpl.reset_index(inplace=True, drop=True)
		self.df_good_gpl = df_good_gpl
		if verbose: print(">>>> LFC table symbols valid and annotated: df_good_gpl", len(df_good), 'x', len(df_good_gpl))

		cols_good_gpl  = ['probe', 'symbol', 'symbol_prev', 'symb_or_syn', '_type', 'lfc', 'abs_lfc', 'pval', 'fdr', 'mean_exp', 't', 'B',
						  'accession', 'ensembl_id',  'geneid', 'cytoband', 'description', 'description_prev',
						  'symbol_pipe', 'seqname', 'start', 'end', 'go_id', 'seq' ]
		df_good_gpl = df_good_gpl[cols_good_gpl]

		if len(df_good_gpl) != len(df_good):
			print("df_good files have diferent sizes: df_gpl did not map correctly")
			return False


		df_empty_gpl = pd.merge(df_empty, self.df_gpl, how='inner', on='probe')
		self.df_empty_gpl = df_empty_gpl
		if verbose: print(">>>> LFC table symbols empty and annotated: df_good_gpl", len(df_empty), 'x', len(df_empty_gpl))
		df_empty_gpl = df_empty_gpl[cols_good_gpl]

		if df_empty_gpl is None or df_empty_gpl.empty:
			print("No data for df_empty_gpl")
			return False

		if len(df_empty_gpl) != len(df_empty):
			print("df_empty files have diferent sizes: df_gpl did not map correctly")
			return False

		''' after merging, some df_empty_gpl are annotated, some not '''
		df_empty_gpl_symb = df_empty_gpl[ ~pd.isnull(df_empty_gpl.symbol)].copy()
		df_empty_gpl_symb.reset_index(inplace=True, drop=True)
		self.df_empty_gpl_symb = df_empty_gpl_symb

		df_empty_gpl_nosymb = df_empty_gpl[ pd.isnull(df_empty_gpl.symbol)].copy()
		df_empty_gpl_nosymb.reset_index(inplace=True, drop=True)
		self.df_empty_gpl_nosymb = df_empty_gpl_nosymb

		good_lncs = [ True if ( isinstance( df_empty_gpl_nosymb.iloc[i].description, str) and \
							  ('lincRNA' in df_empty_gpl_nosymb.iloc[i].description or \
							   'lncRNA'  in df_empty_gpl_nosymb.iloc[i].description) ) else False \
					   for i in range(len(df_empty_gpl_nosymb))]

		df_empty_gpl_lnc = df_empty_gpl_nosymb[good_lncs].copy()
		df_empty_gpl_lnc.reset_index(inplace=True, drop=True)
		df_empty_gpl_lnc.loc[:,'symbol'] = [f'lncRNA{i+1}' for i in range(len(df_empty_gpl_lnc))]
		df_empty_gpl_lnc['_type'] = 'LNC'
		self.df_empty_gpl_lnc = df_empty_gpl_lnc

		df_empty_gpl_not_lnc = df_empty_gpl_nosymb[np.invert(good_lncs)].copy()
		df_empty_gpl_not_lnc.reset_index(inplace=True, drop=True)
		self.df_empty_gpl_not_lnc = df_empty_gpl_not_lnc

		'''
		df_empty_gpl = df_empty_gpl_symb + df_empty_gpl_nosymb
		df_empty_gpl_nosymb = df_empty_gpl_lnc + df_empty_gpl_not_lnc

		df_empty_gpl_lnc_new  = df_empty_gpl_symb + df_empty_gpl_lnc + df_empty_gpl_not_lnc
		'''
		df_empty_gpl_lnc_new  = pd.concat([df_empty_gpl_symb, df_empty_gpl_lnc, df_empty_gpl_not_lnc])
		df_empty_gpl_lnc_new.reset_index(inplace=True, drop=True)
		self.df_empty_gpl_lnc_new = df_empty_gpl_lnc_new

		''' ----------- GFF3 annoation -----------------'''
		dfgff = self.gene.prepare_final_gff(force=False, verbose=False)
		# & ( ~pd.isnull(dfgff.ensembl_id) )
		dfgff_symb = dfgff[ ~pd.isnull(dfgff.symbol) ].copy()
		dfgff_symb.reset_index(inplace=True, drop=True)

		cols_gff = ['symbol', 'ensembl_id', 'biotype', 'description', ]
		dfgff_symb = dfgff_symb[cols_gff]

		dfgff_symb = dfgff_symb.sort_values(['symbol', 'biotype'])
		if verbose: print(">>> GFF3 annotation: before removing = ", len(dfgff_symb))

		'''
		There are repeated gff annotations for one gene with different ensembl_id

			 	symbol 	ensembl_id 			biotype						   description
		1556 	ZPLD2P 	ENSG00000236155 	transcribed_unitary_pseudogene   zona pellucida like dom...
		1557 	ZPLD2P 	ENSG00000293494 	lncRNA						   zona pellucida like domain ...
		'''
		previous = ''; goods = []
		for i in range(len(dfgff_symb)):

			if dfgff_symb.iloc[i].symbol != previous:
				previous = dfgff_symb.iloc[i].symbol
				goods.append(True)
			else:
				goods.append(False)

		dfgff_symb = dfgff_symb[goods]

		if len(dfgff_symb) == len(dfgff_symb.symbol.unique()):
			if verbose: print(">>> GFF3 annotation: after remove duplicates", len(dfgff_symb), 'Ok.')
		else:
			print("Error: GFF3 annotation: after remove duplicates {len(dfgff_symb)} != {len(dfgff_symb.symbol.unique())}")
		self.dfgff_symb = dfgff_symb

		'''------ reviewing tables ---------------------------------------------'''
		dflfc_new = pd.concat([df_good_gpl, df_empty_gpl_lnc_new])
		dflfc_new = dflfc_new.sort_values('probe')
		dflfc_new.reset_index(inplace=True, drop=True)
		self.dflfc_new = dflfc_new.copy()

		''' some symbols may be fulfilled '''
		df_new_empty = dflfc_new[ pd.isnull(dflfc_new.symbol) ].copy()
		df_new_empty.reset_index(inplace=True, drop=True)
		self.df_new_empty = df_new_empty

		df_new_symbol = dflfc_new[ (~pd.isnull(dflfc_new.symbol)) & (dflfc_new._type != 'LNC') ].copy()
		df_new_symbol.reset_index(inplace=True, drop=True)
		self.df_new_symbol = df_new_symbol

		df_lnc_new = dflfc_new[ (~pd.isnull(dflfc_new.symbol)) & (dflfc_new._type == 'LNC') ].copy()
		df_lnc_new.reset_index(inplace=True, drop=True)
		self.df_lnc_new2 = df_lnc_new

		df_new_symbol_gff = pd.merge(df_new_symbol, dfgff_symb, how='inner', on=['symbol'])
		''' ensembl_id --> ensembl_transc_id is a transcript '''
		cols_rename = ['probe', 'symbol', 'symbol_prev', 'symb_or_syn', '_type', 'lfc', 'abs_lfc', 'pval', 'fdr', 'mean_exp', 't',
					   'B', 'accession', 'ensembl_transc_id', 'geneid', 'cytoband', 'description', 'description_prev',
					   'symbol_pipe', 'seqname', 'start', 'end', 'go_id', 'seq', 'ensembl_id', 'biotype', 'desc_gff']
		df_new_symbol_gff.columns = cols_rename

		''' final columns to be saved and used '''
		cols_final = ['probe', 'symbol', 'symbol_prev', 'symb_or_syn', 'biotype', '_type', 'lfc', 'abs_lfc',
					  'pval', 'fdr', 'mean_exp', 't', 'B', 'description', 'desc_gff',
					  'description_prev', 'accession', 'ensembl_id', 'ensembl_transc_id',
					  'geneid', 'cytoband', 'symbol_pipe', 'seqname', 'start', 'end', 'go_id', 'seq' ]
		df_new_symbol_gff = df_new_symbol_gff[cols_final]
		self.df_new_symbol_gff = df_new_symbol_gff

		df_new_symbol_not_gff = df_new_symbol[~df_new_symbol.probe.isin(df_new_symbol_gff.probe)].copy()
		''' did not merge, remove last 3 columns '''
		df_new_symbol_not_gff.columns = cols_rename[:-3]
		df_new_symbol_not_gff['ensembl_id'] = None
		df_new_symbol_not_gff['biotype'] = df_new_symbol_not_gff['_type']
		df_new_symbol_not_gff['desc_gff'] = None
		df_new_symbol_not_gff = df_new_symbol_not_gff[cols_final]
		df_new_symbol_not_gff.reset_index(inplace=True, drop=True)
		self.df_new_symbol_not_gff = df_new_symbol_not_gff

		'''--------------------- fix mistakes for LOC --------------------------'''

		''' LNC -> there are descriptions like
			"BROAD lincRNAs version v2 (http://www.broadinstitute.org/genome_bio/human_lincrnas/)"
		'''
		df_new_unk = df_new_symbol_not_gff[ df_new_symbol_not_gff._type == 'UNK' ]
		if not df_new_unk.empty:

			probe_list = [df_new_unk.iloc[i].probe for i in range(len(df_new_unk)) \
						  if  isinstance(df_new_unk.iloc[i].description, str) and \
						  ('lincRNA' in df_new_unk.iloc[i].description or 'lncRNA' in df_new_unk.iloc[i].description)]
			if len(probe_list) > 0:
				df_new_symbol_not_gff.loc[df_new_symbol_not_gff.probe.isin(probe_list), ('biotype', '_type')] = [('LNC', 'LNC')]*len(probe_list)

			''' LOC.???? - uncharacterized RNA - ncRNA '''
			probe_list = [df_new_unk.iloc[i].probe for i in range(len(df_new_unk)) if df_new_unk.iloc[i].symbol.startswith('LOC')]
			if len(probe_list) > 0:
				'''												non-coding, uncharacterized '''
				df_new_symbol_not_gff.loc[df_new_symbol_not_gff.probe.isin(probe_list), ('biotype', '_type')] = [('ncRNA', 'UNC')]*len(probe_list)

		'''----------------------- end fix mistakes ----------------------------'''


		''' did not merge, remove last 3 columns '''
		df_lnc_new.columns = cols_rename[:-3]
		df_lnc_new['ensembl_id'] = None
		df_lnc_new['biotype'] = df_lnc_new['_type']
		df_lnc_new['desc_gff'] = None
		df_lnc_new = df_lnc_new[cols_final]
		self.df_lnc_new = df_lnc_new

		df_empty_new = dflfc_new[ pd.isnull(dflfc_new.symbol) ].copy()
		df_empty_new.reset_index(inplace=True, drop=True)
		''' did not merge, remove last 3 columns '''
		df_empty_new.columns = cols_rename[:-3]
		df_empty_new['ensembl_id'] = None
		df_empty_new['biotype'] = df_empty_new['_type']
		df_empty_new['desc_gff'] = None
		df_empty_new = df_empty_new[cols_final]
		self.df_empty_new = df_empty_new

		''' concat: good_gpl, lnc, not_lnc '''
		dflfc_new_final = pd.concat([df_new_symbol_gff, df_new_symbol_not_gff, df_lnc_new, df_empty_new])
		dflfc_new_final = dflfc_new_final.sort_values('symbol')
		dflfc_new_final.reset_index(inplace=True, drop=True)

		previous = ''; goods = []
		for i in range(len(dflfc_new_final)):

			if not isinstance(dflfc_new_final.iloc[i].symbol, str):
				goods.append(False)
			elif dflfc_new_final.iloc[i].symbol != previous:
				previous = dflfc_new_final.iloc[i].symbol
				goods.append(True)
			else:
				goods.append(False)

		dflfc_new_final = dflfc_new_final[goods]
		dflfc_new_final = dflfc_new_final.sort_values('probe')
		dflfc_new_final.reset_index(inplace=True, drop=True)
		self.dflfc_new_final = dflfc_new_final

		'''------ again ---  calc_interm_tables --------------------------------'''
		if not exists(filename) or force:
			pdwritecsv(dflfc_new_final, fname_final_lfc_ori, self.root_data, verbose=verbose)

		print(f"There are {len(self.dflfc_all)} unique RNAs/probes in the microarray expression table")
		print(f"And the processed data len(dflfc_new_final) = {len(self.dflfc_new_final)}")

		if len(self.dflfc_all) == len(self.dflfc_new_final):
			print(f"It is ok")
		else:
			print(f"Warning: lost probes: {len(self.dflfc_all)} != {len(self.dflfc_new_final)}")

		dflfc_lncRNA = self.dflfc_new_final[self.dflfc_new_final.biotype == 'lncRNA']
		len(dflfc_lncRNA)

		print(f"Related to the processed RNAs ({len(self.dflfc_new_final)})")
		print(f"\tSymbols	 annotated in GFF = {len(self.df_new_symbol_gff)}")
		print(f"\tSymbols not annotated in GFF = {len(self.df_new_symbol_not_gff)}")
		print(f"\tEmpty symbols = {len(self.df_empty_new)}")
		print(f"\tNot annotated LNC = {len(self.df_lnc_new)}")
		print(f"\tAnnotated	 LNC = {len(dflfc_lncRNA)}")
		print("")

		self.dflfc_ori = dflfc_new_final
		print(f"\tdflfc_ori: contains {len(dflfc_new_final)}")

		self.valid_genes = [x for x in dflfc_new_final.symbol if isinstance(x, str)]
		self.n_total_genes = len(self.valid_genes)

		print(f"\tdflfc_ori: has {self.n_total_genes} valid symbols.")

		return True

	def get_dflfc_biotypes(self):
		dfg = self.dflfc.groupby('biotype').count().reset_index().iloc[:, :2]
		dfg.columns = ['biotype', 'n']
		return dfg

	def get_dflfc_ori_biotypes(self):
		dfg = self.dflfc_ori.groupby('biotype').count().reset_index().iloc[:, :2]
		dfg.columns = ['biotype', 'n']
		return dfg

	def get_df_biotypes(self, df:pd.DataFrame):
		dfg = df.groupby('biotype').count().reset_index().iloc[:, :2]
		dfg.columns = ['biotype', 'n']
		return dfg

	def split_chr_location(self, x:str):
		if x == 'unmapped':
			return None, None, None
		if pd.isnull(x) or x == 'nan':
			return None, None, None

		start, end = -1, -1
		try:
			mat =  x[3:].split(":")
			seqname = mat[0]
			start, end = mat[1].split('-')
			start = int(start) if isint(start) else start
			end   = int(end)   if isint(end)   else end
		except:
			print(x)
			return '', start, end

		return seqname, start, end

	def open_affymetrix_final_table(self, fname_affy:str="GPL14550-9757.tsv", verbose:bool=False) -> pd.DataFrame:

		fname_geo = fname_affy.replace('.tsv', '_location.tsv')
		self.fname_geo = fname_geo
		filename = osjoin(self.root_affy, fname_geo)

		if not exists(filename):
			print(f"Could not find '{filename}'")
			self.df_gpl = pd.DataFrame()
			return self.df_gpl

		df_gpl = pdreadcsv(fname_geo, self.root_affy, verbose=verbose)
		self.df_gpl = df_gpl

		return df_gpl

	def open_affymetrix_table(self, fname_geo:str, verbose:bool=False) -> pd.DataFrame:
		self.fname_geo = fname_geo
		df_gpl = pdreadcsv(fname_geo, self.root_affy, verbose=verbose)

		cols = ['probe', 'probe_id', 'CONTROL_TYPE', 'refseq', 'accession', 'geneid', 'symbol', 'name',
	   'unigene_id', 'ensembl_id', 'TIGR_ID', 'ACCESSION_STRING', 'chr_location', 'cytoband', 'description', 'go_id', 'seq']
		df_gpl.columns = cols
		cols2  = ['probe', 'refseq', 'accession', 'geneid', 'symbol', 'name',
				  'unigene_id', 'ensembl_id', 'chr_location', 'cytoband', 'description', 'go_id', 'seq']

		df_gpl = df_gpl[cols2]
		self.df_gpl = df_gpl

		return df_gpl

	def prepare_gpl(self, fname_affy:str="GPL14550-9757.tsv", force:bool=False, verbose:bool=False) -> bool:

		fname = fname_affy.replace('.tsv', '_location.tsv')
		filename = osjoin(self.root_affy, fname)

		if exists(filename) and not force:
			df_gpl = pdreadcsv(fname, self.root_affy, verbose=verbose)
			self.df_gpl = df_gpl
			return True

		df_gpl = self.open_affymetrix_table(fname_geo=fname_affy, verbose=verbose)

		df_gpl.loc[:, ['seqname', 'start', 'end']] = [self.split_chr_location(x) for x in df_gpl.chr_location]

		cols = ['probe', 'refseq', 'accession', 'geneid', 'symbol', 'name', 'description',
				'unigene_id', 'ensembl_id', 'cytoband', 'seqname', 'start', 'end', 'go_id', 'seq']

		df_gpl = df_gpl[cols]
		self.df_gpl = df_gpl

		ret = pdwritecsv(df_gpl, fname, self.root_affy, verbose=verbose)
		return ret


	def how_many_degs_given_LFC_FDR_cutoffs(self, dflfc:pd.DataFrame, FDR_cutoff:float, LFC_cut:float):
		return len(dflfc[ (dflfc.fdr < FDR_cutoff) & (dflfc.abs_lfc >= LFC_cut) ])

	def calc_norm_cdf_for_bayesian(self, lfc:float, mean_LFC:float=1.0, std_LFC:float=0.2):
		p = norm.cdf(lfc, loc=mean_LFC, scale=std_LFC)
		if p > 0.5: p = 1 - p
		return p


	lista_bayes = list(np.round(np.arange(0.05, 0.80, 0.05), 2))

	def calc_bayesian_cutoffs(self, case:str, mean_LFC:float=1.0, std_LFC:float=0.2, ndraws:int=1000, 
					 		  fdr_lista:List=lista_bayes, perc_delta:float=0.01,
						      xaxis_title:str='LFC', yaxis_title:str='p',
					 		  width:int=1100, height:int=600, 
							  horizontal_spacing:float=0.1, vertical_spacing:float=0.1,
							  plot_bgcolor:str='lightgray', verbose:bool=False) -> dict:


		_, _, _, _ = self.open_case(case, verbose=False)
		dflfc = self.dflfc_ori
		N = len(dflfc)

		lfc_samples = np.random.normal(loc=mean_LFC, scale=std_LFC, size=ndraws)

		# minimum LFC cutoff = 0.2
		# lfc_samples = [x if x > .2 else 0.2 for x in lfc_samples]
		lfc_samples = [x for x in lfc_samples if x >= 0]

		subtitles = [f'LFC ~ N({mean_LFC:.1f}, {std_LFC:.2f})', 'p(DEG|LFC)', 'nDEG', 'p(LFC) = a prioriri', 'p(LFC|DEG) = a posteriori', None]

		dic = {}; icount = -1
		for FDR in fdr_lista:
			title = f'Bayesian simulation for {case} with {ndraws} draws and FDR = {FDR}'

			ndegs = [ self.how_many_degs_given_LFC_FDR_cutoffs(dflfc, FDR, lfc_sample) for lfc_sample in lfc_samples]
			pis   = [ ndeg/N for ndeg in ndegs ]

			df = pd.DataFrame( {'lfc': lfc_samples, 'p': pis, 'ndegs': ndegs} )
			df = df.sort_values('lfc')
			df.reset_index(inplace=True, drop=True)

			nrows = len(df)

			df['FDR'] = FDR
			df['p_lfc'] = [ self.calc_norm_cdf_for_bayesian(lfc, mean_LFC, std_LFC) for lfc in df.lfc]
			# normalizing
			total = df['p_lfc'].sum()
			df['p_lfc'] = df['p_lfc']/total
			
			#--- Dr. Bayes ------------------
			df['p_posterior'] = df['p'] * df['p_lfc']
			total = df['p_posterior'].sum()
			df['p_posterior'] = df['p_posterior']/total
			
			mu = np.sum( df['p_posterior'] * df['lfc'] )
			delta = perc_delta * mu

			self.mu = mu
			self.delta = delta
			self.df = df

			if total == 0:
				print(f"Warning: total posterior probability is zero for FDR={FDR}")
				continue

			dfm = df[ (df.lfc >= mu-delta) & (df.lfc <= mu+delta) ]
			i_mu = int(np.median(dfm.index.to_list()))
			
			maxi = df.p_posterior.max()
			delta = 0.001 * maxi
			dfm = df[ (df.p_posterior >= maxi-delta) & (df.p_posterior <= maxi+delta) ]
			lfc_max = dfm.lfc.mean()

			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]
			dic2['case'] = case
			dic2['fdr'] = FDR
			dic2['lfc_max'] = lfc_max
			dic2['lfc_mean'] = mu

			stri = f'For FDR={FDR}, the lfc_max = {lfc_max:.3f}'
			dic2['str_stat'] = stri

			for alpha in [0.9, 0.8, 0.7]:
				delta = int(np.round(nrows * alpha/ 2, 0))

				lim_inf = i_mu - delta
				if lim_inf < 0: lim_inf = 0
				lim_sup = i_mu + delta
				if lim_sup > nrows: lim_sup = nrows-1

				try:
					ci = np.round([df.iloc[lim_inf].lfc, df.iloc[lim_sup].lfc], 3)
				except:
					print("Error", nrows, alpha, delta,  lim_inf, lim_sup)
					ci = None

				# stri = f"\tfor alpha={alpha} the CI = {ci}"

				dic2[f'CI_{alpha}'] = ci
			
			fig = make_subplots(rows=2, cols=3, subplot_titles=subtitles, 
					    		horizontal_spacing=horizontal_spacing, vertical_spacing=vertical_spacing)
			
			arrow = None
			for i in range(5):
				xs = df.lfc

				if i == 0:
					nrow=1
					ncol=1
					vals = df.lfc
				elif i == 1:
					nrow=1
					ncol=2
					vals = df.p
				elif i == 2:
					nrow=1
					ncol=3
					vals = df.ndegs
				elif i == 3:
					nrow=2
					ncol=1
					vals = df.p_lfc
				else:
					nrow=2
					ncol=2
					vals = df.p_posterior
					# here is different
					# xs = df.ndegs
			
				name = subtitles[i]
			
				if i == 0:
					fig.add_trace(go.Histogram(x=vals, name=name) , row=nrow, col=ncol)
				else:
					fig.add_trace(go.Scatter(x=xs, y=vals, name=name) , row=nrow, col=ncol)

					fig.update_yaxes(range=[0, None], row=nrow, col=ncol)

					if i == 4:
						expected_lfc = np.sum(df.p_posterior * df.lfc)
						y_max = df.p_posterior.max()/2

						arrow = go.layout.Annotation(dict(
							x=expected_lfc,
							y=0,
							xref="x5", yref="y5",
							text=f"E[LFC] = {expected_lfc:.2f}",
							showarrow=True,
							axref="x5", ayref='y5',
							ax=expected_lfc,
							ay=y_max,
							arrowhead=3,
							arrowwidth=1.5,
							arrowcolor='red')
						)

				fig.update_xaxes(title_text=xaxis_title, row=nrow, col=ncol)
				fig.update_yaxes(title_text=yaxis_title, row=nrow, col=ncol)

			fig.add_annotation(arrow)
			
			fig.update_layout(
						autosize=True,
						title=title,
						width=width,
						height=height,
						showlegend=False,
						legend_title='probs',
						plot_bgcolor=plot_bgcolor,
						margin=dict(t=120),
						font=dict(
							family="Arial",
							size=14,
							color="Black"
						), 
			)
			
			dic2['fig'] = fig

			figname = title_replace(title)
			figname = osjoin(self.root_figure, figname+'.html')

			fig.write_html(figname)
			if verbose: print(">>> HTML and png saved:", figname)
			fig.write_image(figname.replace('.html', '.png'))

		return dic



	def plot_bayesian_cutoff_series(self, case:str, dfc:pd.DataFrame, xaxis_title:str='FDR',
									width:int=1000, height:int=600, plot_bgcolor:str='lightgray', 
									verbose:bool=False):

		title = f'Bayesian Analysis: (max and mean) x FDR for {case}'
			
		subtitles = ['max LFC', 'mean LFC']
		fig = make_subplots(rows=2, cols=1, subplot_titles=subtitles)

		for i in range(2):
		
			for i in range(2):
				if i == 0:
					nrow=1
					vals = dfc.lfc_max
				else:
					nrow=2
					vals = dfc.lfc_mean
			
				name = subtitles[i]
			
				fig.add_trace(go.Scatter(x=dfc.fdr, y=vals, name=name) , row=nrow, col=1)
			
			
		fig.update_layout(
					autosize=True,
					title=title,
					width=width,
					height=height,
					xaxis2_title=xaxis_title,
					yaxis_title ='max LFC',
					yaxis2_title ='mean LFC',
					showlegend=False,
					legend_title='probs',
					plot_bgcolor=plot_bgcolor,
					font=dict(
						family="Arial",
						size=14,
						color="Black"
					) )
		
		figname = title_replace(title)
		figname = osjoin(self.root_figure, figname+'.html')

		fig.write_html(figname)
		if verbose: print(">>> HTML and png saved:", figname)
		fig.write_image(figname.replace('.html', '.png'))		

		return fig

	def calc_sectors(self, dff_case:pd.DataFrame, sel_list:List=[], lfc_threshold:float=1.) -> Tuple[dict, dict]:

		sectors = {}; dic={}

		if sel_list == []:
			print("Error: sel_list is empty")
			return sectors, dic

	
		for pathway in sel_list:
			dfa = dff_case[dff_case.pathway == pathway]

			if dfa.empty:
				continue

			row = dfa.iloc[0]
			
			pathway_id   = row.pathway_id
			gene_up_list = row.mod_up_in_pathway
			gene_up_list = eval(gene_up_list) if isinstance(gene_up_list, str) else gene_up_list
		
			lfc_up_list = row.lfc_up
			lfc_up_list = eval(lfc_up_list) if isinstance(lfc_up_list, str) else lfc_up_list
		
			gene_dw_list = row.mod_dw_in_pathway
			gene_dw_list = eval(gene_dw_list) if isinstance(gene_dw_list, str) else gene_dw_list
		
			lfc_dw_list = row.lfc_dw
			lfc_dw_list = eval(lfc_dw_list) if isinstance(lfc_dw_list, str) else lfc_dw_list
		
			dic_gene = {}
			for i, gene in enumerate(gene_up_list):
				if lfc_up_list[i] > lfc_threshold:
					dic_gene[gene] = lfc_up_list[i]
					
			for i, gene in enumerate(gene_dw_list):
				if lfc_dw_list[i] < -lfc_threshold:
					dic_gene[gene] = lfc_dw_list[i]
					
			n = len(dic_gene)

			sectors[pathway] = n
		
			dic[pathway] = {}
			dic2 = dic[pathway]
			dic2['pathway_id'] = pathway_id
			dic2['n'] = n
			# print(pathway, n)
			dic2['genes'] = list(dic_gene.keys())
			dic2['lfc'] = list(dic_gene.values())

		return sectors, dic


	def calc_mod_and_plot_circle_sel_list(self, sel_list:List, title_head:str, 
					plot_bars:bool=True, plot_lines:bool=False,
					fontsize_circos:int=14, lfc_threshold:float=1., space_factor:float=3, 
					with_pPMI_obs:bool=False, do_case_translate:bool=True, 
					dpi:int=300, format:str='png', facecolor:str='white', 
					show_fig:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, list]:

		_, dff, _, _ = self.calc_pPMI_summary_and_pivot_tables(pathway_list=sel_list, title_head=title_head, with_pPMI_obs=with_pPMI_obs, 
															   do_case_translate=do_case_translate, verbose=verbose)

		figname_circos_list = []

		self.dfpiv_mod, self.dfpiv_inv = pd.DataFrame(), pd.DataFrame()
		self.dfpiv_high = pd.DataFrame()
		self.dfpiv_high_inter, self.dfpiv_high_female, self.dfpiv_high_male = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

		if self.has_gender:
			case_list = self.group_female_list + self.group_male_list
		else:
			case_list = self.case_list

		df_list = []
		for case in case_list:
			dff_case = dff[dff.case == case]
			
			circos, dfmod =  self.plot_circle(dff_case, sel_list=sel_list, plot_bars=plot_bars, plot_lines=plot_lines, lfc_threshold=lfc_threshold, space_factor=space_factor)
			if circos is None:
				continue

			df_list.append(dfmod)
			
			if circos:
				#circos.savefig("example01.png")
				plt.rcParams["font.size"] = fontsize_circos
				_ = circos.plotfig(figsize=(13,13))
				for text in circos.ax.texts:
					text.set_fontsize(fontsize_circos)

				plt.title(f'Case: {case}')

				figname0 = f"circos_{title_head}_{case}"
				figname0 = title_replace(figname0) 

				figname = osjoin(self.root_figure, figname0 + '.' + format)
				if verbose: print(">>> image saved:", figname)
				plt.savefig(figname, dpi=dpi, format=format, facecolor=facecolor)

				figname = osjoin(self.root_fig_md, figname0 + '.' + format)
				figname_circos_list.append(figname)
				plt.savefig(figname, dpi=dpi, format=format, facecolor=facecolor)

				if show_fig: plt.show()
				plt.clf()
				plt.cla()
				plt.close()

		dfa = pd.concat(df_list)
		self.dfa = dfa
		dfa.reset_index(inplace=True, drop=True)

		dfpiv_cases = pd.pivot_table(dfa, values='lfc', index='symbol', columns='case', fill_value=None)
		self.dfpiv_cases = dfpiv_cases

		#-------------- merge data from original source: dflfc_ori -------------------
		df2 = dfpiv_cases.copy()

		for case in self.case_list:
			self.open_case(case)
			case2 = case + '_2'

			# copy deflc_ori lfc to case2 --> dfpiv_cases[case]
			# a gene may exists in one case and not in another
			df2 = pd.merge(df2, self.dflfc_ori[['symbol','lfc']], how="outer", on='symbol')
			df2 = df2[df2.symbol.isin(dfpiv_cases.index.to_list())]

			cols = list(df2.columns)
			cols[-1] = case2
			df2.columns = cols

			'''
			# cannot compare: to plot it limits to [-3, +3]
			for i in range(len(df2)):
				row = df2.iloc[i]
				if not pd.isnull(row[case]):
					if row[case] != row[case+'_2']:
						print(f"Error for {case}", row.symbol, row[case],  row[case2])
						raise Exception('stop')
			'''
			dfpiv_cases[case] = df2[case2].to_list()

		def greater_equal_abs_mod_diff_cutoff(row: pd.Series):
			ret = False

			row = row.abs()
			for i in range(len(self.case_list)):
				if row[i] >= self.abs_mod_diff_cutoff:
					ret = True
					break
			return ret

		lista = [greater_equal_abs_mod_diff_cutoff(dfpiv_cases.iloc[i][self.case_list]) for i in range(len(dfpiv_cases))]
		self.dfpiv_mod = dfpiv_cases[lista].copy()

		def positive_negative(row: pd.Series):
			is_positive = False
			is_negative = False

			for i in range(len(self.case_list)):
				if row[i] >= self.min_mod_diff_cutoff:
					is_positive = True
				elif row[i] <= -self.min_mod_diff_cutoff:
					is_negative = True

			return is_positive and is_negative
		
		lista = [positive_negative(dfpiv_cases.iloc[i]) for i in range(len(dfpiv_cases))]
		self.dfpiv_inv = dfpiv_cases[lista].copy()

		def is_highly_modulated_simple(row: pd.Series):
			
			lista = np.array(row)
			mini = np.min(lista)
			maxi = np.max(lista)

			return (maxi-mini) >= self.highly_mod_cutoff

		def is_highly_modulated_by_comparing_genders(row_female: pd.Series, row_male:pd.Series):
			
			lista = np.array(row_female) - np.array(row_male)
			mini = np.min(lista)
			maxi = np.max(lista)

			return (maxi-mini) >= self.highly_mod_cutoff

		if self.has_gender:
			lista = [is_highly_modulated_by_comparing_genders(dfpiv_cases.iloc[i][self.group_female_list], \
															  dfpiv_cases.iloc[i][self.group_male_list]) \
															  for i in range(len(dfpiv_cases))]
			dfpiv_aux = dfpiv_cases[lista].copy()

			dfpiv_aux_female = dfpiv_aux[self.group_female_list]
			dfpiv_aux_female.columns = self.group_list
			dfpiv_aux_male = dfpiv_aux[self.group_male_list]
			dfpiv_aux_male.columns = self.group_list

			dfpiv_diff= dfpiv_aux_male.infer_objects().fillna(0) - dfpiv_aux_female.infer_objects().fillna(0)
			self.dfpiv_high_inter = dfpiv_diff

			lista = [is_highly_modulated_simple(dfpiv_cases.iloc[i][self.group_female_list]) for i in range(len(dfpiv_cases))]
			self.dfpiv_high_female = dfpiv_cases[lista][self.group_female_list].copy()

			lista = [is_highly_modulated_simple(dfpiv_cases.iloc[i][self.group_male_list]) for i in range(len(dfpiv_cases))]
			self.dfpiv_high_male = dfpiv_cases[lista][self.group_male_list].copy()
		else:
			lista = [is_highly_modulated_simple(dfpiv_cases.iloc[i]) for i in range(len(dfpiv_cases))]
			self.dfpiv_high = dfpiv_cases[lista].copy()

		return dfpiv_cases, figname_circos_list
	

	def plot_circle(self, dff_case:pd.DataFrame, sel_list:List=[], 
					plot_bars:bool=True, plot_lines:bool=False, 
					lfc_threshold:float=1., maxi:float=+3., mini:float=-3., 
					space_factor:float=1.8) -> Tuple[Circos, pd.DataFrame]:

		colors = self.my_colors
		sectors, dic = self.calc_sectors(dff_case, sel_list=sel_list, lfc_threshold=lfc_threshold)
		self.sectors = sectors
		self.dic = dic

		if len(sectors) == 0:
			print("No sectors")
			return None, pd.DataFrame()

		case = dff_case.iloc[0]['case']

		circos = Circos(sectors, space=len(sectors)*space_factor)

		names = list(dic.keys())

		seq_y = np.arange(mini, maxi+1)
		s_seq_y = [ str(int(val)) if val <= 0 else '+'+str(int(val)) for val in seq_y]

		i_sector = 0
		df_list = []
		for sector in circos.sectors:
			# Plot sector name
			i_sector += 1
			sector.text(f"Sector: {i_sector}", r=110, size=18)
			pathway = sector.name
			
			# Create x positions & random y values
			y = dic[pathway]['lfc']

			genes = dic[pathway]['genes']
			genes = eval(genes) if isinstance(genes, str) else genes

			y = np.round(y, 3)
			y = [val if val <= maxi else maxi for val in y]
			y = [val if val >= mini else mini for val in y]

			dfa = pd.DataFrame({'y':y, 'gene': genes})
			dfa = dfa.sort_values(['y', 'gene'], ascending=[False, True])
			dfa.reset_index(inplace=True, drop=True)

			df2 = dfa.copy()
			df2['case'] = case
			df2['pathway'] = pathway
			df2.columns = ['lfc', 'symbol', 'case', 'pathway']
			cols =  ['case', 'pathway', 'symbol', 'lfc']
			df2 = df2[cols]
			df_list.append(df2)

			self.sector = sector
			self.dfa = dfa

			if dfa.empty:
				continue

			x = np.arange(sector.start, sector.end) + 0.5
			y = dfa['y'].to_list()
			genes = dfa['gene'].to_list()

			dic[pathway]['genes'] = genes
			dic[pathway]['lfc'] = y

			offset_ticks = np.arange(0.5, sector.end + 0.5, 1)
			offset_labels =  [genes[i] if len(genes) > i else '' for i in range(sector.end)]

			y_middle = [0]*len(x)
			
			if plot_lines and plot_bars:
				# track1.xticks_by_interval(interval=1, label_formatter=lambda v: f"{genes[v] if v<len(genes) else ''}", label_orientation="vertical",)
				track1 = sector.add_track((25, 50), r_pad_ratio=0.1)
				track1.axis()
				track1.line(x, y, vmin=mini, vmax=maxi)
				track1.line(x, y_middle, color='black', vmin=mini, vmax=maxi)

				track3 = sector.add_track((50, 80), r_pad_ratio=0.1)
				track3.axis()
				track3.bar(x, y, vmin=mini, vmax=maxi )
				track3.line(x, y_middle, color='black', vmin=mini, vmax=maxi)
				track3.xticks(offset_ticks, offset_labels, label_orientation="vertical")
				track3.yticks(seq_y, s_seq_y, vmin=mini, vmax=maxi)

			
			elif plot_lines:
				track1 = sector.add_track((50, 80), r_pad_ratio=0.1)
				track1.axis()
				# track1.xticks_by_interval(interval=1, label_formatter=lambda v: f"{genes[v] if v<len(genes) else ''}", label_orientation="vertical",)
				track1.line(x, y, vmin=mini, vmax=maxi)
				track1.line(x, y_middle, color='black', vmin=mini, vmax=maxi)
				track1.xticks(offset_ticks, offset_labels, label_orientation="vertical",)
				track1.yticks(seq_y, s_seq_y, vmin=mini, vmax=maxi)
				
			else: # plot_bars
				track3 = sector.add_track((50, 80), r_pad_ratio=0.1)
				track3.axis()
				track3.bar(x, y, vmin=mini, vmax=maxi )
				track3.line(x, y_middle, color='black', vmin=mini, vmax=maxi)
				track3.xticks(offset_ticks, offset_labels, label_orientation="vertical")
				track3.yticks(seq_y, s_seq_y, vmin=mini, vmax=maxi)

		#---------------------- Plot links -----------------------------
		iplot = -1

		N = len(names)
		if N > 1:
			for i in range(N-1):
				for j in range(i+1, N ):
					iplot += 1
					color = colors[iplot]
					rgb = [int(x*255) for x in mpl_colors.to_rgb(color)]

					try:
						name_a = names[i]
						name_b = names[j]
					except:
						print("\nError:", names, i, j, '\n\n')
						continue

					genes_a = dic[name_a]['genes']
					genes_b = dic[name_b]['genes']

					set_a = set(genes_a)
					set_b = set(genes_b)

					commons = list( set_a.intersection(set_b) )

					i_link = -1
					for gene in commons:
						i_a, i_b = -1, -1
						
						for igene in range(len(genes_a)):
							if genes_a[igene] == gene:
								i_a = igene
								break

						for igene in range(len(genes_b)):
							if genes_b[igene] == gene:
								i_b = igene
								break
						
						if i_a >= 0 and i_b >= 0:
							i_link += 1
							'''
							if i_link%2 == 1:
								hatch = "//"
							else:
								hatch = None
							'''

							params = rgb + [i_link*3]
							hex_color = inc_rgb_to_hex(*params)
							circos.link((name_a, i_a, i_a+1), (name_b, i_b, i_b+1), ec="black", color=hex_color, lw=1, ) 
							# hatch=hatch), ls="dotted"

		dfmod = pd.concat(df_list)
		dfmod.reset_index(inplace=True, drop=True)
		
		return circos, dfmod
		
	def list_LFC_modulation(self, selected_sets:List, lfc_threshold:float, diff_cutoff:float,
						    force:bool=False, verbose:bool=False) -> str:

		fname = self.fname_lfc_mod_summary%(lfc_threshold, diff_cutoff)
		filename = osjoin(fname, self.root_ptw_modulation)
		
		if exists(filename) and not force:
			text = read_txt(fname, self.root_ptw_modulation, verbose=verbose)
			return text

		text = ''
		for set_pathway in selected_sets:
			text += self.one_pathway_LFC_modulation(set_pathway, lfc_threshold=lfc_threshold, diff_cutoff=diff_cutoff)
			text += "\n"

		text += "------------ end -------------\n"

		write_txt(text, fname, self.root_ptw_modulation, verbose=verbose)

		return text
	

	def one_pathway_LFC_modulation(self, set_pathway:List, lfc_threshold:float, diff_cutoff:float) -> str:
		
		term = 'type of LFC modulation'
		title_head, pathway_list, _ = set_pathway

		text = f"# LFC modulation summary\n\n"
		text += f"**{title_head}** -> pathways: {'; '.join(pathway_list)}\n\n"

		for type_modulation in ['all', 'similar', 'different']:
			for up_down_all in ['all', 'up', 'dw']:
		
				text += f"## Modulation {type_modulation} and up_down_all {up_down_all}\n"
				_, df_anal, s_filtered = self.calc_pathway_modulation(set_pathway, 
											type_modulation=type_modulation, up_down_all=up_down_all,
											lfc_threshold=lfc_threshold, diff_cutoff=diff_cutoff)
		
				if df_anal.empty:
					text += f"Error: s_filtered {s_filtered} -> no data (empty)\n"
				else:
					gene_list = list(df_anal.symbol)
					gene_list = np.unique(gene_list)

					s_filtered = s_filtered.replace('<br>', '\n')
					
					s_filtered = s_filtered.replace(term, ' -> '+term)
					text += f"{s_filtered}\n\n{len(gene_list)} {self.s_deg_dap}s:\n\n{'; '.join(gene_list)}\n"
				
				text += "\n"


		return text
	

	def import_from_GDC(self, prog_id:str = "TCGA", 
					    force:bool = False, verbose:bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:

		fname_lfc, _, _ = self.set_lfc_names()
		fname_lfc_all = fname_lfc.replace('.tsv', '_all_transcripts.tsv')

		filename = self.root_lfc / fname_lfc
		filename_all = self.root_lfc / fname_lfc

		if filename.exists() and filename_all.exists() and not force:
			df_lfc     = pdreadcsv(fname_lfc, self.root_lfc, verbose=verbose)
			df_lfc_all = pdreadcsv(fname_lfc_all, self.root_lfc, verbose=verbose)
			return df_lfc, df_lfc_all

		gdc = GDC(ROOT0=self.root0, ROOT_DATA0=self.root0_data)
		self.gdc = gdc

		_ = gdc.get_primary_sites(prog_id=prog_id, force=False, verbose=verbose)

		print(">>> psi_id or disease:", self.disease)
		psi_id = self.disease
		self.psi_id = psi_id

		gdc.set_primary_site(psi_id=psi_id, verbose=False)

		df_lfc_all, msg = gdc.calc_lfc_table(
			psi_id=psi_id,
			root_src=self.root_src,
			run_conda=True,
			method="deseq2",
			verbose=verbose,
		)

		if df_lfc_all is None or df_lfc_all.empty:
			print(f"Error: No data available for the specified {psi_id}.")
			return pd.DataFrame(), pd.DataFrame()

		df_lfc_all["abs_lfc"] = np.abs(df_lfc_all.lfc)

		df_lfc = (
			df_lfc_all
			.dropna(subset=["symbol"])
			.sort_values("abs_lfc", ascending=False)
			.drop_duplicates(subset="symbol", keep="first")
			.reset_index(drop=True)
		)

		self.df_lfc = df_lfc

		if verbose:
			print(msg)

		_ = pdwritecsv(df_lfc, fname_lfc, self.root_lfc, verbose=verbose)
		_ = pdwritecsv(df_lfc_all, fname_lfc_all, self.root_lfc, verbose=verbose)

		return df_lfc, df_lfc_all

