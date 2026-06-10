#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2024/01/29
# Updated on 2025/10/21
# @author: Flavio Lichtenstein
# @institution: Butantan Institute, Molecular Biology, Bioinformatics, CENTD

# import numpy as np
import os
import numpy as np
from pathlib import Path
import pandas as pd
from typing import Tuple, List # Optional, Iterable, Set, Any

from libs.Basic import pdreadcsv, pdwritecsv, create_dir, title_replace

class Config(object):
	def __init__(self, root0:Path, root_disease:Path, disease:str, case_list:List):

		self.root0 = root0
		self.root_disease = root_disease
		self.root_colab = create_dir(root0, 'colab')
		self.root_config = create_dir(root_disease, 'config')

		self.case_list = case_list

		self.disease = title_replace(disease)

		self.normalization = 'not_normalized'

		self.dfbest_lfc_cutoff = None
		self.dfbest_cutoffs = pd.DataFrame()

		self.LFC_cut = 1
		self.lfc_FDR_cut = 0.05

		self.n_degs = -1
		self.n_degs_up = -1
		self.n_degs_dw = -1

		self.quantile = -1
		self.geneset_num = 0
		self.ptw_FDR_cut = 0.05
		self.ptw_min_num_of_degs_cut = 3

		self.n_genes_annot_ptw = -1
		self.n_degs = -1
		self.n_degs_in_ptw = -1
		self.n_degs_not_in_ptw = -1
		self.degs_in_all_ratio = -1

		self.toi1 = -1
		self.toi2 = -1

		# row.LFC_cut, row.lfc_FDR_cut, row.n_degs,  row.n_degs_up,  row.n_degs_dw
		self.param_lfc_defaults = (1, 0.05, -1, -1, -1)

		'''
		return row.LFC_cut, row.lfc_FDR_cut, row.ptw_FDR_cut,  \
			   row.n_genes_annot_ptw, row.n_degs, row.n_degs_in_ptw, row.n_degs_not_in_ptw, row.degs_in_all_ratio, row.toi1, row.toi2
		'''
		self.param_ptw_defaults = (0.5, 0.9, 1, 0.05, 0.05, 0.05, 3, -1, -1, -1, -1, -1, -1, -1, -1)

		fname_lfc_cutoff = f'all_lfc_cutoffs_{self.disease}.tsv'
		self.fname_lfc_cutoff = title_replace(fname_lfc_cutoff)

		fname_ptw_cutoff = f'best_ptw_cutoffs_{self.disease}.tsv'
		self.fname_ptw_cutoff = title_replace(fname_ptw_cutoff)

	def set_default_best_lfc_cutoff(self, normalization:str, LFC_cut:float=1, lfc_FDR_cut:float=.05):

		self.quantile = -1
		self.normalization = normalization

		self.LFC_cut = LFC_cut
		self.lfc_FDR_cut = lfc_FDR_cut
		self.cutoff = f"{LFC_cut:.3f} - {lfc_FDR_cut:.3f}"

		self.n_degs = -1
		self.n_degs_up = -1
		self.n_degs_dw = -1


	def open_all_lfc_cutoff(self, verbose=True) -> pd.DataFrame:
		filename = self.root_config / self.fname_lfc_cutoff
		if not filename.exists():
			if verbose: print(f"Best parameter file for LFC does not exist {filename}")
			self.dfbest_lfc_cutoff = pd.DataFrame()
			return pd.DataFrame()

		dfbest_lfc_cutoff = pdreadcsv(self.fname_lfc_cutoff, self.root_config)

		cols = np.array(dfbest_lfc_cutoff.columns)
		# fix - columns rename - remove in the future
		if 'abs_lfc_cutoff' in cols:
			cols = ['case', 'normalization', 'cutoff', 'LFC_cut',
					'lfc_FDR_cut', 'degs', 'n_degs', 'degs_ensembl',
					'n_degs_ensembl', 'degs_up', 'n_degs_up', 'degs_up_ensembl',
					'n_degs_up_ensembl', 'degs_dw', 'n_degs_dw', 'degs_dw_ensembl',
					'n_degs_dw_ensembl']
			dfbest_lfc_cutoff.columns = cols
			_ = pdwritecsv(dfbest_lfc_cutoff, self.fname_lfc_cutoff, self.root_config)

		self.dfbest_lfc_cutoff = dfbest_lfc_cutoff

		return dfbest_lfc_cutoff

	def save_best_lfc_cutoff(self, dfi, verbose:bool=False) -> bool:
		dfbest_lfc_cutoff = self.open_all_lfc_cutoff(verbose=False)

		if dfbest_lfc_cutoff is None or dfbest_lfc_cutoff.empty:
			dfbest_lfc_cutoff = dfi
		else:
			dfbest_lfc_cutoff = dfbest_lfc_cutoff[(dfbest_lfc_cutoff.normalization != self.normalization)]

			if dfbest_lfc_cutoff.empty:
				dfbest_lfc_cutoff = dfi
			else:
				dfbest_lfc_cutoff = pd.concat([dfbest_lfc_cutoff, dfi])

		dfbest_lfc_cutoff.reset_index(inplace=True, drop=True)

		_ = pdwritecsv(dfbest_lfc_cutoff, self.fname_lfc_cutoff, self.root_config, verbose=verbose)
		return ret

	def get_best_lfc_cutoff(self, case:str, normalization:str, verbose:bool=False) -> Tuple[float, float, int, int, int]:
		if self.dfbest_lfc_cutoff is None or self.dfbest_lfc_cutoff.empty:
			_ = self.open_all_lfc_cutoff()
			if self.dfbest_lfc_cutoff is None or self.dfbest_lfc_cutoff.empty:
				if verbose:
					print("Houston we have a problem: no CONFIGURATION file was found.")
					print(">>> run: new03_up_down_simulation")
				return self.param_lfc_defaults

		dfa = self.dfbest_lfc_cutoff[(self.dfbest_lfc_cutoff.case == case) & (self.dfbest_lfc_cutoff.normalization == normalization) ]
		if dfa.empty:
			return self.param_lfc_defaults

		row = dfa.iloc[0]

		return row.LFC_cut, row.lfc_FDR_cut, row.n_degs,  row.n_degs_up,  row.n_degs_dw


	def set_default_best_ptw_cutoff(self, normalization:str, geneset_num:int=0, quantile:float=0.5,
									LFC_cut:float=1, lfc_FDR_cut:float=.05, ptw_FDR_cut:float=0.05,
									n_genes_annot_ptw:int=0, n_degs:int=0, n_degs_in_ptw:int=0,  n_degs_not_in_ptw:int=0, degs_in_all_ratio:int=0):

		self.normalization = normalization
		self.geneset_num = geneset_num
		self.quantile = quantile

		self.cutoff = f"{LFC_cut:.3f} - {lfc_FDR_cut:.3f}"
		self.LFC_cut = LFC_cut
		self.lfc_FDR_cut = lfc_FDR_cut
		self.ptw_FDR_cut = ptw_FDR_cut

		self.n_genes_annot_ptw = n_genes_annot_ptw

		self.n_degs = n_degs

		self.n_degs_in_ptw = n_degs_in_ptw
		self.n_degs_not_in_ptw = n_degs_not_in_ptw
		self.degs_in_all_ratio = degs_in_all_ratio

		self.toi1 = -1
		self.toi2 = -1


	def open_best_ptw_cutoff(self, verbose:bool=False) -> pd.DataFrame:
		filename = os.path.join(self.root_config, self.fname_ptw_cutoff)
		if not os.path.exists(filename):
			if verbose: print(f"Best parameter file for Pathways does not exist {filename}")
			self.dfbest_cutoffs = pd.DataFrame()
			return pd.DataFrame()

		dfbest_cutoffs = pdreadcsv(self.fname_ptw_cutoff, self.root_config, verbose=verbose)
		self.dfbest_cutoffs = dfbest_cutoffs

		return dfbest_cutoffs

	def save_best_ptw_cutoff(self, dfi, verbose:bool=False) -> bool:
		'''
		config table is build in
		pubmed_taubate_new06_summary_degs_and_pathways_save_config
		'''
		dfbest_cutoffs = self.open_best_ptw_cutoff(verbose=False)

		if dfbest_cutoffs is None or dfbest_cutoffs.empty:
			dfbest_cutoffs = dfi
		else:
			dfbest_cutoffs = dfbest_cutoffs[(dfbest_cutoffs.normalization != self.normalization) & (dfbest_cutoffs.geneset_num != self.geneset_num)]

			if dfbest_cutoffs.empty:
				dfbest_cutoffs = dfi
			else:
				dfbest_cutoffs = pd.concat([dfbest_cutoffs, dfi])

		dfbest_cutoffs = dfbest_cutoffs.drop_duplicates(['case', 'normalization', 'geneset_num', 'quantile', 'med_max_ptw'])
		dfbest_cutoffs.reset_index(inplace=True, drop=True)

		self.dfbest_cutoffs = dfbest_cutoffs

		ret = pdwritecsv(dfbest_cutoffs, self.fname_ptw_cutoff, self.root_config, verbose=verbose)
		return ret


	def get_cfg_best_ptw_cutoff(self, case:str, normalization:str, geneset_num:int, med_max_ptw:str='median',
								verbose:bool=False) -> Tuple[float, float, float, float, float, float, 
															 int, int, int, int, int, 
															 float, float, float, float]:
		if self.dfbest_cutoffs is None or self.dfbest_cutoffs.empty:
			_ = self.open_best_ptw_cutoff(verbose=verbose)
			if self.dfbest_cutoffs is None or self.dfbest_cutoffs.empty:
				if verbose:
					print(f"Could not find best params for {case} {normalization} {geneset_num}")
					print("Houston we have a problem: No best parameter file for Pathways was found.")
					print(">>> run: pubmed_taubate_new05_best_cutoffs_sim_and_save_config.ipynb")
				return self.param_ptw_defaults

		''' remove: (self.dfbest_cutoffs.geneset_num == geneset_num) '''
		dfa = self.dfbest_cutoffs[(self.dfbest_cutoffs.case == case) &
								  (self.dfbest_cutoffs.normalization == normalization) ]
		if dfa.empty:
			print(f"Could not find best params for {case} {normalization}") #{geneset_num}
			return self.param_ptw_defaults

		dfa = dfa[dfa.med_max_ptw == med_max_ptw]
		if dfa.empty:
			return self.param_ptw_defaults

		'''
		 return values:

		 ['case', 'geneset_num', 'normalization', 'parameter', 'quantile',
		  'quantile_val', 'LFC_cut', 'lfc_FDR_cut',
		 'ptw_pval_cut', 'ptw_FDR_cut', 'ptw_min_num_of_degs_cut',
		 'n_pathways', 'n_degs_in_pathways', 'n_degs_in_pathways_mean', 'n_degs_in_pathways_median',
		 'n_degs_in_pathways_std', 'toi1_median', 'toi2_median', 'toi3_median', 'toi4_median'],
		'''
		row = dfa.iloc[0]

		if verbose:
			print( 'quantile, LFC_cut, lfc_FDR_cut, ' +
				   'ptw_pval_cut, ptw_FDR_cut, ptw_min_num_of_degs_cut,  ' +
				   'n_pathways, n_degs_in_pathways, ' +
				   'n_degs_in_pathways_mean, n_degs_in_pathways_median, n_degs_in_pathways_std, ' +
				   'toi1_median, toi2_median, toi3_median, toi4_median')

			print( row['quantile'], row.LFC_cut, row.lfc_FDR_cut, \
				   row.ptw_pval_cut, row.ptw_FDR_cut, row.ptw_min_num_of_degs_cut, \
				   row.n_pathways, row.n_degs_in_pathways, \
				   row.n_degs_in_pathways_mean, row.n_degs_in_pathways_median, row.n_degs_in_pathways_std, \
				   row.toi1_median, row.toi2_median, row.toi3_median, row.toi4_median)

		return row['quantile'], row.LFC_cut, row.lfc_FDR_cut, \
			   row.ptw_pval_cut, row.ptw_FDR_cut, row.ptw_min_num_of_degs_cut, \
			   row.n_pathways, row.n_degs_in_pathways, \
			   row.n_degs_in_pathways_mean, row.n_degs_in_pathways_median, row.n_degs_in_pathways_std, \
			   row.toi1_median, row.toi2_median, row.toi3_median, row.toi4_median


	def get_any_ptw_cutoff(self, case:str, normalization:str, geneset_num:int,
						   quantile:float, verbose:bool=False) -> Tuple[float, float, float, float, float, float,
													                    int, int, int, int, int, 
																		float, float, float, float]:
		if self.dfbest_cutoffs is None or self.dfbest_cutoffs.empty:
			_ = self.open_best_ptw_cutoff()
			if self.dfbest_cutoffs is None or self.dfbest_cutoffs.empty:
				if verbose:
					print("Houston we have a problem: No best parameter file for Pathways was found.")
					print(">>> run: new06_enricher_statistics_and_save_config_table.ipynb")
				return self.param_ptw_defaults

		dfa = self.dfbest_cutoffs[(self.dfbest_cutoffs.case == case) & (self.dfbest_cutoffs.normalization == normalization) &
								  (self.dfbest_cutoffs.geneset_num == geneset_num) & (self.dfbest_cutoffs['quantile'] == quantile) ]
		if dfa.empty:
			return self.param_ptw_defaults

		row = dfa.iloc[0]

		return row['quantile'], row.LFC_cut, row.lfc_FDR_cut, \
			   row.ptw_pval_cut, row.ptw_FDR_cut, row.ptw_min_num_of_degs_cut, \
			   row.n_pathways, row.n_genes_in_pahtways, \
			   row.n_degs_in_pathways_mean, row.n_degs_in_pathways_median, row.n_degs_in_pathways_std, \
			   row.toi1_median, row.toi2_median, row.toi3_median, row.toi4_median

