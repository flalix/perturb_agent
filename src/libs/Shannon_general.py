#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-

'''
Created: 22025/03/16
Updated: 2025/05/28; 2025/03/16
@author: Flavio Lichtenstein
@email:  flalix@gmail.com, flavio.lichtenstein@butantan.gov.br
@local:  Instituto Butantan / CENTD / Molecular Biology / Bioinformatics & Systems Biology
'''


import copy, os, re, random, sys, time
from os import path as osp
from collections import OrderedDict
import numpy as np
import pandas as pd
from datetime import datetime
from typing import Optional, Iterable, Set, Tuple, Any, List

from scipy.stats import t
import warnings

from markdown2 import Markdown

import matplotlib
# matplotlib.use('Agg') # Use backend agg to access the figure canvas as an RGB string and then convert it to an array and pass it to Pillow for rendering.
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import plotly
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

import BioPythonClass
import Sequence as ms
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Basic import *
from util_general import *

basic = Basic()

from graphic_lib import plotly_colors, plotly_colors_proteins

b = BioPythonClass.Basic()

colors2 = ['navy', 'red', 'darkcyan', 'gold', 'mediumvioletred',
		   'green', 'darkgreen', 'orange', 'olivedrab', 'bisque', 'gray', 'pink',
		   'aquamarine', 'darkgreen', 'darksalmon', 'yellowgreen', 'magenta',
		   'lime', 'yellow']
"""
Shannon class:
	define variables
	metadata
	consensus
	calculates
		Shannon Entropy
		Mutual Information
"""
class Shannon:
	def __init__(self, prjName:str, version:str, which_db:str, isProtein:bool, root0:str,
				 nCount:int=500, year_max:int=None, month_max:int=None,
				 badAA_cutoff_perc:float=0.10, 
				 min_seqs:int=3, min_gaps:float=0.20, cutoff_good:float=0.25,
				 date_mask:str='%Y-%m-%d', file_template:str='table_template.html',
				 seaviewPath:str = '../../tools/seaview'):

		'''
			Shannon class:\n
				- main functions:
					- init_vars(self, country, protein=None, protein_name=None)
					- init_vars_year_month(self, country, protein, protein_name, year, month)
					- get and calculate
						- consensus
						- metadata (human, animal + environment)
						- entropy & polymorphism
						- mutations
					- Data:
						- open_DNA
						- open_CDS
						- open_protein
						- open_nCount_sample_DNA
						- open_nCount_sample_CDS
						- open_nCount_sample_protein
					- pseudo-variants get and calc
					- variantes get and calc
					- suggested peptides

		'''
		yearmonth = version[:6]
		"""General varibles:"""
		self.prjName   = prjName
		self.yearmonth = yearmonth
		self.version   = version


		self.which_db  = which_db

		self.pid = None

		self.seaviewPath = seaviewPath

		self.nCount	= nCount

		self.badAA_cutoff_perc = badAA_cutoff_perc
		self.min_gaps = min_gaps
		self.min_seqs = min_seqs
		self.cutoff_good = cutoff_good

		now = datetime.now()

		if year_max is None:
			year_max  = now.year

		if month_max is None:
			month_max  = now.month

		self.year_max = year_max
		self.month_max = month_max

		self.isProtein = isProtein
		self.define_isProtein(isProtein)

		self.country = None
		self.protein = None; self.protein_name = None
		self.year	= None; self.month = None

		self.root0 = root0

		self.set_filelog('entropy.log')

		#-- Define Gisaid filenames -----
		base_filename = 'msa_%s.fasta'

		self.base_filename = base_filename

		# --- root0 ---
		self.root0			= create_dir(root0)
		self.root_ncbi_ens  = create_dir(self.root0, which_db)
		self.root_templates = create_dir(self.root0, 'templates' )
		self.root_ensembl   = create_dir(self.root0, 'ensembl' )

		self.root_data	 = create_dir(self.root0, 'data')
		self.root_figure = create_dir(self.root0, 'figure')

		self.root_db			= create_dir(self.root_ncbi_ens, "dna/")
		self.root_db_CDS		= create_dir(self.root_ncbi_ens, "cds/")
		self.root_db_protein	= create_dir(self.root_ncbi_ens, "protein/")
		self.root_db_animal_env = create_dir(self.root_ncbi_ens, "animal/")

		self.root_db_fasta   = create_dir(self.root_ncbi_ens, 'fasta')
		self.root_db_entropy = create_dir(self.root_ncbi_ens, 'entropy')
		self.root_db_sampled = create_dir(self.root_ncbi_ens, 'sampled')
		self.root_db_result  = create_dir(self.root_ncbi_ens, 'results')
		self.root_db_html	 = create_dir(self.root_ncbi_ens, 'html')
		self.root_db_figure	 = create_dir(self.root_ncbi_ens, 'figure')

		self.file_template = file_template
		self.date_mask	 = date_mask

	def set_path(self, root, _dir=None):
		"""set / create path (dir)"""

	def define_isProtein(self, isProtein):
		"""Define if is protein or dna"""
		self.isProtein = isProtein
		self.dna_protein = 'Protein' if isProtein else 'DNA'

	def set_filelog(self, fname:str, root:str='./logs'):

		filefull = osp.join(root, fname)

		if not osp.exists(root):
			os.mkdir(root)

		self.filelog = filefull

	def init_vars(self, protein:str, protein_name:str):
		"""Init - country & protein/protein_name"""
		self.seq_nuc_consensus = None
		self.mseq = None


		self.set_filelog(f'pipeline_shannon_{protein}.log')

		''' if only protein changed must not load again mseq_CDS '''
		self.mseq_DNA = None
		self.mseq_CDS = None

		self.protein   = protein
		self.protein_name = protein_name

		base_filename = self.base_filename.replace("msa_", "prot_")

		self.mseq_prot = None

		self.dicPiList		  = []
		self.dicNList		   = []
		self.HShannonList	   = []
		self.SeHShannonList	 = []
		# self.HShannonCorrList   = []
		# self.SeHShannonCorrList = []
		self.df_ids_descriptions = None
		self.dfvar = None


	def open_protein(self, fname_prot:str, ignore_msa:bool=False, is_msa:bool=True, 
					 fulfilled:bool=False, verbose:bool=False):

		self.fname_prot	 = fname_prot
		self.fname_prot_ref = self.fname_prot.replace('.fasta', '_ref.fasta')
		self.fname_prot_entropy = "entropy_" + self.fname_prot.replace(".fasta", ".pickle")

		if fulfilled:
			fname_prot = fname_prot.replace(".fasta", "_fulfilled.fasta")

		filename = osp.join(self.root_data, fname_prot)
		self.filename_prot = filename

		if not osp.exists(filename):
			print(f"File does not exist: '{filename}'")
			return False, None

		mseq_prot = ms.MySequence(self.prjName, root=self.root_data)
		ret = mseq_prot.readFasta(filename, ignore_msa=ignore_msa, is_msa=is_msa, showmessage=verbose)

		if not ret or mseq_prot is None or len(mseq_prot.seq_records) == 0:
			stri = f"Could not read: '{filename}'"
			log_save(stri, filename=self.filelog, verbose=True)
			return False, None

		try:
			mseq_prot.seqs = np.array(mseq_prot.seqs)
		except:
			mseq_prot.seqs = []

		self.mseq_prot = mseq_prot

		return ret, mseq_prot


	def fulfill_split_good_bad_seqrec(self, fname_prot:str, 
									  ignore_msa:bool=False, is_msa:bool=True,
									  force:bool=False, verbose:bool=False):

		fname = fname_prot.replace('.fasta', '_fulfilled.fasta')
		filename = osp.join(self.root_data, fname)
		self.filename_prot_alig = filename

		if osp.exists(filename) and not force:
			mseq_prot = ms.MySequence(self.prjName, root=self.root_data)	
			ret = mseq_prot.read_fasta(filename=filename, ignore_msa=False, is_msa=True, verbose=verbose)
			if ret:
				return mseq_prot

		ret, mseq = self.open_protein(fname_prot, fulfilled=False, verbose=verbose)

		if ret==False or mseq is None:
			stri = "Could not read fasta for {fname_prot}"
			log_save(stri, filename=self.filelog, verbose=False)
			return mseq

		n = len(mseq.seqs)
		L = len(mseq.seqs[0])

		if n < self.min_seqs or L == 0:
			stri = "fulfill split_good_bad_seqrec(): mseq insuficient for {self.protein}: nrows={n} and ncols={L}"
			log_save(stri, filename=self.filelog, verbose=False)
			return mseq

		aas = np.array(basic.getSeqAA())

		for i in range(len(mseq.seqs)):
			mseq.seqs[i] = [x if x in aas else '-' for x in mseq.seqs[i]]

		''' sequences should not have more than 25% of gaps, ohterwiser, a bad sequence '''
		good_list = [True if list(mseq.seqs[i]).count('-') / L <= self.cutoff_good else False \
					 for i in range(len(mseq.seqs))]

		''' if equal size or to little '''
		seqs_good		= [mseq.seqs[i]		for i in range(len(good_list)) if good_list[i]]
		seq_records_good = [mseq.seq_records[i] for i in range(len(good_list)) if good_list[i]]

		seqs1 = mseq.fulfill_polymorphic_seqs_with_gaps(seqs_good, isProtein=self.isProtein,
														min_gaps=self.min_gaps, min_seqs=self.min_seqs)

		mseq2 = ms.MySequence(self.prjName, root=self.root_data)
		mseq2.seqs = np.array(seqs1)
		mseq2.seq_records = seq_records_good
		
		for i in range(len(mseq2.seq_records)):
			mseq2.seq_records[i].seq = Seq("".join(mseq2.seqs[i]))

		ret = mseq2.write_fasta(mseq2.seq_records, filename, verbose=True)

		return mseq2


	def entropy_get_calc_protein_cym(self, fname_prot:str, fulfilled:bool=True, ini:int=0, length:int=-1, 
									 h_threshold:int=4, warning:bool=False, return_dic:bool=False, 
									 force:bool=False, verbose:bool=False) -> dict:

		ret, mseq_prot = self.open_protein(fname_prot, fulfilled=fulfilled, verbose=verbose)

		if ret==False or mseq_prot is None:
			stri = "Could not read fasta for {fname_prot}"
			log_save(stri, filename=self.filelog, verbose=False)
			return None


		t1 = datetime.now()
		stime = t1.strftime('%Y-%b-%d, %H h %M min %S sec')

		stri = "%s\tstart\t%s\tNone\tNone\tentropy_get_calc_protein_cym"%(stime, self.protein)
		log_save(stri, filename=self.filelog, verbose=False)

		''' pickle entropy year month '''
		filename = osp.join(self.root_db_entropy, self.fname_prot_entropy)

		if osp.exists(filename) and not force:
			stri = "%s\tend-notfound\t%s\tNone\tNone\tentropy_get_calc_protein_cym"%(stime, self.protein)
			log_save(stri, filename=self.filelog, verbose=False)

			dic = loaddic(self.fname_prot_entropy, path=self.root_db_entropy, verbose=verbose) if return_dic else None
			return dic

		stri = "%s\tend\t%s\tNone\tNone\tentropy_get_calc_protein_cym"%(stime, self.protein)
		log_save(stri, filename=self.filelog, verbose=False) 

		''' sample nCount '''
		stri = "%s\tstart\t%s\tNone\tNone\tcalc_protein_entropy"%(stime, self.protein)
		log_save(stri, filename=self.filelog, verbose=False)
		dic = self.calc_protein_entropy(mseq_prot, ini, length, warning, h_threshold, verbose=verbose)

		stri = "%s\tend\t%s\tNone\tNone\tcalc_protein_entropy"%(stime, self.protein)
		log_save(stri, filename=self.filelog, verbose=False)

		return dic if return_dic else None


	def calc_protein_entropy(self, mseq_prot, ini:int, length:int, 
							 warning:bool=False, h_threshold:int=4, 
							 there_are_gaps:bool=False, verbose=False):
		if length == -1:
			end0 = -1
		else:
			end0 = ini+length

		ret = self.entropy_fast_calcHShannon_1pos(mseq_prot, h_threshold=h_threshold, there_are_gaps=there_are_gaps, verbose=verbose, warning=warning)
		if not ret:
			stri = f"calc_protein_entropy(): error in calcHShannon_1pos() for {self.protein}"
			log_save(stri, filename=self.filelog, verbose=verbose)
			return None

		dicPiList = self.dicPiList
		dicNList  = self.dicNList
		hs		= self.HShannonList
		maxVal	= np.max(hs)

		maxi = len(hs)

		if end0 == -1:
			end = maxi
		elif end0 > maxi:
			end = maxi
		else:
			end = end0

		dicPiList = dicPiList[ini : end]
		dicNList  = dicNList[ ini : end]
		hs		= hs[ini : end]

		labelx = []; _seqx = []; _seqn = []; count = 1
		for jj in range(len(dicPiList)):
			#---- percent
			dic2 = dicPiList[jj]
			val = "; ".join(["%s %.2f%%"%(k, dic2[k]*100) for k in dic2.keys() if dic2[k] != 0])
			_seqx.append(val)

			valks = []
			for k in dic2.keys():
				if dic2[k] == 0: continue

				if dic2[k] == 1:
					valks = [k]
					break

				valks.append("%s %.2f%%"%(k, dic2[k]*100) )

			labelx.append(str(count) + '-' + "; ".join(valks))

			#---- n or #
			dic3 = dicNList[jj]
			val = "; ".join(["%s %d"%(k, dic3[k]) for k in dic3.keys() if dic3[k] != 0])
			_seqn.append(val)

			count += 1

		seqx = np.arange(ini, end)
		n		= len(hs)
		nSamples = len(mseq_prot.seqs)

		df = pd.DataFrame({'x': seqx, 'y': hs, 'aas': _seqx, 'nns': _seqn})

		dic={}
		dic[self.protein] = {"df": df, "nSamples": nSamples, "maxVal": maxVal,
							 "dicPiList": dicPiList, "dicNList": dicNList, "df_ids_descriptions": self.df_ids_descriptions}

		ret = pddumpdic(dic, self.fname_prot_entropy, path=self.root_db_entropy, verbose=verbose)
		return dic


	def open_protein_entropy(self, verbose:bool=False):
		filename = osp.join(self.root_db_entropy, self.fname_prot_entropy)

		if not osp.exists(filename):
			stri = f"entropy file could not find '{filename}'"
			log_save(stri, filename=self.filelog, verbose=verbose)
			return None

		return loaddic(self.fname_prot_entropy, path=self.root_db_entropy, verbose=verbose)



	# add_gaps --> there_are_gaps
	def entropy_fast_calcHShannon_1pos(self, mseq, h_threshold=4, there_are_gaps:bool=True, 
									   verbose=False, warning=False):

		# if not fulfill the gaps, there may be bad chars
		self.dicPiList	= []; self.dicNList = []
		self.HShannonList = []; self.SeHShannonList = []
		self.HShannonCorrList = []; self.SeHShannonCorrList = []
		self.df_ids_descriptions  = None

		if self.dna_protein == "DNA":
			valid = b.getDnaNucleotides()
		elif self.dna_protein == "Protein" or self.dna_protein == "Amino acid":
			valid = b.getAA()
		else:
			stri = "Error: only DNA and Protein tables are accepted (type = 'DNA' or 'Protein') for %s %s %d/%d - %s"%(self.country, self.protein_list, self.year, self.month, self.variant_subvar)
			log_save(stri, filename=self.filelog, verbose=True)
			return False

		''' there can be gaps !!! '''
		if there_are_gaps: valid.append('-')
		'''--- add stop codons ---'''
		valid.append('*')

		#---------------- metadata for filtering - dic['metadata'] = dfn ----------------------
		descs = [seqrec.description  for seqrec in mseq.seq_records]
		ids   = [desc.split('||')[0] for desc   in descs]

		if len(ids) == 0:
			stri = "Error?? No records found in fast calcHShannon_1pos for %s %s %d/%d - %s"%(self.country, self.protein_list, self.year, self.month, self.variant_subvar)
			log_save(stri, filename=self.filelog, verbose=verbose)
			return False

		dfa = pd.DataFrame(data={'id': ids, 'desc': descs})
		dfa['i'] = dfa.index
		dfa.index = dfa['id']
		self.df_ids_descriptions = dfa

		seqs = np.array(mseq.seqs)
		nrow, ncol = seqs.shape

		''' remove stop codon if last char '''
		if seqs[0][-1] == '*':
			ncol -= 1

		t1 = datetime.now()
		if verbose: print("Start calc Shannon 1 letter")

		# --- looping each column ---------
		msa_error_msg = False
		for jj in range(ncol):
			dic = char_frequency(seqs[:,jj])

			badkeys = [c for c in dic.keys() if c not in valid]
			for badk in badkeys:
				del(dic[badk])

			''' gaps are allowed
				'-' must be the major count
			'''
			if len(dic) > 1 and '-' in dic.keys():
				num_c = dic['-']
				total = np.sum(list(dic.values()))

				''' if gap is less than 5% '''
				if num_c/total < 0.05:
					del(dic['-'])

			# if self.protein == 'S': print(jj, dic)

			counts = np.array(list(dic.values()))
			dicPi  = OrderedDict(); dicN=OrderedDict()

			if len(counts) == 1:
				k = list(dic.keys())[0]
				dicPi[k] = 1
				dicN[k]  = counts[0]

				self.dicPiList.append(dicPi)
				self.dicNList.append(dicN)
				self.HShannonList.append(0)
				self.SeHShannonList.append(0)
				# print(">>> jj %d - len %d"%(jj, len(self.dicPiList) ))
				continue

			n	  = np.sum(counts)
			percs  = counts / n
			hShan = 0; varPi = 0

			kcount = 0;
			for k in dic.keys():
				pi = percs[kcount]
				dicPi[k] = pi
				dicN[k]  = counts[kcount]
				if (pi != 0.) and (pi != 1.):
					hShan -= pi * np.log2(pi)
					varPi += (1. + np.log2(pi))**2  * pi * (1.-pi)

				kcount += 1

				''' warnings '''
				if warning:
					if hShan > h_threshold:
						print('>>> jj', jj, hShan)
						print("percs", percs, '\n')

			''' Roulston Eq. 13 '''
			'''
			B = len(dic.keys())
			if hShan > 0:
				hShanCorr = hShan + ((B-1)/ (2*n))
			else:
				hShanCorr = 0
			if hShanCorr > 4:
				print('>>> jj1', jj, hShan)self.dic_meta_country.keys()
				print('>>> jj2', jj, "B", B, "n", n, "hShanCorr", hShanCorr)
				print("percs", percs, '\n')

			VarCorr = 0
			# Roulston Eq. 40
			for key in dicPi.keys():
				pi = dicPi[key]
				if (pi != 0.) and (pi != 1.):
					VarCorr += (np.log2(pi)+hShan)**2 * pi * (1-pi)
			'''

			SeHShannon = 0 if hShan == 0 else hShan * np.sqrt(varPi/n)
			# SeHShannonCorr = np.sqrt(VarCorr/n)

			self.dicPiList.append(dicPi)
			self.dicNList.append(dicN)
			self.HShannonList.append(hShan)
			self.SeHShannonList.append(SeHShannon)
			# self.HShannonCorrList.append(hShanCorr)
			# self.SeHShannonCorrList.append(SeHShannonCorr)
			# print(">>> jj %d - len %d"%(jj, len(self.dicPiList) ))

		if verbose:
			t2 = datetime.now()
			print("Calc fast calcHSannon 1pos:", round( (t2-t1).microseconds / 1000), "ms")

		return True


	def entropy_plot(self, hs:List, nrow:int, factor:int=1000, title:str='Entropy', title_region:str='',
					xlabel:str="protein residues", color:str='navy',
					fontsize:int=14, fontsize_title:int=15, show_grid:bool=False, ggplot_like:bool=True,
					zoomRegion:List=[], figsize:tuple=(8,6), onlycalc:bool=False,
					save_plot:bool=True, filetype:str='png', dpi:int=300, 
					orientation:str='landscape', facecolor:str='w', edgecolor:str='k', verbose:bool=False):

		# zoomRegion = [min, max]
		if isinstance(zoomRegion, list) and len(zoomRegion) == 2:
			mini = zoomRegion[0]; maxi=zoomRegion[1]
		else:
			mini = 0; maxi = len(hs)

		seqx = np.arange(mini+1, maxi+1)
		hs = hs[mini:maxi]
		L = len(hs)

		totalH = np.sum(hs)
		meanH  = totalH * factor / L

		if onlycalc:
			return None, totalH, meanH, L

		fig = plt.figure(figsize=figsize, dpi=dpi, facecolor=facecolor, edgecolor=edgecolor)
		plt.plot(seqx, hs, color=color)

		if factor==1:
			unit = 'bits'
			ylabel:str='H (bits)'
		elif factor==1000:
			unit = 'mbits'
			ylabel:str='H (mbits)'
		else:
			unit = "Factor is 1 or 1000 ???"
			ylabel:str='H (bits) ???'


		if title_region is not None or title_region == '':
			title += f' - {title_region}'
		title += f"\nTotal H = {totalH:.2f} bits, mean H = {meanH:.2f} {unit}, {nrow} sequences, length = {L}"
		if zoomRegion != []:
			title += f'\nzoom between {mini} and {maxi}'

		if ggplot_like: plt.style.use('ggplot')
		if show_grid: plt.grid()
		plt.xlabel(xlabel, fontsize=fontsize)
		plt.ylabel(ylabel, fontsize=fontsize)
		plt.title(title, fontsize=fontsize_title)

		if save_plot:
			filefig = title_replace(title) + '.' + filetype
			filefig = osp.join(self.root_db_figure, filefig)
			plt.savefig(filefig, format=filetype, dpi=dpi, facecolor='w', edgecolor='w', orientation=orientation);
			if verbose: print("File saved: '%s'"%(filefig))

		return plt, totalH, meanH, L

	def calc_all_entropy(self, verbose:bool=False):

		fname0 = copy.deepcopy(self.fname_prot0)

		for species_name in self.species_names_good:
			fname_spec = fname0.replace('.fasta', f'_species_{species_name}.fasta')
			print(">>>", species_name, end=': ')
			
			self.entropy_get_calc_protein_cym(fname_prot=fname_spec, fulfilled=False, ini=0, length=-1, 
											  h_threshold=4, warning=False, return_dic=False,
											  force=False, verbose=verbose)
			

	def map_polimorphic_regions(self, zoomRegion:list, title_region:str, verbose:bool=False):

		fname0 = copy.deepcopy(self.fname_prot0)

		df_mut_list = []
		for species_name in self.species_names_good:
			print(">>", species_name)
			fname_spec = fname0.replace('.fasta', f'_species_{species_name}.fasta')
			
			dic = self.entropy_get_calc_protein_cym(fname_prot=fname_spec, fulfilled=False, ini=0, length=-1, 
												h_threshold=4, warning=False, return_dic=True,
												force=False, verbose=verbose)
			
			self.dic = dic
			
			key = list(dic.keys())[0]
			dfgb = dic[key]['df']
			nrow = dic[key]['nSamples']

			if zoomRegion == []:
				df2 = dfgb
			else:
				mini = zoomRegion[0]
				maxi = zoomRegion[1]
				lista = np.arange(mini, maxi)
				
				df2 = dfgb[dfgb.index.isin(lista)].copy()
				df2.index = np.arange(len(df2))
			
			df_list=[]
			dic={}; icount=-1
			for i in range(len(df2)):
				row = df2.iloc[i]
				if row.y > 0.:
					dfa = pd.DataFrame(row).T
					df_list.append(dfa)
			
					ass = row.aas
					mut_list = ass.split('; ')
					for i_mut, mut in enumerate(mut_list):
						mat = mut.split(' ')
						aa = mat[0]
						perc = mat[1]
						
						icount += 1
						dic[icount] = {}
						dic2 = dic[icount]
			
						dic2['species'] = species_name
						dic2['protein'] = 'GBA1'
						dic2['posi'] = row['x'] + 1
						dic2['i_mut'] = i_mut+1
						dic2['aa_perc'] = mut
						dic2['aa'] = aa
						dic2['perc'] = perc
		
			if df_list == []:
				print(f"Nothing found for {species_name}")
				continue
				
			dfpoli = pd.concat(df_list)

			if title_region is None:
				fname = f'polymorphism_{species_name}_for_all_protein.tsv'
			else:
				fname = f'polymorphism_{species_name}_for_{title_region}.tsv'

			pdwritecsv(dfpoli, fname, self.root_db_result, verbose=verbose)
		
			dfmut = pd.DataFrame(dic).T
			df_mut_list.append(dfmut)

		dfmut = pd.concat(df_mut_list)

		df_piv = pd.pivot_table(dfmut, values=['aa_perc'], index=['species', 'i_mut'], columns='posi', aggfunc='sum')
		df_piv = df_piv.fillna('')
		df_piv = df_piv.reset_index()
		cols = list(df_piv.columns)
		cols[0] = 'species'
		cols[1] = 'poly'

		for i_col in np.arange(2, len(cols)):
			cols[i_col] = cols[i_col][1]
		
		df_piv.columns = cols

		if title_region is None:
			fname = f'pivot_table_for_all_protein.tsv'
		else:
			fname = f'pivot_table_for_{title_region}.tsv'
		pdwritecsv(df_piv, fname, self.root_db_result, verbose=verbose)

		return df_piv, dfmut


	def calc_and_plot_entropy(self, species_name:str, title_region:str=None, zoomRegion:list=[], factor:int=1000,
								xlabel:str="protein residues", ylabel:str='H (bits)', color:str='navy',
								fontsize:int=14, fontsize_title:int=15,
								figsize:tuple=(8,6), onlycalc:bool=False,
								save_plot:bool=True, filetype:str='png', dpi:int=300, 
								orientation:str='landscape', facecolor:str='w', edgecolor:str='k', verbose:bool=False):


		fname0 = copy.deepcopy(self.fname_prot0)
		fname_spec = fname0.replace('.fasta', f'_species_{species_name}.fasta')

		dic = self.entropy_get_calc_protein_cym(fname_prot=fname_spec, fulfilled=False, ini=0, length=-1, 
											h_threshold=4, warning=False, return_dic=True,
											force=False, verbose=verbose)

		key = list(dic.keys())[0]
		dfgb = dic[key]['df']
		nrow = dic[key]['nSamples']

		hs = dfgb.y.to_list()
		posi_h = [i for i in range(len(hs)) if hs[i] != 0]

		title=f'Entropy for {species_name}'
	
		plt, totalH, meanH, L = self.entropy_plot(hs, nrow, zoomRegion=zoomRegion, factor=factor, 
												title=title, title_region=title_region,
												xlabel=xlabel, ylabel=ylabel,
												fontsize=fontsize, fontsize_title=fontsize_title,
												figsize=figsize, onlycalc=onlycalc,
												save_plot=save_plot, filetype=filetype, dpi=dpi,
												orientation=orientation, facecolor=facecolor, verbose=verbose)
	
		return plt, totalH, meanH, L, hs, nrow, posi_h



	def plot_all_entropy(self, zoomRegion:list=[], title_region:str=None, factor:int=1000,
						xlabel:str="protein residues", ylabel:str='H (bits)', color:str='navy',
						fontsize:int=14, fontsize_title:int=15,
						figsize:tuple=(8,6), onlycalc:bool=False,
						save_plot:bool=True, filetype:str='png', dpi:int=300, 
						orientation:str='landscape', facecolor:str='w', edgecolor:str='k', verbose:bool=False):

		dic={}; icount=-1
		for species_name in self.species_names_good:
			print(">>>", species_name, end=': ')
			
			plt, totalH, meanH, L, hs, nrow, posi_h = \
				self.calc_and_plot_entropy(species_name, zoomRegion=zoomRegion, title_region=title_region, 
											factor=factor, xlabel=xlabel, ylabel=ylabel,
											fontsize=fontsize, fontsize_title=fontsize_title,
											figsize=figsize, onlycalc=onlycalc,
											save_plot=save_plot, filetype=filetype, dpi=dpi,
											orientation=orientation, facecolor=facecolor, verbose=verbose)

						
	
			icount += 1
			dic[icount] = {}
			dic2 = dic[icount]

			dic2['species'] = species_name
			dic2['protein'] = self.protein
			dic2['zoomRegion'] = zoomRegion
			dic2['title_region'] = title_region
			dic2['totalH'] = totalH
			dic2['meanH'] = meanH
			dic2['length'] = L
			dic2['hs'] = hs
			dic2['n_samples'] = nrow
			dic2['posi_h'] = posi_h

	
		dfh = pd.DataFrame(dic).T
	
		if title_region is None:
			fname = f'entropy_summary_for_{title_region}.tsv'
		else:
			fname = f'entropy_summary_for_all_protein.tsv'

		pdwritecsv(dfh, fname, self.root_db_result, verbose=verbose)

