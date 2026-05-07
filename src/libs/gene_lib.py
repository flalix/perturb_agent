#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2023/06/25
# Udated  on 2023/08/28 - 2023/08/16
# @author: Flavio Lichtenstein
# @local: Bioinformatics: CENTD/Molecular Biology; Instituto Butatan

import os
from pathlib import Path
import numpy as np
import pandas as pd
from typing import List, Tuple #  Optional, Iterable, Set, Any

# import mygene
from libs.Basic import create_dir, pdwritecsv, pdreadcsv, to_roman_numeral

class Gene(object):
	def __init__(self, root0:Path, root_project:Path):

		self.root0 = root0
		self.root_colab = create_dir(root0, 'colab')
		self.root_project = root_project

		self.root_refseq = create_dir(self.root_colab, 'refseq')

		files = os.listdir(self.root_refseq)
		if len(files)==0:
			print(f"No files found in the refseq directory {self.root_refseq}.")
			raise Exception("\n------- stop -------\n")
		
		'''
		fileMeta	= "gene_de_para_metacore_clarivate_python.csv"
		fileBadUniprotSymbs = "bad_uniprot_symbols.txt"

		Entrez.email = email
		econv = Entrez_conversions(email, root_refseq, fileBadUniprotSymbs)
		econv.badunis
		'''

		self.fname_refseq_ncbi = 'refseq_gene_ncbi.tsv'
		''' or,  'refseq_homo_sapiens_gene_info_edited.tsv' '''
		self.fname_ncbi_syn	= 'refseq_synonyms_ncbi.tsv'
		self.fname_df_ref_seq  = 'refseq_df_dict.tsv'

		self.fname_conv_uniprot = "conversion_symbol_uniprot_entrez.tsv"
		self.fname_uniprot = 'uniprot_symbol.tsv'
		self.fname_my_gene = "refseq_my_gene.tsv"
		self.fname_my_gene_synonym = "my_gene_synonym.tsv"
		self.fname_dic_change	  = "synonym_table_syn_plus_change_tuple.tsv"

		self.fname_refseq_gene_info  = 'gene_info_ncbi_refseq.tsv'

		self.fname_biomart_entrez = 'biomart_entrez.tsv'
		self.fname_biomart_protein = 'biomart_swissprot.tsv'

		self.fname_hgnc_set	= 'hgnc_complete_set.tsv'
		self.fname_refseq_hgnc = 'biomart_hgnc.tsv'
		self.fname_refseq_omim = 'biomart_omim_mutation.tsv'

		'''
			gene annotation

			NCBI: Genome assembly GRCh38
				https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
				https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/
			
			Release 90:
				https://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/

			Ensemble: https://www.ensembl.org/Homo_sapiens/Info/Index 
					  https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/


		'''
		self.fname_gff_ensembl0 = 'Homo_sapiens.GRCh38.111.chr.gff3'
		self.fname_gff_ensembl  = 'Homo_sapiens_GRCh38_111_chr_gff3_final.tsv'
		self.fname_gff_ensembl0_tsv = 'Homo_sapiens.GRCh38.111.chr.gff3.tsv'
		'''
		  uniprot API queries: 
			https://www.uniprot.org/help/programmatic_access
			https://www.uniprot.org/help/api_queries
			https://www.uniprot.org/help/query-fields
		  ensembl FTP
			https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
			https://www.ensembl.org/info/data/ftp/index.html

		  tools:
			https://pypi.org/project/eutils/
			https://pypi.org/project/pyensembl/

		  tutorials:
			https://www.ncbi.nlm.nih.gov/books/NBK25500/
			https://www.ensembl.org/info/website/tutorials/sequence.html


		'''
		self.fname_gff_selected  = 'Homo_sapiens_GRCh38_111_chr_gff3_selected.tsv'
		self.df_gfff, self.df_gff_full = None, None

		self.gene_term_list = [('IL8', 'CXCL8'), ('von Willebrand factor', 'VWF'),
							   ('factor 2', 'F2'), ('factor 3', 'F3'), ('factor 4', 'F4'), ('factor 5', 'F5'),
							   ('factor 6', 'F6'), ('factor 7', 'F7'), ('factor 8', 'F8'), ('factor 9', 'F9'),
							   ('factor 10', 'F10'), ('factor 11', 'F11'), ('factor 12', 'F12'), ('factor 13', 'F13')]  # F13 = F13A1, F13B

		self.gene_term_list_rom = [('factor I', 'F1'), ('factor II', 'F2'), ('factor III', 'F3'), ('factor IV', 'F4'),  ('factor V', 'F5'),
								   ('factor VI', 'F6'), ('factor VII', 'F7'),  ('factor VIII', 'F8'), ('factor IX', 'F9'),
								   ('factor X', 'F10'), ('factor XI', 'F11'), ('factor XII', 'F12'), ('factor XIII', 'F13'), ]


		self.factor_nums = [str(x+1) for x in range(99)]
		self.factor_roman_nums = [to_roman_numeral(int(x)) for x in self.factor_nums]
		self.dic_roman = {self.factor_roman_nums[i]: self.factor_nums[i] for i in range(len(self.factor_nums))}

		self.dic_gene_term_list = {}
		for gene_desc, gene in self.gene_term_list:
			self.dic_gene_term_list[gene] = gene_desc


		self.dic_gene_term_list_rom = {}
		for gene_desc, gene in self.gene_term_list_rom:
			self.dic_gene_term_list_rom[gene] = gene_desc

		self.dic_gene_mnem = {}
		self.dic_gene_mnem['F'] = 'F'
		self.dic_gene_mnem['factor'] = 'F'
		self.dic_gene_mnem['Factor'] = 'F'
		self.dic_gene_mnem['G'] = 'G'
		self.dic_gene_mnem['group'] = 'G'
		self.dic_gene_mnem['Group'] = 'G'
		self.dic_gene_mnem['caspase'] = 'CASP'
		self.dic_gene_mnem['Caspase'] = 'CASP'
		self.dic_gene_mnem['metaloproteinase'] = 'MMP'
		self.dic_gene_mnem['Metaloproteinase'] = 'MMP'

		self.dic_gene_mnem['groups'] = 'G'
		self.dic_gene_mnem['Groups'] = 'G'
		self.dic_gene_mnem['factors'] = 'F'
		self.dic_gene_mnem['Factors'] = 'F'
		self.dic_gene_mnem['caspases'] = 'CASP'
		self.dic_gene_mnem['Caspases'] = 'CASP'
		self.dic_gene_mnem['metaloproteinases'] = 'MMP'
		self.dic_gene_mnem['Metaloproteinases'] = 'MMP'


		self.lista_gene_many = ['MMP', 'IL', 'CCL','CXCL', 'KLK' ]

		'''self.genes_dash = ['CXCL', 'CCL', 'IL', 'IGF', 'STAT', 'SOCS', 'MMP', 'IFN', 'ESR', 'CSF', 'RIG', 'CTLA', 'AKT', 'IFNB', 'IFNA',
							  'SERPINE', 'PAI', 'F', 'KRT', 'NOD']'''

		drop_symbs = ['Yes', 'COPD', 'C02', 'FBS', 'NHS', 'P1', 'P2', 'P3', 'CTPA', 'JAMA', 'ANA', 'AIM', 'EMA', 'NLR',
					  'JH', 'CPAP', 'NOD', 'COX', 'Cox', 'MAP', 'PRR', 'PCI', 'MRC', 'ICH', 'LDL', 'HDL', 'VLHDL', 'TG',
					  'TH1', 'TH2', 'ABS', 'ASP', 'CAP',  'CRA', 'HCP', 'INI', 'IRR', 'MRSA', 'NPI', 'SMR', 'SAS', 'Sara',
					  'MB', 'NB', 'NHL', 'ALT', 'AST', 'AT', 'RT', 'RTX', 'RTx', 'CRT', 'DIC', 'PCT', 'DOHH', 'CDR', 'CNR',
					  'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'not', 'NOT', 'polymerase', 'NBCCS', 'BCC',
					  'ANOVA', 'SD', 'SE', 'SEM', 'CI', 'IV', 'EV', 'SI', 'FISH', 'NTCP', 'RPC', 'IHG',
					  'SARS', 'Tech', 'envelope', 'kallikrein', 'peptidase', 'homeobox', 'B7', 'H3', 'ASA',
					  'mol', 'Delta', 'Protease', 'MIR', 'P', 'Dicer', 'TO', 'T', 'L1',
					  'POST', 'TOP', 'IgA', 'IgG', 'BT', 'EWS', 'II', 'KD', 'UK', 'KILLER', 'SERA', 'MRI', 'ABMR',
					  'EFS', 'OS', 'CR', 'PR',  'ML', 'ROS', 'OA', 'AA', 'TS', 'ADC', 'MV', 'FSH', 'LH', 'MS',
					  '1A',  'AR', 'AS', 'AU', 'B1', 'BG', 'BK', 'BS',  'CD', 'CF', 'CG', 'CH', 'CK', 'CL', 'CN',
					  'CO', 'CP', 'CT', 'D1', 'D2', 'D3',  'D4', 'D9', 'DD', 'DF', 'DI', 'DL', 'DS', 'DT',
					  'E1', 'E2', 'E3', 'E5', 'EL', 'EP', 'ER', 'ES', 'ET', 'El', 'FA',  'FR', 'FS', 'FU',
					  'GC', 'GD', 'GO',  'GS', 'GT', 'Gi',
					  'H', 'H1', 'H4', 'H5', 'HC', 'HD', 'HH', 'HI', 'HL', 'HR', 'HS', 'IB',
					  'ID', 'IF', 'IP', 'J4', 'K4', 'K5', 'KO', 'Ki', 'L2', 'L3', 'LI', 'M1', 'MA', 'ME', 'MG',
					  'MI', 'MK', 'MR', 'MT', 'NA', 'ND', 'NE', 'NF', 'NM', 'NN', 'NS', 'NT', 'O1', 'O2', 'O3',
					  'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'OB', 'OC',
					  'PA', 'PB', 'PC', 'PH', 'PI', 'PK', 'PL', 'PN', 'Ps', 'QC', 'R1', 'R2', 'RB', 'RC', 'RD',
					  'RL', 'RN', 'RS', 'S2', 'S3', 'S6', 'SA', 'SC', 'SF', 'SK', 'SN', 'SP', 'SR', 'SS', 'T1',
					  'T2', 'T3', 'T4', 'T6', 'TC',  'TG', 'TK', 'VA', 'VP', 'WD', 'WS', 'Y1', 'ZF',
					  'D54', 'DLT', 'DR4', 'DR5', 'Era', 'E11', 'E12', 'E13', 'E14',
					  '1D', 'A1', 'A2', 'A5', 'AF', 'A3', 'A4', 'A6', 'A7', 'A8', 'A9', 'a1', 'a3',
					  'AC', 'AH', 'AK', 'AM', 'AN', 'AO', 'B2', 'B3', 'B4',
					  'B5', 'B8', 'B9', 'BF', 'BH', 'BI', 'BN', 'BP', 'BR',  'D5', 'D6', 'D8',
					  'DC', 'DK', 'DM', 'DP', 'E4', 'EE', 'EG', 'EI', 'EK', 'FB', 'FD', 'FE',
					  'FH', 'FI', 'FL', 'FM', 'FP', 'FX', 'FY', 'GA', 'GE', 'GK', 'GM',
					  'GO', 'GP', 'GR', 'GW', 'H2', 'H6', 'HE', 'HF', 'HK', 'Hc', 'I1', 'I3', 'IH', 'IN',
					  'IO', 'J2', 'JN', 'K1', 'K2', 'K3', 'K8', 'KC', 'KF', 'KL', 'KN', 'KS', 'KY', 'L4', 'L5',
					  'L8', 'LD', 'LE', 'LF', 'LP', 'LS', 'LT', 'LU', 'LW', 'La', 'M8', 'MN', 'MP', 'MU', 'MX',
					  'NC', 'NP', 'NU', 'OF', 'OK', 'OM', 'ON', 'OT', 'PE', 'PG', 'PM', 'PP', 'PT', 'PZ',
					  'Po', 'Pr', 'QM', 'QV', 'R3', 'R4', 'RG', 'RH', 'RK', 'RM', 'RO', 'RP', 'RY', 'S1', 'S4',
					  'S5', 'S7', 'S8', 'S9', 'SW', 'T9', 'TH', 'TM', 'TN', 'TP', 'TR', 'TX', 'UI', 'UP', 'UR',
					  'V1', 'V2', 'VF', 'VH', 'VL', 'VN', 'VT', 'WI', 'WR', 'WW', 'XB', 'XO', 'XR', 'XS',
					  'hang', 'helicase', 'hole', 'homodimer', 'mu2', 'pol', 'quo', 'tTG', 'vRNA', 'zeta',
					  'Yes', 'PBS', 'MHC', 'PDB', 'NTD', 'ICD', 'Rho', 'MIS', 'NET', 'LPS', 'ESR',
					  'HEK', 'TCR', 'MFI', 'CRL', 'PCA', 'RBM', 'PLT',  'CAR', 'Paul', 'CRS',
					  'TEM', 'BAL', '8B', 'S12', 'IgA', 'IgG', 'IgE', 'S5A', 'MMR', 'UTR', 'ALP', 'VLP', 'CTA', 'MAS',
					  'JP', 'OG', 'X3', '10C', '16A', '19A', '42A', 'A10', 'env', 'mucin', 'myosin', 'neu', 'pk', 'syndecan',
					  '105A', '105B', '1A1', '1C7', '2B4', '2F1', '42C', '4EBP1', '51A', '5A3', '5T4', '7SL', '9F',
					  'A11', 'A12', 'A121', 'A13', 'A135', 'A14', 'A15', 'A16', 'A17', 'A18', 'A19', 'A1B',
					  'A20', 'A21', 'A22', 'A23', 'A24', 'A25', 'A26', 'A27', 'A28', 'A29', 'A2AP', 'A2a',
					  'A30', 'A33', 'A43', 'A49', 'cell1', 'cell3', 'chymase', 'cob', 'cofilin', 'dUTPase', 'emb',
					  'endocan', 'feat', 'fractalkine', 'fru', 'gC1qR', 'glypican', 'hAPN', 'hexokinase', 'irisin',
					  'leB', 'lnc', 'lymphokine', 'memB', 'morpheus', 'mover', 'nNIF', 'nNOS', 'nephrin', 'omentin',
					  'pRb', 'pT7', 'rols', 'slim', 'snRNP', '104p', '11B6', '170A',
					  '1C5', '1F5', '1G5', '24p3', '24p3R', '28H', '2A9', '2DD', '2P1', '3D3', '3DL2', '3G2', '3OST',
					  '4F2', '4F5', '4F9', '53BP1', '5B10', '5HTT', '60H', '6A1', '6PGD', '7A5', '7B2', '7B4', '7SK',
					  '8D6', '9G8', 'A170', 'A190', 'A2aR', 'A2c', 'A3A', 'A3AR', 'A3B', 'mDia1', 'mER', 'mL37', 'mL39',
					  'mL50', 'mL54', 'mago', 'miPEP155', 'mili', 'nc886', 'neuroserpin', 'nrs', 'onzin', 'opiorphin',
					  'pG2', 'pS2', 'pT1', 'pT3', 'pVHL', 'phosphacan', 'pp62', 'rbc', 'sentan', 'slv', 'sulfatase',
					   'sws', 'syd', 'tryptaseB', 'uORF2',  'wim', 'wts', 'Proteinase']

		'''
			List of RNAs: https://en.wikipedia.org/wiki/List_of_RNAs

			GH growth hormone

		  - DOHH - directly observe hand hygiene or Deoxyhypusine Hydroxylase
		  - Aggressive Behaviour Scale (ABS)
		  - CAP community-acquired pneumonia
		  - COPD - Chronic Obstructive Pulmonary Disease or Pulmonary Disease, Chronic Obstructive, Severe Early-Onset (it is a GENE !)
		  - CDR - Clinical Dementia Rating Scale
		  - CRA - carbapenem-resistant Acinetobacter species
		  - HCP - health care professional
		  - IRR - incident rate ratio
		  - MRSA - methicillin-resistant Staphylococcus aureus
		  - NPI-NH - NPI - Neuropsychiatric Inventory nursing home version
		  - NLR - neutrophil/lymphocyte ratio
		  - SMR - standardized mortality ratios
		  - ASL - Frosinone Bruno Sala // Argininosuccinate Lyase
		  - ASP - Arenare Golgi Redaelli
		  - Sezione Invecchiamento CNR // CNR1 - Cannabinoid Receptor 1
		  - IHG Guidonia RM // IHG1 - Iris Hypoplasia With Glaucoma 1
		  - Yes - Subjects level // YES1 - YES Proto-Oncogene 1, Src Family Tyrosine Kinase
		  - ASA - acetylsalicylic aci
		  - DIC - Disseminated Intravascular Coagulation
		  - PCT - Procalcitonin time

		  - LDL, HDL, VLHDL, Triglicerides (TG) - they are not gene
		  - TH1, TH2 - T helper type
		  - ICH - Intracerebral hemorrhage
		  - MRC - Department MRC
		  - PCI - percutaneous coronary intervention (PCI)
		  - NOD - NOD COVID India Study
		  - MAP - maps
		  - COX - COX regression
		  - PRR - populist radical right // NECTIN1 (official name)
		  - CPAP - continuous positive airway (therapy) // CENPJ (official name)
		'''

		self.drop_symbs = np.unique(drop_symbs)

		self.drop_symbs_dic = {}
		for symb in drop_symbs:
			self.drop_symbs_dic[symb] = 1

		'''------------ ambiguous symbols --------------'''
		self.ambiguous_symbol_dic = {}
		self.ambiguous_symbol_set = []

		for symb, term in self.ambiguous_symbol_set:
			self.ambiguous_symbol_dic[symb] = term


		self.change_tuple = \
			[('C2orf71', 'PCARE'), ('COLCA2', 'POU2AF3'), ('FAM155A', 'NALF1'), ('FAM189A2', 'ENTREP1'), ('FAM198A', 'GASK1A'),  ('FAM198B', 'GASK1B'),
			('C12orf80', 'LINC02874'), ('C2orf40', 'ECRG4'), ('C2orf71', 'PCARE'), ('C8orf31', 'LY6S-AS1'), ('C9orf24', 'SPMIP6'), ('C9orf3', 'AOPEP'),
			('CSTF3-AS1', 'CSTF3-DT'), ('FAM189A1', 'ENTREP2'), ('FAM46A', 'TENT5A'),  ('FAM46C', 'TENT5C'), ('FAM49A', 'CYRIA'), ('FAM84B', 'LRATD2'),
			('GRASP', 'TAMALIN'), ('HIST1H2AC', 'H2AC6'), ('HIST2H2AA3', 'H2AC18'), ('HIST2H2AA4', 'H2AC19'), ('FAM189A2', 'ENTREP1'), ('HIST2H3PS2', 'H3-7'),
			('KIAA1024', 'MINAR1'), ('KIAA1257', 'CFAP92'), ('LINC01185', 'REL-DT'), ('LINC01197', 'LETR1'), ('METTL7A', 'TMT1A'), ('PER4', 'PER3P1'),
			('AAED1', 'PRXL2C'), ('AARS', 'AARS1'), ('ACPP', 'ACP3'), ('ADPRHL2', 'ADPRS'), ('ADSS', 'ADSS2'), ('ADSSL1', 'ADSS1'), ('AES', 'TLE5'),
			('AKR1C8P', 'AKR1C8'), ('ALPPL2', 'ALPG'), ('ALS2CR12', 'FLACC1'), ('ANKRD20A2', 'ENSG00000279381'), ('APOPT1', 'COA8'), ('ARMC4', 'ODAD2'),
			('ARNTL', 'BMAL1'), ('ARNTL2', 'BMAL2'), ('ASNA1', 'GET3'), ('ATP5MD', 'ATP5MK'), ('ATP5MPL', 'ATP5MJ'), ('ATP5S', 'ATP5F1E'),
			('AURKAPS1', 'AURKAP1'), ('AZIN1-AS1', 'MAILR'), ('BHMG1', 'MEIOSIN'), ('BMPR1B-AS1', 'BMPR1B-DT'), ('BREA2','LINC02878'),
			('C10orf111','RPP38-DT'), ('C10orf113','NEBL'), ('C10orf142','LINC02881'), ('C10orf25','ZNF22-AS1'), ('C10orf82','SPMIP5'),
			('C10orf91','LINC02870'), ('C10orf99','GPR15LG'), ('C11orf1','C11orf1'), ('C11orf44','LINC02873'), ('C11orf45','KCNJ5-AS1'),
			('C11orf49','CSTPP1'), ('C11orf53','POU2AF2'), ('C11orf80','TOP6BL'), ('C11orf88','HOATZ'), ('C11orf94', 'FREY1'), ('C11orf1','CFAP68'),
			('C11orf95','ZFTA'), ('C12orf10','MYG1'), ('C12orf29','RLIG1'), ('C12orf45',''), ('C12orf49','NOPCHAP1'), ('C12orf65','MTRFR'), ('C12orf66','KICS2'),
			('C12orf73','UQCC6'), ('C12orf74','PLEKHG7'), ('C12orf77','LINC02909'), ('C14orf144','LINC02691'), ('C14orf177','LINC02914'), ('C15orf41','CDIN1'),
			('C15orf53','LINC02694'), ('C15orf56','PAK6-AS1'), ('C15orf59','INSYN1'), ('C15orf65','PIERCE2'), ('C16orf45','BMERB1'), ('C16orf47','LOC108251797'),
			('C16orf58','RUSF1'), ('C16orf70','PHAF1'), ('C16orf71','DNAAF8'), ('C16orf91','UQCC4'), ('C16orf97','LINC02911'), ('C17orf102','TMEM132E-DT'),
			('C17orf112','LINC02876'), ('C17orf47','SEPTIN4'), ('C17orf51','ENSG00000266466'), ('C17orf53','HROB'), ('C17orf64','CHCT1'), ('C17orf82','LINC02875'),
			('C17orf97', 'LIAT1'), ('C17orf98','ENSG00000259015'), ('FAM69A','DIPK1A'),  ('FLJ46906', 'NHSL1-AS1'), ('TMEM56', 'ENSG00000237249'),
			('C3orf72', 'FOXL2NB'), ('C9orf64', 'QNG1'), ('C3orf72', 'FOXL2NB'), ('GALNTL4', 'GALNT18'), ('C3orf51', 'ERC2-IT1'), ('FAM208B', 'TASOR2'),
			('GALNTL4', 'GALNT18'), ('C3orf51', 'ERC2-IT1'), ('FAM208B', 'TASOR2'), ('SPG20OS', 'SPART-AS1'),('KIAA1324L', 'ELAPOR2'), ('ICK', 'CILK1'),
			('CCDC36', 'IHO1'), ('LOC100288911', 'PPARG'), ('IGJ','JCHAIN'),  ('SEPP1','SELENOP'), ('MCP', 'CD46'), ('CSF', 'CSF2'), ('GMCSF', 'CSF2'),
			('GCSF', 'CSF3'), ('LNP', 'LNPK'), ('K18', 'KRT18'), ('TGFB', 'TGIF1'), ('ARB', 'BEST1'), ('CTPA', 'EPHA2'), ('IgG3', 'IGHG3'), ('BNP', 'NPPB'),
			('MCP1', 'CC2'), ('COX', 'COX5A'), ('P450', 'POR'), ('NOS', 'NOS2') ]

		'''
			'IGHA2' is and A2M marker; however 'A2M' exists
		'''

		print("Start opening tables ....")
		# self.df_refseq, self.df_synonyms, self.dic_ref_seq = self.open_refseq(verbose=False)
		# self.df_omim = self.open_omim(verbose=False)
		# self.df_hgnc = self.open_hgnc(verbose=False)
		self.df_my_gene, self.df_my_gene_syn = self.open_my_gene(force=False, verbose=False)
		# self.dfunip = self.open_uniprot(verbose=False)
		# self.dfc = self.open_uniprot_conversion(verbose=False)

		self.dfunip = pd.DataFrame()
		self.dfc = pd.DataFrame()

		print("Building synonym dictionary ...")
		''' synonym dictionary - good one!!! 2024/04/27 '''
		dic_synonyms={}; multiple_synonym_list = []

		df = self.df_my_gene_syn.sort_values(['synonym', 'symbol'])
		for i in range(len(df)):
			row = df.iloc[i]
			try:
				symbol = dic_synonyms[row.synonym]
				# print(f"{row.synonym} already exists: {symbol} and new {row.symbol}")
				multiple_synonym_list.append(row.synonym)
			except:
				dic_synonyms[row.synonym] = df.iloc[i].symbol

		multiple_synonym_list = np.unique(multiple_synonym_list)
		# print(len(multiple_synonym_list))

		for gene, symbol in self.change_tuple:
			try:
				symb = dic_synonyms[gene]
			except:
				dic_synonyms[gene] = symbol
				# print("Adding", gene, symbol)

		self.dic_synonyms = dic_synonyms
		del(df)
		del(multiple_synonym_list)
		print("")


	def filter_my_gene(self):
		df_mg = self.df_my_gene.copy()

		df_mg = df_mg[ ~ df_mg.symbol.str.startswith('LOC')]
		df_mg = df_mg[ ~ df_mg.symbol.str.startswith('LINC')]
		df_mg = df_mg[ ~ df_mg.symbol.str.startswith('LNC')]
		df_mg = df_mg[ ~ df_mg.symbol.str.startswith('MIR')]
		df_mg = df_mg[ [False if '-' in x else True for x in df_mg.symbol] ]

		df_mg.reset_index(inplace=True, drop=True)

		return df_mg


	def save_change_table(self) -> bool:
		df = pd.DataFrame({'alias':self.dic_synonyms.keys(), 'symbol': self.dic_synonyms.values()})
		ret = pdwritecsv(df, self.fname_dic_change, self.root_refseq, verbose=True)
		return ret


	def open_my_gene(self, force:bool=False, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame]:

		df_my_gene = pdreadcsv(self.fname_my_gene, self.root_refseq, verbose=verbose)
		df_my_gene = df_my_gene[~df_my_gene.geneid.isna()]
		df_my_gene.geneid = df_my_gene.geneid.astype(int)
		if verbose: print(f"There are {len(df_my_gene)} symbols in mygene table (df_my_gene)")

		df_my_gene_syn = self.my_gene_to_synonym(df_my_gene, force=force, verbose=verbose)

		if df_my_gene_syn is None:
			self.dfmyg_index_symbol  = pd.DataFrame()
			self.dfmyg_index_synonym = pd.DataFrame()
		else:
			self.dfmyg_index_symbol  = df_my_gene_syn.set_index('symbol')
			self.dfmyg_index_synonym = df_my_gene_syn.set_index('synonym')

		return df_my_gene, df_my_gene_syn

	def open_biomart_protein(self, verbose:bool=False) -> pd.DataFrame:

		df_uni_prot = pdreadcsv(self.fname_biomart_protein, self.root_refseq, verbose=verbose)

		if df_uni_prot is not None and not df_uni_prot.empty:
			df_uni_prot.columns = ['entrezid', 'swissprot', 'uniprot_accession']
		else:
			df_uni_prot = pd.DataFrame()

		return df_uni_prot


	def open_biomart_entrez(self, verbose:bool=False) -> pd.DataFrame:

		df_etz = pdreadcsv(self.fname_biomart_entrez, self.root_refseq, verbose=verbose)

		if df_etz is not None and not df_etz.empty:
			df_etz.columns = ['entrezid', 'geneid', 'entrez_accession', 'entrez_desc']
		else:
			df_etz = pd.DataFrame()

		return df_etz


	def my_gene_to_synonym(self, df_my_gene:pd.DataFrame, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		filefull = os.path.join(self.root_refseq, self.fname_my_gene_synonym)

		if os.path.exists(filefull) and not force:
			dfd = pdreadcsv(self.fname_my_gene_synonym, self.root_refseq, verbose=verbose)
			return dfd

		icount=-1; dic={}
		for i in range(len(df_my_gene)):
			symbol = df_my_gene.iloc[i].symbol
			synonyms = df_my_gene.iloc[i].synonyms

			if isinstance(synonyms, float):
				continue

			if isinstance(synonyms, str):
				if '[' in synonyms:
					synonyms = eval(synonyms)
				else:
					synonyms = [synonyms]

			if not isinstance(synonyms, list):
				print(f"Error: {i} {symbol} synonyms {type(synonyms)} - {synonyms}")
				continue

			for synonym in synonyms:
				icount += 1
				dic[icount] = {}
				dic2 = dic[icount]

				dic2['symbol'] = symbol
				dic2['synonym'] = synonym

		if dic == {}:
			dfd = pd.DataFrame()
		else:
			dfd = pd.DataFrame(dic).T
			ret = pdwritecsv(dfd, self.fname_my_gene_synonym, self.root_refseq, verbose=verbose)

		return dfd

	def is_mygene_symbol(self, symb:str) -> bool:
		return symb in self.dfmyg_index_symbol.index.to_list()

	def find_mygene_symbol_to_synonym(self, symb:str, only_one:bool=True):
		dfa = self.df_my_gene_syn[self.df_my_gene_syn.symbol == symb]

		if dfa.empty:
			ret = None
		elif len(dfa) == 1:
			ret = dfa.iloc[0].synonym
		else:
			if only_one:
				ret = dfa.iloc[0].synonym
			else:
				ret = dfa.synonym.to_list()

		return ret

	def is_mygene_synonym(self, synonym:str) -> bool:
		return synonym in self.dfmyg_index_synonym.index.to_list()

	def find_mygene_synonym_to_symbol(self, synonym:str, only_one:bool=True):
		dfa = self.df_my_gene_syn[self.df_my_gene_syn.synonym == synonym]

		if dfa.empty:
			ret = None
		elif len(dfa) == 1:
			ret = dfa.iloc[0].symbol
		else:
			if only_one:
				ret = dfa.iloc[0].symbol
			else:
				ret = dfa.symbol.to_list()

		return ret


	def find_mygene_symbol(self, symb:str) -> object:
		dfa = self.df_my_gene[self.df_my_gene.symbol == symb]
		return None if dfa.empty else symb

	def find_mygene_geneid(self, symb:str) -> object:
		dfa = self.df_my_gene[self.df_my_gene.symbol == symb]
		return None if dfa.empty else dfa.iloc[0].geneid

	''' Ensembl Human GFF3 annotation - Ensembl
		https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
		Homo_sapiens.GRCh38.111.chr.gff3.gz  '''

	def open_gff(self, verbose:bool=False) -> pd.DataFrame:
		if self.df_gfff is not None and not self.df_gfff.empty:
			return self.df_gfff

		filefull = os.path.join(self.root_refseq, self.fname_gff_ensembl)

		if not os.path.exists(filefull):
			print('GFF table does not exists: %s'%(filefull))
			return pd.DataFrame()

		df_gfff = pdreadcsv(self.fname_gff_ensembl, self.root_refseq, verbose=verbose)

		self.df_gfff  = df_gfff

		return df_gfff


	'''----------------- RefSeq NCBI -------------------------------------------
		https://ftp.ncbi.nlm.nih.gov/refseq/
		https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/

		--- HGNC - did not finished webservice call (Flavio 16/08/2023) --------

		https://www.genenames.org/help/rest/
		https://rest.genenames.org/search/symbol:ZNF3+OR+ZNF12
	'''
	def open_refseq(self, verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
		
		df_refseq, df_synonyms, dic_ref_seq = pd.DataFrame(), pd.DataFrame(), {}
		self.df_refseq, self.df_synonyms, self.dic_ref_seq = df_refseq, df_synonyms, dic_ref_seq

		filefull = os.path.join(self.root_refseq, self.fname_refseq_ncbi)

		if not os.path.exists(filefull):
			print('RefSeq table does not exists: %s'%(filefull))
			return df_refseq, df_synonyms, dic_ref_seq

		df_refseq = pdreadcsv(self.fname_refseq_ncbi, self.root_refseq, verbose=verbose)
		if verbose: print(f"Found {len(df_refseq)} genes for refseq table")

		''' --------------- synonyms table ----------------------'''
		filefull = os.path.join(self.root_refseq, self.fname_ncbi_syn)

		if not os.path.exists(filefull):
			print(f"Synonyms table does not exists: '{filefull}'")
			return df_refseq, df_synonyms, dic_ref_seq

		df_synonyms = pdreadcsv(self.fname_ncbi_syn, self.root_refseq, verbose=verbose)
		if verbose: print(f"Found {len(df_refseq)} synonyms in refseq synonym table")

		filefull_dict = os.path.join(self.root_refseq, self.fname_df_ref_seq)

		if not os.path.exists(filefull_dict):
			dic_ref_seq, dfa = self.calc_dic_ref_seq(df_synonyms)
			_ = pdwritecsv(dfa, self.fname_df_ref_seq, self.root_refseq, verbose=verbose)
		else:
			dfa = pdreadcsv(self.fname_df_ref_seq, self.root_refseq, verbose=verbose)
			dic_ref_seq = {dfa.iloc[i].gene: list(dfa.iloc[i].iline) for i in range(len(dfa))}

		if verbose: print(f"dic_ref_seq has {len(dic_ref_seq)}")

		self.df_refseq   = df_refseq
		self.df_synonyms = df_synonyms
		self.dic_ref_seq = dic_ref_seq

		return df_refseq, df_synonyms, dic_ref_seq

	def calc_dic_ref_seq(self, df_synonyms) -> Tuple[dict, pd.DataFrame]:
		dic_ref_seq = {}

		for i in range(len(df_synonyms)):
			gene = df_synonyms.iloc[i].symbol
			if not isinstance(gene, str): continue

			if gene in dic_ref_seq.keys():
				dic_ref_seq[gene] += [i]
			else:
				dic_ref_seq[gene] = [i]

			mat = eval(df_synonyms.loc[i, 'synonyms'])
			for gene in mat:
				if not isinstance(gene, str): continue

				if gene in dic_ref_seq.keys():
					dic_ref_seq[gene] += [i]
				else:
					dic_ref_seq[gene] = [i]

		dfa = pd.DataFrame( {'gene': list(dic_ref_seq.keys()), 'iline': [np.unique(x) for x in dic_ref_seq.values()]} )

		return dic_ref_seq, dfa

	def open_hgnc(self, verbose:bool=False) -> pd.DataFrame:
		filefull = os.path.join(self.root_refseq, self.fname_hgnc_set)

		if not os.path.exists(filefull):
			df_hgnc = pd.DataFrame()
			print(f"HGNC table does not exist: '{filefull}'")
			print("Download it in https://www.genenames.org/download/archive/")
		else:
			df_hgnc = pdreadcsv(self.fname_hgnc_set, self.root_refseq, verbose=verbose)

		return df_hgnc

	def open_omim(self, verbose:bool=False) -> pd.DataFrame:
		filefull = os.path.join(self.root_refseq, self.fname_refseq_omim)

		if not os.path.exists(filefull):
			print('Table does not exists: %s'%(filefull))
			df_omim = pd.DataFrame()
		else:
			df_omim = pdreadcsv(self.fname_refseq_omim, self.root_refseq, verbose=verbose)

		return df_omim

	'''
	def parse_lista_to_mygene(self, symbols:List, N_symbs:int=500) -> bool:

		dfmyg3 = self.df_my_gene.copy()
		dfmyg3 = dfmyg3.set_index('symbol')

		symbols = [self.replace_synonym(x) for x in symbols]
		symbols = list(np.unique(symbols))

		symbols_not = [x for x in symbols if x not in dfmyg3.index]

		print(f"Was sended {len(symbols)} symbols, after replace_synonym(), and will be processed new {len(symbols_not)} new symbols.")

		if len(symbols_not) == 0:
			return False

		nLoop = int(np.ceil(len(symbols_not)/ N_symbs))

		mg = mygene.MyGeneInfo()
		symbols_not = np.array(symbols_not)

		for iloop in range(nLoop):
			print(f"{iloop}/{nLoop}", sep=' ')
			symbol_sel = symbols_not[(N_symbs*iloop): (N_symbs *(iloop+1)) ]

			lista = mg.querymany(symbol_sel, scopes='symbol',  species=9606, fields='all') # fields='geneid,name,entrezgene,refseq'

			dfq = self.mygene_lista_to_df(lista)
			# dfq = dfq[~dfq.geneid.isna()]
			dfq = dfq.drop_duplicates(['symbol', 'geneid'])
			print(f"...found {len(dfq)} symbos/geneid")

			dfq = pd.concat([self.df_my_gene, dfq])
			dfq.reset_index(inplace=True, drop=True)

			ret = pdwritecsv(dfq, self.fname_my_gene, self.root_refseq, verbose=True)
			if ret:
				self.df_my_gene = dfq
			else:
				break

		return ret
	'''


	def mygene_lista_to_df(self, lista: List) -> pd.DataFrame:
		i = -1
		dicg = {}

		for dic in lista:
			i += 1
			dicg[i] = {}
			dic2 = dicg[i]

			dic2['symbol'] = dic['query']

			try:
				geneid = dic['entrezgene']
			except:
				geneid = None

			dic2['geneid'] = geneid

			try:
				name = dic['name']
			except:
				name = None

			dic2['name'] = name


			try:
				synonyms = dic['alias']
			except:
				synonyms = None

			dic2['synonyms'] = synonyms


			try:
				other_names = dic['other_names']
			except:
				other_names = None

			dic2['other_names'] = other_names


			try:
				ec_enzyme = dic['ec']
			except:
				ec_enzyme = None

			dic2['ec_enzyme'] = ec_enzyme


			try:
				dic_refseq = dic['refseq']
			except:
				dic_refseq = None

			# dic2['dic_refseq'] = dic_refseq

			if dic_refseq is None:
				dic2['refseq_gen'], dic2['refseq_prot'], dic2['refseq_rna'], dic2['refseq_trans'] = None, None, None, None
			else:
				try:
					dic2['refseq_gen'] = dic_refseq['genomic']
				except:
					dic2['refseq_gen'] = None

				try:
					dic2['refseq_prot'] = dic_refseq['protein']
				except:
					dic2['refseq_prot'] = None


				try:
					dic2['refseq_rna'] = dic_refseq['rna']
				except:
					dic2['refseq_rna'] = None

				try:
					dic2['refseq_trans'] = dic_refseq['translation']
				except:
					dic2['refseq_trans'] = None

			try:
				hgnc = dic['HGNC']
			except:
				hgnc = None

			dic2['hgnc'] = hgnc

			try:
				omim = dic['MIM']
			except:
				omim = None

			dic2['omim'] = omim

			try:
				dic_acessions = dic['accession']
			except:
				dic_acessions = None

			if dic_acessions is None:
				acessions_gen, acessions_prot, acessions_rna, acessions_trans = None, None, None, None
			else:
				try:
					acessions_gen = dic_acessions['genomic']
				except:
					acessions_gen = None

				try:
					acessions_prot= dic_acessions['protein']
				except:
					acessions_prot = None

				try:
					acessions_rna= dic_acessions['rna']
				except:
					acessions_rna = None

				try:
					acessions_trans = dic_acessions['translation']
				except:
					acessions_trans = None


			dic2['acessions_gen'] = acessions_gen
			dic2['acessions_prot'] = acessions_prot
			dic2['acessions_rna'] = acessions_rna
			dic2['acessions_trans'] = acessions_trans

			try:
				synonyms = dic['alias']
			except:
				synonyms = None

			dic2['synonyms'] = synonyms

			try:
				map_location = dic['map_location']
			except:
				map_location = None

			dic2['map_location'] = map_location

			try:
				dic_panther = dic['pantherdb']
			except:
				dic_panther = None

			if dic_panther is None:
				dic2['dic_panther'], dic2['ortholog'] = None, None
			else:
				dic2['dic_panther'] = dic_panther
				dic2['ortholog'] = dic_panther['ortholog']

			try:
				dic_uniprot = dic['uniprot']
			except:
				dic_uniprot = None

			if dic_uniprot is None:
				dic2['dic_uniprot'], dic2['swissprot'] = None, None
			else:
				dic2['dic_uniprot'] = dic_uniprot

				try:
					dic2['swissprot'] = dic_uniprot['Swiss-Prot']
				except:
					dic2['swissprot'] = None

			try:
				wikipedia = dic['wikipedia']
			except:
				wikipedia = None

			dic2['wikipedia'] = wikipedia

			# print(i, symbol, geneid, name, refseq_gen, acessions, synonyms, ensemblid)#

		df_my_gene = pd.DataFrame(dicg).T

		return df_my_gene


	def open_uniprot(self, verbose:bool=False) -> pd.DataFrame:
		dfunip = pdreadcsv(self.fname_uniprot, self.root_project, verbose=verbose)
		self.dfunip = dfunip
		return dfunip

	def open_uniprot_conversion(self, verbose:bool=False) -> pd.DataFrame:
		dfunip = pdreadcsv(self.fname_conv_uniprot, self.root_refseq, verbose=verbose)
		self.dfunip = dfunip
		return dfunip

	def uniprot_add_symbol(self, df, fname) -> pd.DataFrame:
		fname_backup = fname.replace('.tsv', '_backup.tsv')
		ret = pdwritecsv(df, fname_backup, self.root_project)
		if not ret:
			return pd.DataFrame()

		if self.dfunip is None or self.dfunip.empty:
			return pd.DataFrame()

		cols = ['uniID', 'symbol', 'uniprot_id', 'description']
		dfn = pd.merge(df, self.dfunip[cols], how='outer', on='uniID')
		dfn = dfn[~pd.isnull(dfn.lfc)].copy()
		dfn = dfn.sort_values('lfc', ascending=False)
		dfn.reset_index(inplace=True, drop=True)

		ret = pdwritecsv(dfn, fname, self.root_project)
		return dfn


	def replace_symbol_to_synonym(self, symbol_syn) -> object:
		if not isinstance(symbol_syn, str) or len(symbol_syn) == 0:
			return None

		try:
			ret = self.dic_synonyms[symbol_syn]
		except:
			ret = symbol_syn

		return ret

	'''
	HUGO - HGNC - https://www.genenames.org/download/statistics-and-files/
	# gene = 'Yes'
	# gene = 'ACE2'
	'''

	'''
	def find_hugo_symbol(self, gene):
		try:
			i = self.dic_my_genes[gene]
			if isinstance(i, list):
				i = i[0]

			mat = self.df_synonyms.iloc[i]['synonyms']
			# print(">>>", gene, i, mat, type(mat))
			if isinstance(mat, str):
				mat = eval(mat)

			gene0 = mat[0]
		except:
			gene0 = gene

		return gene0
	'''

	def merge_refseq_synonyms(self, no_symbols:List=[]) -> bool:
		print('Searching in synonyms in refseq')

		if not isinstance(no_symbols, list):
			return False
		
		if no_symbols == []:
			return False
		
		df_list = []
		no_symbols = np.unique(no_symbols)

		for symbol in no_symbols:
			# print('>>', symbol)
			dfnew = self.df_synonyms[self.df_synonyms.symbol == symbol]
			if dfnew.empty: continue

			dfnew = pd.DataFrame(dfnew.iloc[0]).T.copy()
			df_list.append(dfnew)

		if len(df_list) > 1:
			dfnew = pd.concat(df_list)
			dfref = pd.concat([self.df_refseq, dfnew])

			ret = pdwritecsv(dfref, self.fname_refseq_ncbi, self.root_refseq, verbose=True)
			print("Merging %d new genes into refseq table"%(len(dfnew)))

			if ret: 
				self.df_refseq = dfref
		else:
			print("Nothing to merge")

		return True


	def open_refseq_ncbi_info(self):

		filefull = os.path.join(self.root_refseq, self.fname_refseq_gene_info)

		if not os.path.exists(filefull):
			print('Table does not exists: %s'%(filefull))
			return None

		df = pdreadcsv(self.fname_refseq_gene_info, self.root_refseq, verbose=True)
		cols = ['tax_id', 'geneid', 'symbol', 'locus_tag', 'synonyms', 'dbxrefs',
			   'chromosome', 'map_location', 'description', 'type_of_gene',
			   'symbol_from_nomenclature_authority', 'full_name_from_nomenclature_authority',
			   'nomenclature_status', 'other_designations', 'modification_date', 'feature_type']
		df.columns = cols
		selcols = ['tax_id', 'geneid', 'symbol', 'locus_tag', 'synonyms',  'chromosome', 'map_location', 'description', 'type_of_gene']
		df = df[selcols]
		df.synonyms = [x.split('|') for x in df.synonyms]

		print("Found %d genes for refseq gene info"%(len(df)))

		return df

	def search_symbol_and_synonyms(self, symbol:str) -> List:

		if symbol in self.dic_gene_term_list.keys():
			symbol_new = self.dic_gene_term_list[symbol]
			symbol_list = [symbol_new]
		else:
			symbol_list = []


		if symbol in self.dic_gene_term_list_rom.keys():
			symbol_new = self.dic_gene_term_list_rom[symbol]
			symbol_list += [symbol_new]

		synonyms = self.find_mygene_symbol_to_synonym(symbol, only_one=False)

		if synonyms is not None:
			if isinstance(synonyms, str):
				symbol_list += [synonyms]
			else:
				symbol_list += synonyms

		return symbol_list

	def prepare_full_gff(self, force:bool=False, verbose:bool=True) -> bool:

		fullname = os.path.join(self.root_refseq, self.fname_gff_ensembl0)
		fullname_tsv = os.path.join(self.root_refseq, self.fname_gff_ensembl0_tsv)

		if not os.path.exists(fullname):
			print(f"Error: gff must be downloaded '{fullname}'")
			return False

		if os.path.exists(fullname_tsv) and not force:
			print(f"Already done: file '{fullname_tsv}'")
			return True

		fh, fhw = None, None
		try:
			fh = open(fullname, mode="r")
			fhw = open(fullname_tsv, mode="w")

			mat = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
			line = "\t".join(mat) + '\n'
			fhw.write(line)

			i = -1
			icount = -1
			while True:
				i += 1
				line = fh.readline()
				if not line: break
				if line.startswith("#"): continue

				icount += 1
				if icount%1000 == 0:
					print(i, len(mat), line[:-1])

				fhw.write(line)

				# if icount == 1000:break

		except:
			stri = f"Could not open {fullname} / {fullname_tsv}"
			print(stri)
			# log_save(stri, filename=self.filelog, verbose=True)
			return False
		finally:
			if fh:
				fh.close()
			if fhw:
				fhw.close()

		if verbose:
			print("------------- end ------------")
			print(f"saved {icount} lines")

		return True

	def open_gff_full(self, verbose:bool=False) -> pd.DataFrame:
		if self.df_gff_full is not None and not self.df_gff_full.empty:
			return self.df_gff_full

		fullname = os.path.join(self.root_refseq, self.fname_gff_ensembl0_tsv)
		if not os.path.exists(fullname):
			print(f"Could not find {fullname}")
			self.df_gff_full = pd.DataFrame()
			return self.df_gff_full

		self.df_gff_full = pdreadcsv(self.fname_gff_ensembl0_tsv, self.root_refseq)
		return self.df_gff_full

	def group_full_gff(self, verbose:bool=False) -> pd.DataFrame:

		df_gff_full = self.open_gff_full(verbose=verbose)
		if df_gff_full is None:
			return pd.DataFrame()

		dfcount = df_gff_full.groupby("feature").count().reset_index().iloc[:,:2]
		return dfcount

	def split_gff_attribute(self, x:str):

		mat = x[3:].split(';')

		'''
		gene, gene_id, transcript, transcript_id, name, logic_name, biotype, description, version, \
		tag, tag_basic, parent, transcript_support_level = \
		'gene', 'gene_id', 'transcript', 'transcript_id', 'name', 'logic_name', 'biotype', 'description', 'version', \
		'tag', 'tag basic', 'parent_gene', 'transcript_support_level'
		'''

		gene_val, gene_id_val, transcript_val, transcript_id_val, name_val, logic_name_val, biotype_val, \
		description_val, version_val, tag_val, tag_basic_val, parent_val, transcript_support_level_val = \
		None, None, None, None, None, None, None, None, None, None, None, None, None

		for stri in mat:
			try:
				if stri.startswith('description'):
					_type = 'description'
					val = stri[ (len('description')+1): ]
				elif stri.startswith('transcript_id'):
					_type = 'transcript_id'
					val = stri[ (len('transcript_id')+1): ]
				elif stri.startswith('tag=basic'):
					_type = 'tag=basic'
					val = stri[ (len('tag=basic')+1): ]
				elif ":" in stri:
					_type, val = stri.split(':')
				elif "=" in stri:
					_type, val = stri.split('=')
				elif "," in stri:
					_type, val = stri.split(',')
				else:
					print("????", stri)
					continue
			except:
				print("????", stri)
				continue

			try:
				if _type == 'gene':
					gene_val = val
				elif _type == 'gene_id':
					gene_id_val = val
				elif _type == 'transcript':
					transcript_val = val
				elif _type == 'transcript_id':
					transcript_id_val = val
				elif _type == 'Name':
					name_val = val
				elif _type == 'logic_name':
					logic_name_val = val
				elif _type == 'biotype':
					biotype_val = val
				elif _type == 'description':
					description_val = val
				elif _type == 'version':
					version_val = val
				elif _type == 'Parent=gene':
					parent_val = val
				elif _type == 'tag':
					tag_val = val
				elif _type == 'tag=basic':
					tag_basic_val = val
				elif _type == 'tag=basic':
					tag_basic_val = val
				elif _type == 'transcript_support_level':
					transcript_support_level_val = val

				else:
					print("Include!", stri, "->", _type, val)

			except:
				print("Error:", x)

		# print(gene, gene_val, transcript, transcript_val, name, name_val, biotype, biotype_val, description, description_val )
		# raise Exception('stop')
		return gene_val, gene_id_val, transcript_val, transcript_id_val, name_val, logic_name_val, biotype_val, \
			   description_val, version_val, tag_val, tag_basic_val, parent_val, transcript_support_level_val


	'''
	def split_gff_attribute_new_old(self, x:str, old_version:bool=False):

		x = x.strip()
		if old_version:
			mat = x.split(';')
		else:
			mat = x[3:].split('; ')

		gene_val, gene_id_val, version_val, source_val, transcript_val, transcript_id_val, \
		name_val, logic_name_val, biotype_val, description_val, exon_number_val, exon_id_val, \
		tag, tag_basic_val, parent_val, transcript_support_level_val, protein_id_val, ccds_id_val = [None]*18

		for stri in mat:
			stri = stri.strip()
			if stri == '': continue

			if old_version:
				# gene_name "DDX11L1" --> gene_name:DDX11L1 
				if ' "' in stri:
					stri = stri.replace(' "', ':')
					stri = stri[:-1]

			mat=None
			try:
				if stri.startswith('description'):
					_type = 'description'
					val = stri[ (len('description')+1): ]
		
				elif stri.startswith('transcript_id'):
					_type = 'transcript_id'
					val = stri[ (len('transcript_id')+1): ]
				elif stri.startswith('tag=basic'):
					_type = 'tag=basic'
					val = stri[ (len('tag=basic')+1): ]
				elif ":" in stri:
					mat = stri.split(':')
					char = ':'
				elif "=" in stri:
					mat = stri.split(':')
					char = '='
				elif "," in stri:
					mat = stri.split(':')
					char = ','
				else:
					print("Error1:", stri)
					raise Exception('stop')

				if mat is not None:
					if len(mat) == 2:
						_type, val = mat
					else:
						_type = mat[0]
						val   = char.join(mat)
						
			except:
				print("Error2:", stri)
				raise Exception('stop')

			_type = _type.strip()
			val  = val.strip()
		
			# print("???", _type, val)

			if _type == 'gene':
				gene_val = val
			elif _type == 'gene_id':
				gene_id_val = val
			elif _type == 'exon_id':
				exon_id_val = val
			elif _type == 'protein_id':
				protein_id_val = val
			elif _type == 'ccds_id':
				ccds_id_val = val
			elif 'version' in _type:
				version_val = val
			elif 'source' in _type:
				source_val = val
			elif _type == 'transcript':
				transcript_val = val
			elif _type == 'transcript_id':
				transcript_id_val = val
			elif _type == 'Name' or _type == 'gene_name' or _type=='transcript_name':
				name_val = val
			elif _type == 'logic_name':
				logic_name = val
			elif 'biotype' in _type:
				biotype_val = val
			elif _type == 'description':
				description_val = val
			elif _type == 'exon_number':
				exon_number_val = val
			elif _type == 'Parent=gene':
				parent_val = val
			elif _type == 'tag':
				tag_val = val
			elif _type == 'tag=basic':
				tag_basic_val = val
			elif _type == 'tag=basic':
				tag_basic_val = val
			elif _type == 'transcript_support_level':
				transcript_support_level_val = val
			else:
				print("Include!", stri, "->", _type, val)


		# print(gene, gene_val, transcript, transcript_val, name, name_val, biotype, biotype_val, description, description_val )
		# raise Exception('stop')
		return gene_val, gene_id_val, version_val, source_val, transcript_val, transcript_id_val, \
			   name_val, logic_name_val, biotype_val, description_val, exon_id_val, exon_number_val, \
			   tag, tag_basic_val, parent_val, transcript_support_level_val, protein_id_val, ccds_id_val
	'''


	def prepare_final_gff(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		fullname = os.path.join(self.root_refseq, self.fname_gff_ensembl)
		if os.path.exists(fullname) and not force:
			df_gfff = pdreadcsv(self.fname_gff_ensembl, self.root_refseq)
			self.df_gfff = df_gfff
			return df_gfff

		df_gff_full = self.open_gff_full(verbose)

		if df_gff_full is None:
			return pd.DataFrame()

		feature_list = ['gene', 'lnc_RNA', 'miRNA', 'ncRNA_gene',
						'pseudogene', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'tRNA', ]

		df_gffp = df_gff_full[df_gff_full.feature.isin(feature_list)]
		self.df_gffp = df_gffp

		fname = self.fname_gff_selected
		ret = pdwritecsv(df_gffp, fname, self.root_refseq)

		if not ret:
			return pd.DataFrame()

		df_gfff = df_gffp.copy()
		cols_ret =  ['gene', 'gene_id', 'transcript', 'transcript_id', 'name', 'logic_name', 'biotype', \
					 'description', 'version', 'tag', 'tag_basic', 'parent', 'transcript_support_level']

		df_gfff.loc[:, cols_ret] =  [self.split_gff_attribute(x) for x in df_gfff.attribute]

		cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
				'frame', 'name', 'gene_id', 'biotype', 'description',
				'version',  'tag_basic', 'parent',
				'transcript', 'transcript_id','transcript_support_level'] # 'attribute',

		df_gfff = df_gfff[cols]
		cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
			   'frame', 'symbol', 'ensembl_id', 'biotype', 'description',
			   'version',  'tag_basic', 'parent',
			   'transcript', 'transcript_id','transcript_support_level']

		df_gfff.columns = cols
		self.df_gfff = df_gfff

		ret = pdwritecsv(df_gfff, self.fname_gff_ensembl, self.root_refseq)
		return df_gfff


	def prepare_any_gff_gtd(self, fname_gft90:str, root_gff:str, skiprows:int=5,
							header:int=0, force:bool=False, verbose:bool=False) -> pd.DataFrame:

		fname_tsv = fname_gft90 + '.tsv'

		fullname	 = os.path.join(root_gff, fname_gft90)
		fullname_tsv = os.path.join(root_gff, fname_tsv)


		if os.path.exists(fullname_tsv) and not force:
			df = pdreadcsv(fullname_tsv, root_gff, verbose=verbose)
			return df

		if not os.path.exists(fullname):
			print(f"Error: gff/gft must be downloaded '{fullname}'")
			return pd.DataFrame()

		df90 = pdreadcsv(fname_gft90, root_gff, skiprows=skiprows, header=header)

		cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
		dfa = df90.copy()
		dfa.columns = cols

		cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
				'frame', 'name', 'gene_id', 'biotype', 'description',
				'version',  'tag_basic', 'parent',
				'transcript', 'transcript_id','transcript_support_level'] # 'attribute',

		dfa = dfa[cols]
		cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
				'frame', 'symbol', 'ensembl_id', 'biotype', 'description',
				'version',  'tag_basic', 'parent',
				'transcript', 'transcript_id','transcript_support_level']

		dfa.columns = cols

		_ = pdwritecsv(dfa, fname_tsv, root_gff, verbose=True)
		return dfa




