#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2024/06/10
# Udated  on 2024/06/11
# @self.author: Flavio Lichtenstein

import os
from pathlib import Path
import numpy as np
from typing import Optional, Iterable, Set, Tuple, Any, List

# !pip install biomart
from biomart import BiomartServer

from libs.Basic import create_dir, pdwritecsv, pdreadcsv

class Biomart(object):
    def __init__(self, root_colab:Path=Path('../colab/'), dataset='hsapiens_gene_ensembl'):

        self.biomart_web_service_url = "https://www.ensembl.org/biomart/martservice"
        self.biomart_url = "http://www.ensembl.org/biomart"
        self.dataset = 'hsapiens_gene_ensembl'

        self.attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 
                           'hgnc_symbol', 'description', 'gene_biotype', 
                           'chromosome_name', 'start_position', 'end_position', 'strand', 'go_id', 'version']

        self.root_biomart = create_dir(root_colab, 'biomart')

    def download_biomart_hsapiens(self, attributes):
        lista = []
        
        lista.append("wget -O result.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query>")
        lista.append("<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >")
        lista.append(f"<Dataset name = \"{self.dataset}\" interface = \"default\" >")
        
        for attribute in attributes:
            lista.append(f"<Attribute name=\"{attribute}\"/>")
        
        lista.append("</Dataset></Query>'")

        wget_command = "".join(lista)

        os.system("".join(lista))
        print("exports to ./biomart.txt")


    def open_biomart_hsapiens(self):
        fname = 'biomart_hsapiens_no_symbol_repeated.tsv'
        dfbm = pdreadcsv(fname, self.root_biomart)

        self.dfbm = dfbm

        return dfbm

    def prepare_downloaded_biomart_table(self):
        '''
        problems COVID-19 proteomics:
            IGL --> IGLC1  
            KIAA2026 --> BRD10  
            LILRA3 OK  
        '''

        fname = 'biomart_hsapiens_complete.tsv'
        dfbm = pdreadcsv(fname, self.root_biomart, header=None)
        dfbm.columns = ['ensembl_transcript_id',  'ensembl_gene_id', 'symbol',  'description', 'gene_biotype', 
                        'chromosome', 'start_position',  'end_position', 'strand']

        pdwritecsv(dfbm, fname, self.root_biomart)

        fname = 'biomart_hsapiens_no_symbol_repeated.tsv'
        dfbm1 = dfbm.sort_values(['symbol', 'ensembl_gene_id', 'ensembl_transcript_id'])
        dfbm1.index = np.arange(0, len(dfbm1))

        dfbm1 = dfbm1.drop_duplicates('symbol').copy()
        dfbm1.index = np.arange(0, len(dfbm1))

        pdwritecsv(dfbm1, fname, self.root_biomart)
        return dfbm1


    def open_biomart(self, verbose:bool=True):
        server = BiomartServer(self.biomart_url)
        server.verbose = verbose
        mart = server.datasets[self.dataset]

        self.server = server
        self.mart = mart


 