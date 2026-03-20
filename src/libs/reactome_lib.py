#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2024/04/11
# Udated  on 2024/04/11
# @author: Flavio Lichtenstein
# @local: Bioinformatics: CENTD/Molecular Biology; Instituto Butatan

import numpy as np
import os, json # , gzip, pickle, sys
import pandas as pd
# from collections import defaultdict, Counter, OrderedDict
# from datetime import datetime

from typing import List, Tuple #Optional, Iterable, Set, Any,

from libs.Basic import *

class Reactome(object):
    def __init__(self, root_data:str='../data/reactome/'):

        self.root_data = root_data
        self.root_json = os.path.join(root_data, 'pathway_json')

        self.fname_reactome          = 'reactome_pathways_human.tsv'
        self.fname_reactome_2024     = 'Reactome_Pathways_2024.tsv'
        self.fname_reactome_abstract = 'reactome_pathways_abstract.tsv'
        self.fname_reactome_gmt      = 'reactome_pathways_human_gmt.tsv'

        self.dfr, self.df_gmt, self.dfr_merge = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        self.pathway_list = []
        self.n_pathways = 0

        self.script_curl = "curl -X GET --header 'Accept: application/json' 'https://reactome.org/ContentService/data/query/%s' --silent > %s"

    def refresh_reactome_table(self, pathw_list:List, force:bool=False, verbose:bool=True) -> bool:

        dfr = self.open_reactome_abstract(verbose=False)
        self.dfr = dfr

        if dfr is None:
            redo = True
        else:
            pathw_list = [x for x in pathw_list if x not in dfr.pathway_id]
            redo = len(pathw_list) > 0

        if redo:
            self.get_many_reactome_content_service_pathway(pathw_list, verbose=verbose)
            ret = self.build_reactome_table(dfr, force=force, verbose=verbose)

            if verbose: print("reactome table was updated.")
        else:
            if verbose: print("reactome table is Ok.")
            ret = False

        # return  if has reviewed
        return ret

    def get_many_reactome_content_service_pathway(self, pathw_list:List, force:bool=False, verbose:bool=False):

        for pathway_id in pathw_list:
            # if not verbose: print(".", end='')
            self.get_reactome_content_service_pathway(pathway_id, force=force, verbose=verbose)


    def get_reactome_content_service_pathway(self, pathway_id:str, force:bool=False, verbose:bool=False) -> bool:
        pathway_id_json = f'{pathway_id}.json'
        fullname = os.path.join(self.root_json, pathway_id_json)

        if os.path.exists(fullname) and not force:
            if verbose: print(f"File already exists: {fullname}")
            return True

        command = self.script_curl%(pathway_id, fullname)
        ret = os.system(command)
        if verbose:
            if ret == 0:
                 print(f"File downloaded {ret}: {pathway_id} -> {fullname}")
            else:
                 print(f"Error: could not download {ret}: {pathway_id} -> {fullname}")

        if os.path.exists(fullname) and os.path.getsize(fullname) == 0:
            os.remove(fullname)
            print(f"Error: size zero json file: {pathway_id} -> {fullname}")

        return ret==0
    

    def open_reactome(self, verbose:bool=False) -> pd.DataFrame:

        fullname = os.path.join(self.root_data, self.fname_reactome)

        if not os.path.exists(fullname):
            self.pathway_list = []
            self.dfr = pd.DataFrame()
            self.n_pathways = 0

            print(f"Could not find: {fullname}")
            return pd.DataFrame()

        '''
            bpx.dfr['pathway'] = [x.strip() for x in bpx.dfr.pathway]
            pdwritecsv(bpx.dfr, bpx.reactome.fname_reactome, bpx.reactome.root_data, verbose=True)
        '''
        dfr = pdreadcsv(self.fname_reactome, self.root_data, verbose=verbose)

        self.dfr = dfr
        self.pathway_list = list(np.unique(dfr.pathway))
        self.n_pathways = len(self.pathway_list)

        return dfr

    def open_reactome_abstract(self, verbose:bool=False) -> pd.DataFrame:

        fullname = os.path.join(self.root_data, self.fname_reactome_abstract)

        if not os.path.exists(fullname):
            print(f"Could not find: {fullname}")
            return None

        df = pdreadcsv(self.fname_reactome_abstract, self.root_data, verbose=verbose)
        return df



    def open_reactome_gmt(self, verbose:bool=False) -> pd.DataFrame:
        fullname = os.path.join(self.root_data, self.fname_reactome_gmt)

        if not os.path.exists(fullname):
            print(f"There is no reactome gst file: '{fullname}'")
            return pd.DataFrame()

        df_gmt = pdreadcsv(self.fname_reactome_gmt, self.root_data, verbose=verbose)
        self.df_gmt = df_gmt

        return df_gmt

    def merge_reactome_table_gmt(self, verbose:bool=False)-> bool:
        self.dfr_merge = pd.DataFrame()

        dft = self.open_reactome_abstract(verbose=verbose)
        dfg = self.open_reactome_gmt(verbose=verbose)

        if dft is None or dft.empty:
            return False
        if dfg is None or dfg.empty:
            return False

        cols = ['pathway_id', 'genes', 'n']
        dfg = dfg[cols].copy()
        dfg.columns =  ['pathway_id', 'genes_pathway', 'ngenes_pathway']

        self.dfr_merge = pd.merge(dft, dfg, how="inner", on='pathway_id')

        return True

    def build_reactome_table(self, dfr:pd.DataFrame, force:bool=False, verbose:bool=False) -> bool:

        if force:
            dfr = None
            files = [x for x in os.listdir(self.root_json) if x.endswith('.json')]

            if len(files) == 0:
                print(f"There are no json files in '{self.root_json}'")
                return False

        else:
            pathway_id_list = dfr.pathway_id.to_list()
            files = [x for x in os.listdir(self.root_json) if x.endswith('.json') if x+'.json' not in pathway_id_list]

            if len(files) == 0:
                if verbose: print("There are no json files to update")
                return False

        if verbose: print(f"Updating {len(files)} json files.")

        dicf = {}; i = -1

        for fname in files:
            fullname = os.path.join(self.root_json, fname)

            if not os.path.exists(fullname):
                print(f"Could not find: {fullname}")
                continue

            try:
                f = open(fullname)
                dic = json.load(f)
            except:
                dic = None
                print(f"Could not open: {fullname}")
                continue

            try:
                ret_error = dic["code"]
                os.remove(fullname)
            except:
                ret_error = None
                pass

            if ret_error is not None:
                try:
                    ret_msg = dic["messages"]
                except:
                    ret_msg = '???'
                print(f"Error in json download: {fullname} - ret_error={ret_error} {ret_msg}")

                continue

            i += 1
            dicf[i] = {}
            dic2 = dicf[i]

            for key, val in dic.items():
                dic2[key] = val

                if key == 'summation':
                    all_text = ''
                    summaries = val
                    if len(summaries) > 0:
                        for dic3 in summaries:
                            if all_text != '':
                                all_text += '\n\n'
                            all_text += dic3['text']

                    dic2['abstract'] = all_text

        if len(dicf) == 0:
            print("Nothing updated.")
            return False

        df = pd.DataFrame(dicf).T
        cols = list(df.columns)

        # manually add R-HSA-71406
        special_pathwya_id = 'R-HSA-71406'

        if dfr is None or dfr.empty:
            import_special = True
        else:
            import_special = False if special_pathwya_id in dfr.pathway_id else True

        if import_special:
            i += 1
            dicf[i] = {}
            dic2 = dicf[i]

            dic2[cols[1]] = 'Pyruvate metabolism and Citric Acid (TCA) cycle'
            dic2[cols[2]] = special_pathwya_id
            dic2['abstract'] = 'In the citric acid or tricarboxylic acid (TCA) cycle, the acetyl group of acetyl CoA (derived primarily from oxidative decarboxylation of pyruvate, beta-oxidation of long-chain fatty acids, and catabolism of ketone bodies and several amino acids) can be completely oxidized to CO2 in reactions that also yield one high-energy phosphate bond (as GTP or ATP) and four reducing equivalents (three NADH + H+, and one FADH2). Then, the electron transport chain oxidizes NADH and FADH2 to yield nine more high-energy phosphate bonds (as ATP). All reactions of the citric acid cycle take place in the mitochondrion.'

            # again
            df = pd.DataFrame(dicf).T

        cols[1] = 'pathway'
        cols[2] = 'pathway_id'
        df.columns = cols

        if dfr is not None and not dfr.empty:
            if verbose: print(f"There are {len(dfr)} dfr regs, incrementing {len(df)}.")

            cols = list(dfr.columns)
            # df with the same columns
            df = df[cols]

            df = pd.concat( [dfr, df] )

        goods = [True if isfloat(df.iloc[i]['dbId']) or pd.isnull(df.iloc[i]['dbId']) else False for i in range(len(df))]
        df = df[goods]

        df = df.drop_duplicates("pathway_id")
        df.index = np.arange(0, len(df))

        pdwritecsv(df, self.fname_reactome_abstract, self.root_data, verbose=verbose)
        self.dfr_merge = df
        return True
