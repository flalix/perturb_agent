#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

import numpy as np
import pandas as pd
import requests
# from fileinput import filename
import json
import os
import re
import time
import math
import warnings
from collections import Counter
from pathlib import Path
from tabnanny import verbose
from typing import Any, Iterable, List, Optional, Tuple


from scipy.stats import hypergeom, ttest_ind, zscore
from sklearn.cluster import KMeans
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

import umap
import hdbscan

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns


from libs.Basic import create_dir, pdreadcsv, pdwritecsv, read_txt, title_replace, write_txt
from libs.calc_degs_lib import CALC_DEGS
from libs.stat_lib import fdr
from project_context_GDC import load_project_context


class GDC(object):
    def __init__(self, root0: Path, root0_data: Path, memory_restriction: bool=True):

        self.memory_restriction = memory_restriction

        self.url_gdc_project = "https://api.gdc.cancer.gov/projects"
        self.url_gdc_cases = "https://api.gdc.cancer.gov/cases"
        self.url_gdc_files = "https://api.gdc.cancer.gov/files"
        self.url_gdc_data = "https://api.gdc.cancer.gov/data"

        self.url_cbioportal = "https://www.cbioportal.org/api"

        self.prog_id, self.psi_id = "", ""

        self.root0 = Path(root0)
        self.root_src   =  create_dir(self.root0, 'src')

        self.root0_data = Path(root0_data)
        self.root_colab = create_dir(self.root0_data, 'colab')
        self.root_gtex  = create_dir(self.root_colab, "GTEx")

        self.fname_gtex_table = "tcga_primary_site_to_gtex_ids.tsv"
        self.df_gtex_to_tcga = pd.DataFrame()
        self.gtex_id = ""
        self.fname_gtex_exp_counts = "gtex_expression_counts_%s.tsv"

        self.fname_tpm_exp = "gtex_TPM_%s.tsv"
        self.fname_GTEx_counts = "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.gct.gz"
        self.fname_GTEx_meta = "GTEx_Analysis_v11_Annotations_SampleAttributesDS.tsv"
        self.fname_GTEx_pheno = "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.tsv"

        self.GTEX_API = "https://gtexportal.org/api/v2"

        self.df_gtex_counts = pd.DataFrame()
        self.df_gtex_pheno = pd.DataFrame()
        self.df_meta = pd.DataFrame()
  
        self.fname_exp_tumor = 'expression_tumor_for_%s.tsv'
        self.fname_exp_normal = 'expression_normal_for_%s.tsv'
        self.fname_exp_gtex = 'expression_gtex_for_%s.tsv'

        # self.get_gtex_control(Nsamples=15, force=False, verbose=verbose)
        # GTEx control - per tissue
        self.df_gtex_ctrl = pd.DataFrame()
        self.df_meta_prep = pd.DataFrame()

        self.df_gtex_normal = pd.DataFrame()
        self.df_normal = pd.DataFrame()
        self.df_tumor = pd.DataFrame()

        self.root_project = Path()
        self.root_summary = Path()
        self.root_disease = Path()
        self.root_samples = Path()
        self.root_lfc = Path()
        self.root_mutations = Path()

        self.GENE_COLS = ["geneid", "symbol", "biotype"]

        self.clean_gdc_files()

        self.fname_all_cases = "%s_summ_cases.tsv"
        self.fname_all_samples = "%s_summ_samples.tsv"
        self.fname_all_mutations = "%s_summ_mutations.tsv"

        self.fname_lfc = "lfc_%s.tsv"
        self.fname_degs = "degs_%s.tsv"
        self.fname_degs_txt = "degs_%s.txt"
        self.fname_sample_txt = "samples_msg_%s.txt"

        
        ctx = load_project_context( )

        self.SUBTYPE_GENES = ctx.SUBTYPE_GENES
        self.HISTOLOGY_GENES = ctx.HISTOLOGY_GENES
        self.TUMOR_CLASS = ctx.TUMOR_CLASS
        self.GLOBAL_SUBTYPE = ctx.GLOBAL_SUBTYPE
        self.HISTOLOGY = ctx.HISTOLOGY
        self.SITE_MAP = ctx.SITE_MAP
        self.colors = ctx.colors

    def clean_gdc_files(self):

        self.gdc_file_name = ""
        self.gdc_file_id = ""
        self.gdc_data_type = ""

        self.s_case = ""

        self.fname_programs = "gdc_programs.txt"

        # primary_site
        self.fname_prim_site_tcga = "primary_site_program_TCGA.tsv"
        self.fname_prim_site_cbio = "gdc_to_cbioportal_study_mapping.tsv"
        self.fname_cases0 = "cases_for_%s.tsv"
        self.fname_subtype0 = "subtype_for_%s.tsv"
        self.fname_samples0 = "samples_for_%s.tsv"
        self.fname_vcf_files0 = "vcf_files_for_case_%s_sample_%s_to_%s.tsv"
        self.fname_rnaseq_exp = "rnaseq_exp_files_for_PS_%s_Subtype_%s_Stage_%s.tsv"

        self.fname_cases_deprecated = "cases_for_PS_%s_Subtype_%s_Stage_%s.tsv"

        self.fname_case_file = "%s_%s_for_%s_case_%s_file_%s.%s"
        self.fname_mut_anal0 = "mutations_anal_for_study_%s.tsv"
        self.fname_mut_summ0 = "mutations_summ_for_study_%s.tsv"

        self.gdc_fname = ""
        self.gdc_filename = ""
        self.gdc_ouptut_fname = ""
        self.gdc_ouptut_filename = ""

        self.exp_unit = ""
        self.value_col = ""

        # program, primary site, subtype, stage, case_id, samples
        self.df_psi, self.df_subt, self.df_cases = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        self.df_stage, self.df_samples, self.df_files = (
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

 
    def text_normalization(self, x):
        if pd.isna(x):
            return ""

        x = x.lower().strip()
        x = re.sub(r"\bnos\b", "", x)
        x = re.sub(r"[,;()]", " ", x)
        x = re.sub(r"\s+", " ", x)
        return x.strip()

    def map_global_subtype(self, text: str) -> str:
        for k, patterns in self.GLOBAL_SUBTYPE.items():
            if any(p in text for p in patterns):
                return k
        return "other"

    def map_tumor_class(self, text: str) -> str:
        for k, patterns in self.TUMOR_CLASS.items():
            if any(p in text for p in patterns):
                return k
        return "other"

    def map_histology(self, subtype: str) -> str:
        for h, patterns in self.HISTOLOGY.items():
            if any(p in subtype for p in patterns):
                return h
        return "other"

    def map_tissue_subtype(self, global_subtype: str) -> str:
        if self.psi_id in self.SITE_MAP:
            return self.SITE_MAP[self.psi_id].get(global_subtype, global_subtype)

        return global_subtype

    def validate_consistency(self, global_subtype: str, disease_type: str) -> str:
        disease_type = self.text_normalization(disease_type)

        if "squamous" in disease_type and global_subtype != "squamous":
            return "conflict"

        return "ok"

    def simplify_stage(self, stage: str) -> Any:
        if not isinstance(stage, str):
            return None

        stage = stage.strip().upper()

        if stage.startswith("STAGE X"):
            return "missing"
        elif stage.startswith("STAGE IV"):
            return "IV"
        elif stage.startswith("STAGE III"):
            return "III"
        elif stage.startswith("STAGE II"):
            return "II"
        elif stage.startswith("STAGE I"):
            return "I"
        elif stage.startswith("UNKNOWN"):
            return "missing"

        return None

    def set_program(self, prog_id: str):
        self.prog_id = prog_id

        self.root_project = create_dir(self.root0_data, prog_id)
        self.root_summary = create_dir(self.root_project, "summary")

        self.clean_gdc_files()

    def get_primary_sites(
        self, prog_id: str = "TCGA", force: bool = False, verbose: bool = False
    ) -> pd.DataFrame:
        '''
        A primary site (like TCGA-BRCA) is a 'disease'
        root_disease = root0_data / psi_id

        input: project or prog_id, force, verbose
        output: df_psi (dataframe)
        '''

        self.set_program(prog_id)

        if prog_id == 'TCGA':
            fname = self.fname_prim_site_tcga
        else:
            fname = self.fname_prim_site_cbio

        filename = self.root0_data / fname

        if filename.exists() and not force:
            df_psi = pdreadcsv(fname, self.root0_data, verbose=verbose)
            self.df_psi = df_psi

            if prog_id != 'TCGA':
                df_psi = df_psi[df_psi.prog_id == prog_id].copy()
                df_psi.reset_index(drop=True, inplace=True)

            self.df_psi = df_psi

            return df_psi
        
        self.df_psi = pd.DataFrame()
        print("Could not find primary site information for:", prog_id)
        return self.df_psi

    def get_priamry_site_TCGA(self, force: bool = False, verbose: bool = False):

        prog_id = "TCGA"
        self.set_program(prog_id)
        self.df_psi = pd.DataFrame()

        fname = self.fname_prim_site_tcga
        filename = self.root0_data / fname

        if filename.exists() and not force:
            df_psi = pdreadcsv(fname, self.root_project, verbose=verbose)
            self.df_psi = df_psi

            return df_psi

        filters = {"op": "in", "content": {"field": "program.name", "value": [self.prog_id]}}

        params = {
            "filters": json.dumps(filters),
            "fields": "project_id,name,primary_site,disease_type",
            "format": "JSON",
            "size": 1000,
        }

        response = None
        try:
            res = requests.get(self.url_gdc_project, params=params)
            response = res.json()

            if "data" not in response.keys():
                print(f"No data found while searching for '{self.prog_id}'")
                print(">>> response", response)
                self.df_psi = pd.DataFrame()
                return self.df_psi

            hits = response["data"]["hits"]
            print(">>> hits", len(hits))

            df_psi = pd.DataFrame(hits)
            # fix list columns
            for col in df_psi.columns:
                df_psi[col] = df_psi[col].apply(
                    lambda x: ", ".join(x) if isinstance(x, list) else x
                )

            df_psi = df_psi.rename(columns={"id": "psi_id"})

            df_psi = df_psi.sort_values(["primary_site", "disease_type"])

            _ = pdwritecsv(df_psi, fname, self.root_project, verbose=verbose)

        except Exception as e:
            print(f"Error searching for '{self.prog_id}': {e}")
            print(">>> response", response)
            return self.df_psi

        self.df_psi = df_psi

        return df_psi

    def set_primary_site(
        self, psi_id: Any = None, primary_site: Any = None, disease_id: Any = None, verbose: bool = False
    ) -> bool:
        '''
        primary site, here, is a disease
        given psi_id --> root_disease and disease

        input: psi_id (primary site identifier, like TCGA-BRCA) - OR - primary_site (its name/description)
        output: bool (success or failure)
        '''

        self.psi_id = ""
        self.primary_site, self.disease_type, self.disease_name = "", "", ""

        if isinstance(psi_id, str) and psi_id != "":
            dfa = self.df_psi[self.df_psi.psi_id == psi_id]
            if dfa.empty:
                print("No primary site information found for:", psi_id)
                return False
        elif isinstance(primary_site, str) and primary_site != "":
            dfa = self.df_psi[self.df_psi.primary_site == primary_site]
            if dfa.empty:
                print("No primary site information found for:", primary_site)
                return False
        elif isinstance(disease_id, str) and disease_id != "":
            dfa = self.df_psi[self.df_psi.disease_id == disease_id]
            if dfa.empty:
                print("No primary site information found for:", disease_id)
                return False
        else:
            print("No primary site information provided.")
            return False

        row = dfa.iloc[0]

        self.many_cbioportal = True if len(dfa) > 1 else False

        self.psi_id = row.psi_id
        self.primary_site = row.primary_site

        if self.prog_id == 'TCGA':
            self.disease_id = None
            self.disease_type = row.disease_type
            self.disease_name = row['name']
            s_name = 'name'
        else:
            self.disease_id = row.disease_id
            self.disease_type = row.disease_id
            self.disease_name = row.disease_context
            s_name = 'context'

        self.root_disease = create_dir(self.root_project, self.psi_id)
        self.root_samples = create_dir(self.root_disease, 'samples')
        self.root_lfc = create_dir(self.root_disease, 'lfc')
        self.root_mutations = create_dir(self.root_disease, 'mutations')

        if verbose:
            print("\n-----------------------------")
            print(">> psi_id:", self.psi_id)
            print(">> primary_site:", self.primary_site)
            print(">> disease_id:", self.disease_id)
            print(">> disease_type:", self.disease_type)
            print(f">> disease_{s_name}:", self.disease_name)
            print("\n-----------------------------")
            print(">> root disease:", self.root_disease)
            print(">> root samples:", self.root_samples)
            print(">> root lfc:", self.root_lfc)
            print(">> root mutations:", self.root_mutations)
            print("-----------------------------\n")

        self.set_filenames()

        return True

    def get_gdc_progams(self, force: bool = False, verbose: bool = False) -> List:

        filename = os.path.join(self.root0_data, self.fname_programs)

        if os.path.exists(filename) and not force:
            txt = read_txt(filename, verbose=verbose)
            prog_list = eval(txt)
            return prog_list

        params = {"facets": "program.name", "size": 0}

        try:
            res = requests.get(self.url_gdc_project, params=params)
            buckets = res.json()["data"]["aggregations"]["program.name"]["buckets"]

            prog_list = [b["key"] for b in buckets]

            write_txt(str(prog_list), filename, verbose=verbose)

        except Exception as e:
            print(f"No programs found. Error: {e}")
            return []

        return prog_list

    def list_disease_types(self, psi_id: str) -> List:

        self.psi_id = psi_id

        try:
            row = self.df_psi[self.df_psi.psi_id == psi_id].iloc[0]
            deas_type_list = row.disease_type

            if isinstance(deas_type_list, str):
                deas_type_list = eval(deas_type_list)
        except ValueError:
            print("No disease types were found.")
            deas_type_list = []

        self.deas_type_list = deas_type_list

        return deas_type_list

    def set_filenames(self):
        self.fname_cases = self.fname_cases0 % (self.psi_id)
        self.filename_cases = self.root_disease / self.fname_cases

        self.fname_subt = self.fname_subtype0 % (self.psi_id)
        self.filename_subt = self.root_disease / self.fname_subt

    def apply_filter_cases(self, df_cases: pd.DataFrame) -> pd.DataFrame:
        df_cases = df_cases[df_cases.validity == "valid"].copy()
        df_cases = df_cases[df_cases["consistency"] == "ok"]
        df_cases.reset_index(drop=True, inplace=True)

        # frac_threshold:float=0.01,
        # df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()
        # df_cases = df_cases[df_cases["frac"] > frac_threshold]
        # df_cases.reset_index(drop=True, inplace=True)'
        # df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()

        return df_cases

    def get_cases_and_subtypes(
        self,
        batch_size: int = 200,
        do_filter: bool = True,
        debug: bool = False,
        force: bool = False,
        verbose: bool = False,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        calc all subtypes, given and psi_id --> df_cases
        group by ["subtype_global", "subtype_tissue", "stage"] --> df_subt

        filter: NOS (Not Otherwise Specified) → the pathologist could not (or did not) assign a more specific subtype.
        e.g.: "Yes, it's an adenocarcinoma — but we don’t have finer classification"

        input: psi_id = primary site ID
        output: df_cases, df_subt, df_prof
        """

        self.set_filenames()

        if self.filename_cases.exists() and self.filename_subt.exists() and not force:
            df_cases = pdreadcsv(self.fname_cases, self.root_disease, verbose=verbose)
            if "pid" in df_cases.columns:
                df_cases = df_cases.rename(columns={"pid": "psi_id"})
                pdwritecsv(df_cases, self.fname_cases, self.root_disease)

            self.df_cases = df_cases

            if do_filter:
                df_cases = self.apply_filter_cases(df_cases)

            df_subt = self.groupby_case_by_subtypes(df_cases)
            df_prof = self.build_profile(df_cases)

            self.df_cases = df_cases
            self.df_subt = df_subt
            self.df_prof = df_prof

            return df_cases, df_subt, df_prof

        def build_tcga_ontology(df):
            """
            'id', 'primary_site', 'disease_type', 'case_id', 'diagnoses',
            'project.project_id', 'subtype_global', 'stage_ajcc', 'tumor_grade',
            'stage_clin', 'figo_stage', 'tumor_stage', 'stage'],


            print("-------- build -------------")
            print(df.columns)
            print("---------------------------")
            """

            df["primary_site_norm"] = df["primary_site"].apply(self.text_normalization)
            df["disease_type_norm"] = df["disease_type"].apply(self.text_normalization)
            df["diagnosis_norm"] = df["primary_diagnosis"].apply(self.text_normalization)

            df["tumor_class"] = df["diagnosis_norm"].apply(self.map_tumor_class)
            
            """
			for psi_id=='TCGA-ACC' if subtype_global = other --> change to adrenal_cortical_carcinoma
			"""
            if self.psi_id=='TCGA-ACC':
                df["subtype_global"] = [
                    df.iloc[i]["tumor_class"]
                    if (df.iloc[i]["psi_id"] == "TCGA-ACC" and df.iloc[i]["subtype_global"] == "other")
                    else df.iloc[i]["subtype_global"]
                    for i in range(len(df))
                ]
            else:
                df["subtype_global"] = df["diagnosis_norm"].apply(self.map_global_subtype)

            # histology
            df["histology"] = df["subtype_global"].apply(self.map_histology)

            # tissue-specific subtype
            df["subtype_tissue"] = df.apply(
                lambda r: self.map_tissue_subtype(r["subtype_global"]), axis=1
            )

            # consistency check
            df["consistency"] = df.apply(
                lambda r: self.validate_consistency(r["subtype_global"], r["disease_type"]), axis=1
            )

            return df

        # nos -> removes valid dominant classes, like "Endometrioid adenocarcinoma, NOS"
        def classify_validity(row) -> str:

            diag = row["diagnosis_norm"]

            if diag in ["", "unknown", "not reported"]:
                return "invalid"

            if row["subtype_global"] == "other":
                return "ambiguous"

            return "valid"

        def extract_any(x, key, debug: bool = False):
            if debug:
                print("#", x)

            if isinstance(x, list):
                for d in x:
                    if isinstance(d, dict) and key in d and d[key] is not None:
                        return d[key]
            elif isinstance(x, dict):
                return x.get(key)
            return None

        def extract_all(x, key, main_diag):
            # print(">>> extract_all")
            # print(f">>> main_diag '{main_diag}' - key '{key}' - '{x}'")

            dicf = {"diagnosis": "unknown", "ajcc": "unknown"}

            i = 0
            if isinstance(x, list):
                for d in x:
                    # print("int loop", i)
                    if isinstance(d, dict):
                        diag_aux = d.get("primary_diagnosis")

                        if diag_aux == main_diag:
                            if key in d.keys():
                                dicf = {"diagnosis": main_diag, "ajcc": d.get(key)}
                                # print(">> exit1:", i, dicf, '\n\n')
                                return dicf
                            else:
                                dicf = {"diagnosis": main_diag, "ajcc": "unknown"}
                    i += 1

            elif isinstance(x, dict):
                dicf = {"diagnosis": x.get("primary_diagnosis"), "ajcc": x.get(key)}

            # print(">> exit2:", i, dicf, '\n\n')
            return dicf

        def calc_main_diagnosis(df_cases: pd.DataFrame) -> str:
            diag_list = df_cases["diagnoses"].map(lambda x: extract_any(x, "primary_diagnosis"))
            diag_list = [x for x in diag_list if isinstance(x, str) and x.strip()]
            dic = Counter(diag_list)

            dfa = pd.DataFrame({"diag": list(dic.keys()), "n": list(dic.values())})

            dfa = dfa.sort_values("n", ascending=False)

            return dfa.iloc[0].diag

        def unpack_diagnoses(df, main_diag):

            # series = df["diagnoses"].map(lambda x: extract_all(x, "ajcc_pathologic_stage", main_diag))
            series = [
                extract_all(diags, "ajcc_pathologic_stage", main_diag)
                for diags in df_cases["diagnoses"]
            ]

            df["subtype_global"] = [dic["diagnosis"] for dic in series]
            df["stage_ajcc"] = [dic["ajcc"] for dic in series]

            df["primary_diagnosis"] = df["diagnoses"].map(
                lambda x: extract_any(x, "primary_diagnosis")
            )
            df["tumor_grade"] = df["diagnoses"].map(lambda x: extract_any(x, "tumor_grade"))
            df["stage_clin"] = df["diagnoses"].map(lambda x: extract_any(x, "ajcc_clinical_stage"))
            df["figo_stage"] = df["diagnoses"].map(lambda x: extract_any(x, "figo_stage"))
            df["tumor_stage"] = df["diagnoses"].map(lambda x: extract_any(x, "tumor_stage"))

            df["stage"] = (
                df["stage_ajcc"]
                .fillna(df["stage_clin"])
                .fillna(df["figo_stage"])
                .fillna(df["tumor_stage"])
            )

            # df["stage"] = df["stage"].fillna('unknown')
            return df

        # -------------------------- batch loop ---------------------------
        filters = {
            "op": "in",
            "content": {"field": "cases.project.project_id", "value": [self.psi_id]},
        }

        all_hits = []
        from_ = 0
        size_ = batch_size
        total = None
        df_cases = pd.DataFrame()

        try:
            while True:
                print(".", end="")

                params = {
                    "filters": json.dumps(filters),
                    "fields": ",".join(
                        [
                            "case_id",
                            "project.project_id",
                            "primary_site",
                            "disease_type",
                            "diagnoses.primary_diagnosis",
                            "diagnoses.tumor_grade",
                            "diagnoses.ajcc_pathologic_stage",
                            "diagnoses.ajcc_clinical_stagediagnoses.figo_stage",
                            "diagnoses.tumor_stage",
                        ]
                    ),
                    "format": "JSON",
                    "size": size_,
                    "from": from_,
                }

                res = requests.get(self.url_gdc_cases, params=params)
                response = res.json()

                if "data" not in response.keys():
                    print(f"No data found while searching for '{self.psi_id}'")
                    print(">>> response", response)
                    self.df_cases = pd.DataFrame()
                    self.df_subt = pd.DataFrame()
                    self.df_prof = pd.DataFrame()
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
                print(f"No subtypes found for {self.psi_id} ")
                self.df_cases = response
                self.df_subt = pd.DataFrame()
                self.df_prof = pd.DataFrame()
                return self.df_cases, self.df_subt, self.df_prof

            # ------------ lost data? ------------------
            N = len(all_hits)

            if N < total:
                print(
                    f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} "
                )
            else:
                print(f"👉 Returned {N} / Total paginated {total}")

            # ------------ having all hits -------------

            df_cases = pd.json_normalize(all_hits)
            self.df_cases = df_cases

            """
			print("> 1")
			print("----------- 1 ---------------")
			print('rows', len(df_cases), '\ncolumns', df_cases.columns)
			print("---------------------------")


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

            # rename for sanity
            self.df_cases = df_cases

            # ------------------- main_diag -------------------------------------------------------
            main_diag = calc_main_diagnosis(df_cases)
            self.main_diag = main_diag

            df_cases = unpack_diagnoses(df_cases, main_diag)
            """
			'id', 'primary_site', 'disease_type', 'case_id', 'diagnoses',
			'project.project_id', 'subtype_global', 'stage_ajcc', 'tumor_grade',
			'stage_clin', 'figo_stage', 'tumor_stage', 'stage'],
			"""

            df_cases = df_cases.rename(columns={"project.project_id": "psi_id"})

            if debug:
                print("----------- 2 ---------------")
                print(df_cases.head(3).T)
                print("---------------------------")

            self.df_cases2 = df_cases

            if debug:
                print("----------- 3 ---------------")
                print(df_cases.head(3).T)
                print("---------------------------")

            df_cases = build_tcga_ontology(df_cases)

            df_cases["validity"] = df_cases.apply(classify_validity, axis=1)

            df_cases["n"] = 1
            df_cases["frac"] = df_cases["n"] / df_cases["n"].sum()
            df_cases = df_cases.sort_values("n", ascending=False).reset_index(drop=True)
            df_cases.reset_index(drop=True, inplace=True)

            df_cases = df_cases.drop(columns=["id"])
            df_subt = self.groupby_case_by_subtypes(df_cases)

            _ = pdwritecsv(df_cases, self.fname_cases, self.root_disease, verbose=verbose)
            _ = pdwritecsv(df_subt, self.fname_subt, self.root_disease, verbose=verbose)

        except Exception as e:
            print(f"Error for searching diags for '{self.psi_id}'. error: {e}")
            self.df_cases = df_cases
            self.df_subt = pd.DataFrame()
            self.df_prof = pd.DataFrame()
            return self.df_cases, self.df_subt, self.df_prof

        if do_filter:
            df_cases = self.apply_filter_cases(df_cases)

        df_prof = self.build_profile(df_cases)

        self.df_cases = df_cases
        self.df_subt = df_subt
        self.df_prof = df_prof

        return df_cases, df_subt, df_prof

    def group_file_types(self, df_samples: pd.DataFrame) -> pd.DataFrame:
        dic = Counter(df_samples.data_type)
        dfu = pd.DataFrame(dic.items(), columns=["data_type", "n"])
        dfu = dfu.sort_values("n", ascending=False).reset_index(drop=True)
        return dfu

    def groupby_sstate(self, df_cases: pd.DataFrame):

        df_cases["sstage"] = df_cases["stage"].map(lambda x: self.simplify_stage(x))
        df_subt = (
            df_cases.groupby(
                ["psi_id", "subtype_global", "tumor_class", "subtype_tissue", "sstage"],
                dropna=False,
            )
            .size()
            .reset_index(name="n")
        )
        df_subt = df_subt.sort_values("n", ascending=False).reset_index(drop=True)

        self.df_subt = df_subt

        return df_subt

    def groupby_case_by_subtypes(self, df_cases: pd.DataFrame):

        # df_subt = df_cases[cols].copy().drop_duplicates()
        # df_subt = df_subt.sort_values(cols).reset_index(drop=True)

        cols = ["psi_id", "subtype_global", "tumor_class", "subtype_tissue", "stage"]
        df_subt = df_cases.groupby(cols, dropna=False).size().reset_index(name="n")
        df_subt = df_subt.sort_values("n", ascending=False).reset_index(drop=True)

        self.df_subt = df_subt

        return df_subt

    def build_profile(self, df_cases: pd.DataFrame) -> pd.DataFrame:

        def clean_case_profile(diagnoses) -> dict:

            result = {
                "primary_diagnosis": None,
                "stage": None,
                "tumor_grade": None,
                "diagnosis_conflict": False,
                "n_diagnoses": 0,
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
                        d.get("ajcc_pathologic_stage")
                        or d.get("ajcc_clinical_stage")
                        or d.get("figo_stage")
                        or d.get("tumor_stage")
                    )

                    result["stage"] = stage

                    # tumor grade
                    grade = d.get("tumor_grade")
                    if grade and result["tumor_grade"] is None:
                        result["tumor_grade"] = grade

            return result

        profiles = df_cases["diagnoses"].apply(clean_case_profile)
        df_prof = pd.DataFrame(profiles.tolist())
        self.df_prof = df_prof

        return df_prof

    def get_stage_from_cases(self, case_id_list: list):

        df2 = self.df_cases[~pd.isnull(self.df_cases.stage)]
        lista = []
        for case_id in case_id_list:
            dfa = df2[df2.case_id == case_id]

            if dfa.empty:
                lista.append(None)
            else:
                lista.append(dfa.iloc[0].stage)
        return lista

    def set_s_case(self, subtype_global: str, tumor_class: str, subtype_tissue: str):

        self.s_case = f"{self.psi_id}_{self.primary_site}_subtype_{subtype_global}_tumor_{tumor_class}_subtype_tissue_{subtype_tissue}"

        if len(self.s_case) > 180:
            self.s_case = f"{self.psi_id}_{self.primary_site[:40]}_subtype_{subtype_global[:40]}_tumor_{tumor_class[:40]}_tissue_{subtype_tissue[:40]}"

        self.s_case = title_replace(self.s_case)

    def get_samples_for_subtypes(
        self,
        subtype_global: str,
        tumor_class: str,
        subtype_tissue: str,
        batch_cases: int = 50,
        batch_size: int = 200,
        force: bool = False,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """
        return all samples given a list of cases
        for psi_id, subtype_global, tumor_class, subtype_tissue

        input: psi_id, subtype_global, tumor_class, subtype_tissue
        output: df_samples
        """
        self.df_samples = pd.DataFrame()

        self.set_s_case(subtype_global, tumor_class, subtype_tissue)

        df_cases, _, _ = self.get_cases_and_subtypes(
            batch_size=batch_size, do_filter=False, force=False, verbose=verbose
        )

        self.df_cases = df_cases

        if df_cases is None or df_cases.empty:
            print(f"No cases found while searching for '{self.s_case}'")
            return self.df_samples

        """
		lista=[]
		if isinstance(sstage, str):
			if sstage.startswith('I'):
				stage = 'Stage ' + sstage
			elif sstage ==  'missing':
				lista = ['unknown', 'X']
		stage = sstage

		if len(lista) > 0:
			df_cases = df_cases[(df_cases.subtype_global == subtype_global) & 
								(df_cases.tumor_class == tumor_class) &
								(df_cases.subtype_tissue == subtype_tissue) &
								(df_cases.stage.isin(lista))]
		else:
			df_cases = df_cases[(df_cases.subtype_global == subtype_global) & 
								(df_cases.tumor_class == tumor_class) &
								(df_cases.subtype_tissue == subtype_tissue) &
								(df_cases.stage == stage)]
		"""

        df_cases = df_cases[
            (df_cases.subtype_global == subtype_global) & 
            (df_cases.tumor_class == tumor_class) & 
            (df_cases.subtype_tissue == subtype_tissue)
        ]

        df_cases = df_cases.copy().reset_index(drop=True)
        self.df_cases = df_cases

        if df_cases.empty:
            print(f"No cases found for {self.s_case}")
            return self.df_samples

        fname = self.fname_samples0 % (self.s_case)
        fname = title_replace(fname)
        filename = self.root_samples / fname

        if filename.exists() and not force:
            df_samples = pdreadcsv(fname, self.root_samples, verbose=verbose)
            self.df_samples = df_samples

            return df_samples

        case_id_list = list(df_cases.case_id)
        case_id_list.sort()

        # s_case_id_list3 = f"[{','.join(case_id_list[:3])}]"

        N_cases = len(case_id_list)
        print(f">>> {N_cases} cases")

        # -------------------------- batch loop ---------------------------
        all_hits = []
        from_ = 0
        size_ = batch_size
        total = None
        df_samples = pd.DataFrame()

        # print("Searching: ", end='')

        ini = -batch_cases
        end = 0

        while True:
            ini += batch_cases
            end += batch_cases

            if ini >= N_cases:
                break

            if end > N_cases:
                end = N_cases

            print(f"{ini}-{end} ", end="")

            lista = case_id_list[ini:end]
            # print("\n>>>", len(lista), lista)

            filters = {"op": "in", "content": {"field": "cases.case_id", "value": lista}}

            try:
                while True:
                    print(".", end="")

                    params = {
                        "filters": json.dumps(filters),
                        "fields": ",".join(
                            [
                                "file_id",
                                "file_name",
                                "data_type",
                                "data_format",
                                "experimental_strategy",
                                "cases.case_id",
                                "cases.submitter_id",
                                "cases.samples.sample_id",
                                "cases.samples.submitter_id",
                                "cases.samples.sample_type",
                            ]
                        ),
                        "format": "JSON",
                        "size": size_,
                        "from": from_,
                    }

                    res = requests.get(self.url_gdc_files, params=params)
                    response = res.json()

                    if "data" not in response.keys():
                        print(f"No data found while searching for '{self.psi_id}' cases {case_id_list}")
                        print(">>> response", response)
                        self.df_samples = pd.DataFrame()
                        return self.df_samples

                    hits = response.get("data", {}).get("hits", [])

                    if total is None:
                        total = response["data"]["pagination"]["total"]

                    if not hits:
                        break

                    all_hits.extend(hits)
                    from_ += size_

                # print("\n")

                if all_hits == []:
                    print(f"No files were found for {self.psi_id} cases {case_id_list}")
                    self.df_samples = pd.DataFrame()
                    return self.df_samples

                # ------------ lost data? ------------------
                N = len(all_hits)

                if N < total:
                    print(
                        f"⚠️ Warning: results truncated — consider pagination - all hits = {N};  Total paginated {total} "
                    )
                else:
                    if verbose:
                        print(f"👉 Returned {N} / Total paginated {total}")

                # ------------ having all hits -------------
                records = []

                for hit in all_hits:
                    for case in hit.get("cases", []):
                        for sample in case.get("samples", []):
                            records.append(
                                {
                                    "case_id": case["case_id"],
                                    "submitter_id": case["submitter_id"],
                                    "sample_id": sample["sample_id"],
                                    "sample_type": sample["sample_type"],
                                    "barcode_sample": sample["submitter_id"],
                                    "file_id": hit["file_id"],
                                    "file_name": hit["file_name"],
                                    "data_type": hit["data_type"],
                                    "data_format": hit["data_format"],
                                }
                            )

                df_samples = pd.DataFrame(records)
                self.df_samples = df_samples
                cols = list(df_samples.columns)

                # 🔹 Metadata
                df_samples["psi_id"] = self.psi_id
                df_samples["subtype_global"] = subtype_global
                df_samples["tumor_class"] = tumor_class
                df_samples["subtype_tissue"] = subtype_tissue
                df_samples["stage"] = self.get_stage_from_cases(df_samples.case_id.tolist())

                cols = ["psi_id", "subtype_global", "tumor_class", "subtype_tissue", "stage"] + cols

                df_samples = df_samples.sort_values(
                    ["case_id", "sample_type"], ascending=[False, False]
                ).reset_index(drop=True)
                df_samples.reset_index(drop=True, inplace=True)

                _ = pdwritecsv(df_samples, fname, self.root_samples, verbose=verbose)

            except Exception as e:
                print(f"Error for searching files for {self.s_case}'. error: {e}")
                self.df_samples = df_samples
                return df_samples

        self.df_samples = df_samples

        return df_samples

    def get_table_given_fileID(
        self,
        case_id: str,
        file_id: str,
        sample_type: str,
        file_type: str,
        timeout: int = 120,
        force: bool = False,
        debug: bool = False,
        verbose: bool = False,
    ) -> tuple[Any, Any]:
        """
        Retrieve any kind of table like: RNA or Proteomic expression
        input: case_id and file_id
        output: the desired file
        """

        if debug:
            print("case_id:", case_id)
            print("file_id:", file_id)
            print("sample_type:", sample_type)
            print("file_type:", file_type)

        if not file_id and not isinstance(file_id, str):
            print("No file_id defined.")

            self.df_table = pd.DataFrame()
            return self.df_table

        file_type = file_type.strip()

        if file_type == "Gene Expression Quantification":
            is_expression = True
            type_of_file = "tsv"
            root = self.root_lfc
        elif file_type == "Raw Simple Somatic Mutation":
            is_expression = False
            type_of_file = "tar.gz"
            root = self.root_mutations
        else:
            print(f"Develope the method for this file type {file_type}")
            raise Exception("\n------------ stop ---------------\n")

        # self.fname_case_file = "%s_%s_for_%s_case_%s_file_%s.%s"
        fname = self.fname_case_file % (
            file_type,
            sample_type,
            self.psi_id,
            case_id,
            file_id,
            type_of_file,
        )
        fname = title_replace(fname)
        filename = root / fname

        if os.path.exists(filename) and not force:
            if is_expression:
                df_table = pdreadcsv(fname, root, verbose=verbose)

                changed = False
                if "gene_id" in df_table.columns:
                    df_table = df_table.rename(columns={"gene_id": "geneid"})
                    changed = True
                if "gene_type" in df_table.columns:
                    df_table = df_table.rename(columns={"gene_type": "biotype"})
                    changed = True

                if changed:
                    _ = pdwritecsv(df_table, fname, root, verbose=verbose)


                self.df_table = df_table
                return df_table, filename
            else:
                return filename, filename

        if verbose:
            print("Downloading: ", end="")
        try:
            url_file = f"{self.url_gdc_data}/{file_id}"

            with requests.get(url_file, stream=True, timeout=timeout) as r:
                if r.status_code != 200:
                    print("Error:", r.status_code)
                    print("URL:", url_file)
                    print("Content-Type:", r.headers.get("Content-Type"))
                    try:
                        print("Preview:", r.text[:500])
                    except Exception:
                        print("Could not decode error body as text for file:", file_id)
                    return None, filename

                with open(filename, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

            if not is_expression:
                return filename, filename

            df_table = pd.read_csv(filename, sep="\t", comment="#")
            df_table = self.clean_expression_table(df_table)

            cols = ["gene_id", "symbol", "gene_type", "unstranded", "counts", "stranded_second",
                    "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded",]
            df_table = df_table[cols]

            cols = ["geneid", "symbol", "biotype", "unstranded", "counts", "stranded_second",
                    "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded",]
            df_table.columns = cols

            _ = pdwritecsv(df_table, fname, self.root_lfc, verbose=verbose)
               

        except Exception as e:
            s_error = f"Download error for '{self.psi_id}', '{file_type}', case {case_id} and {file_id}: {e}"
            if verbose:
                print(s_error)

            self.df_table = pd.DataFrame()
            return self.df_table, filename

        self.df_table = df_table
        return df_table, filename


    def clean_expression_table(self, df: pd.DataFrame) -> pd.DataFrame:

        # Remove summary rows (N_*)
        df = df[~df["gene_id"].str.startswith("N_")]

        # Keep only valid Ensembl genes
        df = df[df["gene_id"].str.startswith("ENSG")]

        # Remove version from gene_id (ENSG... -> ENSG...)
        df["gene_id"] = df["gene_id"].str.split(".").str[0]

        df = df.rename(columns={"gene_name": "symbol", "stranded_first": "counts"})

        return df

    def get_table_searching_for_fileID(
        self, data_type: str, sample_type: str, file_id: str, verbose: bool = False
    ) -> pd.DataFrame:

        data_type2 = title_replace(data_type)
        sample_type2 = title_replace(sample_type)

        files = [
            x
            for x in os.listdir(self.root_disease)
            if file_id in x and data_type2 in x and sample_type2 in x
        ]

        if len(files) == 0:
            print(f"No files found for {file_id}.")
            self.df_table = pd.DataFrame()
            return self.df_table

        if len(files) > 1:
            print(f"Multiple files found for {file_id}. Using the first one.")

        fname = files[0]
        df_table = pdreadcsv(fname, self.root_disease, verbose=verbose)
        self.df_table = df_table

        return df_table

    def calc_file_expression_tumor_normal_gtex(self, imax_samples: int = 200, force: bool = False,
                                               verbose: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        '''

        flow:
            calc_file_expression_tumor_normal_gtex
                get_dic_expression_tumor_and_normal

                    _, df_tumor_samples, _, _ = self.get_filtered_tables( sample_type_term="Primary Tumor", verbose=verbose )
                    _, df_normal_samples, _, _ = self.get_filtered_tables( sample_type_term="Solid Tissue Normal", verbose=verbose )

                    get expression filenames:
                    dff_normal = df_normal_samples[df_normal_samples.data_type == "Gene Expression Quantification"]
                    dff_tumor = df_tumor_samples[df_tumor_samples.data_type == "Gene Expression Quantification"]

                    build 2 dicts having the expression tables
                    output:  dic_tumor, dic_normal

                # from dict to df -> df tumor and normal
                df_tumor, df_normal = self.prepare_normal_tumor_tables()

                # get GTEx normal control tissue data
                df_gtex_ctrl, _ = self.get_gtex_control()


            output: df_tumor, df_normal, df_gtex_ctrl
        
        '''

        fname_exp_tumor = self.fname_exp_tumor%(self.psi_id)
        filename_tumor = self.root_lfc / fname_exp_tumor

        fname_exp_normal = self.fname_exp_normal%(self.psi_id)
        filename_normal = self.root_lfc / fname_exp_normal

        fname_exp_gtex = self.fname_exp_gtex%(self.psi_id)
        filename_gtex = self.root_lfc / fname_exp_gtex

        if filename_tumor.exists() and filename_normal.exists() and filename_gtex.exists() and not force:
            df_tumor = pdreadcsv(fname_exp_tumor, self.root_lfc, verbose=verbose)
            df_normal = pdreadcsv(fname_exp_normal, self.root_lfc, verbose=verbose)
            df_gtex_ctrl = pdreadcsv(fname_exp_gtex, self.root_lfc, verbose=verbose)

            self.df_tumor = df_tumor
            self.df_normal = df_normal
            self.df_gtex_ctrl = df_gtex_ctrl

            return df_tumor, df_normal, df_gtex_ctrl

        dic_tumor, dic_normal = self.get_dic_expression_tumor_and_normal(verbose=verbose)
        self.dic_tumor = dic_tumor
        self.dic_normal = dic_normal

        if len(dic_tumor) == 0 and len(dic_normal) == 0:
            print(f"Insufficient expression data for {self.psi_id}.")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

        df_tumor, df_normal = self.prepare_normal_tumor_tables(
            dic_tumor, dic_normal, imax_tumor=imax_samples, imax_normal=imax_samples, verbose=verbose
        )

        df_gtex_ctrl, _ = self.get_gtex_control(Nsamples=15, force=False, verbose=verbose)

        if "gene_id" in df_tumor.columns:
            df_tumor = df_tumor.rename(columns={"gene_id": "geneid"})
        if "gene_type" in df_tumor.columns:
            df_tumor = df_tumor.rename(columns={"gene_type": "biotype"})

        if "gene_id" in df_normal.columns:
            df_normal = df_normal.rename(columns={"gene_id": "geneid"})
        if "gene_type" in df_normal.columns:
            df_normal = df_normal.rename(columns={"gene_type": "biotype"})

        self.df_tumor = df_tumor
        self.df_normal = df_normal
        self.df_gtex_ctrl = df_gtex_ctrl

        _ = pdwritecsv(df_tumor, fname_exp_tumor, self.root_lfc)
        _ = pdwritecsv(df_normal, fname_exp_normal, self.root_lfc)
        _ = pdwritecsv(df_gtex_ctrl, fname_exp_gtex, self.root_lfc)

        return df_tumor, df_normal, df_gtex_ctrl

    def get_dic_expression_tumor_and_normal(self, verbose: bool = False) -> Tuple[dict, dict]:

        #----------- tumor --------------------------------------------------
        _, df_tumor_samples, _, _ = self.get_filtered_tables(
            sample_type_term="Primary Tumor", verbose=verbose
        )
        #----------- normal --------------------------------------------------
        _, df_normal_samples, _, _ = self.get_filtered_tables(
            sample_type_term="Solid Tissue Normal", verbose=verbose
        )

        if df_tumor_samples is None or df_tumor_samples.empty:
            if verbose:
                print(f"No tumor expression data found for {self.psi_id}.")
            return {}, {}

        if df_normal_samples is None or df_normal_samples.empty:
            if verbose:
                print(f"No normal expression data found for {self.psi_id}.")
            return {}, {}

        self.file_type_list = np.unique(df_tumor_samples.data_type)

        dff_normal = df_normal_samples[df_normal_samples.data_type == "Gene Expression Quantification"]
        dff_normal.reset_index(drop=True, inplace=True)
        self.dff_normal = dff_normal

        dff_tumor = df_tumor_samples[df_tumor_samples.data_type == "Gene Expression Quantification"]
        dff_tumor.reset_index(drop=True, inplace=True)
        self.dff_tumor = dff_tumor

        # raise ValueError("\n\n------------ stop --------------\n\n")

        if verbose:
            print(
                f"There are {len(dff_tumor)} tumor and {len(dff_normal)} normal Gene Expression tables"
            )

        if len(dff_tumor) == 0 and len(dff_normal) == 0:
            print(f"No valid expression data found for {self.psi_id}.")
            return {}, {}

        print("Dowloading normal files:", end=" ")
        cols = ["geneid", "symbol", "biotype", "counts"]

        dic_normal = {}
        for i, row in dff_normal.iterrows():
            if i % 10 == 0:
                print(i, end="")
                  
            case_id = row.case_id
            file_id = row.file_id

            dfexp, filename_normal = self.get_table_given_fileID(
                                            case_id=case_id,
                                            file_id=file_id,
                                            sample_type="normal",
                                            file_type="Gene Expression Quantification",
                                            force=False,
                                            verbose=verbose,
                                        )
            if dfexp is None or dfexp.empty:
                print("x", end="")
                continue
            print(".", end="")

            self.dfexp_normal = dfexp

            try:
                dfexp = dfexp[cols]
            except Exception as e:
                print(f"Error occurred while processing normal file {filename_normal}: {e}")
                print(dfexp.shape)
                print(dfexp.columns)
                continue
            dic_normal[f"normal_{file_id}"] = dfexp

        print("")
        if verbose:
            print(f" -> {len(dff_normal)}")

        print("Dowloading tumor files:", end=" ")
        dic_tumor = {}
        for i, row in dff_tumor.iterrows():
            if i % 10 == 0:
                print(i, end="")

            case_id = row.case_id
            file_id = row.file_id

            dfexp, filename_tumor = self.get_table_given_fileID(
                                            case_id=case_id,
                                            file_id=file_id,
                                            sample_type="tumor",
                                            file_type="Gene Expression Quantification",
                                            force=False,
                                            verbose=verbose,
                                        )
            if dfexp is None or dfexp.empty:
                print("x", end="")
                continue
            print(".", end="")

            try:
                dfexp = dfexp[cols]
            except Exception as e:
                print(f"Error occurred while processing tumor file {filename_tumor}: {e}")
                print(dfexp.shape)
                print(dfexp.columns)
                continue
            dic_tumor[f"tumor_{file_id}"] = dfexp

        print("")
        if verbose:
            print(f" -> {len(dff_tumor)}")

        self.dic_tumor = dic_tumor
        self.dic_normal = dic_normal

        return dic_tumor, dic_normal

    def get_case_id(self, barcode: str) -> str:

        self.clean_gdc_files()

        if not isinstance(barcode, str) or len(barcode) < 3:
            print(f"Barcode bad formated {barcode}.")
            return ""

        try:
            filters = {"op": "in", "content": {"field": "submitter_id", "value": [barcode]}}

            params = {
                "filters": json.dumps(filters),
                "fields": "case_id,submitter_id",
                "format": "JSON",
                "size": 1,
            }

            response = requests.get(self.url_gdc_cases, params=params)
            data = response.json()

            hits = data.get("data", {}).get("hits", [])
            if not hits:
                raise ValueError(f"No case found for barcode {barcode}")

        except Exception as e:
            print(f"No data found for {barcode}. error: {e}")
            return ""

        return hits[0]["case_id"]

    def get_representative_geneids(self, 
        dfs: list[pd.DataFrame],
        min_fraction: float = 0.75,
    ) -> pd.DataFrame:
        """
        Return genes present in more than min_fraction of dataframes.

        For 10 dataframes and min_fraction=0.75:
        strict >75% means present in at least 8 dataframes.

        Presence is counted once per dataframe, even if duplicated inside a dataframe.
        """

        gene_cols: list[str] = ["geneid", "symbol"]

        n = len(dfs)
        if n == 0:
            return pd.DataFrame(columns=gene_cols + ["n_dfs", "fraction"])

        min_count = math.floor(n * min_fraction) + 1  # strict > min_fraction

        counter = Counter()

        for df in dfs:
            missing = [c for c in gene_cols if c not in df.columns]
            if missing:
                raise ValueError(f"Missing columns in dataframe: {missing}")

            genes_in_df = (
                df[gene_cols]
                .dropna(subset=gene_cols)
                .astype(str)
                .drop_duplicates()
            )

            counter.update(map(tuple, genes_in_df.to_numpy()))

        result = (
            pd.DataFrame(
                [(geneid, symbol, count) for (geneid, symbol), count in counter.items()],
                columns=gene_cols + ["n_dfs"],
            )
            .assign(fraction=lambda x: x["n_dfs"] / n)
            .query("n_dfs >= @min_count")
            .sort_values(["n_dfs"] + gene_cols, ascending=[False] + [True] * len(gene_cols))
            .reset_index(drop=True)
        )

        return result
    
    def get_common_gene_list(self, dic_tumor: dict, min_fraction: float = 0.75) -> np.ndarray:
        df_list = []

        cols = ["geneid", "symbol", "biotype", "counts"]

        for _, dfa in dic_tumor.items():
            if dfa is None or dfa.empty:
                continue

            if "gene_id" in dfa.columns:
                dfa = dfa.rename(columns={"gene_id": "geneid"})
            if "gene_type" in dfa.columns:
                dfa = dfa.rename(columns={"gene_type": "biotype"})

            dfa = dfa[cols]
            df_list.append(dfa)

        dfq = self.get_representative_geneids(df_list, min_fraction=min_fraction)
        lista = np.unique(dfq.geneid.to_list())
        return lista


    def prepare_normal_tumor_tables(
        self,
        dic_tumor: dict,
        dic_normal: dict,
        imax_tumor: int = 12,
        imax_normal: int = 12,
        verbose: bool = False,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process tumor and normal dictionary expression data.

        input:
                dic_tumor
                dic_normal
                    # get the most common geneids to merge all tumor tables
                    lista = self.get_common_gene_list(dic_tumor, min_fraction=0.75)
                verbose: bool, whether to print verbose messages
        output:
                df_tumor and df_normal tables
        """

        cols = ["geneid", "symbol", "counts"]
        common_cols = ["geneid", "symbol"]        

        # ----------- Normal tissue ----------------
        df_normal = pd.DataFrame()
        if verbose:
            print(">>> Processing normal data:", len(dic_normal))
        i = 0
        for _, dfa in dic_normal.items():
            if dfa is None or dfa.empty:
                continue

            i += 1
            # print(i, end=' ')
            if "gene_id" in dfa.columns:
                dfa = dfa.rename(columns={"gene_id": "geneid"})
            if "gene_type" in dfa.columns:
                dfa = dfa.rename(columns={"gene_type": "biotype"})

            dfa = dfa[cols]
            dfa = dfa.rename(columns={"counts": f"normal_{i}"})

            if df_normal.empty:
                df_normal = dfa
            else:
                if i <= imax_normal:
                    df_normal = df_normal.merge(dfa, on=common_cols, how="outer")
                else:
                    if verbose:
                        print(">>> dfa", len(dfa), ",".join(dfa.symbol[:30]))
                    break

        # ----------- tumor ----------------
        lista = self.get_common_gene_list(dic_tumor, min_fraction=0.75)
        df_tumor = pd.DataFrame()

        if len(lista) == 0:
            if verbose:
                print(">>> No common genes found.")
            return df_tumor, df_normal

        if verbose:
            print(">>> Processing tumor data:", len(dic_tumor))

        i = 0
        for _, dfa in dic_tumor.items():
            if dfa is None or dfa.empty:
                continue

            i += 1
            # print(i, end=' ')
            if "gene_id" in dfa.columns:
                dfa = dfa.rename(columns={"gene_id": "geneid"})
            if "gene_type" in dfa.columns:
                dfa = dfa.rename(columns={"gene_type": "biotype"})

            dfa = dfa[cols]

            dfa = (
                dfa.dropna(subset=['geneid', 'symbol'])
                .drop_duplicates(['geneid', 'symbol'])
            )
            if dfa.empty:
                continue

            dfa = dfa.rename(columns={"counts": f"tumor_{i}"})

            dfa = dfa[dfa.geneid.isin(lista)]
            if dfa.empty:
                continue
            
            dfa.reset_index(drop=True, inplace=True)
    
            if df_tumor.empty:
                df_tumor = dfa
            else:
                if i <= imax_tumor:
                    df_tumor = df_tumor.merge(dfa, on=common_cols, how="outer")
                else:
                    if verbose:
                        print(">>> dfa", len(dfa), ",".join(dfa.symbol[:30]))
                    break

        return df_tumor, df_normal

    def resolve_mutation_profile(self, study_id: str, timeout: int = 60) -> str:
        '''
        here is the cBioPortl endpoint

        there is no PSI_ID 
        one must know the correct study_id

        input:  study_id
        output: molecular_profile_id
        '''

        candidates = [
            f"{study_id}_mutations",
            f"{study_id}_mutations_extended",
        ]

        for mp in candidates:
            url = f"{self.url_cbioportal}/molecular-profiles/{mp}"
            if requests.get(url, timeout=timeout).ok:
                return mp

        raise ValueError(f"\n\n------------ No mutation profile found for {study_id} ------------- \n\n")

    def get_cBioportal_mutations_from_samples(
        self,
        barcode_sample_list: Iterable[str],
        study_id: str,
        session: Optional[requests.Session] = None,
        timeout: int = 60,
    ) -> pd.DataFrame:
        """
        Fetch mutation records from cBioPortal for a list of sample IDs.

        Parameters
        ----------
        barcode_sample_list : Iterable[str]
                cBioPortal sample IDs, e.g. ["TCGA-GC-A3BM-01", "TCGA-XF-A9SY-01"].
        study_id : str
                cBioPortal study ID, e.g. "blca_tcga".
        molecular_profile_id : str | None
                Mutation profile ID. If None, defaults to "{study_id}_mutations".
        base_url : str
                cBioPortal API base URL.
        session : requests.Session | None
                Optional requests session.
        timeout : int
                Request timeout in seconds.

        Returns
        -------
        pd.DataFrame
                Mutation table. Empty DataFrame if nothing is returned.

        Notes
        -----
        - This function assumes all sample_ids belong to the same study.
        - If your samples come from multiple studies, call the function per study.
        """

        # molecular_profile_id = f"{study_id}_mutations"
        molecular_profile_id = self.resolve_mutation_profile(study_id)

        http = session or requests.Session()

        url = f"{self.url_cbioportal}/molecular-profiles/{molecular_profile_id}/mutations/fetch"

        payload = {"sampleIds": barcode_sample_list}

        headers = {
            "Accept": "application/json",
            "Content-Type": "application/json",
        }

        resp = http.post(url, json=payload, headers=headers, timeout=timeout)

        # Helpful error message from cBioPortal
        if not resp.ok:
            msg = ""
            try:
                msg = resp.json()
            except Exception:
                msg = resp.text

            raise RuntimeError(
                f"cBioPortal request failed: HTTP {resp.status_code} | "
                f"profile={molecular_profile_id} | details={msg}"
            )

        data = resp.json()
        if not data:
            if verbose:
                print(f"Error: cBioPortal URL: {url}")
                print(
                    f"No mutations found for molecular profile '{molecular_profile_id}' barcodes: {barcode_sample_list}."
                )
            return pd.DataFrame()

        df = pd.DataFrame(data)

        cols_ori = list(df.columns)

        if "tumorAltCount" in df.columns:
            cols_ori.remove("tumorAltCount")
        else:
            df["tumorAltCount"] = None

        cols = cols_ori + ["tumorAltCount"]
        df = df[cols]

        # self.df = df
        # raise Exception('stop3')

        # the selected cols + others not listed
        # cols = [c for c in df.columns if c in df.columns.to_list()] + [c for c in df.columns if c not in cols]
        # self.df = df

        df["keyword"] = [x.split(" ")[0] if isinstance(x, str) else x for x in df["keyword"]]

        dic_rename = {
            "uniqueSampleKey": "unique_sample_key",
            "uniquePatientKey": "unique_patient_key",
            "molecularProfileId": "molecular_profile_id",
            "sampleId": "barcode_sample",
            "patientId": "barcode",
            "entrezGeneId": "entrez_gene_id",
            "studyId": "psi_id",
            "center": "center",
            "mutationStatus": "mutation_status",
            "validationStatus": "validation_status",
            "tumorRefCount": "tumor_ref_count",
            "normalRefCount": "normal_ref_count",
            "startPosition": "start",
            "endPosition": "end",
            "referenceAllele": "ref_allele",
            "proteinChange": "protein_mut",
            "mutationType": "mutation_type",
            "ncbiBuild": "ncbi_build",
            "variantType": "variant_type",
            "keyword": "symbol",
            "chr": "chr",
            "variantAllele": "variant_allele",
            "refseqMrnaId": "refseq_mrna_id",
            "proteinPosStart": "protein_pos_start",
            "proteinPosEnd": "protein_pos_end",
            "tumorAltCount": "tumor_alt_count",
        }

        rename_cols = [dic_rename.get(col, col) for col in df.columns]

        df.columns = rename_cols

        df["sample"] = [x.split("-")[-1] for x in df["barcode_sample"]]

        order_cols = [
            "psi_id",
            "molecular_profile_id",
            "barcode",
            "sample",
            "barcode_sample",
            "symbol",
            "refseq_mrna_id",
            "entrez_gene_id",
            "protein_mut",
            "mutation_type",
            "mutation_status",
            "ref_allele",
            "variant_allele",
            "variant_type",
            "chr",
            "start",
            "end",
            "validation_status",
            "protein_pos_start",
            "protein_pos_end",
            "tumor_alt_count",
            "ncbi_build",
            "center",
            "tumor_ref_count",
            "unique_sample_key",
            "unique_patient_key",
        ]

        df = df[order_cols]
        # self.df = df

        return df

    def change_cbioportal_studyid(self, study_id: str) -> str:
        """
        Normalize TCGA study IDs to cBioPortal PanCancer Atlas studies.

        In cBioPortal:
                COAD = colon adenocarcinoma
                READ = rectum adenocarcinoma

                👉 In PanCancer Atlas they are merged into one cohort:
        """

        dic = {
            "acc_tcga": "acc_tcga_pan_can_atlas_2018",
            "luad_tcga": "luad_tcga_pan_can_atlas_2018",
            "lusc_tcga": "lusc_tcga_pan_can_atlas_2018",
            "coad_tcga": "coadread_tcga_pan_can_atlas_2018",
            "read_tcga": "coadread_tcga_pan_can_atlas_2018",
            "brca_tcga": "brca_tcga_pan_can_atlas_2018",
            "gbm_tcga": "gbm_tcga_pan_can_atlas_2018",
            "ov_tcga": "ov_tcga_pan_can_atlas_2018",
            "skcm_tcga": "skcm_tcga_pan_can_atlas_2018",
            "ucec_tcga": "ucec_tcga_pan_can_atlas_2018",
            "stad_tcga": "stad_tcga_pan_can_atlas_2018",
            "blca_tcga": "blca_tcga_pan_can_atlas_2018",
            "hnsc_tcga": "hnsc_tcga_pan_can_atlas_2018",
            "kirc_tcga": "kirc_tcga_pan_can_atlas_2018",
            "kirp_tcga": "kirp_tcga_pan_can_atlas_2018",
            "lihc_tcga": "lihc_tcga_pan_can_atlas_2018",
            "prad_tcga": "prad_tcga_pan_can_atlas_2018",
            "thca_tcga": "thca_tcga_pan_can_atlas_2018",
            "esca_tcga": "esca_tcga_pan_can_atlas_2018",
            "paad_tcga": "paad_tcga_pan_can_atlas_2018",
            "kich_tcga": "kich_tcga_pan_can_atlas_2018",  # kidney chromophobe
            "sarc_tcga": "sarc_tcga_pan_can_atlas_2018",
            "pcpg_tcga": "pcpg_tcga_pan_can_atlas_2018",
            "tgct_tcga": "tgct_tcga_pan_can_atlas_2018",
            "thym_tcga": "thym_tcga_pan_can_atlas_2018",
            "meso_tcga": "meso_tcga_pan_can_atlas_2018",
            "ucs_tcga": "ucs_tcga_pan_can_atlas_2018",
            "uvm_tcga": "uvm_tcga_pan_can_atlas_2018",
            "chol_tcga": "chol_tcga_pan_can_atlas_2018",
            "dlbc_tcga": "dlbc_tcga_pan_can_atlas_2018",
        }

        return dic.get(study_id, study_id)

    def set_mutation_filenames(self):
        self.fname_mut_anal = self.fname_mut_anal0 % (self.s_case)
        self.fname_mut_anal = title_replace(self.fname_mut_anal)
        self.filename_mutanal = self.root_mutations / self.fname_mut_anal

        self.fname_mut_summ = self.fname_mut_summ0 % (self.s_case)
        self.fname_mut_summ = title_replace(self.fname_mut_summ)
        self.filename_mutsumm = self.root_mutations / self.fname_mut_summ

    def get_df_mut_transform_mutation_table(
        self,
        psi_id: str,
        barcode_sample_list: List[str],
        session: Optional[requests.Session] = None,
        timeout: int = 60,
        force: bool = False,
        verbose: bool = False,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:

        """
		if TCGA remove the last characters if len > 2
		TCGA-OR-A5J2-01A -> TCGA-OR-A5J2-01
		"""
        barcode_list = self.prepare_barcode_sample_list(barcode_sample_list)

        if self.prog_id == 'TCGA':
            study_id = psi_id

            if study_id[0].isupper():
                mat = study_id.lower().split("-")
                # cBioPortal disease - tcga
                study_id = mat[1] + "_" + mat[0]

            study_id = self.change_cbioportal_studyid(study_id)
        elif self.prog_id == 'CPTAC':
            if psi_id == 'CPTAC-3':
                study_id = 'paad_cptac_2021'
        
        self.study_id = study_id

        print(f"\n>>> {study_id} --> {self.s_case} len = {len(barcode_list)} - {barcode_list[:5]}...")

        self.set_mutation_filenames()

        if (
            os.path.exists(self.filename_mutanal)
            and os.path.exists(self.filename_mutsumm)
            and not force
        ):
            dff = pdreadcsv(self.fname_mut_summ, self.root_disease, verbose=verbose)
            df_mut = pdreadcsv(self.fname_mut_anal, self.root_disease, verbose=verbose)
            return dff, df_mut

        """
			df_mut cols: ["sample_id", "barcode_sample", "psi_id", "mol_profile_id","gene",
			"entrez_gene_id", "protein_mut", "mutation_type", "mutation_status",
			"variant_type", "chr", "start", "end",
			"ref_allele", "tumor_seq_allele"]		
		"""
        df_mut = self.get_cBioportal_mutations_from_samples(
            barcode_sample_list=barcode_sample_list,
            study_id=study_id,
            session=session,
            timeout=timeout,
        )

        self.df_mut = df_mut

        if df_mut.empty:
            print("No mutations found for these samples.")
            return pd.DataFrame(), pd.DataFrame()

        # --------------- map main cols from df_mut ------------------------
        """
		order_cols = ['barcode_sample', 'barcode_sample', 'psi_id', 'mol_profile_id',
			'symbol', 'refseq_mrna_id', 'entrez_gene_id', 
			'protein_mut', 'mutation_type', 'mutation_status',
			'ref_allele', 'variant_allele', 'variant_type', 
			'chr', 'start', 'end', 
			'validation_status', 'protein_pos_start', 'protein_pos_end', 'tumor_alt_count',
			'ncbi_build', 'center', 'tumorRefCount', 'unique_sample_key',
			'unique_patient_key']
		"""
        dff = (
            df_mut.groupby(
                [
                    "psi_id",
                    "barcode",
                    "barcode_sample",
                    "symbol",
                    "refseq_mrna_id",
                    "entrez_gene_id",
                    "protein_mut",
                    "mutation_type",
                    "variant_type",
                    "chr",
                ]
            )
            .size()
            .reset_index(name="n_mutations")
        )

        dff = dff[dff.barcode.notna()]
        dff = dff[dff.barcode_sample.notna()]
        dff = dff[dff.entrez_gene_id.notna()]

        dff = dff.sort_values(["barcode", "symbol", "protein_mut"])
        dff = dff.reset_index(drop=True)

        self.dff = dff

        _ = pdwritecsv(dff, self.fname_mut_summ, self.root_mutations, verbose=False)
        _ = pdwritecsv(df_mut, self.fname_mut_anal, self.root_mutations, verbose=False)

        return dff, df_mut

    def cbioportal_studies(self):
        url = "https://www.cbioportal.org/api/studies"

        res = requests.get(url, headers={"Accept": "application/json"})
        res.raise_for_status()

        studies = res.json()

        study_ids = [s["studyId"] for s in studies]

        return study_ids

    def loop_program_psi_samples(
        self, prog_id: str = "TCGA", force: bool = False, verbose: bool = True
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

        df_psi = self.get_primary_sites(prog_id=prog_id, force=force, verbose=verbose)

        df_cases, df_subt = pd.DataFrame(), pd.DataFrame()

        fname_all_cases = self.fname_all_cases % (self.prog_id)
        filename_cases = os.path.join(self.root_summary, fname_all_cases)

        fname_all_samples = self.fname_all_samples % (self.prog_id)
        filename_samples = os.path.join(self.root_summary, fname_all_samples)

        fname_all_mutations = self.fname_all_mutations % (self.prog_id)
        filename_mutations = os.path.join(self.root_summary, fname_all_mutations)

        if (
            os.path.exists(filename_cases)
            and os.path.exists(filename_samples)
            and os.path.exists(filename_mutations)
            and not force
        ):
            df_all_cases = pdreadcsv(fname_all_cases, self.root_summary)
            df_all_samples = pdreadcsv(fname_all_samples, self.root_summary)
            df_all_mutations = pdreadcsv(fname_all_mutations, self.root_summary)

            self.df_all_cases = df_all_cases
            self.df_all_samples = df_all_samples
            self.df_all_mutations = df_all_mutations

            return df_all_cases, df_all_samples, df_all_mutations

        df_list_cases, df_list_samples, df_list_mutations = [], [], []

        for ipsi, row in df_psi.iterrows():

            psi_id = row.psi_id
            primary_site = row.primary_site

            self.set_primary_site(psi_id)

            print(f"{ipsi}) {psi_id} -{primary_site}", end=" - ")

            df_cases, df_subt, _ = self.get_cases_and_subtypes(
                batch_size=200, do_filter=False, force=force, verbose=verbose
            )

            if df_cases.empty:
                print("No cases found for PSI_ID:", psi_id)
                continue

            if isinstance(df_cases, pd.DataFrame):
                df_list_cases.append(df_cases)
            else:
                print("Unexpected type for df_cases:", type(df_cases))
                raise Exception("Stope: unexpected type for df_cases")

            for isubt, row in df_subt.iterrows():
                subtype_global = row.subtype_global
                tumor_class = row.tumor_class
                subtype_tissue = row.subtype_tissue

                df_samples = self.get_samples_for_subtypes(
                    subtype_global=subtype_global,
                    tumor_class=tumor_class,
                    subtype_tissue=subtype_tissue,
                    batch_size=200,
                    force=force,
                    verbose=verbose,
                )
                print(f"{isubt}) {self.s_case}")

                if df_samples.empty:
                    if verbose:
                        print(
                            f"No samples found for PSI_ID: {psi_id} subtype: {subtype_global} tumor_class: {tumor_class} subtype_tissue: {subtype_tissue}"
                        )
                    continue

                df_list_samples.append(df_samples)

                df2 = df_samples[
                    ~df_samples.sample_type.str.contains("Blood", case=False, na=False)
                ]

                if df2.empty:
                    print("No samples having non-blood types for PSI_ID:", psi_id)
                    continue

                barcode_sample_list = list(np.unique(df2.barcode_sample))
                self.barcode_sample_list = barcode_sample_list

                print("Getting mutations", end=" ")
                dff, _ = self.get_df_mut_transform_mutation_table(
                    psi_id=psi_id,
                    barcode_sample_list=barcode_sample_list,
                    force=force,
                    verbose=verbose,
                )

                if dff.empty:
                    print("Could not find mutations for :", self.s_case)
                    continue

                df_list_mutations.append(dff)

        if len(df_list_cases) > 0:
            df_all_cases = pd.concat(df_list_cases, ignore_index=True)
            df_all_cases = df_all_cases.drop_duplicates()
            df_all_cases = df_all_cases.reset_index(drop=True)
        else:
            df_all_cases = pd.DataFrame()

        if len(df_list_samples) > 0:
            df_all_samples = pd.concat(df_list_samples, ignore_index=True)
            df_all_samples = df_all_samples.drop_duplicates()
            df_all_samples = df_all_samples.reset_index(drop=True)
        else:
            df_all_samples = pd.DataFrame()

        if len(df_list_mutations) > 0:
            df_all_mutations = pd.concat(df_list_mutations, ignore_index=True)
            df_all_mutations = df_all_mutations.drop_duplicates()
            df_all_mutations = df_all_mutations.reset_index(drop=True)
        else:
            df_all_mutations = pd.DataFrame()

        _ = pdwritecsv(df_all_cases, fname_all_cases, self.root_summary)
        _ = pdwritecsv(df_all_samples, fname_all_samples, self.root_summary)
        _ = pdwritecsv(df_all_mutations, fname_all_mutations, self.root_summary)

        self.df_all_cases = df_all_cases
        self.df_all_samples = df_all_samples
        self.df_all_mutations = df_all_mutations

        return df_all_cases, df_all_samples, df_all_mutations

    def get_filtered_tables(
        self, 
        sample_type_term: str = "Primary Tumor", 
        verbose: bool = False
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:

        df_cases, df_all_samples, df_all_mut, all_barcode_list = self.get_filtered_tables_subtypes(
            sample_type_term=sample_type_term, do_filter=True, verbose=verbose
        )

        if verbose:
            if df_cases.empty:
                print("No cases found for primary site:", self.primary_site)

            if df_all_samples.empty:
                print("No samples found for primary site:", self.primary_site)

            if df_all_mut.empty:
                print("No mutations found for primary site:", self.primary_site)

            if len(all_barcode_list) == 0:
                print("No barcodes found for primary site:", self.primary_site)

        return df_cases, df_all_samples, df_all_mut, all_barcode_list

    def get_filtered_tables_subtypes(
        self, sample_type_term: str = "Primary Tumor", do_filter: bool = True, verbose: bool = True
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:

        self.df_cases = pd.DataFrame()
        self.df_all_samples = pd.DataFrame()
        self.df_all_mut = pd.DataFrame()
        self.all_barcode_list = []

        self.set_filenames()

        if not os.path.exists(self.filename_cases):
            print("Error: could not find cases file:", self.filename_cases)
            return self.df_cases, self.df_all_samples, self.df_all_mut, self.all_barcode_list

        df_cases = pdreadcsv(self.fname_cases, self.root_disease, verbose=verbose)
        if "pid" in df_cases.columns:
            df_cases = df_cases.rename(columns={"pid": "psi_id"})
            pdwritecsv(df_cases, self.fname_cases, self.root_disease, verbose=verbose)

        if do_filter:
            df_cases = self.apply_filter_cases(df_cases)

        self.df_cases = df_cases

        df_subt = self.groupby_case_by_subtypes(df_cases)

        df_list_samples = []
        df_list_mut = []
        list_all_barcodes = []

        for _, row in df_subt.iterrows():
            subtype_global = row.subtype_global
            tumor_class = row.tumor_class
            subtype_tissue = row.subtype_tissue

            self.set_s_case(subtype_global, tumor_class, subtype_tissue)

            df3 = df_cases[
                (df_cases.subtype_global == subtype_global)
                & (df_cases.tumor_class == tumor_class)
                & (df_cases.subtype_tissue == subtype_tissue)
            ]

            if df3.empty:
                print(f"No cases found for {subtype_global} {tumor_class} {subtype_tissue}")
                continue

            df3 = df3.copy().reset_index(drop=True)

            case_id_list = np.unique(df3.case_id)

            fname = self.fname_samples0 % (self.s_case)
            fname = title_replace(fname)
            filename = self.root_samples / fname

            if not filename.exists():
                print("Error: could not find samples file:", filename)
                continue

            df_samples = pdreadcsv(fname, self.root_samples, verbose=verbose)

            if df_samples.empty:
                print("Error: could not read samples file:", filename)
                continue

            df_samples = df_samples[
                (df_samples.case_id.isin(case_id_list))
                & (df_samples.sample_type.str.contains(sample_type_term, case=False))
            ]
            self.df_samples = df_samples

            if df_samples.empty:
                if verbose:
                    print(f"Warning: could not filter df_samples for {subtype_global} {tumor_class} {subtype_tissue}")
                continue

            df_samples = df_samples.copy().reset_index(drop=True)
            self.barcode_list = self.prepare_barcode_sample_list(df_samples.barcode_sample.tolist())

            self.set_mutation_filenames()

            if os.path.exists(self.filename_mutsumm):
                df_mut = pdreadcsv(self.fname_mut_anal, self.root_mutations, verbose=verbose)
            else:
                if verbose:
                    print("No mutation analysis file found for:", self.s_case)
                df_mut = pd.DataFrame()

            df_list_samples.append(df_samples)
            if not df_mut.empty:
                df_list_mut.append(df_mut)
            list_all_barcodes += list(self.barcode_list)

        df_all_samples = (
            pd.concat(df_list_samples, ignore_index=True) if df_list_samples else pd.DataFrame()
        )

        if self.memory_restriction:
            df_all_samples = df_all_samples.head(200).copy()
        
        df_all_samples.reset_index(drop=True, inplace=True)

        cases_ids = np.unique(df_all_samples.case_id)
        df_cases = df_cases[df_cases["case_id"].isin(cases_ids)].copy()
        df_cases.reset_index(drop=True, inplace=True)

        self.df_all_samples = df_all_samples
        self.df_cases = df_cases
                
        if len(df_list_mut) == 0:
            df_all_mut = pd.DataFrame()
            list_all_barcodes = []
        else:
            df_all_mut = pd.concat(df_list_mut, ignore_index=True) if df_list_mut else pd.DataFrame()
            list_all_barcodes = list(np.unique(list_all_barcodes))

        return self.df_cases, df_all_samples, df_all_mut, list_all_barcodes

    def prepare_barcode_sample_list(self, barcode_sample_list: list[str]) -> list[str]:
        """
        if TCGA remove the last characters if len > 2
        TCGA-OR-A5J2-01A -> TCGA-OR-A5J2-01
        """

        barcode_sample_list = [self.to_cbioportal_barcode_sample(x) for x in barcode_sample_list]
        if not barcode_sample_list:
            raise ValueError("barcode_sample_list is empty.")

        # 01A, 01B, 01Z → all collapse to 01
        barcode_sample_list = list(np.unique(barcode_sample_list))

        return barcode_sample_list

    def to_cbioportal_barcode_sample(self, x: str) -> str:
        parts = x.split("-")

        if len(parts) >= 4 and parts[0] == "TCGA":
            sample_code = parts[3][:2]  # 01A -> 01, 11A -> 11
            return "-".join([parts[0], parts[1], parts[2], sample_code])

        return x

    def build_pivot_table(
        self, df_all_mut: pd.DataFrame, min_barcodes: int = 2, min_genes: int = 2
    ) -> pd.DataFrame:
        """
        Build a barcode x gene binary mutation matrix (0/1).

        Rows represent barcodes (samples), columns represent gene symbols.
        A value of 1 indicates that at least one mutation was observed for that
        barcode-gene pair.

        Parameters
        ----------
        df_all_mut : pd.DataFrame
                Input mutation table. Must contain at least:
                - 'barcode': sample identifier
                - 'symbol': gene symbol

        Returns
        -------
        pd.DataFrame
                Binary mutation matrix with:
                - index   = barcode
                - columns = gene symbol
                - values  = 0 or 1

                Returns an empty DataFrame if the input is empty or required columns
                are missing.
        """
        if df_all_mut is None or df_all_mut.empty:
            return pd.DataFrame()

        required_cols = {"barcode", "symbol"}
        missing_cols = required_cols - set(df_all_mut.columns)
        if missing_cols:
            raise ValueError(
                f"build_pivot_table requires columns {sorted(required_cols)}, "
                f"but is missing {sorted(missing_cols)}."
            )

        # Keep only the columns needed for the mutation matrix
        dfa = df_all_mut.loc[:, ["barcode", "symbol"]].copy()
        dfa["barcode"] = dfa["barcode"].astype(str).str.strip()
        dfa["symbol"] = dfa["symbol"].astype(str).str.strip()

        dfa = dfa[
            (dfa["barcode"] != "")
            & (dfa["symbol"] != "")
            & (dfa["barcode"].str.lower() != "nan")
            & (dfa["symbol"].str.lower() != "nan")
        ]

        if dfa.empty:
            return pd.DataFrame()

        # Mark presence of at least one mutation per barcode-gene pair
        dfa["mutated"] = 1

        dfpiv = dfa.pivot_table(
            index="barcode",
            columns="symbol",
            values="mutated",
            aggfunc="max",
            fill_value=0,
        ).astype(np.uint8)

        # Remove empty samples and genes
        dfpiv = dfpiv.loc[dfpiv.sum(axis=1) >= min_genes, :]

        dfpiv = dfpiv.loc[:, dfpiv.sum(axis=0) >= min_barcodes]

        if dfpiv.shape[0] < 3:
            print("dfpiv has less than 3 samples.")
            return pd.DataFrame()

        if dfpiv.shape[1] < 3:
            print("dfpiv has less than 3 genes.")
            return pd.DataFrame()

        """
		It sorts your matrix in a consistent order:

		axis=0 → sort rows (barcodes)
		axis=1 → sort columns (genes)

		So after this:
			barcodes are alphabetically (or lexicographically) ordered
			genes are alphabetically ordered
		"""
        dfpiv = dfpiv.sort_index(axis=0).sort_index(axis=1)

        return dfpiv

    def calc_HDBSCAN(
        self, dfpiv: pd.DataFrame, min_cluster_size: int = 8, min_samples: int = 3
    ) -> tuple[List, List, Any]:
        """
        Cluster with HDBSCAN, not KMeans
        pairwise_distances with jaccard
        Multidimensional Scaling (MDS)
        If there are a few dense groups plus many ambiguous samples, HDBSCAN can work better.

        In HDBSCAN:
                min_cluster_size → minimum size of a cluster
                min_samples → minimum local neighborhood density
                A point is considered “core” if it has at least min_samples neighbors nearby.

        input: dfpiv, min_cluster_size (minimum number of samples in a cluster)
        output: embedding and labels
        """

        X = dfpiv.to_numpy(dtype=bool)

        n_samples = X.shape[0]
        n_genes = X.shape[1]

        if n_samples < 3:
            print("Need at least 3 non-empty samples to compute HDBSCAN + clustering.")
            return [], [], None

        if min_cluster_size > n_samples:
            print(
                f"min_cluster_size={min_cluster_size} is larger than number of samples ({X.shape[0]}). Using min_cluster_size={X.shape[0]}."
            )
            min_cluster_size = max(2, n_samples - 1)

        if n_genes < 3:
            print("Need at least 3 non-empty genes to compute HDBSCAN + clustering.")
            return [], [], None

        min_cluster_size = min(min_cluster_size, n_samples)

        # print(">>> calc_HDBSCAN MIN_CLUSTER_SIZE", min_cluster_size)

        # where distance_matrix is already your precomputed Jaccard dissimilarity matrix.
        dist_matrix = pairwise_distances(X, metric="jaccard")

        # metric is not exactly the same concept as dissimilarity="precomputed" in older versions.
        # In current MDS, the distance matrix behavior is usually controlled by normalized_stress/input handling 
        # and the API changes around dissimilarity can be confusing.
        embedding = MDS(
            n_components=2,
            metric='precomputed',
            n_init=1,
            random_state=42,
        ).fit_transform(dist_matrix)

        if isinstance(embedding, tuple):
            print("embedding return as a tuple")
            embedding = embedding[0]

        embedding = np.asarray(embedding)

        if embedding.shape[0] < 3:
            print("Too few valid embedded samples after filtering.")
            return [], [], None

        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples if min_samples is not None else min_cluster_size,
            metric="euclidean",
        )
        labels = clusterer.fit_predict(embedding)

        return embedding.tolist(), labels.tolist(), dist_matrix

    def calc_UMAP(self, dfpiv: pd.DataFrame, k: int = 8) -> tuple[List, List]:
        # Binary mutation matrix for Jaccard
        # Force numeric/binary and remove bad values

        X = dfpiv.to_numpy(dtype=np.uint8)

        n_samples = X.shape[0]
        n_genes = X.shape[1]

        if n_samples < 3:
            print("Need at least 3 non-empty samples to compute UMAP + clustering.")
            return [], []

        if k > n_samples:
            print(f"k={k} is larger than number of samples ({X.shape[0]}). Using k={X.shape[0]}.")
            k = max(2, n_samples - 1)

        if n_genes < 3:
            print("Need at least 3 non-empty genes to compute UMAP + clustering.")
            return [], []

        k = min(k, n_samples)

        n_neighbors = min(15, n_samples - 1)
        n_neighbors = max(2, n_neighbors)

        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=0.1,
            metric="jaccard",
            random_state=42,
            init="random",  # important
            output_dens=False,  # avoid tuple output
        )

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="umap")
            embedding = reducer.fit_transform(X)

        if isinstance(embedding, tuple):
            print("embedding return as a tuple")
            embedding = embedding[0]

        embedding = np.asarray(embedding)

        """
		good = np.isfinite(embedding).all(axis=1)
		embedding = embedding[good]
		"""
        if not np.isfinite(embedding).all():
            print("UMAP embedding contains NaN or infinite values.")
            return [], []

        if embedding.shape[0] < 3:
            print("Too few valid embedded samples after filtering.")
            return [], []

        labels = KMeans(
            n_clusters=min(k, embedding.shape[0] - 1),
            random_state=42,
            n_init=10,
        ).fit_predict(embedding)

        return embedding.tolist(), labels.tolist()

    def plot_UMAP(
        self, dfpiv: pd.DataFrame, k: int = 8, figsize: tuple = (14, 10)
    ) -> Tuple[Any, Any, Any]:

        n_samples = dfpiv.shape[0]
        n_genes = dfpiv.shape[1]

        embedding, labels = self.calc_UMAP(dfpiv, k)
        embedding = np.array(embedding)

        if len(embedding) == 0 or len(labels) == 0:
            print("No valid UMAP embedding or labels.")
            return None, None, None

        fig, ax = plt.subplots(figsize=figsize)

        # cmap = plt.cm.get_cmap("tab10", k)
        plt.scatter(
            embedding[:, 0], embedding[:, 1], c=[self.colors[label] for label in labels], s=20
        )

        ax.set_title(
            f"Clustering using UMAP mutation profiles: (k={k})\nPrimary Site: '{self.primary_site}' #{n_samples} samples and #{n_genes} genes"
        )
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

        counts = Counter(labels)

        legend_handles = []

        for cluster_id in sorted(counts.keys()):
            color = self.colors[cluster_id]
            n = counts[cluster_id]

            patch = mpatches.Patch(color=color, label=f"Cluster {cluster_id} (n={n})")
            legend_handles.append(patch)

        ax.legend(handles=legend_handles, title="Groups", loc="best")

        return fig, embedding, labels

    def plot_HDBSCAN(
        self,
        dfpiv: pd.DataFrame,
        min_cluster_size: int = 8,
        min_samples: int = 3,
        figsize: tuple = (14, 10),
    ) -> Tuple[Any, Any, Any, Any]:

        n_samples = dfpiv.shape[0]
        n_genes = dfpiv.shape[1]

        embedding, labels, d = self.calc_HDBSCAN(
            dfpiv=dfpiv, min_cluster_size=min_cluster_size, min_samples=min_samples
        )
        embedding = np.array(embedding)

        if len(embedding) == 0 or len(labels) == 0:
            print("No valid HDBSCAN embedding or labels.")
            return None, None, None, None

        fig, ax = plt.subplots(figsize=figsize)

        plt.scatter(
            embedding[:, 0], embedding[:, 1], c=[self.colors[label] for label in labels], s=20
        )

        stri = "Clustering using HDBSCAN mutation profiles"
        stri += f"\nPrimary Site: {self.psi_id} - '{self.primary_site}' #{n_samples} samples and #{n_genes} genes"
        stri += f"\nmin_cluster_size={min_cluster_size} and min_samples={min_samples}"
        ax.set_title(stri)
        ax.set_xlabel("embedding1")
        ax.set_ylabel("embedding2")

        counts = Counter(labels)

        legend_handles = []

        for cluster_id in sorted(counts.keys()):
            color = self.colors[cluster_id]
            n = counts[cluster_id]

            patch = mpatches.Patch(color=color, label=f"Cluster {cluster_id} (n={n})")
            legend_handles.append(patch)

        ax.legend(handles=legend_handles, title="Groups", loc="best")

        return fig, embedding, labels, d

    def cluster_mutation_table(
        self, dfpiv: pd.DataFrame, labels, cluster: int = 1, min_barcodes: int = 2
    ) -> pd.DataFrame:

        if len(labels) != dfpiv.shape[0]:
            stri = "Error: build_pivot_table filter empty lines."
            stri += "\nNumber of labels does not match number of samples. "
            stri += "\n----------- stop execution -----------\n"
            raise Exception(stri)

        labels = pd.Series(labels, index=dfpiv.index)
        sel_barcodes = labels[labels == cluster].index
        dff = dfpiv.loc[sel_barcodes].copy()
        dff = dff.loc[:, dff.sum(axis=0) >= min_barcodes]

        return dff

    def calc_shannon_entropy_from_dfstat(self, dfstat: pd.DataFrame) -> pd.DataFrame:
        rows = []

        for (k, cluster), dfsub in dfstat.groupby(["k", "cluster"]):
            deg_list = dfsub["degree"].to_list()

            if len(deg_list) == 0:
                H = np.nan
                Hmax = np.nan
                Hnorm = np.nan
                n_genes = 0
            else:
                w = np.array(deg_list, dtype=float)
                p = w / w.sum()
                H = -np.sum(p * np.log2(p))
                n_genes = len(p)
                Hmax = np.log2(n_genes) if n_genes > 1 else 0.0
                Hnorm = H / Hmax if Hmax > 0 else 0.0

            rows.append(
                {
                    "k": k,
                    "cluster": cluster,
                    "n_genes": n_genes,
                    "cluster_size": dfsub["cluster_size"].iloc[0],
                    "entropy": H,
                    "entropy_max": Hmax,
                    "entropy_norm": Hnorm,
                }
            )

        if rows == []:
            return pd.DataFrame()

        dfh = pd.DataFrame(rows)
        dfh = dfh.sort_values("entropy_norm", ascending=True)

        return dfh

    def score_k_from_entropy_table(self, dfh: pd.DataFrame) -> pd.DataFrame:
        rows = []

        for k, sub in dfh.groupby("k"):
            total_n = sub["cluster_size"].sum()

            weighted_mean_entropy = (
                (sub["entropy_norm"] * sub["cluster_size"]).sum() / total_n
                if total_n > 0
                else np.nan
            )

            mean_entropy = sub["entropy_norm"].mean()
            std_entropy = sub["entropy_norm"].std()
            min_cluster_size = sub["cluster_size"].min()
            max_cluster_size = sub["cluster_size"].max()
            n_clusters = sub.shape[0]
            n_small_clusters = (sub["cluster_size"] < 3).sum()

            rows.append(
                {
                    "k": k,
                    "n_clusters": n_clusters,
                    "weighted_mean_hnorm": weighted_mean_entropy,
                    "mean_hnorm": mean_entropy,
                    "std_hnorm": std_entropy,
                    "min_cluster_size": min_cluster_size,
                    "max_cluster_size": max_cluster_size,
                    "n_small_clusters_lt3": n_small_clusters,
                }
            )

        return pd.DataFrame(rows).sort_values("weighted_mean_hnorm", ascending=True)

    def entropy_analysis_for_primary_site(
        self,
        cluster_type: str,
        sample_type_term: str = "Primary Tumor",
        Kmin: int = 2,
        Kmax: int = 10,
        min_barcodes: int = 2,
        min_genes: int = 2,
        verbose: bool = False,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:

        _, _, df_all_mut, _ = self.get_filtered_tables(
            sample_type_term=sample_type_term, verbose=verbose
        )

        dfempty = pd.DataFrame()

        if df_all_mut.empty:
            return dfempty, dfempty, dfempty, dfempty, dfempty

        dfpiv = self.build_pivot_table(df_all_mut, min_barcodes=min_barcodes, min_genes=min_genes)
        self.dfpiv = dfpiv

        if dfpiv.shape[0] < 3 or dfpiv.shape[1] < 3:
            return dfempty, dfempty, dfempty, dfpiv, df_all_mut

        if Kmax >= dfpiv.shape[0]:
            Kmax = dfpiv.shape[0] - 1

        if Kmax <= 3:
            return dfempty, dfempty, dfempty, dfpiv, df_all_mut

        df_list = []
        for k in range(Kmin, Kmax + 1):
            if cluster_type == "UMAP":
                _, labels = self.calc_UMAP(dfpiv, k)
            elif cluster_type == "HDBSCAN":
                _, labels, _ = self.calc_HDBSCAN(dfpiv, k)
            else:
                raise Exception(
                    f"\n---------- Define the cluster_type like UMAP or HDBSCAN, got: {cluster_type}"
                )

            if labels is None or len(labels) == 0:
                continue

            for cluster in np.unique(labels):
                dfc = self.cluster_mutation_table(
                    dfpiv=dfpiv, labels=labels, cluster=cluster, min_barcodes=min_barcodes
                )

                n_cluster = dfc.shape[0]

                if n_cluster < 3 or dfc.shape[1] < 3:
                    continue

                gene_degree = dfc.sum(axis=0).sort_values(ascending=False)
                gene_freq = (gene_degree / n_cluster).sort_values(ascending=False)

                df_cluster_stat = pd.DataFrame(
                    {
                        "k": k,
                        "cluster": cluster,
                        "gene": gene_degree.index,
                        "degree": gene_degree.values,
                        "cluster_size": n_cluster,
                        "freq": gene_freq.values,
                    }
                )
                df_list.append(df_cluster_stat)

        if len(df_list) == 0:
            return dfempty, dfempty, dfempty, dfpiv, df_all_mut

        dfstat = pd.concat(df_list, ignore_index=True)

        dfh = self.calc_shannon_entropy_from_dfstat(dfstat)

        dfw = self.score_k_from_entropy_table(dfh)

        cols = dfw.columns.to_list()

        dfw["psi_id"] = self.psi_id
        dfw["primary_site"] = self.primary_site
        dfw["min_barcodes"] = min_barcodes
        dfw["min_genes"] = min_genes

        cols = ["psi_id", "primary_site", "min_barcodes", "min_genes"] + cols
        dfw = dfw[cols]

        return dfw, dfh, dfstat, dfpiv, df_all_mut

    def buid_purity_table(self, dfpiv: pd.DataFrame, labels: list) -> pd.DataFrame:
        lab_list = np.unique(labels)
        pu_list = []
        n_list = []
        pairs = []

        for idx, label in enumerate(lab_list):
            # idx = labels == label
            X = dfpiv[idx].to_numpy(dtype=bool)

            n_bardodes = len(X)
            n_list.append(n_bardodes)
            pairs.append(n_bardodes * (n_bardodes - 1) / 2)

            if len(X) > 1:
                dist = pairwise_distances(X, metric="jaccard")
                similarity = 1 - dist
                purity = similarity[np.triu_indices_from(similarity, k=1)].mean()
            else:
                purity = 0

            pu_list.append(np.round(purity, 3))

        dfa = pd.DataFrame(
            {"label": lab_list, "n_barcodes": n_list, "purity": pu_list, "n_pairs": pairs}
        )

        max_pairs = dfa["n_pairs"].max()
        dfa["purity_norm"] = dfa["purity"] * (dfa["n_pairs"] / max_pairs)

        dfa = dfa.sort_values(by="purity_norm", ascending=False)

        return dfa

    def plot_purity(
        self, dfpur: pd.DataFrame, dfclu: pd.DataFrame, good_clusters: list, min_perc: float = 0.10
    ):

        # top_n_genes = 30

        ngood = len(good_clusters)
        nrows = int(np.ceil(ngood / 2))
        height = 6 * nrows

        fig, axes = plt.subplots(nrows, 2, figsize=(12, height))
        axes = axes.flatten()

        for ax, label in zip(axes, good_clusters):
            if label == -1:
                continue

            dfb = dfclu[label]

            dfb = dfb[dfb.values > min_perc]
            dfb = dfb.sort_values(ascending=False)
            # dfb = dfb.sort_values(ascending=False).head(top_n_genes)

            ax.bar(dfb.index, dfb.values)

            stri = f"Label {label} | purity_norm={dfpur.loc[dfpur['label'] == label].iloc[0].purity_norm:.3f}"

            ax.set_title(stri)
            ax.set_ylabel("Representative percentage")
            ax.set_xlabel("Genes")
            ax.tick_params(axis="x", rotation=70)

            print(stri, ", ".join(dfb.index.to_list()))

        plt.tight_layout()

        return fig

    def enrichment_test(self, sample_genes, subtype_genes, background_genes):
        N = len(background_genes)
        K = len(subtype_genes)
        n = len(sample_genes)
        overlap_genes = set(sample_genes) & set(subtype_genes)
        overlap = len(overlap_genes)

        # P(X >= k)
        pval = hypergeom.sf(overlap - 1, N, K, n)

        return pval, overlap_genes

    def cluster_analysis(
        self,
        cluster_type: str,
        sample_type_term: str,
        k: int = 5,
        Kmin: int = 2,
        Kmax: int = 10,
        min_barcodes: int = 3,
        min_genes: int = 5,
        pur_threshold: float = 0.05,
        min_represent_perc=0.10,
        verbose: bool = False,
    ) -> Tuple[
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
    ]:
        

        dfw, dfh, dfstat, dfpiv, df_all_mut = self.entropy_analysis_for_primary_site(
            cluster_type, sample_type_term, Kmin, Kmax, min_barcodes, min_genes, verbose
        )

        dfempty = pd.DataFrame()

        if dfpiv.empty:
            print("Did not define the pivot table")
            return dfempty, dfempty, dfempty, dfempty, dfempty, dfempty, dfempty, dfempty

        if cluster_type == "UMAP":
            print(f"Chose {cluster_type} with k={k}")
            _, labels = self.calc_UMAP(dfpiv, k)

        elif cluster_type == "HDBSCAN":
            min_cluster_size = 5
            min_samples = 3
            print(
                f"Chose {cluster_type} with min_cluster_size={min_cluster_size} and min_samples={min_samples}"
            )

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _, labels, _ = self.calc_HDBSCAN(
                    dfpiv, min_cluster_size=min_cluster_size, min_samples=min_samples
                )
        else:
            print("Did not defined the clustering method")
            return dfempty, dfempty, dfempty, dfw, dfh, dfstat, dfpiv, df_all_mut

        # ----------- cluster ----------------
        dfpiv2 = dfpiv.copy()
        dfpiv2["cluster"] = labels
        dfclu = dfpiv2.groupby("cluster").mean().T

        dfpur = self.buid_purity_table(dfpiv, labels)
        if dfpur.empty:
            return dfempty, dfempty, dfempty, dfw, dfh, dfstat, dfpiv, df_all_mut

        good_clusters = dfpur.loc[dfpur["purity_norm"] >= pur_threshold, "label"]

        # ------------ hypergeometric statistics ---------------------
        background_genes = np.unique(df_all_mut.symbol.to_list())

        dic = self.SUBTYPE_GENES.get(self.psi_id, {})

        if dic == {}:
            print(f"No subtype genes found for PSI ID: {self.psi_id}")
            return dfempty, dfpur, dfclu, dfw, dfh, dfstat, dfpiv, df_all_mut

        lista = []
        for label in good_clusters:
            purity_norm = dfpur.loc[dfpur["label"] == label].iloc[0].purity_norm

            dfb = dfclu[label]

            n_barcodes = len(dfpiv2[dfpiv2["cluster"] == label])

            dfb = dfb[dfb.values > min_represent_perc]
            # dfb = dfb.sort_values(ascending=False)

            sample_genes = dfb.index.to_list()

            for subtype, annotated_genes in dic.items():
                pval, overlap_genes = self.enrichment_test(
                    sample_genes, annotated_genes, background_genes
                )
                # print(f"Subtype: {subtype}, overlap: {overlap}, p-value: {pval}")

                overlap = len(overlap_genes)

                if overlap >= 2:
                    mat = [
                        cluster_type,
                        k,
                        n_barcodes,
                        label,
                        purity_norm,
                        subtype,
                        overlap,
                        len(sample_genes),
                        len(annotated_genes),
                        len(background_genes),
                        pval,
                        overlap_genes,
                    ]
                    lista.append(mat)

        if lista == []:
            df = pd.DataFrame()
        else:
            df = pd.DataFrame(
                lista,
                columns=[
                    "cluster_type",
                    "k",
                    "n_barcodes",
                    "label",
                    "purity_norm",
                    "subtype",
                    "overlap",
                    "sample_genes",
                    "annotated_genes",
                    "background_genes",
                    "pval",
                    "overlap_genes",
                ],
            )
            df["fdr"] = fdr(df["pval"])

        return df, dfpur, dfclu, dfw, dfh, dfstat, dfpiv, df_all_mut

    def get_VCF_file(self, case_id: str, sample_id_list: list, subtype_global: str, tumor_class: str, subtype_tissue: str, 
                     timeout: int = 100, force: bool = False, verbose: bool = False,) -> pd.DataFrame:

        df_vcf = pd.DataFrame()
        self.df_vcf = df_vcf


        self.set_s_case(subtype_global, tumor_class, subtype_tissue)

        fname = self.fname_vcf_files0 % (case_id, sample_id_list[0], sample_id_list[-1]) 
        fname = title_replace(fname)
        filename = self.root_mutations / fname

        if filename.exists() and not force:
            df_vcf = pdreadcsv(fname, self.root_mutations, verbose=verbose)
            return df_vcf

        df_samples = self.get_samples_for_subtypes(
                        subtype_global=subtype_global,
                        tumor_class=tumor_class,
                        subtype_tissue=subtype_tissue,
                        batch_size=200,
                        force=False,
                        verbose=verbose,
                    )

        if df_samples.empty:
            print(f"No cases found while searching for '{self.s_case}'")
            return df_vcf

        df_samples = df_samples[df_samples.case_id == case_id]

        if df_samples.empty:
            print(f"No cases found while searching for case_id = '{case_id}'")
            return df_vcf

        sample_id_list = df_samples.sample_id.to_list()
        N = len(df_samples)

        # -------------------------- batch loop ---------------------------

        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "cases.samples.sample_id",
                        "value": sample_id_list,
                    },
                },
                {
                    "op": "=",
                    "content": {
                        "field": "files.data_format",
                        "value": "VCF",
                    },
                },
                {
                    "op": "in",
                    "content": {
                        "field": "files.data_type",
                        "value": [
                            "Raw Simple Somatic Mutation",
                            "Masked Somatic Mutation",
                        ],
                    },
                },
            ],
        }

        fields = [
            "file_id",
            "file_name",
            "data_type",
            "data_format",
            "access",
            "analysis.workflow_type",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
        ]

        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(fields),
            "format": "JSON",
            "size": N,
        }


        res = requests.get(self.url_gdc_files, params=params, timeout=timeout)
        res.raise_for_status()

        response = res.json()

        if "data" not in response:
            print("No data returned from GDC.")
            print(">>> response", response)
            return pd.DataFrame()

        hits = response.get("data", {}).get("hits", [])

        if not hits:
            print(f"No VCF files found for {len(sample_id_list)} sample_id(s) for case_id {case_id}.")
            return pd.DataFrame()


        records = []

        for hit in hits:
            file_id = hit.get("file_id")
            file_name = hit.get("file_name")
            data_type = hit.get("data_type")
            data_format = hit.get("data_format")

            analysis = hit.get("analysis") or {}
            workflow_type = analysis.get("workflow_type")

            for case in hit.get("cases", []):
                case_submitter_id = case.get("submitter_id")

                for sample in case.get("samples", []):
                    records.append(
                        {
                            # case metadata
                            "case_id": case_id,
                            "case_submitter_id": case_submitter_id,

                            # sample metadata
                            "sample_id": sample.get("sample_id"),
                            "barcode_sample": sample.get("submitter_id"),
                            "sample_type": sample.get("sample_type"),

                            # file metadata
                            "file_id": file_id,
                            "file_name": file_name,
                            "data_type": data_type,
                            "data_format": data_format,
                            "workflow_type": workflow_type,

                            # project/subtype metadata
                            "psi_id": self.psi_id,
                            "subtype_global": subtype_global,
                            "tumor_class": tumor_class,
                            "subtype_tissue": subtype_tissue,
                        }
                    )

        df_vcf = pd.DataFrame(records)

        if df_vcf.empty:
            return df_vcf

        cols = [
            "psi_id",
            "case_id",
            "case_submitter_id",
            "sample_id",
            "barcode_sample",
            "sample_type",
            "subtype_global",
            "tumor_class",
            "subtype_tissue",
            "file_id",
            "file_name",
            "data_type",
            "data_format",
            "workflow_type",
        ]

        cols = [c for c in cols if c in df_vcf.columns]
        df_vcf = df_vcf[cols]
        self.df_vcf = df_vcf

        return df_vcf



    def open_df_VCF(self, case_id, sample_id_list, verbose: bool = False) -> pd.DataFrame:

        fname = self.fname_vcf_files0 % (case_id, sample_id_list[0], sample_id_list[-1])
        fname = title_replace(fname)
        filename = self.root_mutations / fname

        if filename.exists():
            df_vcf = pdreadcsv(fname, self.root_mutations, verbose=verbose)
        else:
            print(f"VCF file not found for case_id: {case_id}, sample_id: {sample_id_list[0]} to {sample_id_list[-1]}")
            df_vcf = pd.DataFrame()

        self.df_vcf = df_vcf
        return df_vcf

    def prepare_gtex(self, df_gtex_ctrl: pd.DataFrame) -> pd.DataFrame:

        df_gtex_ctrl.rename(columns={"ensemblid": "geneid"}, inplace=True)
        df_gtex_ctrl['geneid'] = [gene.split('.')[0] for gene in df_gtex_ctrl.geneid]

        cols = list(df_gtex_ctrl.columns)[2:]

        df_gtex_ctrl["biotype"] = "protein_coding"
        df_gtex_ctrl = df_gtex_ctrl[self.GENE_COLS + cols]

        _ = [
            df_gtex_ctrl.rename(columns={cols[i]: f"normal_{i + 1}"}, inplace=True)
            for i in range(len(cols))
        ]

        return df_gtex_ctrl


    def get_tumor_normal_tables(self, imax_samples: int = 200, force: bool = False,
                                verbose: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, str]:
        
        df_tumor, df_normal, df_gtex_ctrl = \
            self.calc_file_expression_tumor_normal_gtex(imax_samples=imax_samples, force=force, verbose=verbose)

        if df_tumor.empty:
            msg = f"No tumor expression data found for {self.psi_id}"
            if verbose:
                print(msg)
            return pd.DataFrame(),pd.DataFrame(), msg

        # geneid, symbol, biotype, samples
        min_N_cols = 3 + 3
        # df_gtex_ctrl - has no biotype
        if df_normal.shape[1] < min_N_cols and df_gtex_ctrl.shape[1] < (min_N_cols-1):
            msg = "Error: Normal samples and GTEx control do not have enough samples."
            if verbose:
                print(msg)
            return pd.DataFrame(),pd.DataFrame(), msg

        df_normal = self.cdegs.deduplicate_by_max_reads(df_normal)

        if df_normal.empty or df_normal.shape[1] < min_N_cols:
            df_normal2 = self.prepare_gtex(df_gtex_ctrl)
            df_normal2 = self.cdegs.deduplicate_by_max_reads(df_normal2)
            msg = f"not enough normal samples --> substituting with GTEx control {df_normal2.shape[1] - 3}."
            if verbose:
                print(msg)
        else:
            msg = f"enough normal samples {df_normal.shape[1] - 3}."
            if verbose:
                print(msg)
            df_normal2 = df_normal

        self.df_normal2 = df_normal2

        msg += f"\nThere are {df_tumor.shape[1] - 3} tumor samples; {msg}"

        df_tumor = self.cdegs.deduplicate_by_max_reads(df_tumor)
        self.df_tumor = df_tumor

        if df_tumor.empty or df_tumor.shape[1] < min_N_cols:
            msg = "Error: Tumor expression data has fewer than 3 samples."
            return pd.DataFrame(),pd.DataFrame(), msg

        if df_normal2.empty or df_normal2.shape[1] < min_N_cols:
            msg = "Error: Normal expression data has fewer than 3 samples."
            return pd.DataFrame(),pd.DataFrame(), msg
        
        return df_tumor, df_normal2, msg


    def calc_lfc_table(
        self,
        psi_id: str,
        run_conda: bool = False,
        method: str = "edger",
        verbose: bool = False,
    ) -> tuple[pd.DataFrame, str]:

        self.set_primary_site(psi_id=psi_id)

        cdegs = CALC_DEGS(root_src=self.root_src, run_conda=run_conda)
        self.cdegs = cdegs

        df_tumor, df_normal, msg = self.get_tumor_normal_tables(force=False, verbose=verbose)

        if df_tumor.empty:
            return pd.DataFrame(), msg

        if df_normal.empty:
            return pd.DataFrame(), msg

        df_lfc = cdegs.run_deg_rscript(
            df_tumor=df_tumor,
            df_normal=df_normal,
            method=method,
            manual_dispersion=0.1,
            min_total_count=10,
            merge_how="inner",
            keep_temp=False,
        )

        cols = df_lfc.columns.to_list()

        commons =  ["geneid"]
        tum_cols = ["geneid", "symbol", "biotype"]

        self.df_lfc5 = df_lfc.copy()
        self.df_tumor = df_tumor

        # biotype can be loose: biotypes come from df_tumor
        df_lfc = pd.merge(df_lfc, df_tumor[tum_cols], on=commons, how="inner")

        cols2 = ["geneid", "symbol", "biotype"] + cols[1:]
        df_lfc = df_lfc[cols2]

        df_lfc = df_lfc.rename(columns={"geneid": "ensembl_id", "log2FoldChange": "lfc", "pvalue": "pval", "padj": "fdr"})

        return df_lfc, msg


    def calc_degs(
        self,
        psi_id: str,
        root_src: Path = Path('.'),
        run_conda: bool = False,
        lfc_cutoff: float = 1.0,
        fdr_cutoff: float = 0.05,
        method: str = "edger",
        imax_samples: int = 200,
        force: bool = False,
        verbose: bool = False,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, str, str]:

        self.set_primary_site(psi_id=psi_id)

        df_tumor, df_normal, df_gtex_ctrl = self.calc_file_expression_tumor_normal_gtex(imax_samples=imax_samples, verbose=verbose)

        if df_tumor.empty:
            if verbose:
                print(f"No tumor expression data found for {self.psi_id}")
            return pd.DataFrame(), pd.DataFrame(), "", ""

        fname_degs = self.fname_degs % self.psi_id
        fname_lfc = self.fname_lfc % self.psi_id
        fname_degs_txt = self.fname_degs_txt % self.psi_id
        fname_sample_txt = self.fname_sample_txt % self.psi_id

        filename_degs = self.root_lfc / fname_degs
        filename_lfc = self.root_lfc / fname_lfc
        # filename_degs_txt = self.root_lfc / fname_degs_txt
        # filename_sample_txt = self.root_lfc / fname_sample_txt

        if filename_degs.exists() and filename_lfc.exists() and not force:
            df_degs = pdreadcsv(fname_degs, self.root_lfc, verbose=verbose)
            df_lfc = pdreadcsv(fname_lfc, self.root_lfc, verbose=verbose)
            
            try:
                degs_txt = read_txt(fname_degs_txt, self.root_lfc, verbose=verbose)
            except ValueError:
                degs_txt = ''

            try:
                sample_txt = read_txt(fname_sample_txt, self.root_lfc, verbose=verbose)
            except ValueError:
                sample_txt = ''

            return df_degs, df_lfc, degs_txt, sample_txt


        # geneid, symbol, biotype, samples
        min_N_cols = 3 + 3
        if df_normal.shape[1] < min_N_cols and df_gtex_ctrl.shape[1] < min_N_cols:
            msg = "Error: Normal samples and GTEx control do not have enough samples."
            return pd.DataFrame(), pd.DataFrame(), "", msg

        cdegs = CALC_DEGS(root_src=root_src, run_conda=run_conda)

        df_normal = cdegs.deduplicate_by_max_reads(df_normal)

        if df_normal.empty or df_normal.shape[1] < min_N_cols:
            df_normal2 = self.prepare_gtex(df_gtex_ctrl)
            df_normal2 = cdegs.deduplicate_by_max_reads(df_normal2)
            msg = f"not enough normal samples --> substituting with GTEx control {df_normal2.shape[1] - 3}."
        else:
            msg = f"enough normal samples {df_normal.shape[1] - 3}."
            df_normal2 = df_normal

        msg += f"\nThere are {df_tumor.shape[1] - 3} tumor samples; {msg}"

        df_tumor = cdegs.deduplicate_by_max_reads(df_tumor)
        self.df_tumor = df_tumor

        # if method == "deseq2":
        # at least 3 columns

        if df_tumor.empty or df_tumor.shape[1] < min_N_cols:
            print("Error: Tumor expression data has fewer than 3 samples.")
            return pd.DataFrame(), pd.DataFrame(), "", ""

        if df_normal2.empty or df_normal2.shape[1] < min_N_cols:
            print("Error: Normal expression data has fewer than 3 samples.")
            return pd.DataFrame(), pd.DataFrame(), "", ""

        df_lfc = cdegs.run_deg_rscript(
            df_tumor=df_tumor,
            df_normal=df_normal2,
            method=method,
            manual_dispersion=0.1,
            min_total_count=10,
            merge_how="inner",
            keep_temp=False,
        )
        self.df_lfc2 = df_lfc.copy()

        print(">>> columns:", df_lfc.columns.tolist())

        df_lfc = df_lfc.rename(columns={"log2FoldChange": "lfc", "padj": "fdr"})


        df_degs = df_lfc[(df_lfc.lfc >= lfc_cutoff) & (df_lfc.fdr < fdr_cutoff)].copy()
        df_degs.reset_index(drop=True, inplace=True)

        _ = pdwritecsv(df_lfc, fname_lfc, self.root_lfc)
        _ = pdwritecsv(df_degs, fname_degs, self.root_lfc)

        degs_txt = "\n".join(df_degs.symbol)
        _ = write_txt(degs_txt, fname_degs_txt, self.root_lfc)

        _ = write_txt(msg, fname_sample_txt, self.root_lfc)

        return df_degs, df_lfc, degs_txt, msg


    def read_GTEx_table(self, verbose: bool = False) -> pd.DataFrame:
        """
        read self.fname_gtex_table = 'tcga_primary_site_to_gtex_ids.tsv'
        output: df_gtex_to_tcga
        """

        filename = self.root_gtex / self.fname_gtex_table

        if not filename.exists():
            print(f"GTEx to TCGA table not found in {self.root_gtex}.")
            return pd.DataFrame()

        self.df_gtex_to_tcga = pdreadcsv(self.fname_gtex_table, self.root_gtex, verbose=verbose)

        return self.df_gtex_to_tcga

    def find_GTEx(self, verbose: bool = False) -> Tuple[str, str]:

        self.gtex_id = ""
        self.gtex_tissue_ids = ""

        if self.df_gtex_to_tcga.empty:
            df_gtex_to_tcga = self.read_GTEx_table(verbose=verbose)

            if df_gtex_to_tcga.empty:
                    print("GTEx to TCGA table is empty.")
                    return "", ""

        dfa = self.df_gtex_to_tcga[self.df_gtex_to_tcga.tcga_project_id == self.psi_id]

        if len(dfa) == 1:
            row = dfa.iloc[0]

            self.gtex_id = row.preferred_gtex_id
            self.gtex_tissue_ids = row.gtex_tissue_site_detail_ids

            if verbose:
                print(f"Found '{self.gtex_id}' tissue '{self.gtex_tissue_ids}'")

        elif len(dfa) == 0:
            if verbose:
                print("GTEx metada was not found")
        else:
            print("Multiple matches found")
            for i, row in dfa.iterrows():
                gtex_id = row.preferred_gtex_id
                gtex_tissue_ids = row.gtex_tissue_site_detail_ids

                print(f"{row.tcga_project_id} -> '{gtex_id}' tissue '{gtex_tissue_ids}'")

        return self.gtex_id, self.gtex_tissue_ids

    def cluster_data(self, df_tumor: pd.DataFrame, perc_min_samples: float = 0.25, 
                     top_n: int = 5_000) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray]:
        
        df_sel, df_cpm,  dfg_filt = self.calc_cpm_and_filter_data(df_tumor, perc_min_samples, top_n)
        
        # Scale genes
        df_scaled = StandardScaler().fit_transform(df_sel)

        return df_sel, df_cpm,  dfg_filt, df_scaled

    def calc_PCA(self, df_scaled: np.ndarray, n_components: int = 10, verbose: bool = False) -> pd.DataFrame:
        pca = PCA(n_components=n_components, random_state=42)

        df_pca = pca.fit_transform(df_scaled)

        df_pca = pd.DataFrame(
            df_pca[:, :3],
            index=self.df_sel.index,
            columns=["PC1", "PC2", "PC3"]
        )

        if verbose:
            print(pca.explained_variance_ratio_[:5])   

        return df_pca
    
    def calc_best_cluster(self, df_pca: pd.DataFrame, min_clusters: int = 3, max_clusters: int = 8) -> tuple[pd.DataFrame, pd.DataFrame]:

        cluster_results = []

        for k in range(min_clusters, max_clusters + 1):
            model = KMeans(n_clusters=k, random_state=42, n_init="auto")
            labels = model.fit_predict(df_pca)

            sil = silhouette_score(df_pca, labels)

            cluster_results.append({
                "k": k,
                "silhouette": sil,
                "labels": labels
            })

        df_eval = pd.DataFrame([
            {"k": r["k"], "silhouette": r["silhouette"]}
            for r in cluster_results
        ])

        # Choose best k
        best = max(cluster_results, key=lambda x: x["silhouette"])

        df_samp_clusters = pd.DataFrame({
            "sample": self.df_sel.index,
            "cluster": best["labels"] + 1
        })

        return df_eval, df_samp_clusters


    def plot_PCA(self, df_pca: pd.DataFrame, figsize : tuple = (6, 5)):
        plt.figure(figsize=figsize)
        plt.scatter(df_pca["PC1"], df_pca["PC2"], s=80)

        for sample in df_pca.index:
            plt.text(x=df_pca.loc[sample, "PC1"], y=df_pca.loc[sample, "PC2"], s=sample, fontsize=8)

        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("PCA of tumor samples")
        plt.tight_layout()
        plt.show()

             
    def calc_PCA_UMAP(self, df_pca: pd.DataFrame, df_samp_clusters: pd.DataFrame, n_neighbors: int = 5, 
                      min_dist: float = 0.2, metric: str = "euclidean") -> pd.DataFrame:

        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=42
        )

        X_umap = reducer.fit_transform(df_pca)

        df_umap = pd.DataFrame(
            X_umap,
            index=self.df_sel.index,
            columns=["UMAP1", "UMAP2"]
        )

        df_umap = df_umap.merge(
            df_samp_clusters,
            left_index=True,
            right_on="sample",
            how="left"
        )

        return df_umap

    def plot_PCA_UMAP(self, df_umap: pd.DataFrame, n_neighbors: int, min_dist: float, figsize : tuple = (6, 5)):
        plt.figure(figsize=figsize)
        plt.scatter(df_umap["UMAP1"], df_umap["UMAP2"], s=80)

        for sample in df_umap.index:
            plt.text(x=df_umap.loc[sample, "UMAP1"], y=df_umap.loc[sample, "UMAP2"], s=sample, fontsize=8)

        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title(f"UMAP of tumor samples (n_neighbors={n_neighbors}, min_dist={min_dist})")
        plt.tight_layout()
        plt.show()

    def plot_HCA_PCA(self, df_pca: pd.DataFrame,  method: str = "ward", figsize : tuple = (6, 5)):
        Z = linkage(df_pca, method=method)

        plt.figure(figsize=figsize)
        dendrogram(Z, labels=self.df_sel.index.tolist(), leaf_rotation=90)
        plt.title("PCA Hierarchical clustering of tumor samples")
        plt.tight_layout()
        plt.show()

    def cut_HCA_PCA(self, df_pca: pd.DataFrame, n_clusters: int = 3, 
                    method: str = "ward", criterion="maxclust", verbose:bool = True) -> pd.DataFrame:

        Z = linkage(df_pca, method=method)
        hc_labels = fcluster(Z, t=n_clusters, criterion=criterion)

        df_samp_clust_hc = pd.DataFrame({
            "sample": self.df_sel.index,
            "cluster": hc_labels
        })

        if verbose:
            print( df_samp_clust_hc.groupby("cluster").size() )

        return df_samp_clust_hc


    def plot_HCA_PCA_UMAP(self, df_umap: pd.DataFrame, method: str = "ward", figsize : tuple = (6, 5)):
        df2 = df_umap[ ['sample', 'UMAP1', 'UMAP2'] ]
        df2.set_index('sample', inplace=True)

        Z = linkage(df2, method=method)

        plt.figure(figsize=figsize)
        dendrogram(Z, labels=df2.index.tolist(), leaf_rotation=90)
        plt.title("PCA-UMAP Hierarchical clustering of tumor samples")
        plt.tight_layout()
        plt.show()


    def cut_HCA_PCA_UMAP(self, df_umap: pd.DataFrame, n_clusters: int = 3, 
                    method: str = "ward", criterion="maxclust", verbose:bool = True) -> pd.DataFrame:
        
        df2 = df_umap[ ['sample', 'UMAP1', 'UMAP2'] ]
        df2.set_index('sample', inplace=True)

        Z = linkage(df2, method=method)

        hc_labels = fcluster(Z, t=n_clusters, criterion=criterion)

        df_samp_clust_hc = pd.DataFrame({
            "sample": self.df_sel.index,
            "cluster": hc_labels
        })

        if verbose:
            print( df_samp_clust_hc.groupby("cluster").size() )

        return df_samp_clust_hc
    



    def calc_cpm_and_filter_data(self, df_tumor: pd.DataFrame, perc_min_samples: float = 0.25, 
                                 top_n: int = 5_000) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        ### Data treatment

        1. get raw dfc
        2. filter low-expression genes
        3. normalize for library size
        4. variance-stabilizing transformation
        5. select most variable genes
        6. cluster samples into k = 3..8 groups
        7. evaluate clusters
        8. find gene dfsig for each cluster

        A low-expression gene can be biologically important and even differentially expressed, especially if it is a transcription factor, cytokine, receptor, lncRNA, or rare-cell marker.

        But for unsupervised tumor clustering, we usually do not want thousands of genes with mostly zero/very low counts because they add noise and unstable distances.        
        """

        gene_cols = ["geneid", "symbol"]
        sample_cols = [c for c in df_tumor.columns if c not in gene_cols]

        dfc = (
            df_tumor[sample_cols]
            .apply(pd.to_numeric, errors="coerce")  # non-numeric -> NaN
            .fillna(0)                              # NaN -> 0
        ).copy()

        dfg = df_tumor[gene_cols].copy()

        dfc.index = df_tumor["geneid"]

        # filter low-count genes
        min_samples = int(perc_min_samples * len(sample_cols))

        print(f"sample_cols {len(sample_cols)} and min_samples")

        keep = list ((dfc >= 10).sum(axis=1) >= min_samples)

        dfc_filt = dfc.loc[keep]
        dfg_filt = dfg.loc[keep]

        # normalize by library size

        library_sizes = dfc_filt.sum(axis=0)

        df_cpm = dfc_filt.div(library_sizes, axis=1) * 1_000_000
        self.df_cpm = df_cpm

        dfc_log = np.log2(df_cpm + 1)

        # Select most variable genes
        gene_var = dfc_log.var(axis=1)

        top_genes = (
            gene_var
            .sort_values(ascending=False)
            .head(top_n)
            .index
        )

        df_sel = dfc_log.loc[top_genes].T.copy()
        self.df_sel = df_sel

        return df_sel, df_cpm,  dfg_filt


    def find_cluster_signature_genes(self, 
        df_logcpm: pd.DataFrame,
        df_samp_clusters: pd.DataFrame,
        gene_annot: pd.DataFrame,
        sample_col: str = "sample",
        cluster_col: str = "cluster",
        lfc_cutoff: float = 1.0,
        fdr_cutoff=0.05,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Find marker/signature genes for each cluster.

        df_logcpm:
            genes x samples matrix, log2(CPM + 1)

        df_samp_clusters:
            dataframe with columns: sample, cluster

        gene_annot:
            optional dataframe with geneid, symbol
        """

        results = []

        for cluster_id in sorted(df_samp_clusters[cluster_col].unique()):

            in_samples  = df_samp_clusters.loc[df_samp_clusters[cluster_col] == cluster_id, sample_col].tolist()
            out_samples = df_samp_clusters.loc[df_samp_clusters[cluster_col] != cluster_id, sample_col].tolist()

            # keep only samples present in expression matrix
            in_samples  = [s for s in in_samples  if s in df_logcpm.columns]
            out_samples = [s for s in out_samples if s in df_logcpm.columns]

            if len(in_samples) < 2 or len(out_samples) < 2:
                print(f"Skipping cluster {cluster_id}: too few samples")
                continue

            df_mean_in = df_logcpm[in_samples].mean(axis=1)
            df_mean_out = df_logcpm[out_samples].mean(axis=1)

            df_lfc = df_mean_in - df_mean_out

            pvals = []

            for geneid in df_logcpm.index:
                stat, p = ttest_ind(
                    df_logcpm.loc[geneid, in_samples],
                    df_logcpm.loc[geneid, out_samples],
                    equal_var=False,
                    nan_policy="omit",
                )
                pvals.append(p)

            fdr = multipletests(pvals, method="fdr_bh")[1]

            res = pd.DataFrame({
                "geneid": df_logcpm.index,
                "cluster": cluster_id,
                "n_in": len(in_samples),
                "n_out": len(out_samples),
                "mean_in": df_mean_in.values,
                "mean_out": df_mean_out.values,
                "lfc": df_lfc.values,
                "pvalue": pvals,
                "fdr": fdr,
            })

            if gene_annot is not None:
                res = res.merge(gene_annot, on="geneid", how="left")

            res = res.sort_values(
                ["lfc", "fdr"],
                ascending=[False, True]
            )

            results.append(res)

        dfall = pd.concat(results, ignore_index=True)

        dfsig = (
            dfall
            .query("lfc >= @lfc_cutoff and fdr <= @fdr_cutoff")
            .sort_values(["cluster", "lfc", "fdr"], ascending=[True, False, True])
            .reset_index(drop=True)
        )

        return dfall, dfsig


    def write_clusters(self, dfall: pd.DataFrame, dfsig: pd.DataFrame,
                       LFC_cutoff: float = 1, FDR_cutoff: float = 0.05, verbose: bool = True) -> pd.DataFrame:
    
        lista = np.unique(dfall.cluster)
        dic = {}; icount=-1

        for ncluster in lista:
            df2 = dfsig[dfsig.cluster == ncluster]
            df2 = df2[ (df2['lfc'].abs() > LFC_cutoff) & (df2['fdr'] < FDR_cutoff) ]
            
            s_genes = '\n'.join(df2.symbol)

            icount += 1
            dic[icount] = {}
            dic2 = dic[icount]
            dic2['ncluster'] = ncluster
            dic2['ngenes'] = len(df2)
            dic2['genes'] = df2.symbol.to_list()

            fname = f"cluster_{ncluster}_{self.psi_id}_signature_genes.txt"
            write_txt(s_genes, fname, self.root_lfc)

            if verbose:
                print(f"Cluster {ncluster} -> {len(df2)} signatures: {s_genes}")

        df = pd.DataFrame(dic).T

        fname = f"clusters_signatures_for_{self.psi_id}.txt"
        _ = pdwritecsv(df, fname, self.root_lfc, verbose=verbose)

        return df

    def add_entropy(self, df, read_limit: int = 50, min_read: int = 200, n_quantiles: int = 10):
        # select tumor columns
        cols = [c for c in df.columns if c.startswith("tumor_")]

        def row_entropy(values):
            values = np.asarray(values, dtype=float)

            # remove NaNs
            values = values[np.isfinite(values)]

            if len(values) == 0:
                return np.nan

            # all equal → entropy = 0
            if np.allclose(values, values[0]):
                return 0.0

            # shift to positive (important if negative values exist)
            values = values - values.min()

            # avoid all zeros
            if np.allclose(values.sum(), 0):
                return 0.0

            probs = values / values.sum()

            # avoid log(0)
            probs = probs[probs > 0]

            entropy = -np.sum(probs * np.log2(probs))
            return entropy

        min_cols = int(len(cols) * 0.4)
        good = (df[cols] < read_limit).sum(axis=1) <= min_cols
        df = df[good].copy()

        df["total_sum"] = df[cols].sum(axis=1)

        df = df[df.total_sum > len(cols) * min_read]

        df["entropy"] = df[cols].apply(lambda row: row_entropy(row.values), axis=1)

        df["entropy_q"] = pd.qcut(df["entropy"], q=n_quantiles, labels=False, duplicates="drop")

        return df

    def select_top_entropy(self, df: pd.DataFrame, q: float = 0.1):

        df = df[df.entropy_q <= q].copy()
        df = df.sort_values("entropy", ascending=True)
        df.reset_index(drop=True, inplace=True)

        return df

    def resolve_gtex_gencode_id(
        self, gene_id: str, datasetId: str = "gtex_v8", timeout: int = 60
    ) -> str:
        gene_base = str(gene_id).split(".")[0]

        url = f"{self.GTEX_API}/reference/gene"

        params = {
            "geneId": gene_base,
            "datasetId": datasetId,
        }

        r = requests.get(url, params=params, timeout=timeout)
        r.raise_for_status()

        data = r.json().get("data", [])

        if not data:
            return ""

        return data[0].get("gencodeId")

    def get_gtex_TPM_expression_for_geneid_list(
        self,
        geneid_list: List[str],
        datasetId: str = "gtex_v8",
        page_size: int = 1000,
        sleep: float = 0.1,
        timeout: int = 60,
        force: bool = False,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """
        Download GTEx normal tissue expression for selected genes.
        Returns long-format dataframe:
        geneSymbol | gencodeId | tissueSiteDetailId | sampleId | tpm | log2Tpm
        """

        fname = self.fname_tpm_exp % (self.gtex_id)
        filename = self.root_gtex / fname

        if filename.exists() and not force:
            return pdreadcsv(fname, self.root_gtex)

        all_rows = []

        print(">>>", len(geneid_list))
        for i, geneid in enumerate(geneid_list):
            page = 0

            while True:
                gencode_id = self.resolve_gtex_gencode_id(geneid)
                print(">>>", gencode_id, end=" ")

                if gencode_id is None:
                    print("?")
                    if verbose:
                        print(f"Nothing found for {geneid}")
                    break

                url = f"{self.GTEX_API}/expression/geneExpression"

                params = {
                    "datasetId": datasetId,
                    "gencodeId": gencode_id,
                    "tissueSiteDetailId": self.gtex_id,
                    "page": page,
                    "itemsPerPage": page_size,
                }

                r = requests.get(url, params=params, timeout=timeout)
                r.raise_for_status()

                js = r.json()
                data = js.get("data", [])

                if not data:
                    print("x")
                    if verbose:
                        print(
                            f"Error: could not retrieve data for {self.gtex_id} for geneid '{geneid}'"
                        )
                    break

                print("Ok")

                all_rows.extend(data)

                paging = js.get("paging_info", {})
                n_pages = paging.get("numberOfPages", 1)

                page += 1
                if page >= n_pages:
                    break

                time.sleep(sleep)

        print("")
        df_gtex = pd.DataFrame(all_rows)

        df_gtex.reset_index(drop=True, inplace=True)

        _ = pdwritecsv(df_gtex, fname, self.root_gtex, verbose=verbose)

        return df_gtex

    def read_GTEx_counts_pheno_meta(self, verbose: bool = False):
        # load matrix'1

        if self.df_gtex_counts.empty:
            print("Reading GTEx count table... be patient.")
            df_gtex_counts = pdreadcsv(self.fname_GTEx_counts, self.root_gtex, skiprows=2, verbose=verbose)
            self.df_gtex_counts = df_gtex_counts

        if self.df_gtex_pheno.empty:
            _ = self.read_GTEx_table_pheno(verbose=verbose)

        if self.df_meta.empty:
            _ = self.read_GTEx_table_metadata(verbose=verbose)

        return self.df_gtex_counts

    def read_GTEx_table_pheno(self, verbose: bool = False):
        # load phenotype
        df_gtex_pheno = pdreadcsv(self.fname_GTEx_pheno, self.root_gtex, verbose=verbose)
        self.df_gtex_pheno = df_gtex_pheno
        return df_gtex_pheno

    def read_GTEx_table_metadata(self, verbose: bool = False):
        # load metadata
        df_meta = pdreadcsv(self.fname_GTEx_meta, self.root_gtex, verbose=verbose)
        self.df_meta = df_meta
        return df_meta

    def get_gtex_control(
        self, Nsamples: int = 15, force: bool = False, verbose: bool = False
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:

        self.df_gtex_ctrl = pd.DataFrame()
        self.df_meta_prep = pd.DataFrame()

        gtex_id, _ = self.find_GTEx(verbose=verbose)
        self.gtex_id = gtex_id

        if gtex_id == "":
            print(f"Error: could not find GTEx ID for {self.prog_id} ID '{self.psi_id}'")
            return pd.DataFrame(), pd.DataFrame()

        # prepare metadata
        df_meta_prep = self.prepare_gtex_df_meta()

        # get all GTEx samples for self.psi_id
        df_gtex_ctrl = self.prepare_gtex_count_table(Nsamples=Nsamples, force=force, verbose=verbose)

        self.df_gtex_ctrl = df_gtex_ctrl
        self.df_meta_prep = df_meta_prep

        return df_gtex_ctrl, df_meta_prep

    def prepare_gtex_df_meta(self) -> pd.DataFrame:
        """
        1. Filter for Tissue
        2. Filter for high quality (Hardy Scale 1 or 2 are 'fast' deaths, less stress)
          . DTHHRDY: 1 = Ventilator, 2 = Fast death of natural causes
        3. Sort by RIN score (SMRIN) to get the best preserved RNA - Then take the top 15

        Test:
            gtex_id = 'Skin_Sun_Exposed_Lower_leg'
            gtex_id = 'Muscle_Skeletal'
            gtex_id = 'Colon_Transverse'
            gtex_id = 'Brain'
            gtex_id = 'Whole Blood'

            print(gtex_id)
            term = " - ".join(gtex_id.split('_')[:2])
            print("term", term)
            df_meta = mtd.gdc.df_meta

            df2 = df_meta[df_meta["SMTSD"].str.startswith(term)]
            print(len(df2))
            df2

        output: df_meta_prep
        """

        print("Preparing GTEx metadata...")

        if self.df_meta.empty:
            self.read_GTEx_table_metadata()

            if self.df_meta.empty:
                print("Metadata DataFrame is empty.")
                return pd.DataFrame()

        if self.df_gtex_pheno.empty:
            self.read_GTEx_table_pheno()

            if self.df_gtex_pheno.empty:
                print("Phenotype DataFrame is empty.")
                return pd.DataFrame()

        # 1. Filter for Tissue
        # df_meta_prep is df_meta filtered by gtex_id
        # 'Colon_Transverse' --> 'Colon - Transverse'

        if self.gtex_id == 'Adrenal_Gland' or self.gtex_id == 'Whole_Blood':
            term = self.gtex_id.replace('_', ' ')
        else:
            term = " - ".join(self.gtex_id.split('_')[:2])
        df_meta_prep = self.df_meta[self.df_meta["SMTSD"].str.startswith(term)].copy()
        df_meta_prep.reset_index(drop=True, inplace=True)

        self.df_meta_prep = df_meta_prep

        if df_meta_prep.empty:
            print(f"Error: could not find metadata for gtex_id '{self.gtex_id}' to term ''")
            return pd.DataFrame()

        if len(df_meta_prep) < 3:
            print(f"Warning: could not find enough metadata for gtex_id '{self.gtex_id}' to term '{term}'")


        lista = df_meta_prep["SAMPID"].str.split("-").str[:2].str.join("-")
        df_meta_prep.loc[:, "SUBJID"] = lista
        df_meta_prep = df_meta_prep.merge(self.df_gtex_pheno, on="SUBJID", how="left")

        # 2. Filter for high quality (Hardy Scale 1 or 2 are 'fast' deaths, less stress)
        # DTHHRDY: 1 = Ventilator, 2 = Fast death of natural causes
        df_meta_prep = df_meta_prep[df_meta_prep["DTHHRDY"].isin([1, 2])]

        # 3. Sort by RIN score (SMRIN) to get the best preserved RNA
        # Then take the top 15
        df_meta_prep = df_meta_prep.sort_values("SMRIN", ascending=False)
        df_meta_prep.reset_index(drop=True, inplace=True)

        self.df_meta_prep = df_meta_prep
        print(f"GTEx metadata prepared on df_meta_prep length: {len(df_meta_prep)}")

        return df_meta_prep

    def running_on_render(self) -> bool:
        return "RENDER_SERVICE_ID" in os.environ

    def prepare_gtex_count_table(
        self, Nsamples=15, force: bool = False, verbose: bool = False
    ) -> pd.DataFrame:

        fname = self.fname_gtex_exp_counts % (self.gtex_id)
        filename = self.root_gtex / fname

        if filename.exists() and not force:
            df_gtex_normal = pdreadcsv(fname, self.root_gtex, verbose=verbose)
            self.df_gtex_normal = df_gtex_normal            
            return df_gtex_normal

        if self.df_gtex_counts.empty:
            # Reading GTEx count super-file
            # to big, do not store in Render
            if not self.running_on_render():
                _ = self.read_GTEx_counts_pheno_meta(verbose=verbose)

            if self.df_gtex_counts.empty:
                print("Count DataFrame is empty.")
                return pd.DataFrame()

        self.df_meta_prep = self.df_meta_prep.sort_values("SMRIN", ascending=False)
        samples = self.df_meta_prep.head(Nsamples * 3)["SAMPID"].to_list()

        cols = list(samples)
        good_cols = [x for x in cols if x in self.df_gtex_counts.columns]
        if len(good_cols) > Nsamples:
            good_cols = good_cols[:Nsamples]

        good_cols = ["Name", "Description"] + good_cols

        df_gtex_normal = self.df_gtex_counts.loc[:, good_cols].copy()
        df_gtex_normal.rename(columns={"Name": "ensemblid", "Description": "symbol"}, inplace=True)
        df_gtex_normal.reset_index(drop=True, inplace=True)

        self.df_gtex_normal = df_gtex_normal

        _ = pdwritecsv(df_gtex_normal, fname, self.root_gtex, verbose=verbose)

        return df_gtex_normal

    def download_file(self, url: str, filename_out: str):
        with requests.get(url, stream=True) as r:
            r.raise_for_status()

            try:
                with open(filename_out, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            except Exception as e:
                print(f"Error: as writing file {filename_out}: {e}")


    def build_df_exp_and_filter(self,
        df_counts: pd.DataFrame,
        df_meta: pd.DataFrame,
        gene_col: str = "geneid",
        condition_col: str = "condition",
        sample_col: str = "sample",
        tumor_label: str = "tumor",
        normal_label: str = "normal",
        equal_var: bool = False,   # False = Welch t-test, safer when n differs
    ) -> tuple[pd.DataFrame, list, list]:
        
        df = df_counts.copy()

        # Samples by condition
        normal_samples = df_meta.loc[
            df_meta[condition_col] == normal_label, sample_col
        ].tolist()

        tumor_samples = df_meta.loc[
            df_meta[condition_col] == tumor_label, sample_col
        ].tolist()

        # Keep only samples present in df_counts
        normal_samples = [s for s in normal_samples if s in df.columns]
        tumor_samples = [s for s in tumor_samples if s in df.columns]

        sample_cols = normal_samples + tumor_samples

        ncols_normal = len(normal_samples)
        ncols_tumor  = len(tumor_samples)

        nmin_cols = min(ncols_normal, ncols_tumor)

        df["total"] = df[sample_cols].sum(axis=1)

        df = df.loc[
            df["total"] > nmin_cols * 25
        ].reset_index(drop=True, inplace=False)



        # Convert counts to numeric
        df[normal_samples + tumor_samples] = df[normal_samples + tumor_samples].apply(
            pd.to_numeric, errors="coerce"
        )

        # Optional but recommended for RNA-seq counts:
        # log-transform before t-test
        normal_mat = np.log2(df[normal_samples] + 1)
        tumor_mat = np.log2(df[tumor_samples] + 1)

        # Row-wise t-test: tumor vs normal
        t_stat, pval = ttest_ind(
            tumor_mat,
            normal_mat,
            axis=1,
            equal_var=equal_var,
            nan_policy="omit",
        )

        df["t_stat"] = t_stat
        df["pval"] = pval

        # Useful summaries
        df["mean_normal"] = normal_mat.mean(axis=1)
        df["mean_tumor"] = tumor_mat.mean(axis=1)
        df["lfc"] = df["mean_tumor"] - df["mean_normal"]
        df["abs_lfc"] = df["lfc"].abs()

        # Order by p-value
        df = df.sort_values("pval", ascending=True)

        # Keep the 40% lowest p-values
        df = df[df.lfc < 0.01]
        df.reset_index(drop=True, inplace=True)

        return df, normal_samples, tumor_samples


    def plot_heatmap_expression(self, dff: pd.DataFrame, normal_samples: list, tumor_samples: list, 
                                title0: str = "", figsize: tuple = (14, 10)):
        cols = ['geneid'] + normal_samples + tumor_samples

        dff2 = dff[cols].copy()
        dff2.set_index('geneid', inplace=True)

        title = 'Hierarchical Clustering of Expression Data'
        if title0 != '':
            title += '\n' + title0

        # numeric matrix
        dff2 = dff2.apply(pd.to_numeric, errors="coerce").fillna(0)
        mat = np.log2(dff2 + 1)

        # gene-wise z-score using pandas/numpy
        row_mean = mat.mean(axis=1)
        row_std = mat.std(axis=1)

        mat_z = mat.sub(row_mean, axis=0).div(row_std.replace(0, np.nan), axis=0)
        mat_z = mat_z.replace([np.inf, -np.inf], np.nan).fillna(0)

        cg = sns.clustermap(
            mat_z,
            metric="correlation",
            method="average",
            figsize=figsize,
            cmap="viridis",
            cbar=True,
        )

        title = "Hierarchical Clustering of Expression Data"
        cg.figure.suptitle(title, y=1.02)

        return cg


    def plot_umap_expression(
        self,
        dff: pd.DataFrame,
        samples: list,
        which_samples: str = 'Tumor',
        title: str = "",
        figsize: tuple = (8, 6),
        n_neighbors: int = 10,
        min_dist: float = 0.2,
        n_clusters: int = 3,
        random_state: int = 42,
    ):
        # ------------------------------------------------------------
        # 1. Use only tumor samples
        # ------------------------------------------------------------
        cols = ["geneid"] + samples

        dff2 = dff[cols].copy()
        dff2.set_index("geneid", inplace=True)

        # numeric matrix
        dff2 = dff2.apply(pd.to_numeric, errors="coerce").fillna(0)

        # log transform
        mat = np.log2(dff2 + 1)

        # ------------------------------------------------------------
        # 2. Gene-wise z-score across tumor samples only
        # ------------------------------------------------------------
        row_mean = mat.mean(axis=1)
        row_std = mat.std(axis=1)

        mat_z = mat.sub(row_mean, axis=0).div(row_std.replace(0, np.nan), axis=0)
        mat_z = mat_z.replace([np.inf, -np.inf], np.nan).fillna(0)

        # ------------------------------------------------------------
        # 3. UMAP uses samples x genes
        # ------------------------------------------------------------
        X_umap = mat_z.T.copy()

        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric="correlation",
            random_state=random_state,
        )

        emb = reducer.fit_transform(X_umap)

        # ------------------------------------------------------------
        # 4. Cluster tumor samples
        # ------------------------------------------------------------
        kmeans = KMeans(
            n_clusters=n_clusters,
            random_state=random_state,
            n_init="auto",
        )

        clusters = kmeans.fit_predict(emb)

        df_umap = pd.DataFrame(
            {
                "sample": X_umap.index,
                "UMAP1": emb[:, 0],
                "UMAP2": emb[:, 1],
                "cluster": clusters.astype(str),
            }
        )

        # ------------------------------------------------------------
        # 5. Plot UMAP colored by tumor cluster
        # ------------------------------------------------------------
        if not title:
            title = f"UMAP expression clustering using {which_samples} samples"
        else:
            title += f" using {which_samples} samples"

        fig_umap, ax_umap = plt.subplots(figsize=figsize)

        sns.scatterplot(
            data=df_umap,
            x="UMAP1",
            y="UMAP2",
            hue="cluster",
            s=80,
            ax=ax_umap,
        )

        ax_umap.set_title(title)
        ax_umap.set_xlabel("UMAP1")
        ax_umap.set_ylabel("UMAP2")
        ax_umap.legend(title="Tumor cluster")

        fig_umap.tight_layout()

        return fig_umap, ax_umap, df_umap
