#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/05/08
# Udated  on 2026/05/08
# @author: Flavio Lichtenstein
# @local: Home sweet home

from dataclasses import dataclass
from pathlib import Path
import os
import numpy as np

@dataclass
class ProjectContext:
    email: str | None

    i_project: int
    project_list: list
    project: str
    s_project_list: list
    s_project: str

    PSI_ID: str
    disease: str

    gene_protein: str
    s_omics: str

    has_age: bool
    has_gender: bool

    exp_normalization: bool
    normalization: str

    LFC_cut_inf: float
    s_pathw_enrichm_method: str
    ptw_min_num_of_degs_cut: int

    tolerance_pPMI: float
    type_sat_ptw_index: str
    saturation_lfc_param: float

    pval_pathway_cutoff: float
    fdr_pathway_cutoff: float
    num_of_genes_cutoff: int
    enr_db_list: list

    case_list: list
    dic_case_list: dict

    std_filename: str
    std_filename_list: list

    min_lfc_modulation: float
    num_of_genes_list: list
    pPMI_normalized: bool

    s_len_case: int

    n_sentences: int
    run_list: list
    chosen_model_list: list
    i_dfp_list: list
    chosen_model_sampling: str

    fdr_ptw_cutoff_list: np.ndarray
    lfc_list: np.ndarray
    fdr_list: np.ndarray

    colors: list
  

def load_project_context(dic_yml: dict, PSI_ID:str = "TCGA-BRCA", i_project: int = 0) -> ProjectContext:
    
    email = os.getenv("email")

    project_list = dic_yml["project_list"]
    n = len(project_list)

    project = project_list[i_project]

    s_project_list = dic_yml["s_project_list"]
    s_project = s_project_list[i_project]

    assert n == len(s_project_list), (
        f"Error: project_list has {n} projects, "
        f"but s_project_list has {len(s_project_list)} projects"
    )

    disease = PSI_ID

    gene_protein = dic_yml["gene_protein"]
    s_omics = dic_yml["s_omics"]

    has_age = dic_yml["has_age"]
    has_gender = dic_yml["has_gender"]

    exp_normalization = dic_yml["exp_normalization"]
    normalization = "quantile_norm" if exp_normalization is True else "not_normalized"

    LFC_cut_inf = dic_yml["LFC_cut_inf"]
    s_pathw_enrichm_method = dic_yml["s_pathw_enrichm_method"]
    ptw_min_num_of_degs_cut = dic_yml["ptw_min_num_of_degs_cut"]

    tolerance_pPMI = dic_yml["tolerance_pPMI"]
    type_sat_ptw_index = dic_yml["type_sat_ptw_index"]
    saturation_lfc_param = dic_yml["saturation_lfc_param"]

    pval_pathway_cutoff = dic_yml["pval_pathway_cutoff"]
    fdr_pathway_cutoff = dic_yml["fdr_pathway_cutoff"]
    num_of_genes_cutoff = dic_yml["num_of_genes_cutoff"]
    enr_db_list = dic_yml["enr_db_list"]

    case_list = dic_yml["case_list"]
    dic_case_list = dic_yml["dic_case_list"]

    std_filename = dic_yml["std_filename"]
    std_filename_list = dic_yml["std_filename_list"]

    min_lfc_modulation = dic_yml["min_lfc_modulation"]
    num_of_genes_list = dic_yml["num_of_genes_list"]
    pPMI_normalized = dic_yml["pPMI_normalized"]

    s_len_case = dic_yml["s_len_case"]

    n_sentences = dic_yml["n_sentences"]
    run_list = dic_yml["run_list"]
    chosen_model_list = dic_yml["chosen_model_list"]
    i_dfp_list = dic_yml["i_dfp_list"]
    chosen_model_sampling = dic_yml["chosen_model_sampling"]

    fdr_ptw_cutoff_list = np.arange(0.05, 0.80, 0.05)
    lfc_list = np.round(np.arange(1.0, -0.01, -0.025), 3)
    fdr_list = np.arange(0.05, 0.76, 0.01)

    colors = ["red", "green", "blue", "orange", "pink", "purple", "black", "cyan", "tomato", "lime", 
              "magenta", "yellow", "gray", "brown", "olive", "navy", "teal", "maroon", "silver",]
    
    return ProjectContext(
        email=email,

        i_project=i_project,
        project_list=project_list,
        project=project,
        s_project_list=s_project_list,
        s_project=s_project,

        PSI_ID=PSI_ID,
        disease=disease,

        gene_protein=gene_protein,
        s_omics=s_omics,

        has_age=has_age,
        has_gender=has_gender,

        exp_normalization=exp_normalization,
        normalization=normalization,

        LFC_cut_inf=LFC_cut_inf,
        s_pathw_enrichm_method=s_pathw_enrichm_method,
        ptw_min_num_of_degs_cut=ptw_min_num_of_degs_cut,

        tolerance_pPMI=tolerance_pPMI,
        type_sat_ptw_index=type_sat_ptw_index,
        saturation_lfc_param=saturation_lfc_param,

        pval_pathway_cutoff=pval_pathway_cutoff,
        fdr_pathway_cutoff=fdr_pathway_cutoff,
        num_of_genes_cutoff=num_of_genes_cutoff,
        enr_db_list=enr_db_list,

        case_list=case_list,
        dic_case_list=dic_case_list,

        std_filename=std_filename,
        std_filename_list=std_filename_list,

        min_lfc_modulation=min_lfc_modulation,
        num_of_genes_list=num_of_genes_list,
        pPMI_normalized=pPMI_normalized,

        s_len_case=s_len_case,

        n_sentences=n_sentences,
        run_list=run_list,
        chosen_model_list=chosen_model_list,
        i_dfp_list=i_dfp_list,
        chosen_model_sampling=chosen_model_sampling,

        fdr_ptw_cutoff_list=fdr_ptw_cutoff_list,
        lfc_list=lfc_list,
        fdr_list=fdr_list,
        colors=colors,
        )
