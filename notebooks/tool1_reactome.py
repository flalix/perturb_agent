#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
'''
# @Updated on 2026/03/20
# @Created on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home
'''

#import requests
#import json
import os, sys
import pandas as pd
from typing import List, Tuple # Optional, Iterable, Set, Any,

import streamlit as st

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

st.write(f"ROOT: {ROOT}", )
st.write(f"SRC added: {SRC}")

from libs.Basic import *
from libs.reactome_lib import *

st.set_page_config(page_title="Target Discovery")
st.title("🧬 Reactome Tool")

root_data='../data/reactome/'
rea = Reactome(root_data=root_data)

def open_tables(verbose:bool=True):
 
    _ = rea.open_reactome(verbose=verbose)
    _ = rea.open_reactome_abstract(verbose=verbose)
    _ = rea.open_reactome_gmt(verbose=verbose)

    return


def tool1(pathway:str) -> Tuple[str, List]:
    '''
    Given a pathway, find it in df_gmt and return pathway_id and genes

    input: pathway
    output: pathway_id, genes (list)
    '''
 
    row = rea.df_gmt[rea.df_gmt.pathway == pathway].iloc[0]

    genes = row.genes
    if isinstance(genes, str):
        genes = eval(genes)

    return row.pathway_id, genes



verbose=True
open_tables(verbose=verbose)

REACTOME_PATHWAYS = rea.pathway_list

n = rea.n_pathways

if len(rea.pathway_list) == 0:
    st.write("No pathways found. Stopping")
    exit(-1)

# --- Dropdown selection ---
selected_pathway = st.selectbox(
    f"Select a pathway ({n} pathways were annotated):",
    REACTOME_PATHWAYS,
    # index=None,
    placeholder="Start typing..."
)

# --- Button to trigger computation ---
if st.button("Run find the pathway"):

    pathway = ''

    with st.spinner(f"Pathway: {selected_pathway}..."):

        # Call your backend tool
        pathway=selected_pathway

        pathway_id, genes = tool1(pathway)

    st.success(f"Pathway found: '{pathway}' ({pathway_id})")

    # --- Display results ---
    st.subheader("Annotated genes")
    st.write(f"{genes}")




    

