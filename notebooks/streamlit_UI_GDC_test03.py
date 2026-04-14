#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

#=============== to run =====================
#
# export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
#
# uv run streamlit run streamlit_UI_GDC_test03.py 
#
#============================================


import os, sys
from pprint import pprint
# from marshmallow import pprint, Schema, fields
import numpy as np
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from collections import defaultdict

#----------- fix incompatibilities ---------------------
import pandas as pd
setattr(pd.Series, "iteritems", pd.Series.items)
setattr(pd.DataFrame, "iteritems", pd.DataFrame.items)

import matplotlib.pyplot as plt
import plotly.express as px

import seaborn as sns
from sklearn.cluster import KMeans
import umap

from pathlib import Path

ROOT = Path().resolve().parent
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

from libs.calc_degs_lib import CALC_DEGS
from libs.tcga_gdc_lib import *
from libs.Basic import *

root0 = ROOT / "data"

print("root:", ROOT)
print("src added:", SRC)
print("data:", root0)

gdc = GDC(root0=root0)

verbose = True
colors = ['red', 'green', 'blue', 'orange', 'pink', 'purple', 'black', 'cyan', 'tomato', 'lime', 'magenta', 'yellow',
          'gray', 'brown', 'olive', 'navy', 'teal', 'maroon', 'silver']


# -----------------------------------------------------------------------------
# PAGE
# -----------------------------------------------------------------------------
st.set_page_config(page_title="GDC / TCGA Explorer", layout="wide")

st.markdown("""
<style>
/* Label (title) */
div[data-baseweb="select"] > label {
    font-size: 22px !important;
    font-weight: 700 !important;
    color: navy !important;
}

/* Selectbox text */
div[data-baseweb="select"] div {
    font-size: 22px !important;
    font-weight: 500 !important;
    color: navy !important;
}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<style>
.block-container {
    max-width: 96%;
    padding-top: 0.5rem;
    padding-bottom: 0.5rem;
    padding-left: 1.5rem;
    padding-right: 1.5rem;
}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<style>
div.stButton > button {
    width: 100%;
}
</style>
""", unsafe_allow_html=True)

st.title("GDC / TCGA Explorer")
st.caption("Explore cases, tumor samples, and mutation matrices by primary site")


force=False
verbose=False

min_barcodes=3
min_genes=5

# --- FOOTER / BOX BELOW ALL TABS ---
def show_profile_box():
    st.markdown("""
    <div style="display:flex; justify-content:left; margin-top:40px;">
        <div style="
            margin-top: 40px;
            background-color: lightblue;
            padding: 25px;
            border-radius: 12px;
            border: 1px solid #A0C4FF;
            max-width: 600px;
        ">
            <h3 style="margin-bottom:5px;">PhD Flavio Lichtenstein</h3>
            <p style="color:#A0A0A0;">
                Bioinformatics, Immunoinformatics, Biostatistics,<br>
                Systems Biology and Artificial Intelligence
            </p>
            <hr>
            <p>📍 Sao Paulo, SP - Brazil</p>
            <p>📞 (+55) 11 96560-1960</p>
            <p>✉️ flalix@gmail.com</p>
            <p>🔗 <a href="https://www.linkedin.com/in/flaviolichtenstein/" target="_blank" style="color:#0A66C2; text-decoration:none;">LinkedIn Profile</a>
            </p>
    </div>
    """, unsafe_allow_html=True)


def make_streamlit_safe(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df

    out = df.copy()
    
    # Make labels plain Python strings
    out.index = out.index.map(str)
    out.columns = [str(c) for c in out.columns]
    
    # Flatten MultiIndex columns if present
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [" | ".join(map(str, col)).strip() for col in out.columns.to_flat_index()]


    # Convert every column to a Streamlit-safe pandas dtype
    for col in out.columns:
        s = out[col]

        # nullable / pyarrow / string extension -> Python objects
        try:
            if pd.api.types.is_string_dtype(s.dtype):
                out[col] = s.astype("object")
            elif str(s.dtype).lower().find("pyarrow") >= 0:
                out[col] = s.astype("object")
            else:
                out[col] = s
        except Exception:
            out[col] = s.astype("object")

        # Replace pandas NA/NaT with None so Arrow does not try fancy typing
        out[col] = out[col].where(pd.notna(out[col]), None)


    return out


def make_aggrid_safe(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df

    out = df.copy()

    # Flatten columns
    if isinstance(out.columns, pd.MultiIndex):
        out.columns = [" | ".join(map(str, c)).strip() for c in out.columns.to_flat_index()]
    else:
        out.columns = [str(c) for c in out.columns]

    # String index
    out.index = [str(i) for i in out.index]

    # Remove duplicated column names
    if pd.Index(out.columns).duplicated().any():
        seen = {}
        new_cols = []
        for c in out.columns:
            if c in seen:
                seen[c] += 1
                new_cols.append(f"{c}_{seen[c]}")
            else:
                seen[c] = 0
                new_cols.append(c)
        out.columns = new_cols

    # Force plain Python scalar values only
    for col in out.columns:
        out[col] = out[col].map(
            lambda x: None if pd.isna(x)
            else str(x) if isinstance(x, (list, dict, set, tuple, np.ndarray))
            else x
        )
        out[col] = out[col].astype(object)

    return out


def show_df_AgGrid2(df, height=800, page_size=25, key="grid"):
    if df is None or df.empty:
        st.info("Empty dataframe")
        return

    df = make_aggrid_safe(df).copy()

    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_default_column(sortable=True, filter=True, resizable=True)

    grid_options = gb.build()
    grid_options["pagination"] = True
    grid_options["paginationPageSize"] = page_size
    grid_options["suppressPaginationPanel"] = False

    AgGrid(
        df,
        gridOptions=grid_options,
        height=height,
        key=key,
        # fit_columns_on_grid_load=True,
    )

def show_df_AgGrid(df, height:int=800, page_size:int=25, key:str="grid"):
    if df is None or df.empty:
        st.info("Empty dataframe")
        return

    df = make_aggrid_safe(df)

    gb = GridOptionsBuilder.from_dataframe(df)

    gb.configure_default_column(
        sortable=True,
        minWidth=150,
        filter=True, 
        resizable=True)
    

    gb.configure_pagination(
        enabled=True,
        paginationAutoPageSize=False,
        paginationPageSize=page_size,
    )

    grid_options = gb.build()
    grid_options["pagination"] = True
    grid_options["paginationPageSize"] = page_size

    if not isinstance(grid_options, dict):
        raise TypeError(f"grid_options must be dict, got {type(grid_options)}")

    AgGrid(
        df,
        gridOptions=grid_options,
        height=height,
        # fit_columns_on_grid_load=True, # nice UX, no horizontal scroll
        allow_unsafe_jscode=False, # safe (good default)
        enable_enterprise_modules=False, # lightweight
        reload_data=False,  # avoids flicker / rerender
        key=key,
    )


def show_df(df, height:int=800, page_size:int=25, key:str="grid"):
    show_df_AgGrid(df, height=height, page_size=page_size, key=key)


def show_df_html(df, height: int = 800):
    if df is None or df.empty:
        st.info("Empty dataframe")
        return

    df = df.copy()

    # flatten columns first
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = [" | ".join(map(str, c)).strip() for c in df.columns.to_flat_index()]
    else:
        df.columns = [str(c) for c in df.columns]

    df.index = [str(x) for x in df.index]

    # force plain Python values
    for col in df.columns:
        df[col] = df[col].map(lambda x: None if pd.isna(x) else x)

    # last-resort HTML rendering: no Arrow, no LargeUtf8
    html = df.to_html(index=False, escape=False)
    st.markdown(
        f"""
        <div style="height:{height}px; overflow:auto; border:1px solid #ddd; padding:0.25rem;">
            {html}
        </div>
        """,
        unsafe_allow_html=True,
    )


def plot_top_mutated_genes(dfpiv: pd.DataFrame, top_n:int=20, figsize=(12,6)):
    if dfpiv is None or dfpiv.empty:
        st.info("No mutation matrix available.")
        return

    if dfpiv.shape[0] == 0:
        st.info("No barcodes available.")
        return

    gene_freq = (dfpiv.sum(axis=0) / dfpiv.shape[0]).sort_values(ascending=False)
    top_genes = gene_freq.head(top_n)

    fig, ax = plt.subplots(figsize=figsize)
    top_genes.plot(kind="bar", ax=ax)

    ax.set_ylabel("Fraction of barcodes mutated")
    ax.set_xlabel("Gene symbol")
    ax.set_title(f"Top {top_n} mutated genes")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()

    st.pyplot(fig)
    plt.close(fig)

def plot_heatmap(dfpiv: pd.DataFrame, title:str="", figsize:tuple=(14, 10)):
    # Ensure numeric + binary (important for Jaccard)
    data = dfpiv.fillna(0).astype(int)

    cg = sns.clustermap(
                data,
                metric="jaccard",
                method="average",
                figsize=figsize,
                cmap="viridis",
                cbar=False
            )

    if title:
        cg.figure.suptitle(title, y=1.02)

    st.pyplot(cg.figure)
    plt.close(cg.figure)

def plot_umap(dfpiv: pd.DataFrame, k:int=8, figsize:tuple=(14, 10)):

    fig, _, _ = gdc.plot_UMAP(dfpiv=dfpiv, k=k, figsize=figsize)

    if fig:
        st.pyplot(fig)
        plt.close(fig)

def plot_hdbscan(dfpiv: pd.DataFrame, min_cluster_size:int=10, min_samples:int=3, figsize:tuple=(14, 10)):

    fig, embedding, labels, d = gdc.plot_HDBSCAN(dfpiv=dfpiv, 
                                                 min_cluster_size=min_cluster_size,
                                                 min_samples=min_samples, figsize=figsize)

    if fig:
        st.pyplot(fig)
        plt.close(fig)



# prog_list = gdc.get_gdc_progams(force=False, verbose=verbose)

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
# hash error: @st.cache(show_spinner=True)

@st.cache_data(show_spinner=False)
def load_primary_site_data( primary_site:str, 
                           verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list]:

    df_cases, df_all_samples, df_all_mut, barcode_list = gdc.get_filtered_tables(primary_site=primary_site, verbose=verbose)

    return (
        make_streamlit_safe(df_cases),
        make_streamlit_safe(df_all_samples),
        make_streamlit_safe(df_all_mut),
        barcode_list,
    )


# hash error: @st.cache(show_spinner=False)
@st.cache_data(show_spinner=False)
def st_build_pivot_table(df_all_mut: pd.DataFrame, min_barcodes:int=3, min_genes:int=5) -> pd.DataFrame:
    """
    Build barcode x gene boolean mutation matrix.
    """
    if df_all_mut is None or df_all_mut.empty:
        return pd.DataFrame()

    dfpiv = gdc.build_pivot_table(df_all_mut, min_barcodes=min_barcodes, min_genes=min_genes)
    dfpiv = make_streamlit_safe(dfpiv)

    return dfpiv.sort_index(axis=0).sort_index(axis=1)


# hash error:@st.cache(show_spinner=False)
def summarize_mutations(df_all_mut: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize mutated genes by number of patients/barcodes.
    """
    if df_all_mut is None or df_all_mut.empty:
        return pd.DataFrame(columns=["symbol", "n_patients_mutated"])

    df = ( df_all_mut.groupby("symbol")["barcode"]
                .nunique()
                .reset_index(name="n_patients_mutated")
                .sort_values("n_patients_mutated", ascending=False)
                .reset_index(drop=True)
        )

    df = make_streamlit_safe(df)
    return df


def safe_unique_sorted(series):
    vals = [x for x in series.dropna().unique().tolist() if str(x).strip() != ""]
    return sorted(vals)


# -----------------------------------------------------------------------------
# SIDEBAR
# -----------------------------------------------------------------------------
if "loaded" not in st.session_state:
    st.session_state.loaded = False

with st.sidebar:
    st.header("Controls")

    # in the future --> dropdown program selector
    prog_id = 'TCGA'

    st.text(f"Program {prog_id}")
    # force = st.checkbox("Force rebuild", value=False)
    # verbose = st.checkbox("Verbose", value=False)

    load_clicked = st.button("Load data", use_container_width=True)

    if load_clicked:
        st.session_state.loaded = True 


# -----------------------------------------------------------------------------
# MAIN LOAD
# -----------------------------------------------------------------------------
if st.session_state.loaded:

    gdc.set_program(prog_id)
    df_psi = gdc.get_primary_sites(prog_id=prog_id, force=False, verbose=verbose)
    df_psi = make_streamlit_safe(df_psi)

    primary_site = df_psi.iloc[0].primary_site
    gdc.set_primary_site(primary_site=primary_site)    

    #---------- primary sites ----------------------------
    primary_sites = safe_unique_sorted(df_psi.primary_site)

    if len(primary_sites) == 0:
        st.warning("No primary sites found.")
        st.stop()

    col1, col2 = st.columns([4,10])  # adjust ratio as you like

    with col1:
        st.markdown("### Choose a primary site")

    with col2:
        selected_primary_site = st.selectbox(
            "",
            options=primary_sites,
            index=0,
            label_visibility="collapsed"  # removes empty label spacing
        )
    # -------------------------------------------------------------------------
    # FILTERED TABLES
    # -------------------------------------------------------------------------
    with st.spinner("Loading primary site data..."):
        df_cases, df_all_samples, df_all_mut, barcode_list = load_primary_site_data(selected_primary_site, verbose=False)

    with st.sidebar:
        st.subheader(f"Primary site: {selected_primary_site}")

        st.text(f"Cases {len(df_cases)}")
        st.text(f"Tumor samples {len(df_all_samples)}")
        st.text(f"Total mutations {len(df_all_mut)}")
        st.text(f"Patients with mutation {len(barcode_list)}")

    # -------------------------------------------------------------------------
    # BUILD MATRIX
    # -------------------------------------------------------------------------
    dfpiv = st_build_pivot_table(df_all_mut, min_barcodes=min_barcodes, min_genes=min_genes)
    # For the mutation matrix tab, I would also make the boolean matrix explicitly integer before display:
    if not dfpiv.empty:
        dfpiv = dfpiv.astype(int)

    df_gene_counts = summarize_mutations(df_all_mut)

    # -------------------------------------------------------------------------
    # TABS
    # -------------------------------------------------------------------------
    tab = st.radio("Main", ['Cases', 'Tumor Samples', 'Mutations', 'Mutation Matrix', 'Downloads'], horizontal=True)


    # -------------------------------------------------------------------------
    # TAB 1 - CASES
    # -------------------------------------------------------------------------
    if tab == "Cases":
        cols = ['case_id', 'psi_id', 'primary_site', 'disease_type',  'diagnoses', 
       'subtype_global', 'stage_ajcc', 'primary_diagnosis', 'tumor_grade',
        'tumor_stage', 'stage', 'tumor_class', 'histology',
       'subtype_tissue'] # 'stage_clin', 'figo_stage',
        
        if len(df_cases) > 200:
            df_cases2 = df_cases.head(200).copy()
            st.write(f"Cases #{len(df_cases)} limited to 200")
        else:
            df_cases2 = df_cases
            st.write(f"Cases #{len(df_cases)}")
        
        show_df(df_cases2[cols], height=800, key=f"cases_{selected_primary_site}")

    # -------------------------------------------------------------------------
    # TAB 2 - TUMOR SAMPLES
    # -------------------------------------------------------------------------
    elif tab == "Tumor Samples":
        
        if len(df_all_samples) > 200:
            df_all_samples2 = df_all_samples.head(200).copy()
            st.write(f"Tumor samples #{len(df_all_samples)} limited to 200")
        else:
            df_all_samples2 = df_all_samples
            st.write(f"Tumor samples #{len(df_all_samples)}")

        show_df(df_all_samples2, height=800, key=f"samples_{selected_primary_site}")

    # -------------------------------------------------------------------------
    # TAB 3 - MUTATIONS
    # -------------------------------------------------------------------------
    elif tab == "Mutations":

        subtab = st.radio("Main", ["Barplot: top Mutated Genes", "Mutated Genes", "Raw Mutation Rows"], horizontal=True)

        if subtab == "Barplot: top Mutated Genes":
            st.write(f"Most frequently mutated genes across filtered barcodes {dfpiv.shape[0]} samples and {dfpiv.shape[1]} genes")

            if dfpiv.shape[0] > 1:
                top_n = st.slider(
                    "Top N genes",
                    min_value=5,
                    max_value=min(100, max(5, dfpiv.shape[1] if not dfpiv.empty else 5)),
                    value=min(20, max(5, dfpiv.shape[1] if not dfpiv.empty else 5)),
                    step=5,
                )

                plot_top_mutated_genes(dfpiv, top_n=top_n)

        elif subtab == "Mutated Genes":
            st.write("Number of patients/barcodes mutated per gene")
            show_df(df_gene_counts, height=800, key=f"gene_counts_{selected_primary_site}")

        elif subtab == "Raw Mutation Rows":
            
            if len(df_all_mut) > 500:
                df_all_mut2 = df_all_mut.head(500).copy()
                st.write(f"Mutation rows #{len(df_all_mut)} limited to 500")
            else:
                df_all_mut2 = df_all_mut
                st.write(f"Mutation rows #{len(df_all_mut)}")

            cols = ['barcode_sample', 'symbol', 'refseq_mrna_id', 'entrez_gene_id', 'protein_mut',
                    'mutation_type', 'ref_allele', 'variant_allele',
                    'variant_type', 'chr', 'start', 'end', 'mutation_status',]

            show_df(df_all_mut2[cols], height=800, key=f"mut_rows_{selected_primary_site}")

    # -------------------------------------------------------------------------
    # TAB 4 - MUTATION MATRIX
    # -------------------------------------------------------------------------
    elif tab == "Mutation Matrix":

        subtab = st.radio("Main", ["Heatmap", "UMAP - cluster", "HDBSCAN - cluster"], horizontal=True)

        if dfpiv.empty:
            st.info("No mutation matrix available for this primary site.")
        else:

            if subtab == "Heatmap":
                n_samples, n_genes = dfpiv.shape
                st.subheader("Mutation matrix heatmap")
                title = f"Primary Site: '{selected_primary_site}' #{n_samples} samples and #{n_genes} genes"
                plot_heatmap(dfpiv, title)

            elif subtab == "UMAP - cluster":
                st.subheader("UMAP Clustering")

                n_samples = dfpiv.shape[0]

                k = st.slider(
                    "K (number of clusters)",
                    2,
                    min(15, n_samples),
                    min(8, n_samples)
                )

                plot_umap(dfpiv, k=k)

            elif subtab == "HDBSCAN - cluster":
                st.subheader("HDBSCAN Clustering")

                n_samples = dfpiv.shape[0]

                min_cluster_size = st.slider(
                    "Minimum cluster size",
                    min_value = 3,
                    max_value = min(15, n_samples),
                    value = 10
                )

                min_samples = st.slider(
                    "Minimum samples",
                    min_value = 3,
                    max_value = min(10, n_samples),
                    value = 3
                )

                plot_hdbscan(dfpiv, min_cluster_size=min_cluster_size, min_samples=min_samples)

    # -------------------------------------------------------------------------
    # TAB 5 - DOWNLOADS
    # -------------------------------------------------------------------------
    elif tab == "Downloads":
        st.write("Export filtered tables")

        st.download_button(
            "Download filtered cases TSV",
            data=df_cases.to_csv(sep='\t',index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_cases.tsv",
            mime="text/tab-separated-values",
        )

        st.download_button(
            "Download tumor samples TSV",
            data=df_all_samples.to_csv(sep='\t',index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_tumor_samples.tsv",
            mime="text/tab-separated-values",
        )

        st.download_button(
            "Download filtered mutations TSV",
            data=df_all_mut.to_csv(sep='\t',index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_mutations.tsv",
            mime="text/tab-separated-values",
        )

        if not dfpiv.empty:
            st.download_button(
                "Download mutation matrix TSV",
                data=dfpiv.reset_index().to_csv(sep='\t',index=False).encode("utf-8"),
                file_name=f"{prog_id}_{selected_primary_site}_mutation_matrix.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
            )
    show_profile_box()

else:
    st.info("Click **Load data** in the sidebar to start.")
    show_profile_box()

