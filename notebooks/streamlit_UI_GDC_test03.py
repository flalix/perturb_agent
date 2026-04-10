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
# streamlit run streamlit_UI_GDC_test03.py
#
# uv run streamlit run streamlit_UI_GDC_test03.py 
#
#============================================


import os, sys
from pprint import pprint
# from marshmallow import pprint, Schema, fields
import numpy as np
import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from collections import defaultdict

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
.block-container {
    max-width: 96%;
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



def show_df_AgGrid(df, height:int=500, page_size:int=25, key:str="grid"):
    if df is None or df.empty:
        st.info("Empty dataframe")
        return

    df = make_aggrid_safe(df)

    # st.write("shape:", df.shape)

    gb = GridOptionsBuilder.from_dataframe(df)

    gb.configure_default_column(sortable=True, filter=True, resizable=True)
    gb.configure_pagination(
        enabled=True,
        paginationAutoPageSize=False,
        paginationPageSize=page_size,
    )

    grid_options = gb.build()
    grid_options["pagination"] = True
    grid_options["paginationPageSize"] = page_size
    grid_options["domLayout"] = "normal"    

    if not isinstance(grid_options, dict):
        raise TypeError(f"grid_options must be dict, got {type(grid_options)}")

    AgGrid(
        df,
        gridOptions=grid_options,
        height=height,
        fit_columns_on_grid_load=True,
        allow_unsafe_jscode=False,
        enable_enterprise_modules=False,
        key=key,
    )


def show_df(df, height:int=500, page_size:int=25, key:str="grid"):
    show_df_AgGrid(df, height=height, page_size=page_size, key=key)


def show_df_html(df, height: int = 450):
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

def plot_heatmap(dfpiv: pd.DataFrame, figsize:tuple=(14, 10)):
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

    st.pyplot(cg.figure)
    plt.close(cg.figure)

def plot_umap(dfpiv: pd.DataFrame, k:int=8, figsize:tuple=(14, 10)):

    # Binary mutation matrix for Jaccard
    # Force numeric/binary and remove bad values
    data = (
        dfpiv
        .replace(["NaN", "nan", ""], np.nan)   # catch fake NaNs
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .astype(bool)
        .astype(int)
    )

    # Drop empty genes and empty samples
    data = data.loc[data.sum(axis=1) > 0, data.sum(axis=0) > 0]

    if data.empty:
        st.warning("No non-empty mutation matrix after filtering.")
        return
    
    X = data.to_numpy(dtype=np.uint8)

    n_samples = X.shape[0]

    if n_samples < 3:
        st.warning("Need at least 3 non-empty samples to compute UMAP + clustering.")
        return

    if k > X.shape[0]:
        st.warning(f"k={k} is larger than number of samples ({X.shape[0]}). Using k={X.shape[0]}.")
        k = X.shape[0]

    k = min(k, n_samples)

    n_neighbors = min(15, n_samples - 1)
    n_neighbors = max(2, n_neighbors)    

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=0.1,
        metric="jaccard",
        random_state=42,
        init="random",      # important
        output_dens=False   # avoid tuple output
    )

    embedding = reducer.fit_transform(X)

    if isinstance(embedding, tuple):
        st.write("embedding return as a tuple")
        embedding = embedding[0]

    embedding = np.array(embedding)

    good = np.isfinite(embedding).all(axis=1)
    embedding = embedding[good]


    labels = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(embedding)

    fig, ax = plt.subplots(figsize=figsize)

    # cmap = plt.cm.get_cmap("tab10", k)

    sc = plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=[colors[label] for label in labels],
        # c=labels,
        # cmap=cmap,
        s=20
    )

    ax.set_title(f"UMAP of Mutation Profiles (k={k})")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")

    handles, _ = sc.legend_elements()
    ax.legend(handles, [f"Cluster {i}" for i in range(k)], title="Groups", loc="best")

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
def build_pivot_table(df_all_mut: pd.DataFrame) -> pd.DataFrame:
    """
    Build barcode x gene boolean mutation matrix.
    """
    if df_all_mut is None or df_all_mut.empty:
        return pd.DataFrame()

    df_tmp = df_all_mut.copy()
    df_tmp["value"] = True

    dfpiv = df_tmp.pivot_table(
        index="barcode",
        columns="symbol",
        values="value",
        aggfunc="max",
        fill_value=False,
    )
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
    force = st.checkbox("Force rebuild", value=False)
    verbose = st.checkbox("Verbose", value=False)

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


    primary_sites = safe_unique_sorted(df_psi.primary_site)

    st.subheader("Primary site selection")

    if len(primary_sites) == 0:
        st.warning("No primary sites found.")
        st.stop()

    selected_primary_site = st.selectbox(
        "Choose a primary site",
        options=primary_sites,
        index=0,
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
    dfpiv = build_pivot_table(df_all_mut)
    # For the mutation matrix tab, I would also make the boolean matrix explicitly integer before display:
    if not dfpiv.empty:
        dfpiv = dfpiv.astype(int)

    df_gene_counts = summarize_mutations(df_all_mut)

    # -------------------------------------------------------------------------
    # TABS
    # -------------------------------------------------------------------------
    tab = st.radio("Main", ['Cases', 'Tumor Samples', 'Mutations', 'Mutation Matrix', 'Downloads'], horizontal=True)

    print(">>>", tab)

    # -------------------------------------------------------------------------
    # TAB 1 - CASES
    # -------------------------------------------------------------------------
    if tab == "Cases":
        st.write(f"Cases {len(df_cases)}")
        show_df(df_cases, height=450, key=f"samples_{selected_primary_site}")

    # -------------------------------------------------------------------------
    # TAB 2 - TUMOR SAMPLES
    # -------------------------------------------------------------------------
    elif tab == "Tumor Samples":
        st.write("Tumor samples linked to the selected primary site")
        show_df(df_all_samples, height=450, key=f"samples_{selected_primary_site}")

    # -------------------------------------------------------------------------
    # TAB 3 - MUTATIONS
    # -------------------------------------------------------------------------
    elif tab == "Mutations":

        subtab = st.radio("Main", ["Barplot: top Mutated Genes", "Mutated Genes", "Raw Mutation Rows"], horizontal=True)

        if subtab == "Barplot: top Mutated Genes":
            st.write("Most frequently mutated genes across filtered barcodes")

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
            show_df(df_gene_counts, height=450, key=f"gene_counts_{selected_primary_site}")

        elif subtab == "Raw Mutation Rows":
            st.write("Mutation rows after barcode filtering")
            show_df(df_all_mut, height=450, key=f"mut_rows_{selected_primary_site}")

    # -------------------------------------------------------------------------
    # TAB 4 - MUTATION MATRIX
    # -------------------------------------------------------------------------
    elif tab == "Mutation Matrix":

        subtab = st.radio("Main", ["Heatmap", "UMAP - cluster"], horizontal=True)

        if dfpiv.empty:
            st.info("No mutation matrix available for this primary site.")
        else:

            if subtab == "Heatmap":
                st.subheader("Mutation matrix heatmap")
                plot_heatmap(dfpiv)

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

    # -------------------------------------------------------------------------
    # TAB 5 - DOWNLOADS
    # -------------------------------------------------------------------------
    elif tab == "Downloads":
        st.write("Export filtered tables")

        st.download_button(
            "Download filtered cases CSV",
            data=df_cases.to_csv(index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_cases.csv",
            mime="text/csv",
        )

        st.download_button(
            "Download tumor samples CSV",
            data=df_all_samples.to_csv(index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_tumor_samples.csv",
            mime="text/csv",
        )

        st.download_button(
            "Download filtered mutations CSV",
            data=df_all_mut.to_csv(index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_mutations.csv",
            mime="text/csv",
        )

        if not dfpiv.empty:
            st.download_button(
                "Download mutation matrix CSV",
                data=dfpiv.reset_index().to_csv(index=False).encode("utf-8"),
                file_name=f"{prog_id}_{selected_primary_site}_mutation_matrix.csv",
                mime="text/csv",
                use_container_width=True,
            )

else:
    st.info("Click **Load data** in the sidebar to start.")