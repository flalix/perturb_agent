# streamlit_gdc_tcga.py

#=============== to run =====================
#
# export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# streamlit run streamlit_UI_GDC_test03.py
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
import tempfile

# Project root (works locally + Render)
# ROOT = Path(__file__).resolve().parent
ROOT = Path('/opt/render/project/src/')

# Add src to path
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

print("ROOT:", ROOT)
print("SRC added:", SRC)

from libs.tcga_gdc_lib import *
from libs.Basic import *
from libs.calc_degs_lib import CALC_DEGS


# GDRIVE_FOLDER_ID = '1Tp4GONa9Qu1gySZaxEK2izwEZJ4E_xKr'

# root_data = os.path.join(ROOT, "data/tcga")


root_data = Path('/opt/render/project/src/storage/')
print("root_data:", root_data)

gdc = GDC(root_data=root_data)

pid = 'TCGA'
verbose = True
colors = ['red', 'green', 'blue', 'orange', 'pink', 'purple', 'black', 'cyan', 'tomato', 'lime', 'magenta', 'yellow',
          'gray', 'brown', 'olive', 'navy', 'teal', 'maroon', 'silver']


# Optional
import plotly.express as px

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

def show_df(df, height:int=500, page_size:int=50):
    df = make_streamlit_safe(df)

    gb = GridOptionsBuilder.from_dataframe(df)

    gb.configure_default_column(
        sortable=True,
        filter=True,
        resizable=True,
    )

    gb.configure_pagination(
        enabled=True,
        paginationAutoPageSize=False,
        paginationPageSize=page_size,
    )

    grid_options = gb.build()

    AgGrid(df,
        gridOptions=grid_options,
        height=height,
        fit_columns_on_grid_load=True,
        allow_unsafe_jscode=False)

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

def show_df_old(df, height:int=450):
    df = make_streamlit_safe(df)

    try:
        st.dataframe(df, use_container_width=True, height=height)
    except TypeError:
        st.dataframe(df, height=height)


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



# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
# hash error: @st.cache(show_spinner=True)
def load_program_data(pid:str, force:bool=False, verbose:bool=False):
    df_all_cases, df_all_samples, df_all_mutations = gdc.loop_program_psi_samples(
        program=pid,
        force=force,
        verbose=verbose,
    )

    df_all_cases = make_streamlit_safe(df_all_cases)
    df_all_samples = make_streamlit_safe(df_all_samples)
    df_all_mutations = make_streamlit_safe(df_all_mutations)

    return df_all_cases, df_all_samples, df_all_mutations


# hash error: @st.cache(show_spinner=False)
def build_mutation_matrix(df_mut: pd.DataFrame) -> pd.DataFrame:
    """
    Build barcode x gene boolean mutation matrix.
    """
    if df_mut is None or df_mut.empty:
        return pd.DataFrame()

    df_tmp = df_mut.copy()
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
def summarize_mutations(df_mut: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize mutated genes by number of patients/barcodes.
    """
    if df_mut is None or df_mut.empty:
        return pd.DataFrame(columns=["symbol", "n_patients_mutated"])

    df = ( df_mut.groupby("symbol")["barcode"]
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
with st.sidebar:
    st.header("Controls")

    pid = st.text_input("Program", value="TCGA")
    force = st.checkbox("Force rebuild", value=False)
    verbose = st.checkbox("Verbose", value=False)

    run = st.button("Load data")


if "loaded" not in st.session_state:
    st.session_state.loaded = False

if run:
    st.session_state.loaded = True


# -----------------------------------------------------------------------------
# SESSION STATE
# -----------------------------------------------------------------------------
if "loaded" not in st.session_state:
    st.session_state.loaded = False

if run:
    st.session_state.loaded = True


# -----------------------------------------------------------------------------
# MAIN LOAD
# -----------------------------------------------------------------------------
if st.session_state.loaded:
    try:
        with st.spinner(f"Loading all data for program: {pid}"):
            df_all_cases, df_all_samples, df_all_mutations = load_program_data(
                pid=pid, force=False, verbose=False
            )

    except Exception as e:
        df_all_cases = pd.DataFrame()
        df_all_samples = pd.DataFrame()
        df_all_mutations = pd.DataFrame()


        st.error(f"Failed to load program data: {e}")
        st.stop()

    # -------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # -------------------------------------------------------------------------
    primary_sites = safe_unique_sorted(df_all_cases["primary_site"])
    symbols = safe_unique_sorted(df_all_mutations["symbol"]) if not df_all_mutations.empty else []

    st.subheader("Program summary")

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Primary sites", len(primary_sites))
    c2.metric("Cases", len(df_all_cases))
    c3.metric("Samples", len(df_all_samples))
    c4.metric("Annotated mutations", len(df_all_mutations))
    c5.metric("Different genes", len(symbols))

    with st.expander("Raw summary text"):
        st.code(
            "\n".join([
                f"Interfacing GDC {pid} data, one gathered:",
                f"\t- {len(primary_sites)} primary sites.",
                f"\t- {len(df_all_cases)} cases.",
                f"\t- {len(df_all_samples)} samples.",
                f"\t- {len(df_all_mutations)} annotated mutations.",
                f"\t- {len(symbols)} different genes.",
            ])
        )

    # -------------------------------------------------------------------------
    # PRIMARY SITE SELECTION
    # -------------------------------------------------------------------------
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
    df_cases, df_samples, df_mut, barcode_list = gdc.get_filtered_tables(selected_primary_site)

    df_cases = make_streamlit_safe(df_cases)
    df_samples = make_streamlit_safe(df_samples)
    df_mut = make_streamlit_safe(df_mut)    

    # -------------------------------------------------------------------------
    # LOCAL SUMMARY
    # -------------------------------------------------------------------------
    st.subheader(f"Filtered summary: {selected_primary_site}")

    cc1, cc2, cc3, cc4 = st.columns(4)
    cc1.metric("Filtered cases", len(df_cases))
    cc2.metric("Tumor samples", len(df_samples))
    cc3.metric("Total mutations", len(df_mut))
    cc4.metric(
        "Patients with mutation rows",
        df_mut["barcode"].nunique() if not df_mut.empty and "barcode" in df_mut.columns else 0
    )

    # -------------------------------------------------------------------------
    # BUILD MATRIX
    # -------------------------------------------------------------------------
    dfpiv = build_mutation_matrix(df_mut)
    # For the mutation matrix tab, I would also make the boolean matrix explicitly integer before display:
    if not dfpiv.empty:
        dfpiv = dfpiv.astype(int)

    df_gene_counts = summarize_mutations(df_mut)

    # -------------------------------------------------------------------------
    # TABS
    # -------------------------------------------------------------------------
    tab = st.radio("Main", ['Cases', 'Tumor Samples', 'Mutations', 'Mutation Matrix', 'Downloads'], horizontal=True)

    # -------------------------------------------------------------------------
    # TAB 1 - CASES
    # -------------------------------------------------------------------------
    if tab == "Cases":
        st.write("Filtered case table")
        show_df(df_cases, height=450)

    # -------------------------------------------------------------------------
    # TAB 2 - TUMOR SAMPLES
    # -------------------------------------------------------------------------
    elif tab == "Tumor Samples":
        st.write("Tumor samples linked to the selected primary site")
        show_df(df_samples, height=450)

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
            show_df(df_gene_counts, height=450)

        elif subtab == "Raw Mutation Rows":
            st.write("Mutation rows after barcode filtering")
            show_df(df_mut, height=450)

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
            file_name=f"{pid}_{selected_primary_site}_cases.csv",
            mime="text/csv",
        )

        st.download_button(
            "Download tumor samples CSV",
            data=df_samples.to_csv(index=False).encode("utf-8"),
            file_name=f"{pid}_{selected_primary_site}_tumor_samples.csv",
            mime="text/csv",
        )

        st.download_button(
            "Download filtered mutations CSV",
            data=df_mut.to_csv(index=False).encode("utf-8"),
            file_name=f"{pid}_{selected_primary_site}_mutations.csv",
            mime="text/csv",
        )

        if not dfpiv.empty:
            st.download_button(
                "Download mutation matrix CSV",
                data=dfpiv.reset_index().to_csv(index=False).encode("utf-8"),
                file_name=f"{pid}_{selected_primary_site}_mutation_matrix.csv",
                mime="text/csv",
                use_container_width=True,
            )

else:
    st.info("Click **Load data** in the sidebar to start.")