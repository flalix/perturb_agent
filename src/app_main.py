#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Updated  on 2026/05/23
# @author: Flavio Lichtenstein
# @local: Home sweet home

# =============== to run =====================
#
# export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
#
# uv run streamlit run app.py
#
# ============================================

import sys
import os
from urllib import response
import yaml
import numpy as np

# ----------- fix incompatibilities ---------------------
import pandas as pd

setattr(pd.Series, "iteritems", pd.Series.items)
setattr(pd.DataFrame, "iteritems", pd.DataFrame.items)

import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, JsCode

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import seaborn as sns


"""
/opt/render/project/src/
├── src/
│   ├── app_main.py
│   ├── libs/
│   └── styles/
└── storage/
    └── data/
        ├─── gdc_progams.txt
        ├─── colab
        ├─── TCGA/
        └─── vector_store/

"""

ROOT0 = Path(os.environ.get("ROOT0", "/home/flavio/uv/perturb_agent"))
ROOT_DATA = Path(os.environ.get("ROOT_DATA", "/home/flavio/uv/perturb_agent/data"))
ROOT_SRC = ROOT0 / "src"
ROOT_COLAB = ROOT_DATA / "colab"
ROOT_CSS = ROOT_SRC / "styles"

if str(ROOT_SRC) not in sys.path:
    sys.path.insert(0, str(ROOT_SRC))

print("ROOT0:", ROOT0)
print("ROOT_SRC:", ROOT_SRC)
print("ROOT_DATA:", ROOT_DATA)
print("ROOT_COLAB:", ROOT_COLAB)

if str(ROOT_SRC) not in sys.path:
    sys.path.insert(0, str(ROOT_SRC))

from libs.enricher_lib import enricheR
from libs.tcga_gdc_lib import GDC
from libs.config_lib import Config
from libs.dashcyto_lib import DASH_CYTO
from libs.reactome_lib import Reactome
from libs.calc_degs_lib import CALC_DEGS
from project_context_MTD import load_project_context

try:
    with open('params.yml', 'r') as file:
        dic_yml = yaml.safe_load(file)
except:
    print("Error loading params.yml")
    print(">>>", os.listdir("../"))
    raise Exception("\n------------ stop --------------\n")

PROG_ID = 'TCGA'
PSI_ID = 'TCGA-BRCA'

ctx = load_project_context(
    dic_yml=dic_yml,
    PSI_ID=PSI_ID,
    i_project=0,
)

colors = ctx.colors

gdc = GDC(root0=ROOT0, root0_data=ROOT_DATA)

ROOT_OWL = ROOT_COLAB / 'owl'
ROOT_REACTOME = ROOT_COLAB / 'reactome'

rea = Reactome(root_owl=ROOT_OWL, root_reactome=ROOT_REACTOME)

verbose = False

# -----------------------------------------------------------------------------
# PAGE
# -----------------------------------------------------------------------------
st.set_page_config(page_title="GDC / TCGA Explorer", layout="wide")


def load_css(fname: str):
    filename = ROOT_CSS / fname

    with open(filename) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


load_css("main.css")

st.title("GDC / TCGA Explorer")
st.caption("Explore cases, tumor samples, and mutation matrices by primary site")


force = False
verbose = False

min_barcodes = 3
min_genes = 5


# --- FOOTER / BOX BELOW ALL TABS ---
def show_profile_box():
    st.markdown(
        """
        <div style="display:flex; justify-content:left; margin-top:40px; margin-bottom:40px;">
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
                <p>
                    🔗 <a href="https://www.linkedin.com/in/flaviolichtenstein/"
                    target="_blank"
                    style="color:#0A66C2; text-decoration:none;">
                    LinkedIn Profile</a>
                </p>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def make_df_streamlit_safe(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df

    df_safe = df.copy()

    # Make labels plain Python strings
    df_safe.index = df_safe.index.map(str)
    df_safe.columns = [str(c) for c in df_safe.columns]

    # Flatten MultiIndex columns if present
    if isinstance(df_safe.columns, pd.MultiIndex):
        df_safe.columns = [" | ".join(map(str, col)).strip() for col in df_safe.columns.to_flat_index()]

    # Convert every column to a Streamlit-safe pandas dtype
    for col in df_safe.columns:
        s = df_safe[col]

        # nullable / pyarrow / string extension -> Python objects
        try:
            if pd.api.types.is_string_dtype(s.dtype):
                df_safe[col] = s.astype("object")
            elif str(s.dtype).lower().find("pyarrow") >= 0:
                df_safe[col] = s.astype("object")
            else:
                df_safe[col] = s
        except Exception:
            df_safe[col] = s.astype("object")

        # Replace pandas NA/NaT with None so Arrow does not try fancy typing
        df_safe[col] = df_safe[col].where(pd.notna(df_safe[col]), None)

    return df_safe

def show_df( df, height: int | None = None, page_size: int = 25, key: str = "grid", selectable: bool = False):
    return show_df_AgGrid(df, height=height, page_size=page_size, key=key, selectable=selectable)


def show_df_AgGrid( df, height: int | None = None, page_size: int = 25, key: str = "grid", selectable: bool = False):

    if df is None or df.empty:
        st.info("Empty dataframe")
        return

    df = make_df_streamlit_safe(df)

    for col in df.columns:
        converted = pd.to_numeric(df[col], errors="coerce")

        # Only replace if conversion produced at least some real numbers
        if converted.notna().sum() > 0 and converted.notna().sum() >= df[col].notna().sum() * 0.8:
            df[col] = converted

    float_cols = df.select_dtypes(include="number").columns

    if height is None:
        height = 400

    gb = GridOptionsBuilder.from_dataframe(df)

    gb.configure_default_column(
        sortable=True,
        minWidth=150,
        filter=True,
        resizable=True,
    )

    pval_keywords = ["pvalue", "p_value", "pval", "p_val", "padj", "fdr", "qvalue", "q_value"]

    float_cols = df.select_dtypes(include=["float", "float64", "float32"]).columns

    pval_formatter = JsCode("""
    function(params) {
        if (params.value === null || params.value === undefined || isNaN(params.value)) {
            return "";
        }
        return Number(params.value).toExponential(2);
    }
    """)

    float_formatter = JsCode("""
    function(params) {
        if (params.value === null || params.value === undefined || isNaN(params.value)) {
            return "";
        }
        return Number(params.value).toFixed(3);
    }
    """)

    for col in float_cols:
        col_l = col.lower()
        is_pval_col = any(k in col_l for k in pval_keywords)

        gb.configure_column(
            col,
            type=["numericColumn"],
            valueFormatter=pval_formatter if is_pval_col else float_formatter,
        )


    gb.configure_pagination(
        enabled=True,
        paginationAutoPageSize=False,
        paginationPageSize=page_size,
    )

    gb.configure_grid_options(
        pagination=True,
        paginationPageSize=page_size,
        suppressPaginationPanel=False,
        domLayout="autoHeight",
    )

    if selectable:
        gb.configure_selection(
            selection_mode="single",
            use_checkbox=True,
        )
    
    grid_options = gb.build()

    custom_css = {
        ".ag-root-wrapper": {
            "border-bottom": "1px solid #ddd !important",
        },
        ".ag-paging-panel": {
            "display": "flex !important",
            "visibility": "visible !important",
            "height": "45px !important",
            "min-height": "45px !important",
            "overflow": "visible !important",
        },
    }

    response = AgGrid(
        df,
        gridOptions=grid_options,
        height=height,
        allow_unsafe_jscode=True,   # important for valueFormatter strings
        enable_enterprise_modules=False,
        reload_data=True,
        key=key,
        update_mode=(
            GridUpdateMode.SELECTION_CHANGED
            if selectable
            else GridUpdateMode.NO_UPDATE
        ),
        theme="alpine",  # or "streamlit"
        custom_css=custom_css,
    )

    st.markdown("<div style='height:30px;'></div>", unsafe_allow_html=True)

    if selectable:
        selected = response.get("selected_rows", [])

        if isinstance(selected, list) and len(selected) > 0:
            return selected[0]

        if isinstance(selected, pd.DataFrame) and not selected.empty:
            return selected.iloc[0].to_dict()

        return None

    return response

def show_selectable_df(df, height: float | None = None, key=None):
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_selection(selection_mode="single", use_checkbox=True)
    # gb.configure_grid_options(domLayout="normal")

    response = AgGrid(
        df,
        gridOptions=gb.build(),
        height=height,
        key=key,
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        allow_unsafe_jscode=False,
    )

    selected = response.get("selected_rows", [])
    if selected:
        return selected[0]

    return None

def show_df_html(df, height: int = 600):
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


def plot_top_mutated_genes(dfpiv: pd.DataFrame, top_n: int = 20, figsize=(12, 6)):
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


def plot_heatmap_gen(dfpiv: pd.DataFrame, title: str = "", figsize: tuple = (14, 10)):
    # Ensure numeric + binary (important for Jaccard)
    data = dfpiv.fillna(0).astype(int)

    cg = sns.clustermap(
        data, metric="jaccard", method="average", figsize=figsize, cmap="viridis", cbar=False
    )

    if title:
        cg.figure.suptitle(title, y=1.02)

    st.pyplot(cg.figure)
    plt.close(cg.figure)


def plot_heatmap_exp(dff: pd.DataFrame, normal_samples: list, tumor_samples: list, title: str = "", figsize: tuple = (14, 10)):

    cg = gdc.plot_heatmap_expression(dff, normal_samples, tumor_samples, title)

    st.pyplot(cg.figure)
    plt.close(cg.figure)


def plot_umap_exp(dff: pd.DataFrame, samples: list, which_samples:str, n_clusters: int , n_neighbors: int, title: str = "", figsize: tuple = (14, 10)):

    fig_umap, ax_umap, df_umap = gdc.plot_umap_expression(
        dff=dff,
        samples=samples,
        which_samples=which_samples,
        title=title,
        figsize=figsize,
        n_neighbors=n_neighbors,
        n_clusters=n_clusters
    )
    st.pyplot(fig_umap)
    plt.close(fig_umap)

    return df_umap


def plot_umap_gen(dfpiv: pd.DataFrame, k: int = 8, figsize: tuple = (14, 10)):

    fig, _, _ = gdc.plot_UMAP(dfpiv=dfpiv, k=k, figsize=figsize)

    if fig:
        st.pyplot(fig)
        plt.close(fig)


def plot_hdbscan_gen(dfpiv: pd.DataFrame, min_cluster_size: int = 10, min_samples: int = 3, figsize: tuple = (14, 10)):

    fig, embedding, labels, d = gdc.plot_HDBSCAN(
        dfpiv=dfpiv, min_cluster_size=min_cluster_size, min_samples=min_samples, figsize=figsize
    )

    if fig:
        st.pyplot(fig)
        plt.close(fig)

    return fig, embedding, labels



def load_disease(PSI_ID:str, root_disease:Path, LFC_cut:float=1, lfc_FDR_cut:float=0.05, verbose:bool=False):

    cfg = Config(root0=ROOT0, root_disease=root_disease, disease=ctx.disease, case_list=ctx.case_list)

    mtd = enricheR(disease=PSI_ID, gene_protein=ctx.gene_protein, s_omics=ctx.s_omics, project=ctx.project, s_project=ctx.s_project, 
            root0=ROOT0, root0_data=ROOT_DATA, prog_id=PROG_ID, psi_id=PSI_ID,
            case_list=ctx.case_list, dic_case_list=ctx.dic_case_list, 
            has_age=ctx.has_age, has_gender=ctx.has_gender, exp_normalization=ctx.exp_normalization,
            std_filename=ctx.std_filename, std_filename_list=ctx.std_filename_list,
            geneset_num=0, ptw_min_num_of_degs_cut=ctx.ptw_min_num_of_degs_cut,
            tolerance_pPMI=ctx.tolerance_pPMI, s_pathw_enrichm_method=ctx.s_pathw_enrichm_method,
            LFC_cut_inf=ctx.LFC_cut_inf, fdr_ptw_cutoff_list=ctx.fdr_ptw_cutoff_list,
            num_of_genes_list=ctx.num_of_genes_list, lfc_list=ctx.lfc_list, fdr_list=ctx.fdr_list, 
            min_lfc_modulation=ctx.min_lfc_modulation, type_sat_ptw_index=ctx.type_sat_ptw_index,
            saturation_lfc_param=ctx.saturation_lfc_param, enr_db_list=ctx.enr_db_list, pPMI_normalized=ctx.pPMI_normalized)


    mtd.cfg.set_default_best_lfc_cutoff(mtd.normalization, LFC_cut=LFC_cut, lfc_FDR_cut=lfc_FDR_cut)

    return mtd

def classify_biotype(b):
    if b == "protein_coding":
        return "protein_coding"
    elif "pseudogene" in b:
        return "pseudogene"
    elif b == "lncRNA":
        return "lncRNA"
    elif b.startswith("IG_") or b.startswith("TR_"):
        return "immune_receptor"
    elif b in ["miRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "vault_RNA", "scRNA"]:
        return "small_RNA"
    elif b == "TEC":
        return "TEC"
    elif b.startswith("Mt_"):
        return "mitochondrial_RNA"
    else:
        return "other"

def group_biotypes(dflfc: pd.DataFrame) -> pd.DataFrame:

    df_de = dflfc.copy()

    df_de["biotype_class"] = df_de["biotype"].apply(classify_biotype)
    df_de['direction'] = df_de['lfc'].apply(lambda x: "up" if x > 1 else ("down" if x < -1 else "neutral"))

    df_grouped = (
        df_de
        .groupby(["direction", "biotype_class"])
        .size()
        .reset_index(name="n")
        .sort_values(["direction", "n"], ascending=[False, False])
    )

    return df_grouped
    
def calc_enrichment_analysis(verbose: bool = False, force: bool = False):

    _, _, dflfc = mtd.list_of_degs(save_file=False, force=False, prompt_verbose=False, verbose=verbose)

    df2 = dflfc[dflfc.biotype.isin(mtd.biotype_annot)]
    df2 = df2.sort_values('fdr', ascending=True)
    
    degs_to_reactome = df2.symbol.to_list()[:mtd.MAX_GENES_ENRICHR_SAFE]

    mtd.calc_EA_dataset_symbol(list(degs_to_reactome), calc_many_sig=False, default=False, force=force, verbose=verbose)
    return mtd.df_enr0


# prog_list = gdc.get_gdc_progams(force=False, verbose=verbose)

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
# hash error: @st.cache(show_spinner=True)


@st.cache_data(show_spinner=False)
def load_primary_site_data(verbose: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list]:

    df_cases, df_all_samples, df_all_mut, barcode_list = gdc.get_filtered_tables(
        sample_type_term="tumor", verbose=verbose
    )

    return (
        make_df_streamlit_safe(df_cases),
        make_df_streamlit_safe(df_all_samples),
        make_df_streamlit_safe(df_all_mut),
        barcode_list,
    )


# hash error: @st.cache(show_spinner=False)
@st.cache_data(show_spinner=False)
def st_build_pivot_table(
    df_all_mut: pd.DataFrame, min_barcodes: int = 3, min_genes: int = 5
) -> pd.DataFrame:
    """
    Build barcode x gene boolean mutation matrix.
    """
    if df_all_mut is None or df_all_mut.empty:
        return pd.DataFrame()

    dfpiv = gdc.build_pivot_table(df_all_mut, min_barcodes=min_barcodes, min_genes=min_genes)
    dfpiv = make_df_streamlit_safe(dfpiv)

    return dfpiv.sort_index(axis=0).sort_index(axis=1)


# hash error:@st.cache(show_spinner=False)
def summarize_mutations(df_all_mut: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize mutated genes by number of patients/barcodes.
    """
    if df_all_mut is None or df_all_mut.empty:
        return pd.DataFrame(columns=["symbol", "n_patients_mutated"])

    df = (
        df_all_mut.groupby("symbol")["barcode"]
        .nunique()
        .reset_index(name="n_patients_mutated")
        .sort_values("n_patients_mutated", ascending=False)
        .reset_index(drop=True)
    )

    df = make_df_streamlit_safe(df)
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
    prog_id = "TCGA"

    import streamlit

    st.text(f"Streamlit {streamlit.__version__}")

    st.text(f"Program {prog_id}")
    # force = st.checkbox("Force rebuild", value=False)
    # verbose = st.checkbox("Verbose", value=False)

    load_clicked = st.button("Load data", use_container_width=True)

    st.text(f"ROOT0 {ROOT0}")
    st.text(f"ROOT_DATA {ROOT_DATA}")
    st.text(f"ROOT_SRC {ROOT_SRC}")
    st.text(f"ROOT_COLAB {ROOT_COLAB}")
    st.text(f"ROOT_CSS {ROOT_CSS}")

    if load_clicked:
        st.session_state.loaded = True


# -----------------------------------------------------------------------------
# MAIN LOAD
# -----------------------------------------------------------------------------
if st.session_state.loaded:
    gdc.set_program(prog_id)
    df_psi = gdc.get_primary_sites(prog_id=prog_id, force=False, verbose=verbose)
    df_psi = make_df_streamlit_safe(df_psi)

    # ---------- primary sites ----------------------------
    primary_sites = [row.psi_id + " - " + row.primary_site for i, row in df_psi.iterrows()]
    primary_sites.sort()

    if len(primary_sites) == 0:
        st.warning("No primary sites found.")
        st.stop()

    col1, col2 = st.columns([4, 10])  # adjust ratio as you like

    with col1:
        st.markdown("### Choose a primary site")

    with col2:
        selected_primary_site = st.selectbox(
            "Primary site",
            options=primary_sites,
            index=0,
            label_visibility="collapsed",
        )
    # -------------------------------------------------------------------------
    # FILTERED TABLES
    # -------------------------------------------------------------------------
    psi_id = str(selected_primary_site).split(" - ")[0]
    gdc.set_primary_site(psi_id=psi_id)
    with st.spinner("Loading primary site data..."):
        df_cases, df_all_samples, df_all_mut, barcode_list = load_primary_site_data(verbose=False)

        st.write(">>> cases", len(df_cases))
        st.write(">>> samples", len(df_all_samples))
        st.write(">>> mutations", len(df_all_mut))
        st.write(">>> barcodes", len(barcode_list))

        cases_ids = np.unique(df_all_samples.case_id)
        df_cases = df_cases[df_cases["case_id"].isin(cases_ids)].copy()
        df_cases.reset_index(drop=True, inplace=True)

        method = "deseq2"
        mtd = load_disease(PSI_ID=psi_id, root_disease=gdc.root_disease,  LFC_cut=1, lfc_FDR_cut=0.05, verbose=False)

        icase=0
        case = ctx.case_list[icase]

        ret, degs, _, dflfc = mtd.open_case(case=case, prompt_verbose=False, verbose=False)

  
    # --------------------------------------------------
    # SESSION STATE (RIGHT HERE ✅)
    # --------------------------------------------------
    if "case_idx" not in st.session_state:
        st.session_state.case_idx = 0

    case_ids = df_cases["case_id"].dropna().unique().tolist()

    # reset if dataset changed
    if "case_ids_prev" not in st.session_state:
        st.session_state.case_ids_prev = []

    if st.session_state.case_ids_prev != case_ids:
        st.session_state.case_idx = 0
        st.session_state.case_ids_prev = case_ids

    def set_case_from_id(case_id: str):
        if case_id in case_ids:
            st.session_state.case_idx = case_ids.index(case_id)

    def current_case_id():
        if not case_ids:
            return None
        st.session_state.case_idx = max(0, min(st.session_state.case_idx, len(case_ids) - 1))
        return case_ids[st.session_state.case_idx]

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
        st.write("dfpiv is empty.")
    else:
        dfpiv = dfpiv.astype(int)

    df_gene_counts = summarize_mutations(df_all_mut)

    # -------------------------------------------------------------------------
    # TABS
    # -------------------------------------------------------------------------
    (
        tab_cases,
        tab_samples,
        tab_head_mutations,
        tab_head_cluster,
        tab_head_diff_exp,
        tab_donwload,
    ) = st.tabs( ["Cases", "Tumor Samples", "Mutations", "Clusterization", "Diff.Expression", "Downloads"] )

    # -------------------------------------------------------------------------
    # TAB 1 - CASES xxxx
    # -------------------------------------------------------------------------
    with tab_cases:
        st.info(f"Selected case_id: {current_case_id()}")

        cols = ["case_id", "disease_type", "diagnoses", "subtype_global", "subtype_tissue", 
                "primary_diagnosis", "tumor_grade", "tumor_stage", "stage", "tumor_class", "histology",]

        if len(df_cases) > 200:
            df_cases2 = df_cases.head(200).copy()
            st.write(f"Cases #{len(df_cases)} limited to 200")
        else:
            df_cases2 = df_cases
            st.write(f"Cases #{len(df_cases)}")

        selected_row = show_df(df_cases2[cols], height=None, selectable=True, key="cases")

        if selected_row is not None:
            set_case_from_id(selected_row["case_id"])

    # -------------------------------------------------------------------------
    # TAB 2 - TUMOR SAMPLES
    # -------------------------------------------------------------------------
    with tab_samples:

        if df_all_samples.empty:
            st.warning("No tumor samples found.")
        else:
            tab_one_case, tab_statistics = st.tabs(
                ["One case", "Statistics"]
            )

            with tab_one_case:        
                selected_case_id = current_case_id()

                col1, col2, col3 = st.columns([2, 1, 1])

                with col1:
                    st.write(f"Selected case: {selected_case_id}")

                with col2:
                    if st.button("Previous"):
                        st.session_state.case_idx -= 1
                        st.experimental_rerun()

                with col3:
                    if st.button("Next"):
                        st.session_state.case_idx += 1
                        st.experimental_rerun()

                if selected_case_id is None:
                    st.warning("Select a case first.")

                    if len(df_all_samples) > 200:
                        df_all_samples2 = df_all_samples.head(200).copy()
                        st.write(f"Tumor samples #{len(df_all_samples)} limited to 200")
                    else:
                        df_all_samples2 = df_all_samples
                        st.write(f"Tumor samples #{len(df_all_samples)}")

                    show_df(df_all_samples2, height=None, key=f"samples_{psi_id}")
                else:
                    df_case_samples = df_all_samples[df_all_samples["case_id"] == selected_case_id].copy()

                    df_grouped = (
                        df_case_samples.groupby(
                            ["barcode_sample", "sample_type", "data_type", "data_format"], dropna=False
                        )
                        .size()
                        .reset_index(name="n")
                    )

                    show_df(df_grouped, height=None, key=f"samples_{selected_primary_site}_samples_case",)

            with tab_statistics:
                    df_grouped = (
                        df_all_samples.groupby(
                            ["sample_type", "data_type", "data_format"], dropna=False
                        )
                        .size()
                        .reset_index(name="n")
                    )

                    show_df(df_grouped, height=None, key=f"samples_{selected_primary_site}_samples",)
    # -------------------------------------------------------------------------
    # TAB 3 - MUTATIONS
    # -------------------------------------------------------------------------
    with tab_head_mutations:
        # subtab = st.radio("Main", ["Barplot: top Mutated Genes", "Mutated Genes", "Raw Mutation Rows"], horizontal=True)
        tab_mut_barplot, tab_mut_genes, tab_mut_raw_data = st.tabs(
            ["Barplot", "Mutated Genes", "Raw Mutation Rows"]
        )

        with tab_mut_barplot:
            st.write(
                f"Most frequently mutated genes across filtered barcodes {dfpiv.shape[0]} samples and {dfpiv.shape[1]} genes"
            )

            if dfpiv.shape[0] > 1:
                top_n = st.slider(
                    "Top N genes",
                    min_value=5,
                    max_value=min(100, max(5, dfpiv.shape[1] if not dfpiv.empty else 5)),
                    value=min(20, max(5, dfpiv.shape[1] if not dfpiv.empty else 5)),
                    step=5,
                )

                plot_top_mutated_genes(dfpiv, top_n=top_n)

        with tab_mut_genes:
            st.write("Number of patients/barcodes mutated per gene")
            show_df(df_gene_counts, height=None, key=f"gene_counts_{psi_id}")

        with tab_mut_raw_data:
            if df_all_mut.empty:
                st.write("No mutation data available.")
            else:
                if len(df_all_mut) > 500:
                    df_all_mut2 = df_all_mut.head(500).copy()
                    st.write(f"Mutation rows #{len(df_all_mut)} limited to 500")
                else:
                    df_all_mut2 = df_all_mut
                    st.write(f"Mutation rows #{len(df_all_mut)}")

                cols = [
                    "barcode_sample",
                    "symbol",
                    "refseq_mrna_id",
                    "entrez_gene_id",
                    "protein_mut",
                    "mutation_type",
                    "ref_allele",
                    "variant_allele",
                    "variant_type",
                    "chr",
                    "start",
                    "end",
                    "mutation_status",
                ]

                show_df(df_all_mut2[cols], height=None, key=f"mut_rows_{psi_id}")

    # -------------------------------------------------------------------------
    # TAB 4 - MUTATION MATRIX
    # -------------------------------------------------------------------------
    with tab_head_cluster:
        tab_heatmap_gen, tab_umap_gen, tab_hdbscan_gen = st.tabs(["Heatmap", "UMAP - cluster", "HDBSCAN - cluster"])

        if dfpiv.empty:
            st.info("No mutation matrix available for this primary site.")
        else:
            with tab_heatmap_gen:
                n_samples, n_genes = dfpiv.shape
                st.subheader("Mutation matrix heatmap")
                title = f"Primary Site: '{selected_primary_site}' #{n_samples} samples and #{n_genes} genes"
                plot_heatmap_gen(dfpiv, title)

            with tab_umap_gen:
                st.subheader("UMAP Clustering")

                n_samples = dfpiv.shape[0]

                k = st.slider("K (number of clusters)", 2, min(15, n_samples), min(8, n_samples))
                plot_umap_gen(dfpiv, k=k)

            with tab_hdbscan_gen:
                st.subheader("HDBSCAN Clustering")

                n_samples = dfpiv.shape[0]

                min_cluster_size2 = st.slider(
                    "Minimum cluster size", min_value=3, max_value=min(15, n_samples), value=5
                )

                min_samples2 = st.slider(
                    "Minimum samples", min_value=3, max_value=min(10, n_samples), value=3
                )

                fig, embedding, labels = plot_hdbscan_gen(
                    dfpiv, min_cluster_size=min_cluster_size2, min_samples=min_samples2
                )                

    with tab_head_diff_exp:
        tab_degs, tab_echo, tab_biotypes, tab_lnc, tab_heatmap_exp, tab_umap_exp, tab_hdbscan_exp, tab_enrich = \
        st.tabs(["DEGs", "Echo", "Biotypes", "Non-Coding", "Heatmap", "UMAP - cluster", "HDBSCAN - cluster", "Enrichment Analysis"])

        msg = ''
        if dflfc.empty:
            st.write("No differentially expressed genes found.")
        else:

            cdegs = CALC_DEGS(root_src=ROOT_SRC, run_conda=False)
            gdc.cdegs = cdegs

            df_tumor, df_normal, msg = gdc.get_tumor_normal_tables(verbose=False)

            if df_tumor.empty:
                dff = pd.DataFrame()
                normal_samples, tumor_samples = [], []
            else:
                df_counts, df_meta = cdegs.build_counts_and_metadata(
                            df_tumor=df_tumor,
                            df_normal=df_normal,
                            how="inner"
                        )
                
                dff, normal_samples, tumor_samples = gdc.build_df_exp_and_filter(df_counts=df_counts, df_meta=df_meta,
                            gene_col = "geneid",
                            equal_var = False,)

            with tab_degs:

                lfc_cutoff = st.slider(
                    "Log2 fold change cutoff", min_value=0.1, max_value=10.0, value=1.0, key='lfc_cutoff'
                )

                fdr_cutoff = st.slider(
                    "FDR cutoff", min_value=0.01, max_value=1.0, value=0.05
                )

                mtd.LFC_cut = lfc_cutoff
                mtd.lfc_FDR_cut = fdr_cutoff

                degs, _, dflfc = mtd.list_of_degs(save_file=False, force=False, prompt_verbose=False, verbose=verbose)

                cols = ['ensembl_id','symbol','biotype', 'abs_lfc', 'lfc','pval','fdr', 'baseMean']
                dflfc = dflfc[cols]

                st.write(
                    f"There are {len(degs)} DEGs: params = lfc_cutoff={lfc_cutoff}, fdr_cutoff={fdr_cutoff}, and method={method}"
                )

                grid_key = f"degs_{psi_id}_lfc_{lfc_cutoff}_fdr_{fdr_cutoff}"
                show_df(dflfc, height=None, key=grid_key)

            with tab_echo:
                stri = mtd.echo_parameters(want_echo_default=True, jump_line=True, echo=False)
                stri = stri.replace('\n', '<br>').replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;")
                st.markdown(stri, unsafe_allow_html=True)

            with tab_biotypes:
              
                df_grouped = group_biotypes(dflfc)

                show_df(
                    df_grouped,
                    height=None,
                    key=f"biotypes_{selected_primary_site}",
                )

                explain = "protein_coding, pseudogene, lncRNA, "
                explain += "immune_receptor (IG_, TR_), small_RNA (miRNA, snRNA, snoRNA, scaRNA, misc_RNA, vault_RNA, scRNA), "
                explain += "TEC (To be Experimentally Confirmed), mitochondrial_RNA, and other"

                st.write(explain)

            with tab_lnc:

                lfc_cutoff_nc = st.slider(
                    "Log2 fold change cutoff", min_value=0.1, max_value=10.0, value=1.0, key="lfc_cutoff_nc"
                )

                fdr_cutoff_nc = st.slider(
                    "FDR cutoff", min_value=0.01, max_value=1.0, value=0.05, key="fdr_cutoff_nc"
                )

                dfnc, msg_nc = mtd.filter_non_coding(lfc_cut=lfc_cutoff_nc, fdr_cut=fdr_cutoff_nc, verbose=False, force=False)

                st.write(msg_nc)

                cols = ['ensembl_id','symbol','biotype', 'abs_lfc', 'lfc','pval','fdr', 'baseMean']
                grid_key = f"non-coding_{psi_id}_lfc_{lfc_cutoff_nc}_fdr_{fdr_cutoff_nc}"
                show_df(dfnc[cols], height=None, key=grid_key)

            if not dff.empty and not dfpiv.empty:
                
                n_samples, n_genes = dfpiv.shape
                
                with tab_heatmap_exp:

                    st.subheader("Expression matrix heatmap")
                    title = f"Primary Site: '{selected_primary_site}' #{n_samples} samples and #{n_genes} genes"
                    plot_heatmap_exp(dff, normal_samples, tumor_samples, title)

                    if msg:
                        st.write(msg)

                with tab_umap_exp:
                    n_clusters_umap_exp = st.slider(
                        "Number of clusters", min_value=2, max_value=20, value=3, key="n_clusters_umap_exp"
                    )

                    n_neighbors_umap_exp = st.slider(
                        "Number of neighbors", min_value=2, max_value=100, value=10, key="n_neighbors_umap_exp"
                    )

                    n_samples, n_genes = dfpiv.shape
                    st.subheader("UMAP expression")
                    title = f"UMAP expression matrix for '{selected_primary_site}' #{n_samples} samples and #{n_genes} genes"
                    df_umap_exp = plot_umap_exp(dff, samples=tumor_samples, which_samples='Tumor', title=title, 
                                                n_clusters=n_clusters_umap_exp,  n_neighbors=n_neighbors_umap_exp)

                    if msg:
                        st.write(msg)


            with tab_enrich:

                fdr_ptw_cutoff = st.slider(
                    "FDR cutoff", min_value=0.01, max_value=1.0, value=0.05, key="fdr_ptw_cutoff"
                )

                min_genes = st.slider(
                    "Minimum number of genes", min_value=1, max_value=20, value=3, key="min_genes"
                )

                df_enr0 = calc_enrichment_analysis(verbose=False, force=False)
                df_enr = df_enr0[ (df_enr0.fdr < fdr_ptw_cutoff) & (df_enr0.num_of_genes >= min_genes) ]

                if df_enr.empty:
                    st.write("No enriched pathways found.")
                else:
                    cols_pathways = ["pathway",  "pathway_id", "pval", "fdr", "odds_ratio", "combined_score", "genes", "num_of_genes"]

                    grid_key = f"EA_{psi_id}_fdr_{fdr_ptw_cutoff}_min_genes_{min_genes}"
                    selected_pathway_row = show_df(df_enr, height=None, key=grid_key, selectable=True,)

                    if selected_pathway_row is not None:
                        pathway_id = selected_pathway_row["pathway_id"]
                        pathway = selected_pathway_row["pathway"]
                        print(">>> Selected pathway:", pathway_id, pathway, "...\n", selected_pathway_row)

                        if rea.df_gmt.empty:
                            df_gmt = rea.open_reactome_gmt(verbose=True)

                        dfa = rea.df_gmt.loc[rea.df_gmt['pathway_id'] == pathway_id]
                        pathway_genes = [] if dfa.empty else dfa.iloc[0]['genes']
                        if isinstance(pathway_genes, str):
                            pathway_genes = list(eval(pathway_genes))
                            pathway_genes.sort()

                        st.session_state["selected_pathway_id"] = pathway_id
                        st.session_state["selected_pathway"] = pathway

                        dcy = DASH_CYTO(root0=ROOT0, root0_data=ROOT_DATA, dflfc_ori=mtd.dflfc_ori, 
                                        lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff, 
                                        found_degs=degs, pathway_genes=pathway_genes)

                        ret = dcy.read_owl(pathway_id, pathway, verbose=True)
                        if not ret:
                            fname_owl = f"{pathway_id}_level3.owl"
                            filename = dcy.root_owl / fname_owl
                            st.error(f"Failed to load OWL file {filename}")
                        else:
                            height = "95%"
                            width = "100%"
                            marginTop="20px"

                            dcy.run_app(height=height, width=width, marginTop=marginTop)


    # -------------------------------------------------------------------------
    # TAB 5 - DOWNLOADS
    # -------------------------------------------------------------------------
    with tab_donwload:
        st.write("Export filtered tables")

        st.download_button(
            "Download filtered cases TSV",
            data=df_cases.to_csv(sep="\t", index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_cases.tsv",
            mime="text/tab-separated-values",
        )

        st.download_button(
            "Download tumor samples TSV",
            data=df_all_samples.to_csv(sep="\t", index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_tumor_samples.tsv",
            mime="text/tab-separated-values",
        )

        st.download_button(
            "Download filtered mutations TSV",
            data=df_all_mut.to_csv(sep="\t", index=False).encode("utf-8"),
            file_name=f"{prog_id}_{selected_primary_site}_mutations.tsv",
            mime="text/tab-separated-values",
        )

        if not dfpiv.empty:
            st.download_button(
                "Download mutation matrix TSV",
                data=dfpiv.reset_index().to_csv(sep="\t", index=False).encode("utf-8"),
                file_name=f"{prog_id}_{selected_primary_site}_mutation_matrix.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
            )
    show_profile_box()

else:
    st.info("Click **Load data** in the sidebar to start.")
    show_profile_box()
