#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/26
# Udated  on 2026/03/26
# @author: Flavio Lichtenstein
# @local: Home sweet home

#--------------- init commands --------------------------
# 
# cd ~/uv/perturb_agent$/notebooks/
# source .venv/bin/activate
# mamba activate renv
# 
#------------------------------------------------------

import os, sys
import pandas as pd
import streamlit as st

from pathlib import Path

ROOT = Path().resolve().parent.parent
SRC = os.path.join(ROOT, "src")

if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

print("ROOT:", ROOT)
print("SRC added:", SRC)

from libs.tcga_gdc_lib import *
from libs.Basic import *
from libs.calc_degs_lib import CALC_DEGS

ROOT = Path().resolve().parent
root_data = os.path.join(ROOT, "data/tcga")

gdc = GDC(root_data=root_data)

verbose = True


# ---------- STYLE ----------

st.markdown("""
<style>
    .block-container {
        max-width: 95% !important;
        padding-top: .5rem !important;
        padding-bottom: .5rem;
    }
    .app-title {
        font-size: .85rem;
        font-weight: 500;
        line-height: 1.1;
        margin-top: 20px;
        margin-bottom: 0px;
        color: #666;
    }
    .app-subtitle {
        font-size: 0.95rem;
        color: #666;
        margin-bottom: 1rem;
    }
    div.stButton > button {
        height: 3.2em;
        font-size: 1.08rem;
        font-weight: 600;
        border-radius: 10px;
    }
    div[data-testid="stButton"] button[kind="primary"] {
        height: 3.4em;
        font-size: 1.15rem;
        font-weight: 700;
    }
    .filter-box {
        padding: 0.8rem 0.8rem 0.3rem 0.8rem;
        border: 1px solid rgba(128,128,128,0.25);
        border-radius: 12px;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

@st.cache_data(show_spinner=False)
def load_programs(verbose: bool) -> list[str]:
    prog_list = gdc.get_gdc_progams(force=False, verbose=verbose)
    if not isinstance(prog_list, list):
        raise TypeError("gdc.get_gdc_progams() did not return a list")
    return [str(x) for x in prog_list]


@st.cache_data(show_spinner=False)
def load_primary_sites(program: str, verbose: bool) -> pd.DataFrame:
    df_psi = gdc.get_primary_sites(program=program, force=False, verbose=verbose)
    if not isinstance(df_psi, pd.DataFrame):
        raise TypeError("gdc.get_primary_sites() did not return a DataFrame")

    expected = {"pid", "primary_site", "project_id", "disease_type", "name"}
    missing = expected - set(df_psi.columns)
    if missing:
        raise ValueError(f"Missing columns in df_psi: {sorted(missing)}")

    return df_psi.copy()


@st.cache_data(show_spinner=False)
def load_cases_and_subtypes(pid: str, verbose: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_cases, df_subt, _ = gdc.get_cases_and_subtypes(
        pid=pid,
        batch_size=200,
        do_filter=False,
        force=False,
        verbose=verbose,
    )

    if not isinstance(df_cases, pd.DataFrame):
        raise TypeError("df_cases is not a DataFrame")
    if not isinstance(df_subt, pd.DataFrame):
        raise TypeError("df_subt is not a DataFrame")

    subt_cols = ["pid", "subtype_global", "tumor_class", "subtype_tissue", "stage", "n"]
    missing_subt = [c for c in subt_cols if c not in df_subt.columns]
    if missing_subt:
        raise ValueError(f"df_subt missing columns: {missing_subt}")

    case_cols = ["pid", "case_id", "subtype_global", "tumor_class", "subtype_tissue", "stage"]
    missing_cases = [c for c in case_cols if c not in df_cases.columns]
    if missing_cases:
        raise ValueError(f"df_cases missing columns: {missing_cases}")

    # df_subt  = df_subt[subt_cols].copy()
    df_cases = df_cases[case_cols].copy().reset_index(drop=True)

    return df_subt, df_cases


@st.cache_data(show_spinner=False)
def load_lfc_and_expression_tables(pid: str, df_samples: pd.DataFrame, data_type:str='Gene Expression Quantification', 
                                   lfc_cutoff = 1.0, fdr_cutoff = .05, 
                                   verbose:bool=False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    
    case_id_list = list(np.unique(df_samples.case_id))

    df_normal, df_tumor = gdc.get_tumor_normal_tables(
        df_samples=df_samples, case_id_list=case_id_list, data_type=data_type, verbose=verbose)
    
    if not isinstance(df_normal, pd.DataFrame):
        raise TypeError("df_normal is not a DataFrame")
    
    if not isinstance(df_tumor, pd.DataFrame):
        raise TypeError("df_tumor is not a DataFrame")

    df_normal, df_tumor = gdc.merge_normal_tumor_tables(pid=pid, df_normal=df_normal, df_tumor=df_tumor)

    cdegs = CALC_DEGS(root_data=root_data)

    df_normal = cdegs.deduplicate_by_max_reads(df_normal)
    df_tumor  = cdegs.deduplicate_by_max_reads(df_tumor)

    df_counts, df_meta = cdegs.build_counts_and_metadata(
        df_tumor=df_tumor,
        df_normal=df_normal,
        how="inner"
    )

    df_lfc = cdegs.run_deg_rscript(df_tumor=df_tumor, df_normal=df_normal,
                                   method="edger",  manual_dispersion=0.1, min_total_count=10, 
                                   merge_how="inner", keep_temp=False)
    
    if not isinstance(df_lfc, pd.DataFrame):
        raise TypeError("df_lfc is not a DataFrame")
    

    df_degs = df_lfc[ (df_lfc.lfc >= lfc_cutoff) & (df_lfc.fdr < fdr_cutoff)].copy()    

    return df_degs, df_counts, df_meta


def init_state() -> None:
    defaults = {
        "program": "TCGA",
        "primary_site": None,
        "df_primary_sites": pd.DataFrame(),
        "df_subt": pd.DataFrame(),
        "df_cases": pd.DataFrame(),
        "df_counts": pd.DataFrame(),
        "df_meta": pd.DataFrame(),
        "df_lfc": pd.DataFrame(),
        "df_degs": pd.DataFrame(),
        "pid": "",
        "loaded_program": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def main() -> None:

    # ---------- HEADER ----------
    st.markdown("<div class='app-title'>GDC Cases Explorer</div>", unsafe_allow_html=True)
    st.markdown(
        "<div class='app-subtitle'>Explore programs, primary sites, tumor classes, subtypes, and stages from GDC.</div>",
        unsafe_allow_html=True
    )
    init_state()

    st.title("GDC Cases Explorer")

    try:
        prog_list = load_programs(verbose=verbose)
    except Exception as e:
        st.error(f"Error loading programs: {e}")
        st.stop()

    if not prog_list:
        st.warning("No programs returned.")
        st.stop()

    default_program_index = prog_list.index("TCGA") if "TCGA" in prog_list else 0

    col1, col2 = st.columns(2)

    with col1:
        selected_program = st.selectbox(
            "Program",
            options=prog_list,
            index=default_program_index,
        )

    # Reload primary sites when program changes
    if st.session_state.loaded_program != selected_program:
        st.session_state.loaded_program = selected_program
        st.session_state.program = selected_program
        st.session_state.primary_site = None
        st.session_state.df_subt = pd.DataFrame()
        st.session_state.df_cases = pd.DataFrame()
        st.session_state.df_counts = pd.DataFrame()
        st.session_state.df_lfc = pd.DataFrame()
        st.session_state.pid = ""

        try:
            df_psi = load_primary_sites(program=selected_program, verbose=verbose)
            st.session_state.df_primary_sites = df_psi
        except Exception as e:
            st.session_state.df_primary_sites = pd.DataFrame()
            st.error(f"Error loading primary sites: {e}")

    df_psi = st.session_state.df_primary_sites

    primary_site_options: list[str] = []
    if not df_psi.empty:
        primary_site_options = (
            df_psi["primary_site"]
            .dropna()
            .astype(str)
            .sort_values()
            .unique()
            .tolist()
        )

    with col2:
        if primary_site_options:
            selected_primary_site = st.selectbox(
                "Primary Site",
                options=primary_site_options,
                index=None,
                placeholder="Select primary site",
                key="primary_site",
            )
        else:
            selected_primary_site = None
            st.selectbox("Primary Site", ["No primary sites available"], disabled=True)
    
    st.markdown("""
    <style>
    div[data-testid="stButton"] button[kind="primary"] {
        height: 2em;
        font-size: 1.3rem;
    }
    </style>
""", unsafe_allow_html=True)
    
    _, col_btn, _ = st.columns([1, 2, 1])

    with col_btn:
        run_search = st.button("Find cases, subtypes, tumor class and stages", width="stretch", type="primary")

    if run_search:
        if not selected_primary_site:
            st.warning("Please select a primary site.")
        elif df_psi.empty:
            st.warning("Primary site table is empty.")
        else:
            df_match = df_psi[df_psi["primary_site"].astype(str) == str(selected_primary_site)]

            if df_match.empty:
                st.warning("Could not find pid for the selected primary site.")
            else:
                pid = str(df_match.iloc[0]["pid"])
                st.session_state.pid = pid
                # once a widget owns a key, you should not manually overwrite
                # st.session_state.primary_site = selected_primary_site

                try:
                    with st.spinner(f"Loading cases and subtypes for {pid}..."):
                        df_subt, df_cases = load_cases_and_subtypes(pid=pid, verbose=verbose)

                    st.session_state.df_subt = df_subt.reset_index(drop=True)
                    st.session_state.df_cases = df_cases.reset_index(drop=True)

                    # clear previous downstream state
                    for key in [
                        "selected_case_idx",
                        "selected_case_row",
                        "df_samples",
                        "df_lfc",
                        "df_counts",
                        "df_meta",
                    ]:
                        st.session_state.pop(key, None)

                    st.success(
                        f"Loaded {len(df_subt)} subtype rows and {len(df_cases)} case rows for {pid}."
                    )
                except Exception as e:
                    st.error(f"Error loading cases/subtypes: {e}")

    tab1, tab2, tab3 = st.tabs(["Subtype + Cases", "Counts", "DEGs"])

    with tab1:
        st.subheader("Subtype table")
        if "df_subt" in st.session_state and not st.session_state.df_subt.empty:
            event = st.dataframe(
                st.session_state.df_subt,
                width="stretch",
                hide_index=True,
                on_select="rerun",
                selection_mode="single-row",
                key="df_subt_select",
            )

            selected_rows = event.selection.rows

            if selected_rows:
                st.session_state.selected_subt_idx = selected_rows[0]

            if "selected_subt_idx" in st.session_state:
                row = st.session_state.df_subt.iloc[st.session_state.selected_subt_idx]
                st.session_state.selected_subt_row = row.to_dict()

                pid = str(row["pid"])
                subtype_global = row["subtype_global"]
                tumor_class = row["tumor_class"]
                subtype_tissue = row["subtype_tissue"]
                stage = row["stage"]

                st.markdown("**Selected subtype row**")
                st.dataframe(
                    pd.DataFrame([row]),
                    width="stretch",
                    hide_index=True,
                )

                try:
                    with st.spinner("Loading samples and expression tables..."):
                        st.success(f"PID: {pid}, subtype {subtype_global}, tumor {tumor_class}, tissue: {subtype_tissue}, stage {stage}")
                        df_samples = gdc.get_samples_for_pid_subtypes(
                            pid=pid,
                            subtype_global=subtype_global,
                            tumor_class=tumor_class,
                            subtype_tissue=subtype_tissue,
                            stage=stage,
                            batch_size=200,
                            force=False,
                            verbose=verbose,
                        )

                        st.success(f">> Loaded {len(df_samples)} samples.")

                        df_degs, df_counts, df_meta = load_lfc_and_expression_tables(
                            pid=pid,
                            df_samples=df_samples,
                            data_type="Gene Expression Quantification",
                            verbose=verbose,
                        )

                    st.session_state.df_samples = df_samples
                    st.session_state.df_degs = df_degs
                    st.session_state.df_counts = df_counts
                    st.session_state.df_meta = df_meta

                    st.success(f"Loaded {len(df_samples)} samples.")
                except Exception as e:
                    st.error(f"Error loading downstream Count and DEG tables: {e}")
        else:
            st.info("No subtype table loaded yet.")

        st.subheader("Cases table")
        if "df_cases" in st.session_state and not st.session_state.df_cases.empty:
            st.dataframe(
                st.session_state.df_cases,
                width="stretch",
                hide_index=True,
            )
        else:
            st.info("No cases table loaded yet.")


    with tab2:
        st.subheader("Counts")
        if "df_counts" in st.session_state and not st.session_state.df_counts.empty:
            st.dataframe(
                st.session_state.df_counts,
                width="stretch",
                hide_index=True,
            )
        else:
            st.info("No counts loaded yet.")

    with tab3:
        st.subheader("DEGs / LFC")
        if "df_lfc" in st.session_state and not st.session_state.df_lfc.empty:
            st.dataframe(
                st.session_state.df_lfc,
                width="stretch",
                hide_index=True,
            )
        else:
            st.info("Select one row in the subtype table.") 


if __name__ == "__main__":
    main()

