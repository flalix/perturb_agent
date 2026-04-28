#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/26
# Udated  on 2026/03/26
# @author: Flavio Lichtenstein
# @local: Home sweet home

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


ROOT = Path().resolve().parent
root_data = os.path.join(ROOT, "data/tcga")

gdc = GDC(root_data=root_data)

verbose = True


@st.cache_data(show_spinner=False)
def load_programs(force: bool, verbose: bool) -> list[str]:
    prog_list = gdc.get_gdc_progams(force=force, verbose=verbose)
    if not isinstance(prog_list, list):
        raise TypeError("gdc.get_gdc_progams() did not return a list")
    return [str(x) for x in prog_list]


@st.cache_data(show_spinner=False)
def load_primary_sites(program: str, verbose: bool) -> pd.DataFrame:
    dfc = gdc.get_primary_sites(program=program, force=False, verbose=verbose)
    if not isinstance(dfc, pd.DataFrame):
        raise TypeError("gdc.get_primary_sites() did not return a DataFrame")

    expected = {"pid", "primary_site", "project_id", "disease_type", "name"}
    missing = expected - set(dfc.columns)
    if missing:
        raise ValueError(f"Missing columns in dfc: {sorted(missing)}")

    return dfc.copy()


@st.cache_data(show_spinner=False)
def load_cases_and_subtypes(pid: str, verbose: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
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

    subt_cols = ["pid", "subtype_global", "tumor_class", "subtype_tissue", "sstage", "n"]
    missing_subt = [c for c in subt_cols if c not in df_subt.columns]
    if missing_subt:
        raise ValueError(f"df_subt missing columns: {missing_subt}")

    case_cols = ["pid", "case_id", "subtype_global", "tumor_class", "subtype_tissue", "stage"]
    missing_cases = [c for c in case_cols if c not in df_cases.columns]
    if missing_cases:
        raise ValueError(f"df_cases missing columns: {missing_cases}")

    df_subt2 = df_subt[subt_cols].copy()
    df_cases2 = df_cases[case_cols].copy()

    return df_subt2, df_cases2


def init_state() -> None:
    defaults = {
        "program": "TCGA",
        "primary_site": None,
        "df_primary_sites": pd.DataFrame(),
        "df_subt": pd.DataFrame(),
        "df_cases2": pd.DataFrame(),
        "pid": "",
        "loaded_program": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def main() -> None:

    st.markdown("""
        <style>
            .block-container {
                padding-top: 0.5rem !important;
            }

            h1 {
                margin-top: -20px !important;
            }
        </style>
    """, unsafe_allow_html=True)

    st.markdown("""
<h1 style='font-size: 1.2rem;'>GDC Cases Explorer</h1>
""", unsafe_allow_html=True)

    # st.set_page_config(page_title="GDC Cases Explorer", layout="wide")
    init_state()

    st.title("GDC Cases Explorer")

    try:
        prog_list = load_programs(force=False, verbose=verbose)
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
        st.session_state.df_cases2 = pd.DataFrame()
        st.session_state.pid = ""

        try:
            dfc = load_primary_sites(program=selected_program, verbose=verbose)
            st.session_state.df_primary_sites = dfc
        except Exception as e:
            st.session_state.df_primary_sites = pd.DataFrame()
            st.error(f"Error loading primary sites: {e}")

    dfc = st.session_state.df_primary_sites

    primary_site_options: list[str] = []
    if not dfc.empty:
        primary_site_options = (
            dfc["primary_site"]
            .dropna()
            .astype(str)
            .sort_values()
            .unique()
            .tolist()
        )

    with col2:
        selected_primary_site = st.selectbox(
            "Primary Site",
            options=primary_site_options,
            index=None,
            placeholder="Select primary site",
        )

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
        run_search = st.button("Find cases, subtypes, tumor class and stages", use_container_width=True, type="primary")

    if run_search:
        if not selected_primary_site:
            st.warning("Please select a primary site.")
        elif dfc.empty:
            st.warning("Primary site table is empty.")
        else:
            df_match = dfc[dfc["primary_site"].astype(str) == str(selected_primary_site)]

            if df_match.empty:
                st.warning("Could not find pid for the selected primary site.")
            else:
                pid = str(df_match.iloc[0]["pid"])
                st.session_state.pid = pid
                st.session_state.primary_site = selected_primary_site

                try:
                    with st.spinner(f"Loading cases and subtypes for {pid}..."):
                        df_subt, df_cases2 = load_cases_and_subtypes(pid=pid, verbose=verbose)

                    st.session_state.df_subt = df_subt
                    st.session_state.df_cases2 = df_cases2
                    st.success(
                        f"Loaded {len(df_subt)} subtype rows and {len(df_cases2)} case rows for {pid}."
                    )
                except Exception as e:
                    st.error(f"Error loading cases/subtypes: {e}")

    if not st.session_state.df_subt.empty:
        st.subheader("df_subt — subtype summary")
        st.dataframe(st.session_state.df_subt, use_container_width=True, height=260)

    if not st.session_state.df_cases2.empty:
        st.subheader("df_cases2 — cases")
        st.dataframe(st.session_state.df_cases2, use_container_width=True, height=320)


if __name__ == "__main__":
    main()

