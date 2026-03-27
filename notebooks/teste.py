

import time
import pandas as pd
import streamlit as st

st.set_page_config(page_title="GDC Cases Explorer", layout="wide")

# ---------- STYLE ----------
st.markdown("""
<style>
    .block-container {
        padding-top: 0.7rem !important;
        padding-bottom: 1rem;
    }

    .app-title {
        font-size: 1.35rem;
        font-weight: 700;
        margin-top: -10px;
        margin-bottom: 0.2rem;
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

    .filter-box {
        padding: 0.8rem 0.8rem 0.3rem 0.8rem;
        border: 1px solid rgba(128,128,128,0.25);
        border-radius: 12px;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# ---------- HEADER ----------
st.markdown("<div class='app-title'>GDC Cases Explorer</div>", unsafe_allow_html=True)
st.markdown(
    "<div class='app-subtitle'>Explore programs, primary sites, tumor classes, subtypes, and stages from GDC.</div>",
    unsafe_allow_html=True
)

# ---------- MOCK DATA ----------
# Replace this with your real lookup tables / API-driven options
PROGRAMS = ["TCGA", "TARGET"]
PRIMARY_SITES = {
    "TCGA": ["Lung", "Breast", "Brain", "Kidney"],
    "TARGET": ["Blood", "Kidney", "Nervous System"]
}

# ---------- CACHED FUNCTION ----------
@st.cache_data(show_spinner=False)
def fetch_gdc_cases(program: str, primary_site: str) -> pd.DataFrame:
    """
    Replace this function with your real GDC API code.
    Cached so repeated selections are faster.
    """
    time.sleep(1.5)  # simulate request delay

    rows = [
        {
            "program": program,
            "primary_site": primary_site,
            "tumor_class": "carcinoma",
            "subtype": "adenocarcinoma",
            "stage": "Stage I",
            "n": 42,
        },
        {
            "program": program,
            "primary_site": primary_site,
            "tumor_class": "carcinoma",
            "subtype": "squamous cell carcinoma",
            "stage": "Stage II",
            "n": 31,
        },
        {
            "program": program,
            "primary_site": primary_site,
            "tumor_class": "other",
            "subtype": "other",
            "stage": "Stage III",
            "n": 12,
        },
    ]
    return pd.DataFrame(rows)

# ---------- SESSION STATE ----------
if "results_df" not in st.session_state:
    st.session_state.results_df = None

if "has_searched" not in st.session_state:
    st.session_state.has_searched = False

# ---------- FILTERS ----------
st.markdown("<div class='filter-box'>", unsafe_allow_html=True)

col1, col2 = st.columns(2)

with col1:
    selected_program = st.selectbox(
        "Program",
        PROGRAMS,
        index=0
    )

with col2:
    selected_primary_site = st.selectbox(
        "Primary Site",
        PRIMARY_SITES[selected_program],
        index=0
    )

st.markdown("</div>", unsafe_allow_html=True)

# ---------- ACTIONS ----------
col_btn1, col_btn2, col_btn3 = st.columns([2, 1, 1])

with col_btn1:
    run_search = st.button(
        "Find cases, subtypes, tumor class and stages",
        use_container_width=True,
        type="primary"
    )

with col_btn2:
    clear_results = st.button(
        "Clear results",
        use_container_width=True
    )

with col_btn3:
    show_table = st.checkbox("Wide table", value=True)

if clear_results:
    st.session_state.results_df = None
    st.session_state.has_searched = False

if run_search:
    with st.spinner("Querying GDC data..."):
        df = fetch_gdc_cases(selected_program, selected_primary_site)
        st.session_state.results_df = df
        st.session_state.has_searched = True

# ---------- RESULTS ----------
if st.session_state.has_searched:
    st.divider()

    df = st.session_state.results_df

    if df is None or df.empty:
        st.warning("No results found for the selected filters.")
    else:
        c1, c2, c3 = st.columns(3)
        c1.metric("Rows", len(df))
        c2.metric("Tumor classes", df["tumor_class"].nunique())
        c3.metric("Subtypes", df["subtype"].nunique())

        st.markdown("#### Results")

        if show_table:
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.table(df)

        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download CSV",
            data=csv,
            file_name="gdc_cases_results.csv",
            mime="text/csv",
            use_container_width=False
        )
        