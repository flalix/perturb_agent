import streamlit as st
import pandas as pd

st.set_page_config(page_title="Chained Filters Table", layout="wide")

# ----------------------------
# Sample data
# ----------------------------
df = pd.DataFrame({
    "country": ["Brazil", "Brazil", "Brazil", "USA", "USA", "Germany", "Germany"],
    "state":   ["SP", "SP", "RJ", "CA", "NY", "BW", "BY"],
    "city":    ["Santo Andre", "Sao Paulo", "Rio", "Los Angeles", "New York", "Heidelberg", "Munich"],
    "gene":    ["TP53", "EGFR", "BRCA1", "KRAS", "PIK3CA", "MYC", "PTEN"],
    "value":   [10, 20, 15, 30, 12, 18, 25]
})

st.title("Chained Dropdowns + Pandas Table")

# ----------------------------
# Chained dropdowns
# ----------------------------
col1, col2, col3 = st.columns(3)

with col1:
    countries = ["All"] + sorted(df["country"].dropna().unique().tolist())
    selected_country = st.selectbox("Country", countries)

# filter for state options
df1 = df.copy()
if selected_country != "All":
    df1 = df1[df1["country"] == selected_country]

with col2:
    states = ["All"] + sorted(df1["state"].dropna().unique().tolist())
    selected_state = st.selectbox("State", states)

# filter for city options
df2 = df1.copy()
if selected_state != "All":
    df2 = df2[df2["state"] == selected_state]

with col3:
    cities = ["All"] + sorted(df2["city"].dropna().unique().tolist())
    selected_city = st.selectbox("City", cities)

# ----------------------------
# Final filtered dataframe
# ----------------------------
filtered = df.copy()

if selected_country != "All":
    filtered = filtered[filtered["country"] == selected_country]

if selected_state != "All":
    filtered = filtered[filtered["state"] == selected_state]

if selected_city != "All":
    filtered = filtered[filtered["city"] == selected_city]

st.subheader("Filtered table")

# Main table
st.dataframe(filtered, use_container_width=True)

# ----------------------------
# Bottom filters
# ----------------------------
st.markdown("---")
st.subheader("Extra filters")

bottom_col1, bottom_col2 = st.columns(2)

with bottom_col1:
    gene_filter = st.multiselect(
        "Filter by gene",
        options=sorted(filtered["gene"].dropna().unique().tolist()),
        default=[]
    )

with bottom_col2:
    min_value, max_value = int(filtered["value"].min()), int(filtered["value"].max())
    value_range = st.slider(
        "Value range",
        min_value=min_value,
        max_value=max_value,
        value=(min_value, max_value)
    )

final_df = filtered.copy()

if gene_filter:
    final_df = final_df[final_df["gene"].isin(gene_filter)]

final_df = final_df[
    (final_df["value"] >= value_range[0]) &
    (final_df["value"] <= value_range[1])
]

st.subheader("Final result")
st.dataframe(final_df, use_container_width=True)