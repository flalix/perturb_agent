import streamlit as st
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder

st.title("AgGrid smoke test")

df = pd.DataFrame({
    "A": [1, 2, 3],
    "B": ["x", "y", "z"]
})

st.write("before grid")

gb = GridOptionsBuilder.from_dataframe(df)
grid_options = gb.build()

AgGrid(df, gridOptions=grid_options, height=200)

st.write("after grid")