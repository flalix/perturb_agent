import streamlit as st
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder

import streamlit as st
import pandas as pd
import sys
from importlib.metadata import version
import st_aggrid

st.title("AgGrid smoke test")


st.write("Python:", sys.version)
st.write("streamlit:", st.__version__)
st.write("pandas:", pd.__version__)
st.write("streamlit-aggrid:", version("streamlit-aggrid"))
st.write("st_aggrid file:", st_aggrid.__file__)
st.write("st_aggrid module:", getattr(st_aggrid, "__file__", "no file"))
st.write("protobuf:", version("protobuf"))


df = pd.DataFrame({
    "A": [1, 2, 3],
    "B": ["x", "y", "z"]
})

st.write("before grid")

gb = GridOptionsBuilder.from_dataframe(df)
grid_options = gb.build()

AgGrid(df, gridOptions=grid_options, height=200)

st.write("after grid")