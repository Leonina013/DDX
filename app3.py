import streamlit as st
import pandas as pd
from pyDESeq2 import DESeq2

st.set_page_config(page_title="Differential Expression Analysis", page_icon=":guardsman:", layout="wide")

st.title("Differential Expression Analysis")

# Function to upload file
def upload_file():
    file = st.file_uploader("Upload count data in CSV format", type=["csv"])
    if file is not None:
        data = pd.read_csv(file)
        return data

# Main function
def main():
    data = upload_file()
    if data is not None:
        st.write("Data shape: ", data.shape)
        condition_col1 = st.selectbox("Select the column name of the condition known variable", data.columns)
        condition_col2 = st.selectbox("Select the column name of the condition new variable", data.columns)
        known_samples = data[data[condition_col1] == "known"]
        new_samples = data[data[condition_col2] == "new"]
        
        deseq2_known = DESeq2(known_samples, design_formula="~ 1")
        deseq2_known.run_analysis()
        
        deseq2_new = DESeq2(new_samples, design_formula="~ 1")
        deseq2_new.run_analysis()
        
        comparison = deseq2_known.compare_to(deseq2_new)
        st.write("Differentially Expressed genes:")
        st.write(comparison)

if __name__=="__main__":
    main()