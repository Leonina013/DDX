import streamlit as st
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

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
        condition_col = st.selectbox("Select the column name of the condition variable", data.columns)
        known_samples = data[data[condition_col] == "known"]
        new_samples = data[data[condition_col] == "new"]
        
        robjects.r('''
            library(DESeq2)
            ''')
        
        known_dds = robjects.r['DESeqDataSetFromMatrix'](countData = pandas2ri.py2ri(known_samples), 
                                                         colData = pandas2ri.py2ri(known_samples.loc[:,condition_col]),
                                                         design = robjects.Formula('~1'))
        robjects.r('dds <- estimateSizeFactors(dds)')
        robjects.r('dds <- estimateDispersions(dds)')
        robjects.r('res <- nbinomTest(dds,"known","new")')
        
        results = robjects.r['res']
        
        st.write("Differentially Expressed genes:")
        st.write(results)

if __name__=="__main__":
    main()
