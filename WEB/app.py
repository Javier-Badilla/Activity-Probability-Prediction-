import streamlit as st
import pandas as pd
import io
from seq2smiles import clean_sequence, sequence_to_smiles # Importing your core logic

# Set up the page
st.set_page_config(page_title="Protein Sequence → Canonical SMILES", page_icon="🧬")
st.title("Protein Sequence → Canonical SMILES")
import os

# Create three columns for the logos
col1, col2, col3 = st.columns(3)

# Display the images in their respective columns
# We use a try-except block just in case the files aren't found
try:
    with col1:
        st.image("assets/NBC_LOGO_TRANS.png", width=150)
    with col2:
        st.image("assets/PUCV_100.png", width=150)
    with col3:
        st.image("assets/UTFSM_LOGO_TRANS.png", width=150)
except FileNotFoundError:
    st.warning("Make sure the 'assets' folder and logos are uploaded to GitHub!")
    
st.markdown(
    """
    *This creates SMILES according to python's rdkit. 
    Made for the thesis work of Javier Badilla at the Núcleo Biotecnología Curauma 
    of the Pontificia Universidad Católica de Valparaíso.*
    """
)
st.write("**Format:** column 1 = ID, column 2 = sequence (one-letter code). First row is the header.")

# 1. File Uploader
uploaded_file = st.file_uploader("Upload Input Excel File", type=["xlsx", "xls", "xlsm"])

if uploaded_file is not None:
    # 2. Process the file in memory
    try:
        df = pd.read_excel(uploaded_file, header=0, dtype=str)
        
        if df.shape[1] < 2:
            st.error("The Excel file needs at least two columns (ID, sequence).")
        else:
            ids = df.iloc[:, 0]
            seqs = df.iloc[:, 1]
            
            smiles_out = []
            status_out = []
            
            # Progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()
            total = len(seqs)
            
            for i, raw in enumerate(seqs, start=1):
                try:
                    smiles_out.append(sequence_to_smiles(clean_sequence(raw)))
                    status_out.append("OK")
                except Exception as e:
                    smiles_out.append("")
                    status_out.append(f"ERROR: {e}")
                
                # Update progress
                progress_bar.progress(i / total)
                status_text.text(f"Processing… {i}/{total}")
            
            # 3. Create the result DataFrame
            result_df = pd.DataFrame({
                df.columns[0]: ids,
                df.columns[1]: seqs,
                "Canonical_SMILES": smiles_out,
                "Status": status_out,
            })
            
            status_text.text("Done! Ready for download.")
            st.dataframe(result_df.head()) # Show a preview of the results
            
            # 4. Prepare for download using an in-memory buffer
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                result_df.to_excel(writer, index=False)
            output.seek(0)
            
            st.download_button(
                label="Download SMILES Excel File",
                data=output,
                file_name="sequences_smiles.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
    except Exception as e:
        st.error(f"An error occurred reading the file: {e}")