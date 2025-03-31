#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:21:28 2024

@author: Efe Aydın
"""

import streamlit as st
import pandas as pd
from analyzerFunctions import *
import base64
import os


# -------- Configuration --------
REFERENCE_MANUAL_PATH = os.path.join("static", "reference_manual.pdf")

# streamlit secrets
PASSWORD = st.secrets.get("password", "defaultpassword")  

# -------- Pages --------


def main_page():
    st.title("Micro-C derived CREs in ALL")
    st.write("This database contains whole-genome interaction maps of pediatric ALL derived from MicroC experiments performed on 35 primary samples")

    analysis_option = st.sidebar.selectbox("Choose your analysis option:", ["General (Including all BCP-ALL cases)", "Subtype specific"])
    
    if analysis_option == "Subtype specific":
        subtype = st.sidebar.selectbox("Choose subtype:", ["High Hyperdiploidy", "ETV6::RUNX1", "BCR::ABL1", "DUX4r", "TCF3::PBX1", "KMT2Ar", "iAMP21", "nearHaploid"])
        subtype_index_map = {"High Hyperdiploidy": 1, "ETV6::RUNX1": 2, "BCR::ABL1": 3, "DUX4r": 4, "TCF3::PBX1": 5, 
                             "KMT2Ar": 6, "iAMP21": 7, "nearHaploid": 8}
        subChoice = subtype_index_map[subtype]
        resolution = "10k"  # Only option for subtype specific
        st.sidebar.write(f"Resolution: {resolution} (fixed)")
    else:
        subChoice = 0
        resolution = st.sidebar.selectbox("Choose resolution:", ["1k", "10k"])

    focus_option = st.sidebar.selectbox("Are you interested in a specific gene or a location?", ["Gene", "Location"])
    cre_option = st.sidebar.selectbox("Predictions based on:", ["All promoters", "Canonical promoters"], index=0)
    cre_index = int(["All promoters", "Canonical promoters"].index(cre_option))
    
    if focus_option == "Gene":
        gene_of_interest = st.text_input("Enter gene of interest in HGNC SYMBOL (Example: KRAS):")
        if gene_of_interest:
            gene_of_interest = gene_of_interest.strip().upper()
            geneAnalyzer(subChoice, resolution, gene_of_interest, cre_index)
    elif focus_option == "Location":
        chr_input = st.selectbox("Choose chromosome of interest:", [str(i) for i in range(1, 23)] + ["X", "Y"])
        start_input = st.text_input("Start Position")
        end_input = st.text_input("End Position")
        
        if chr_input and start_input and end_input:
            try:
                start_input = int(start_input)
                end_input = int(end_input)
                if start_input >= end_input:
                    st.error("Start position cannot be larger than or equal to the end position.")
                elif end_input - start_input > 1_000_000:
                    st.error("Maximum range allowed is 1 Mb.")
                else:
                    locAnalyzer(subChoice, resolution, chr_input, start_input, end_input, cre_index)
            except ValueError:
                st.error("Only numeric values are allowed for start and end positions.")
    
    st.sidebar.markdown("---")
    if st.sidebar.button("Get Info"):
        st.session_state.page = "info_page"
        st.rerun()

def info_page():
    st.header("Micro-C derived CRE database for pediatric ALL")
    st.write("""
    This database represents visualization of loop calling results from 35 BCP-ALL samples integrated with the multi-omics data as described in the paper.

    Our app allows both gene-centric and region-centric exploration of regulatory landscapes, with support for general and subtype-specific queries.

    See reference manual for information on usage and content.

    If you use this app in your research, please cite:
    []

    Contact:
    Division of Clinical Genetics, Lund University Faculty of Medicine  
    BMC, C13  
    221 84 Lund, Sweden
    """)

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Go to App"):
            st.session_state.page = "main_page"
            st.rerun()
    with col2:
        if st.button("Reference Manual"):
            st.session_state.page = "pdf_page"
            st.rerun()

def pdf_page():
    st.title("Reference Manual")

    if not os.path.exists(REFERENCE_MANUAL_PATH):
        st.error("Reference manual not found. Please place 'reference_manual.pdf' in the 'static/' folder.")
        return

    with open(REFERENCE_MANUAL_PATH, "rb") as pdf_file:
        base64_pdf = base64.b64encode(pdf_file.read()).decode("utf-8")

    pdf_display = f'''
        <iframe src="data:application/pdf;base64,{base64_pdf}" width="100%" height="900px" style="border: none;"></iframe>
    '''
    st.markdown(pdf_display, unsafe_allow_html=True)

    if st.button("Back to Info Page"):
        st.session_state.page = "info_page"
        st.rerun()

def login():
    st.title("Login")
    password_input = st.text_input("Enter Password:", type="password")
    if st.button("Login"):
        if password_input == PASSWORD:
            st.session_state.authenticated = True
            st.session_state.page = "info_page"
            st.rerun()
        else:
            st.error("Incorrect password. Please try again.")

def main():
    if "authenticated" not in st.session_state:
        st.session_state.authenticated = False

    if "page" not in st.session_state:
        st.session_state.page = "login"

    if not st.session_state.authenticated:
        login()
        return

    if st.session_state.page == "info_page":
        info_page()
    elif st.session_state.page == "main_page":
        main_page()
    elif st.session_state.page == "pdf_page":
        pdf_page()

if __name__ == "__main__":
    main()
