#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refactored on Mar 26, 2025
@author: Efe Aydın
"""

import subprocess
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import tempfile
import pyBigWig

# -------------------- Configuration --------------------

DATA_DIR = os.path.join(os.getcwd(), "data")
TRACKS_DIR = os.path.join(os.getcwd(), "tracks")

# -------------------- Load Required Data --------------------

coding_genes_path = os.path.join(DATA_DIR, "coding_genes2")
all_loops_path = os.path.join(DATA_DIR, "concat_loops.tab")

if not os.path.exists(coding_genes_path) or not os.path.exists(all_loops_path):
    st.error("Required input files not found. Please place 'coding_genes2' and 'concat_loops.tab' in the 'data/' folder.")
    st.stop()

coding_genes = pd.read_csv(coding_genes_path, sep="\t")
all_loops = pd.read_csv(all_loops_path, sep="\t", header=0)
all_loops.iloc[:, [1, 2, 4, 5]] = all_loops.iloc[:, [1, 2, 4, 5]].astype(int)


# -------------------- Functions --------------------


def geneAnalyzer(subChoice, res, gene_of_interest, cre_index):
    gene_index = 9 - cre_index
    myTx = "allTX" if cre_index == 0 else "canonTX"

    if gene_of_interest not in coding_genes.iloc[:, 3].values:
        st.write("No protein coding genes with the given name were found.")
        return

    subtype_map = {
        1: "HeH_10k",
        2: "ER_10k",
        3: "BA_10k",
        4: "DUX4r_10k",
        5: "TP_10k",
        6: "KMT2Ar_10k",
        7: "iAMP_10k",
        8: "nearHaploid_10k",
        0: f"merged_{res}"
    }
    mySubtype = subtype_map.get(subChoice, "unknown")

    if res == "10k" and subChoice == 0:
        df = all_loops[~all_loops['loopSource'].isin(["merged_1k"])].drop_duplicates()
    else:
        df = all_loops[all_loops['loopSource'] == mySubtype]

    df.iloc[:, [1, 2, 4, 5]] = df.iloc[:, [1, 2, 4, 5]].astype(int)

    myList, myList2 = [], []

    for i in range(len(df)):
        gene_list = str(df.iloc[i, gene_index]).split(",")
        if gene_of_interest in gene_list:
            myList.append([f"chr{df.iloc[i, 0]}:{int(df.iloc[i, 1])}-{int(df.iloc[i, 2])}"])
            myList2.append([
                df.iloc[i, 3], df.iloc[i, 4], df.iloc[i, 5],
                df.iloc[i, 0], df.iloc[i, 1], df.iloc[i, 2], gene_of_interest
            ])

    myList = pd.DataFrame(myList, columns=['region'])
    myDf = pd.DataFrame(myList2, columns=['chr', 'start', 'end', 'interChr', 'interStart', 'interEnd', 'target']).drop_duplicates()
    myDf.iloc[:, [1, 2, 4, 5]] = myDf.iloc[:, [1, 2, 4, 5]].astype(int)
    
    if myDf.empty:
        st.write("No regulatory elements were identified for this gene")
        return
    else:
        st.write("Putative enhancers for this gene were found in the following regions:")
        for region in myList['region'].drop_duplicates():
            st.write(region)

    csv_data = myDf.to_csv(index=False)
    st.download_button(
        label="Download Enhancer Results as CSV",
        data=csv_data,
        file_name=f"{gene_of_interest}_{myTx}_{res}_{mySubtype}_enhancer_results.csv",
        mime="text/csv"
    )

    extended_min_start, extended_max_end = get_genomic_range(myDf)
    region = f"chr{myDf.iloc[0, 0]}:{int(extended_min_start)}-{int(extended_max_end)}"

    temp_links_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tab')
    for _, row in myDf.iterrows():
        temp_links_file.write(f"{row['chr']}\t{int(row['start'])}\t{int(row['end'])}\t{row['interChr']}\t{int(row['interStart'])}\t{int(row['interEnd'])}\t1.0\n")
    temp_links_file.close()

    new_bed_file = coding_genes.copy()
    new_bed_file.loc[new_bed_file.iloc[:, 3] == gene_of_interest, new_bed_file.columns[6]] = "255,0,0"
    new_bed_file.insert(6, "thickStart", 1000)
    new_bed_file.insert(7, "thickEnd", 2000)
    new_bed_file = new_bed_file.sort_values(by=[new_bed_file.columns[0], new_bed_file.columns[1]])
    new_bed_file.to_csv(os.path.join(TRACKS_DIR, "tempCodingGenes.bed"), index=False, sep="\t", header=False)

    bb_1 = pyBigWig.open("https://data.cyverse.org/dav-anon/iplant/home/efeaydin/h3k27ac.bigWig")
    h3k27ac_path = tempfile.NamedTemporaryFile(delete=False, suffix=".bigBed").name
    bb_1.close()  # Just closes the remote connection
    
    gene_tracks_path = os.path.join(TRACKS_DIR, "geneTracks.ini")
    temp_tracks_path = os.path.join(TRACKS_DIR, "temp_gene_tracks.ini")
    temp_tracks_file = open(temp_tracks_path, "w")
    with open(gene_tracks_path, "r") as base_tracks:
        lines = base_tracks.readlines()
        temp_tracks_file.write(lines[0])
        temp_tracks_file.write("\n")
        # histone h3k27ac
        temp_tracks_file.write("[encode_bigbed]\n")
        temp_tracks_file.write(f"file = {h3k27ac_path}\n") 
        temp_tracks_file.write("file_type = bed\n")
        temp_tracks_file.write("title = ENCODE bigBed\n")
        temp_tracks_file.write("height = 3\n")
        temp_tracks_file.write("color = black\n\n")
        # enhancer links
        temp_tracks_file.write("[enhancer_links]\n")
        temp_tracks_file.write(f"file = {temp_links_file.name}\n")
        temp_tracks_file.write("file_type = links\n")
        temp_tracks_file.write("line_width = 3\n")
        temp_tracks_file.write("color = red\n")
        temp_tracks_file.write("title = Enhancer Links\n")
        temp_tracks_file.write("height = 5\n\n")
        #  Third: Append rest of the tracks
        temp_tracks_file.writelines(lines[1:])
        temp_tracks_file.close()

    output_file = "output_genome_track.png"
    pyGenomeTracks_command = [
        "pyGenomeTracks",
        "--tracks", temp_tracks_path,
        "--region", region,
        "--dpi", "100",
        "-o", output_file
    ]

    try:
        st.write("📣 Running pyGenomeTracks with command:")
        st.code(" ".join(pyGenomeTracks_command))
        st.write("📄 Here's the contents of temp_gene_tracks.ini:")

        with open(temp_tracks_path, "r") as f:
            st.code(f.read(), language="ini")
        subprocess.run(pyGenomeTracks_command, check=True)
        img = mpimg.imread(output_file)
        plt.figure(figsize=(10, 5))
        plt.imshow(img)
        plt.axis('off')
        st.pyplot(plt.gcf())

        os.remove(output_file)
        os.remove(temp_links_file.name)
        os.remove(temp_tracks_file.name)
        os.remove(h3k27ac_path)
        os.remove(os.path.join(TRACKS_DIR, "tempCodingGenes.bed"))
    except subprocess.CalledProcessError as e:
        st.write(f"Error generating genome track: {e}")


def locAnalyzer(subChoice, res, myChr, myStart, myEnd, cre_index):
    from collections import defaultdict
    
    subtype_labels = ["HeH", "ER", "BA", "DUX4r", "TP", "KMT2Ar", "iAMP21", "nearHaploid"]
    subtype_keys = ["HeH_10k", "ER_10k", "BA_10k", "DUX4r_10k", "TP_10k", "KMT2Ar_10k", "iAMP_10k", "nearHaploid_10k"]
    track_files = defaultdict(dict)

    for i, label in enumerate(subtype_labels):
        track_files[i+1]["allTX"] = os.path.join(TRACKS_DIR, f"tracks_{label.lower()}_all.ini")
        track_files[i+1]["canonTX"] = os.path.join(TRACKS_DIR, f"tracks_{label.lower()}_canon.ini")

    merged_tracks = {
        "1k": {"allTX": "tracks_1k_all.ini", "canonTX": "tracks_1k_canon.ini"},
        "10k": {"allTX": "tracks_10k_all.ini", "canonTX": "tracks_10k_canon.ini"},
    }

    gene_index = 9 - cre_index
    annot_index = 7 - cre_index
    myTx = "allTX" if cre_index == 0 else "canonTX"

    if subChoice == 0:
        if res == "10k":
            df = all_loops[~all_loops['loopSource'].isin(["merged_1k"])].drop_duplicates()
        else:
            df = all_loops[all_loops['loopSource'] == f"merged_{res}"]
        mySubtype = "merged"
        myTrack = os.path.join(TRACKS_DIR, merged_tracks[res][myTx])
    else:
        df = all_loops[all_loops['loopSource'] == subtype_keys[subChoice - 1]]
        mySubtype = subtype_labels[subChoice - 1]
        myTrack = track_files[subChoice][myTx]

    myStart, myEnd = int(myStart), int(myEnd)
    myList = []

    for i in range(len(df)):
        cond1 = str(df.iloc[i, 0]) == str(myChr) and df.iloc[i, 1] >= myStart and df.iloc[i, 2] <= myEnd
        cond2 = str(df.iloc[i, 3]) == str(myChr) and df.iloc[i, 4] >= myStart and df.iloc[i, 5] <= myEnd
        if cond1 or cond2:
            myList.append([
                df.iloc[i, 0], df.iloc[i, 1], df.iloc[i, 2], df.iloc[i, gene_index],
                df.iloc[i, 3], df.iloc[i, 4], df.iloc[i, 5], df.iloc[i, annot_index]
            ])

    myDf_all = pd.DataFrame(myList, columns=['chr', 'start', 'end', 'target', 'interChr', 'interStart', 'interEnd', 'type'])
    myDf = myDf_all[(myDf_all['type'] == "CRE") & (myDf_all['target'] != "no")].iloc[:, :7].drop_duplicates()
    myDf_all = myDf_all[myDf_all['target'] == "no"]

    if myDf.empty:
        st.write("No predicted enhancers were found for this location")
        return

    unique_enhancers = myDf[['chr', 'start', 'end']].drop_duplicates()
    st.write(f"There are {len(unique_enhancers)} enhancers found in this region")

    if "show_enhancers" not in st.session_state:
        st.session_state.show_enhancers = False
    if "show_loops" not in st.session_state:
        st.session_state.show_loops = False
    if "last_region" not in st.session_state:
        st.session_state.last_region = None

    if st.button("Show/Hide Enhancers"):
        st.session_state.show_enhancers = not st.session_state.show_enhancers

    enhancer_df = []
    if st.session_state.show_enhancers:
        st.write("Enhancer details:")
        for _, row in unique_enhancers.iterrows():
            chr_val, start_val, end_val = row
            associated_genes = myDf[
                (myDf['chr'] == chr_val) &
                (myDf['start'] == start_val) &
                (myDf['end'] == end_val)
            ]['target'].unique()

            unique_genes = set()
            for gene_list in associated_genes:
                if isinstance(gene_list, str):
                    unique_genes.update(gene_list.split(", "))
            genes_str = ", ".join(sorted(unique_genes))
            st.write(f"CRE on {chr_val}:{start_val}-{end_val} regulates: {genes_str}")
            enhancer_df.append([chr_val, start_val, end_val, genes_str])

        enhancer_df = pd.DataFrame(enhancer_df, columns=['chr', 'start', 'end', 'target'])
        st.download_button(
            label="Download Enhancer Details as CSV",
            data=enhancer_df.to_csv(index=False),
            file_name=f"{myChr}_{myStart}_{myEnd}_{mySubtype}_{res}_{myTx}_enhancer_results.csv",
            mime="text/csv"
        )

    extended_min_start, extended_max_end = get_genomic_range(myDf)
    region = f"chr{myDf.iloc[0, 0]}:{extended_min_start}-{extended_max_end}"

    if st.session_state.last_region != region:
        st.session_state.last_region = region
        output_file = "output_genome_track.png"
        pyGenomeTracks_command = [
            "pyGenomeTracks", "--tracks", myTrack, "--region", region, "-o", output_file
        ]
        try:
            subprocess.run(pyGenomeTracks_command, check=True)
            st.session_state.track_image = mpimg.imread(output_file)
            if os.path.exists(output_file):
                os.remove(output_file)
        except subprocess.CalledProcessError as e:
            st.write(f"Error generating genome track: {e}")
    else:
        st.session_state.track_image = st.session_state.track_image

    plt.figure(figsize=(10, 5))
    plt.imshow(st.session_state.track_image)
    plt.axis('off')
    st.pyplot(plt.gcf())

    if st.button("Show/Hide Additional Loops"):
        st.session_state.show_loops = not st.session_state.show_loops

    if st.session_state.show_loops:
        if myDf_all.empty:
            st.write("No additional loops were found in this region")
        else:
            st.write("The following additional loop interactions were also found:")
            seen = set()
            for _, row in myDf_all.iterrows():
                anchor1 = f"{row['chr']}:{row['start']}-{row['end']}"
                anchor2 = f"{row['interChr']}:{row['interStart']}-{row['interEnd']}"
                sorted_pair = tuple(sorted([anchor1, anchor2]))
                if sorted_pair not in seen:
                    seen.add(sorted_pair)
                    st.write(f"{sorted_pair[0]}    {sorted_pair[1]}")


def get_genomic_range(df):
    min_start = min(df['start'].min(), df['interStart'].min())
    max_end = max(df['end'].max(), df['interEnd'].max())
    return min_start - 100_000, max_end + 100_000
