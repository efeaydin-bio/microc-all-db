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
import matplotlib
import pyBigWig
import shutil
import requests

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
def fetch_bigwig_locally(url):
    filename = os.path.basename(url)
    cache_dir = "cached_bigwigs"
    os.makedirs(cache_dir, exist_ok=True)
    local_path = os.path.join(cache_dir, filename)

    if not os.path.exists(local_path):
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    return local_path

# bigwig loading
@st.cache_resource
def load_bigwig_tracks():
    bw_files = {
        "h3k4me1": fetch_bigwig_locally("https://data.cyverse.org/dav-anon/iplant/home/efeaydin/h3k4me1.bigWig"),
        "h3k4me3": fetch_bigwig_locally("https://data.cyverse.org/dav-anon/iplant/home/efeaydin/h3k4me3.bigWig"),
        "h3k27ac": fetch_bigwig_locally("https://data.cyverse.org/dav-anon/iplant/home/efeaydin/h3k27ac.bigWig"),
        "h3k27me3": fetch_bigwig_locally("https://data.cyverse.org/dav-anon/iplant/home/efeaydin/h3k27me3.bigWig"),
        "dnase": fetch_bigwig_locally("https://data.cyverse.org/dav-anon/iplant/home/efeaydin/dnase.bigWig"),
    }

    return {key: pyBigWig.open(path) for key, path in bw_files.items()}

bigwig_tracks = load_bigwig_tracks()
h3k4me1 = bigwig_tracks["h3k4me1"]
h3k4me3 = bigwig_tracks["h3k4me3"]
h3k27ac = bigwig_tracks["h3k27ac"]
h3k27me3 = bigwig_tracks["h3k27me3"]
dnase = bigwig_tracks["dnase"]

# function to write tempfiles from bigWigs
def write_temp_bedgraph(track, chrom, start, end):
    intervals = track.intervals(chrom, start, end)
    temp_path = tempfile.NamedTemporaryFile(delete=False, suffix=".bedGraph").name
    with open(temp_path, "w") as out:
        for interval in intervals:
            out.write(f"{chrom}\t{interval[0]}\t{interval[1]}\t{interval[2]}\n")
    return temp_path

# function for gene targeted query
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
        temp_links_file.write(
        f"{row['chr']}\t"
        f"{int(row['start'])}\t"
        f"{int(row['end'])}\t"
        f"{row['interChr']}\t"
        f"{int(row['interStart'])}\t"
        f"{int(row['interEnd'])}\t"
        "1.0\n")

    temp_links_file.close()

    # fix genes
    new_bed_file = coding_genes[coding_genes.iloc[:, 3] != gene_of_interest].copy()
    new_bed_file.loc[new_bed_file.iloc[:, 3] == gene_of_interest, new_bed_file.columns[6]] = "255,0,0"
    new_bed_file.insert(6, "thickStart", 1000)
    new_bed_file.insert(7, "thickEnd", 2000)
    new_bed_file.to_csv(os.path.join(TRACKS_DIR, "tempCodingGenes.bed"), index=False, sep="\t", header=False)

    only_gene = coding_genes[coding_genes.iloc[:, 3] == gene_of_interest].copy()
    only_gene = only_gene.iloc[:, [0, 1, 2, 3, 4, 5]]  # BED6
    only_gene.to_csv("tracks/onlyTargetGene.bed", sep="\t", index=False, header=False)
    # fix bigwig
    # define regions of interest
    chrom = f"chr{myDf.iloc[0, 0]}"
    start = int(extended_min_start)
    end = int(extended_max_end)
    bw_region = f"{chrom}:{start}-{end}"

    temp_h3k4me1_path = write_temp_bedgraph(h3k4me1, chrom, start, end)
    temp_h3k4me3_path = write_temp_bedgraph(h3k4me3, chrom, start, end)
    temp_h3k27ac_path = write_temp_bedgraph(h3k27ac, chrom, start, end)
    temp_h3k27me3_path = write_temp_bedgraph(h3k27me3, chrom, start, end)
    temp_dnase_path   = write_temp_bedgraph(dnase, chrom, start, end)
    
    temp_tracks_path = os.path.join(TRACKS_DIR, "temp_gene_tracks.ini")
    with open(temp_tracks_path, "w") as temp_tracks_file:
        # x-axis
        temp_tracks_file.write("[x-axis]\n")
        temp_tracks_file.write("height = 4\n")
        temp_tracks_file.write("fontsize = 14\n")
        temp_tracks_file.write("title = hg38\n\n")
        # enhancer links
        temp_tracks_file.write("[enhancer_links]\n")
        temp_tracks_file.write(f"file = {temp_links_file.name}\n")
        temp_tracks_file.write("file_type = links\n")
        temp_tracks_file.write("line_width = 3\n")
        temp_tracks_file.write("color = red\n")
        temp_tracks_file.write("title = Enhancer Links\n")
        temp_tracks_file.write("height = 5\n\n")
        # only gene
        temp_tracks_file.write("[genes2]\n")
        temp_tracks_file.write("file = tracks/onlyTargetGene.bed\n")
        temp_tracks_file.write("file_type = bed\n")
        temp_tracks_file.write("color = red\n")
        temp_tracks_file.write("height = 2\n")
        temp_tracks_file.write("title = \n")
        temp_tracks_file.write("fontsize = 12\n")
        temp_tracks_file.write("arrow_interval = 5\n")
        temp_tracks_file.write("gene_rows = 1\n\n")
        #  Genes
        temp_tracks_file.write("[genes]\n")
        temp_tracks_file.write("file = tracks/tempCodingGenes.bed\n")
        temp_tracks_file.write("file_type = bed\n")
        temp_tracks_file.write("color = red\n")
        temp_tracks_file.write("height = 2\n")
        temp_tracks_file.write("max_labels = 20\n")
        temp_tracks_file.write("overlay_previous = no\n")
        temp_tracks_file.write("title = Genes\n")
        temp_tracks_file.write("fontsize = 12\n")
        temp_tracks_file.write("arrow_interval = 5\n")
        #temp_tracks_file.write("gene_rows = 10\n\n")
        # Promoters
        temp_tracks_file.write("[promoters]\n")
        temp_tracks_file.write("file = data/promoters.bed\n")
        temp_tracks_file.write("file_type = bed\n")
        temp_tracks_file.write("color = red\n")
        temp_tracks_file.write("height = 2\n")
        temp_tracks_file.write("merge_overlapping_exons: true\n")
        temp_tracks_file.write("title = Promoters\n")
        temp_tracks_file.write("display = collapsed\n")
        temp_tracks_file.write("labels: false\n")
        temp_tracks_file.write("merge_transcripts: true\n")
        # h3k4me1
        temp_tracks_file.write("[h3k4me1]\n")
        temp_tracks_file.write(f"file = {temp_h3k4me1_path}\n")
        temp_tracks_file.write("file_type = bedgraph\n")
        temp_tracks_file.write("color = green\n")
        temp_tracks_file.write("height = 4\n")
        temp_tracks_file.write("title = H3K4me1\n")
        temp_tracks_file.write("min_value = 0\n")
        temp_tracks_file.write("max_value = 30\n\n")
        # h3k4me3
        temp_tracks_file.write("[h3k4me3]\n")
        temp_tracks_file.write(f"file = {temp_h3k4me3_path}\n")
        temp_tracks_file.write("file_type = bedgraph\n")
        temp_tracks_file.write("color = green\n")
        temp_tracks_file.write("height = 4\n")
        temp_tracks_file.write("title = H3K4me3\n")
        temp_tracks_file.write("min_value = 0\n")
        temp_tracks_file.write("max_value = 30\n\n")
        # h3k27ac
        temp_tracks_file.write("[h3k27ac]\n")
        temp_tracks_file.write(f"file = {temp_h3k27ac_path}\n")
        temp_tracks_file.write("file_type = bedgraph\n")
        temp_tracks_file.write("color = green\n")
        temp_tracks_file.write("height = 4\n")
        temp_tracks_file.write("title = H3K27ac\n")
        temp_tracks_file.write("min_value = 0\n")
        temp_tracks_file.write("max_value = 30\n\n")
        # h3k27me3
        temp_tracks_file.write("[h3k27me3]\n")
        temp_tracks_file.write(f"file = {temp_h3k27me3_path}\n")
        temp_tracks_file.write("file_type = bedgraph\n")
        temp_tracks_file.write("color = green\n")
        temp_tracks_file.write("height = 4\n")
        temp_tracks_file.write("title = H3K27me3\n")
        temp_tracks_file.write("min_value = 0\n")
        temp_tracks_file.write("max_value = 30\n\n")
        # dnase
        temp_tracks_file.write("[dnase]\n")
        temp_tracks_file.write(f"file = {temp_dnase_path}\n")
        temp_tracks_file.write("file_type = bedgraph\n")
        temp_tracks_file.write("color = grey\n")
        temp_tracks_file.write("height = 4\n")
        temp_tracks_file.write("title = DNase\n")
        temp_tracks_file.write("min_value = 0\n")
        temp_tracks_file.write("max_value = 1\n\n")

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
        subprocess.run(pyGenomeTracks_command, check=True)
        img = mpimg.imread(output_file)
        plt.figure(figsize=(10, 5)) 
        plt.imshow(img)
        plt.axis('off')
        st.pyplot(plt.gcf())

        os.remove(output_file)
        os.remove(temp_links_file.name)
        os.remove(temp_tracks_file.name)
        os.remove(os.path.join(TRACKS_DIR, "tempCodingGenes.bed"))
        os.remove(os.path.join(TRACKS_DIR, "onlyTargetGene.bed"))
        os.remove(temp_h3k4me1_path) 
        os.remove(temp_h3k4me3_path) 
        os.remove(temp_h3k27ac_path) 
        os.remove(temp_h3k27me3_path) 
        os.remove(temp_dnase_path) 
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
