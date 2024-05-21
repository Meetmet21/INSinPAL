#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 17.05.2024

From scramble MEIS txt file to bed file with chr, pos - 1, pos, MEI type fields.
"""

# MODULES
import os

import pandas as pd

# VARIABLES
tsv = str(snakemake.input.meis)
df = pd.read_csv(tsv, delimiter="\t")

# MAIN
# Extract from Insertion column chromsome and position (initial format chr:POS)
df[["chr", "pos"]] = df.Insertion.str.split(":", expand=True)
# Remove contig chromsomes names as 1_bla
df = df[-df["chr"].str.contains(".*_")]
# Remove MT
df = df[-df["chr"].str.contains("MT")]

# Write to bed
with open("tmp", 'w') as tmp:
    for index, row in df[["chr","pos","MEI_Family"]].iterrows():
        # Write chr, pos -1 , pos and type
        tmp.write(f"{row.chr}\t{int(row.pos) - 1}\t{row.pos}\t{row.MEI_Family}\n")

# Sort bed to chromosomal order
command = f"cat tmp | sort -V > {str(snakemake.output.bed)}; rm tmp"
os.system(command)