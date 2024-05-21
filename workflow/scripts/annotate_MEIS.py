#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 16.05.2024

Merge records from BED SVs with ones annotated with scramble MEIS via pandas. Returns an annotated SV BED with new
column 'Mei_type'.
"""

# MODULES
from os.path import getsize

import pandas as pd

# VARIABLES
# Contains SVs in BED format
bed_SV = str(snakemake.input.bed_SV)
# Contains only SVs from bed_SV that overlap with MEIs from scramble
bed_w_MEIS = str(snakemake.input.bed_w_MEIS)
# Output file containing both information
out_bed = str(snakemake.output.final_bed)

# MAIN
# BED to DF
SV = pd.read_csv(bed_SV, delimiter="\t", header=None,
                 names=["chr", "start", "stop","sv_type","len"])

# CHeck if any MEIs in annotation file
if getsize(bed_w_MEIS) > 0:
    # MEIs BED to DF
    MEIS = pd.read_csv(bed_w_MEIS, delimiter="\t", header=None,
                       names=["chr","start","stop","type"])
    # Set chr field as str
    MEIS.chr = MEIS.chr.astype(str)
    SV.chr = SV.chr.astype(str)
    # Merge both df, replace no matching rows type with "."
    merged = SV.merge(MEIS, how="left",
                      on=["chr", "start", "stop"]).fillna(".")
    # Drop duplicates but keep SV matching different MEIS sites and vis-versa
    merged_no_dup = merged.drop_duplicates()
    # Write to BED
    merged_no_dup.to_csv(out_bed, sep="\t", header=False, index=False)
# File empty, set "." type because no MEIs detected
else:
    SV["type"] = "."
    SV.to_csv(out_bed, sep="\t", header=False, index=False)