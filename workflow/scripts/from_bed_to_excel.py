#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 17.05.2024

Take a BED file and make it an excel file, returns xlsx format file.
"""
# MODULES
import pandas as pd

# VARIABLES
# Input Bed file
bed = str(snakemake.input[0])
out = str(snakemake.output[0])
# Sample ID
sample = bed.split("/")[0]

# MAIN PROGRAM
# To df
df_bed = pd.read_csv(bed,
                     delimiter="\t")
# To excel
df_bed.to_excel(snakemake.output[0],
                index=False,
                sheet_name=f"{sample}_SV_INS")