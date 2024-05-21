# DESCRIPTION
# From the merged SV BED file, removed duplicates calls from different caller
# in a given threshold
import numpy as np
# Modules
import pandas as pd

# VARIABLES
# Merged SV BED
BED = snakemake.input[0]
# DF of bed
df_bed = pd.read_csv(BED, delimiter="\t", header=None,
                     names=["chr","start","stop","type","len","MEI","source"],
                     dtype={"chr":str, "start":int, "stop":int, "type":str, "len":str, "MEI":str, "source":str})
# Merge within this window
window = snakemake.params.window
# Nodup dataframe
df_nodup = []
# CHeck if next SV overlaps/neighbors the current SV
# Remove duplicates if it's the case
for index, sv in df_bed.iterrows():
    # CHeck for next SV
    if index + 1 < len(df_bed):
        next_sv = df_bed.iloc[index + 1]
    else:
        break

    # If duplicate skip row
    if sv.chr == next_sv.chr and (next_sv.start in list(range(sv.start - 10, sv.start + window)) or
                                  next_sv.start in list(range(sv.start, sv.start + window))):
        continue
    else:
        df_nodup.append(sv)

# Write result to BED file
pd.DataFrame(df_nodup).to_csv(snakemake.output[0],  sep="\t", header=False, index=False)