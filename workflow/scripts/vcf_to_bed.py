# DESCRIPTION
# From a SV vcf file build bed file containing
#   SVTYPE=INS
#   SVLEN>=50 OR NA
#   CHROM (1-22 + XY)
#   START = POS - 1 / STOP = POS

# MODULES
import re

# VARIABLES
# Match digit for svlen
RE_SVLEN = r"SVLEN=(\d+)"
# Match sv type
RE_SVTYPE = r"SVTYPE=(\w+)"
# Match chrom
RE_REF = r"^\S+"
# Match pas
RE_POS = r"^\w+\s+(\d+)"
# Keep only these chrom ref
CHROM = [str(val) for val in range(1, 23)] + ["X", "Y"]
# VCF file from snakemake rule
vcf = snakemake.input[0]
# Bed Out name
bed = snakemake.output[0]

# MAIN
# Read non compressed vcf file
with open(vcf, "r") as in_file, open(bed, "w") as out_file:
    for record in in_file:
        # Ignore header
        if record.startswith("#"):
            continue
        else:
            ref = re.findall(RE_REF, record)
            type = re.findall(RE_SVTYPE, record)
            size = re.findall(RE_SVLEN, record)
            # Keep only INS or BND within CHROM
            if type[0] == "INS" and ref[0] in CHROM:
                # Check if len empty
                size_check = size[0] if size else "."
                # Match pos sv
                start = int(re.findall(RE_POS, record)[0]) - 1
                stop = int(re.findall(RE_POS, record)[0])
                # Write record
                out_file.write(f"{ref[0]}\t{start}\t{stop}\t{type[0]}\t{size_check}\n")
