# DESCRIPTION
# From TP BED compute metrics for benchmarking

# MODULES

# Functions
# Compute recall of caller
def recall(TP, call_set):
    return TP / call_set


# Compute precision of caller
def precision(TP, call_caller):
    return TP / call_caller


# VARIABLES
# Caller name
caller = snakemake.params.caller
# TP SV BED annotated
tp_bed = snakemake.input.tp_bed
# Trueset BED
bench_bed = snakemake.input.bench_bed
# Out summary file
outfile = snakemake.output.summary

# MAIN
# Count calls in trueset
with open(bench_bed, "r") as file:
    call_set = len(file.readlines())

# Generate summary file for benchmark of caller
with open(outfile, "w") as out:
    # Number of TP
    # Count calls in caller
    with open(tp_bed, "r") as file:
        TP = len(file.readlines())

    out.write(f"{caller}:\n"
              f"\tTP:{TP}\n"
              f"\tFP:{call_set - TP}\n"
              f"\tRecall: {round(recall(TP, call_set), 4)}")