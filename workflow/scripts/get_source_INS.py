#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 17.05.2024

From a breakpoint extract discordant reads to identify the SV sequence prigin and returns BED file.

Based on the sample BAM and callers breakpoints, this algorithm looks for discordant reads that potentially participate
to the SV event. By extracting the discordant mate (the one mapping on the SV) REF field, the origin of the SV is
identified.

Note that many different origin could be localized, in the contect of insertions, this is potentially MEI hotspot. To
avoid overannatation in case where many origins are found, the algorithm filters origins by selecting the ones having at
least half of the frequency of the most common origin. The threshold can be modified. The idea is to remove potential
noise masking real origin of SV because of a given read mapping to multiple locations.

source == origin
"""

# MODULES
from collections import Counter

import pysam
import pandas as pd

# VARIABLES
# Discordant reads origin filtering: Half of the max source count frequency
THRESHOLD = snakemake.params.threshold
# Cutoff for mapping quality of read and mate
cutoff = snakemake.params.cutoff
# Window range to extract split reads around breakpoint
window = snakemake.params.window

# FUNCTIONS


def split_read(read_bam):
    """ Check if a read pair is discordant and returns True

    A discordant read pair is a pair that is not properly paired during alignment but both reads are mapped with good
    scores which suggest potential SV event.

    :param read_bam: pysam.AlignedSegment instance
    :return: Boolean
    """
    if (not read_bam.is_secondary
            and not read_bam.is_proper_pair
            and not read_bam.is_unmapped
            and not read_bam.mate_is_unmapped
            and read_bam.mapping_quality >= cutoff):

        return True


def select_source(source_list, THRESHOLD):
    """Normalization of origin frequencies given the most frequent one, returns strings

    Normalized frequencies by dividing all counts for each origin by the most frequent origin count. Then, filtering
    step via a user defined threshold representing the minimal frequency to reach compared to the most frequent. Finally
     returns a qualification value: Multiple (if more than 3 origin), Origin id (if less than 3).

    :param source_list: List if int64 containing each discordant mate origin id.
    :param THRESHOLD: float representing the minimal frequency compared to the most frequent origin. (between 0 and 1)
    :return: Strings
    """
    # Counts of each putative source name of inserted sequences
    counts = Counter(source_list)

    # Normalized frequencies based on the most frequent source name
    max_count = max(counts.values())
    normalized = {key: count / max_count for key, count in counts.items()}

    # Select the most frequent source name
    selected = [key for key, value in normalized.items() if value >= THRESHOLD]

    # Annotation based on selected source names
    # If only one source name
    if len(selected) == 1:
        return "".join(selected)
    # If more than one selected but less than 4, return all putative source names
    elif 1 < len(selected) <= 3:
        return ",".join(selected)
    # If more than 3
    else:
        return "Multiple"


# MAIN PROGRAM
# Sample bam file
bam = pysam.AlignmentFile(str(snakemake.input.bam), "rb")
# Sample SV bed file
bed = pd.read_csv(str(snakemake.input.bed),
                  delimiter="\t",
                  header=None,
                  names=["chr","start","stop","sv_type","len","type"])

# Iter over breakpoints
for index, bp in bed.iterrows():
    # Extract reads around breakpoint in a certain window around
    read_around = bam.fetch(str(bp["chr"]),
                            int(bp["start"]) - window,
                            int(bp["stop"]) + window)
    # Split reads RNEXT in current window
    forward = []
    reverse = []
    # iter over reads
    for read in read_around:
        if split_read(read) and bam.mate(read).mapping_quality >= cutoff:
            if read.is_forward:
                reverse.append(read.next_reference_name)
            else:
                forward.append(read.next_reference_name)

    # Only checking sources if supporting split reads on left and righ breakpoint
    if reverse and forward:
        # Keep source names supported by forward and reverse split reads
        forward_set = set(forward)
        reverse_set = set(reverse)
        both_side_splits = [source for source in
                            (forward + reverse) if source in reverse_set.intersection(forward_set)]
        # Check if any reads supporting in both sides
        if both_side_splits:
            # Annotate source in BED
            bed.loc[index, "source"] = select_source(both_side_splits, THRESHOLD)
        else:
            # Source discordant from both side
            bed.loc[index, "source"] = "Rev-For_Discordant"
    else:
        # Not enough support split reads
        bed.loc[index, "source"] = "Missing"

# Write new BED with annotated source for Insertions
bed.to_csv(str(snakemake.output[0]),
           sep="\t",
           header=False,
           index=False)
