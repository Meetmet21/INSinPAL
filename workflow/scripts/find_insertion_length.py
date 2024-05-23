#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 16.05.2024

From a uniq source of insertionin a given breakpoint, extract discordant reads and compute putative insertion length.
Returns the input BED file with modified "len" field if a len is computed.

The idea is that on each side of the breakpoint, discordant reads mate mapping on source ref delimits the insertion
boundaries. By extracting them and subtracting their position, we have a putative insertion length. To have a consensus
on the mapping ref positions, standard devition is assessed among all reverse and forward discordant mates. A great
standard deviation runs a filtering algorithm to keep relatively close positions.
"""

# MODULES
import numpy as np
import pysam
import pandas as pd

# VARIABLES
# Window around which extract the discordant reads (correspond to a read length (read1 + insert)
window = snakemake.params.window
output_file = str(snakemake.output[0])
MAPQ_CUTOFF = snakemake.params.mapq_cutoff
STD_CUTOFF = snakemake.params.std_cutoff
# Read bam file
bam = pysam.AlignmentFile(str(snakemake.input.bam), "rb")
# Read bed file
bed = pd.read_csv(str(snakemake.input.bed), delimiter="\t", header=None,
                  names=["chr","start","stop","sv_type","len","MEI","source"])
# Replace all equals by chr name
bed.loc[bed["source"] == "=", "source"] = bed.chr

# FUNCTIONS


def discrodant_reads(read_bam):
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
            and read_bam.mapping_quality >= MAPQ_CUTOFF):
        return True


def remove_extreme_values(split_list, STD_CUTOFF):
    """Compute the standard deviation of a list and removes extreme points. Returns a list.

    From a list of position corresponding to discordant mate positions on inserted sequence, seek for sublist positions
    in which the standard deviation is low (meaning reads clustering around the same boundaries), removes the extreme
    points and returns new list.

    :param split_list: list of positions
    :param STD_CUTOFF: Threshold for the standard deviation
    :return: List
    """
    # Compute list std
    std = np.std(split_list)

    # Initial list already clustering on the same position
    if std < STD_CUTOFF:
        # No extreme values
        return split_list

    # No consensus clustering in initial list
    else:
        # Sort list to nail extreme values to list edges
        split_list.sort()
        # Calculate distances between each neighbour pairs from left to right
        # Remove one because last index already taken in account with index - 1
        distances = [split_list[i + 1] - split_list[i] for i in range(len(split_list) - 1)]
        # Check for how many extreme values in list
        # Distance between close values should be in the range of STD_CUTOFF otherwise they are divergent
        count_extreme = sum(distance > STD_CUTOFF for distance in distances)
        # Check for how many close values in list
        count_close = len(split_list) - count_extreme

        # Remove extreme values if at least one pair close positions in list
        if count_close >= 2:
            close_positions = []
            # Pairs of positions so don't consider the last index
            for index in range(len(distances) - 1):
                # Keep pairs with distance in STD_CUTOFF
                if distances[index] < STD_CUTOFF:
                    # Add each pair members
                    # Add duplicates but remove with set
                    close_positions.append(split_list[index])
                    close_positions.append(split_list[index + 1])
            # Remove duplicates
            return set(close_positions)

        # No clustering in whole list
        else:
            return False


# MAIN PROGRAM
# Iter over breakpoints
for index, bp in bed.iterrows():
    # Many source of inserted sequence in breakpoint, potential MEI site, so skip
    if (bp.source in ["Multiple", "Missing_Splits", "Rev-For_Discordant"]
            or len(str(bp.source).split(",")) > 1):
        continue
    # One potential source of inserted sequence in breakpoint
    else:
        # Extract reads around breakpoint in a certain window around
        read_around = bam.fetch(str(bp["chr"]),
                                int(bp["start"]) - window,
                                int(bp["stop"]) + window, )
        # Keep RNEXT track of discordant mate
        forward = []
        reverse = []
        for read in read_around:
            if discrodant_reads(read) and bam.mate(read).mapping_quality >= MAPQ_CUTOFF:
                # Check for discordant reads with mates mapping on the source of inserted sequence
                if read.next_reference_name == bp.source:
                    # Add mate position to reverse or forward depending on orientation
                    if read.is_forward:
                        reverse.append(int(read.next_reference_start))
                    else:
                        forward.append(int(read.next_reference_start))

        # Compute insertion size if enough supporting discordant mate positions
        if reverse and forward:
            # Seek for clustering around one position for both orientation
            cleared_reverse = remove_extreme_values(reverse, STD_CUTOFF)
            cleared_forward = remove_extreme_values(forward, STD_CUTOFF)
            # Check for data presence
            if cleared_forward and cleared_reverse:
                # Size of insertion is equal to the most left RNEXT discordant mate - the most right discordant mate
                # which corresponds to the boundaries of the insertion close to breakpoint.
                # Add + 100 because Forward mate position is based on the left most nucleotide, to take in account, the
                # discordant mate length is around 100 bp.
                bed.loc[index, "len"] = abs(max(cleared_forward) - min(cleared_reverse) + 100)
            # Not enough supporting data
            else:
                continue

# Write to Output
bed.to_csv(output_file, index=False, header=False, sep="\t")
