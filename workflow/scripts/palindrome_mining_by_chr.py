#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 17.05.2024

Mine all palindromic sequences in a fasta file and extract interesting data for downstream analysis, returns a BED file.

The algorithm read each base in a fasta file and tries to identify a palindromic sequence seed from this point. The seed
length is user defined. Once a seed is identified, the algorithm tries to etend this palindromic sequence as it can,
until reaching ending criteria:
    - Spacer length > 20 bp
    - Missmatch rate > 15% between left and right sequence
    - End of palindrome

Then, the start and end positions are extracted with other features related to the current palindrome as:
    - Palindrome type: Near or Perfect
    - Palindrome total length (with spacer)
    - Palindrome spacer length
    - Missmatch rate
    - At% of palindromic sequence
"""

# MODULES
import time
import sys
import re


# FUNCTIONS


def complement_dna(dna_seq):
    """ Complementary DNA sequence, returns a string.

    :param dna_seq: String
    :return: String
    """
    # The constant variable translation_table
    global TRANS_TAB

    # translate sequence via the table to complementary dna sequence
    complement = dna_seq.translate(TRANS_TAB)
    return complement


def isPalindrome(dna_seq, midpoint):
    """Check if a DNA sequence in palindromic and return boolean.

    From the midpoint of the sequence, split DNA sequence in two and check if reverse complementary of left arm is equal
    to the right arm sequence.

    :param dna_seq: String
    :param midpoint: Int64
    :return: Boolean
    """
    if dna_seq[:midpoint] == complement_dna(dna_seq[midpoint:])[::-1]:
        return True
    else:
        return False


# From a core of N nuc palindrom, we expand to each side of the core until reaching the
# longest palindromic sequence.
# Perfect and semi-perfect palindrome are mined from SNA sequence but not with spacer.
def LPS_brute(mysequence: str, start_pal: int, end_pal: int, miss_rate: float):
    """Mine the longest palindromic sequence given a palindromic seed sequence. Returns the start and end positions of
    the longest palindrome with its missmatch ratio.

    Extend a palindromic seed sequence to the most left and right until the mined sequence is no more palindromic or the
    missmatch ratio is reached.

    :param mysequence: Reference genome in a string format
    :param start_pal: Seed starting position relative to reference genome
    :param end_pal: Seed ending position relative to reference genome
    :param miss_rate: Threshold for missmatch rate (between 0 and 1)
    :return: A tuple as following (most_left_pos; most_right_pos, miss_rate)
    """
    # Mismatch count
    mismatch = 0
    # Palindromic seed sequence length without spacer
    global MIN_PAL_NUC
    pal_len = (MIN_PAL_NUC * 2) - 2

    # Border conditions to stop the while loop
    # If the beginning of the pal is reaching the end of the main sequence
    # if the end of pal is reaching the end of the main sequence
    # if left and right bases are not palindromic
    while start_pal >= 0 and end_pal <= len(mysequence):
        # in the case of a match
        if mysequence[start_pal] == complement_dna(mysequence[end_pal - 1]):
            # the left border goes to the left
            start_pal -= 1
            # the right border goes to the right
            end_pal += 1
            # Increase len pal by 2 (1 more nuc on each side)
            pal_len += 2
            # Current miss_rate of palindrome
            try:
                # Mismatch multiplied by two because divided by the whole sequence length
                curr_miss_rate = round((mismatch * 2) / pal_len, 4)
            # No mismatch yet
            except ZeroDivisionError:
                curr_miss_rate = 0

            last_match = (start_pal, end_pal, curr_miss_rate)

        # in the case of a mismatch
        else:
            mismatch += 1
            # Check for missmatch condition and look for following bases if under.
            if (mismatch * 2) / pal_len <= miss_rate:
                # the left border goes to the left
                start_pal -= 1
                # the right border goes to the right
                end_pal += 1
            # Condition met so end
            else:
                return ((last_match[0] + 1), (last_match[1] - 1), last_match[2])

        # Case in which the right interval reaches end main sequence
        if end_pal - 1 == len(mysequence) or start_pal + 1 == 0:
            return ((last_match[0] + 1), (last_match[1] - 1), last_match[2])


# Returns the full sequence in a fasta file into a string
def read_fasta(file):
    """Transform to string a TXT or FASTA file. Returns a string

    :param file: FASTA or TXT file
    :return: String
    """
    with open(file, 'r') as fasta:
        sequence = "".join(line.strip() for line in fasta if not line.startswith(">"))
    return sequence.strip()


# N*logN
def unique_palindromes(list_pal: list):
    """Removes palindromes interval that are included within another bigger palindrome interval. Returns a cleared list.

    Algorithm from https://www.geeksforgeeks.org/check-interval-completely-overlaps/

    :param list_pal:
    :return:
    """
    # Sort the list of palindromes to make easy the comparison between each interval edges
    # Sort acccording to starting position pal = [(start,end),(strat,...)...]
    list_pal.sort(key=lambda tup: (tup[0], tup[3]))

    # Keep track of the last index in modified array
    modified = [list_pal[0]]
    index = 0

    # As the list is sorted the right border is always smaller or equal in modified list
    # compared to current list. Just have to compare the left border to know if the current interval is
    # contained in modified interval
    for item in list_pal[1:]:
        # Contained so don't keep it
        if item[1] <= modified[index][1]:
            continue
        # add the current interval as not contained and update index to its index
        # for futur comparisons
        else:
            index += 1
            modified.append(item)
    return modified


# Compute the AT-richness of a sequence in percent
def AT_richness(dna_seq: str):
    return ((dna_seq.count("A") + dna_seq.count("T")) / len(dna_seq)) * 100


if __name__ == "__main__":

    # Inputs from CL
    # Chromosome fasta sequence to mine
    data_file = str(snakemake.input.chr_fa)
    # Chromosome numer as chr22
    output_file = str(snakemake.output[0])

    # Read Fasta file input from CL
    MYSEQUENCE = read_fasta(data_file)
    # index check 4 nucleotides around
    MIN_PAL_NUC = snakemake.params.seed_pal_len
    # Container for perfect palindroms
    palindromes = []
    # Threshold for lack of identity
    MISS_RATE = snakemake.params.miss_rate
    # The maximum spacer length
    MAX_SP = snakemake.params.max_sp
    # translation table using hashtable
    TRANS_TAB = str.maketrans("atgcATGCNn", "tacgTACGNn")

    # Main loop to go through the genomic sequence all possible midpoints for a palindrome (only even number)
    # for now -> only checking for perfect palindromes
    for midpoint in range((MIN_PAL_NUC), (len(MYSEQUENCE) - MIN_PAL_NUC - 1)):
        # Ignore Ns in the fatsa file as they are unsequenced nucleotides within the gennome
        if MYSEQUENCE[midpoint] == "N":
            continue

        # fixing the window in which we will look for a palindrome of 8 nucleotides minimum.
        # end and start sequence represent the index of each bornes within the main sequence
        end_pal = midpoint + MIN_PAL_NUC
        start_pal = midpoint - MIN_PAL_NUC
        # retrieve the sub-sequence from the main sequence
        sub_seq = MYSEQUENCE[start_pal:end_pal]

        # check if it is a palindrom. MIN_PAL_NUC -> the length of one arm which is also the midpoint of the palindrome
        if isPalindrome(sub_seq, MIN_PAL_NUC):
            # if the subsew is a pal check if it is a core for a longer pal
            # append the result to the list of palindromes
            # Here, LPS_brute returns a tuple (start, end) in which we add the soacer length (0) in this case (0 because perfect or nearpal without spacer)
            # So we have in final tuple (start, end, curr_miss_rate, spacer)
            tup = LPS_brute(MYSEQUENCE, start_pal, end_pal, MISS_RATE) + (0,)
            palindromes.append(tup)
        # look for potiential inverted repeats (Palindrome or Nearpalindrome with spacer)
        # Minimum len for spacer is 1 and maximum s max_sp as range doesn't include right interval
        # Algo check from the smallest to the maximum spacer length
        else:
            for sp in range(1, MAX_SP + 1):
                # Remove the midpoint as putative spacer and construct the new 8 core putative palindrome
                if end_pal + sp <= len(MYSEQUENCE):
                    sp_seq = MYSEQUENCE[start_pal:midpoint] + MYSEQUENCE[midpoint + sp: end_pal + sp]

                # check if this new 8 nuc core is a palidrome and extand it if it is the case
                if isPalindrome(sp_seq, MIN_PAL_NUC):
                    tup = LPS_brute(MYSEQUENCE, start_pal, end_pal + sp, MISS_RATE) + (sp,)
                    palindromes.append(tup)
                    break

    # Merge overlapping palindromes
    cleared_palindromes = unique_palindromes(palindromes)

    with open({output_file}, "w") as file:
        # add palindrome type information
        type = ""
        # First pal border is inclusive and 0 based
        # Last pal border is not inclusive
        for it, pal in enumerate(cleared_palindromes):
            midpoint = len(MYSEQUENCE[pal[0]:pal[1]]) // 2
            sub_seq = MYSEQUENCE[pal[0]:pal[1]]

            if isPalindrome(sub_seq, midpoint):
                type = "Perfect"
            else:
                type = "Near"

            file.write(f"{snakemake.wildcards.chr}\t{pal[0]}\t{pal[1]}\t{type}"
                       f"\t{len(sub_seq)}\t{pal[3]}\t{pal[2]*100:.2f}%\t{AT_richness(sub_seq):.2f}%\n")