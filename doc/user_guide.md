INSinPAL User Guide
===================

## Table of content

[//]: #

* [Introduction](#Introduction)
* [Palindromic fragile sites](#Palindromic-fragile-sites)
  * [Palindrome mining algorithm](#Palindrome-mining-algorithm)
  * [Filtering step](#Filtering-step)
* [SV callers](#SV-callers)
  * [Benchmarking](#Benchmarking)
  * [INSurVeyor](#INSurVeyor)
  * [Manta](#Manta)
  * [Basil](#Basil)
* [Annotations](#Annotations)
  * [MEIs annotation](#MEIs-annotation)
  * [Source of inserted sequence](#Source-of-inserted-sequence)
  * [Length of inserted sequence](#Length-of-inserted-sequence)
  * [Functionnal annotations](#Functionnal-annotations)
* [Installation](#Installation)  
* [Execution](#Execution)

[//]: #


## Introduction

When working with Whole Genome Sequencing (WGS) data for diagnosis purposes, it is essential to implement filtering steps to retrain only relevant mutations/calls made by variant callers. Otherwise, thousands of them would need to be filtered by experts. INSinPAL is a Snakemake workflow for analyzing, filtering and formatting large insertion calls in mapped paired-reads data (BAM). INSinPAL utilizes two structural variant (SV) callers ([Manta](https://github.com/Illumina/manta) and [Basil](https://github.com/seqan/anise_basil)) and one insertion-specialized caller ([INSurVeyor](https://github.com/kensung-lab/INSurVeyor)) to form a metacaller for germline large insertion variant detection. The main idea behind INSinPAL is to focus variant selection on palindromic fragile sites, as they are prone to forming DNA secondary strucures and undergoing chromosomal rearrangment. Two modules enable the detection of source and length of the inserted sequence at a given breakpoint, providing more infnromation for variant analysis by experts.


## Palindromic fragile sites

In a genomic context, a palindrome consists of two inverted repeated DNA sequences, where the reverse complement of one sequence matches the sequence of the other. These regions are prone to forming DNA secondary structures through intrastrand base pairing. The resulting structures, such as hairpins on single-strand DNA or cruciforms on double-strand DNA, create hotspots for chromosomal rearrangements. This occurs because they can lead to double-strand breaks (DSBs) through replication fork stalling ([*Mikleni and al.*](https://www.mdpi.com/1422-0067/22/6/2840)). If DSBs are not repaired or are misrepaired, it can lead to deletions. translocations, inversions and large insertions. These alterations can distrupt the genetic functions of nearby genes, potentially leading to various genetic disorders and rare diseases ([*J. Bissler*](https://www.imrpress.com/journal/FBL/3/4/10.2741/A284)).

Thus, understanding the mechanisms leading to palindrome instability and the stability of their secondary structures is crucial for detecting palindromic fragile sites. Indeed, many palindromes do not threaten genome stability. Instead, they contribute to the dynamics of gene regulation. In the literature, three essential characteristics of palindromes are commonly related to their instability (forming secondary strucutures): the length of one repeat, or a stem, the presence and length of a spacer (sequence between both stems) and the sequance mismatching between both stems. The longer the stem, the shorter the spacer, and the higher the identity between stems, the more likely a palindrome is th be unstable and form a secondary structure. These three features also similarly contribute to the stability of the secondary structure ([*Wang and al.*]((https://www.sciencedirect.com/science/article/pii/S0014579306000986))).


### Palindrome mining algorithm

The first step of the process is to find all palindromic sites in a reference genome. This is done by the rules in `mine_palindromes_in_genome.smk` and the main algorithm is in ths script `palindrome_mining_by_chr.py`. Basically the algorithm follows these steps:

```
Reference sequence in string: ref
minimum stem length: min_stem
maximum spacer length: max_spacer
FOR each central_position_of_a_potential_palindrome in all_positions_in_ref:
  IF central_position_nucleotide IS NOT sequenced:
    THEN skip
  ENDIF

  start of seed palindorme sequence in ref: start = midpoint - min_stem
  idem end: end = mindpoint + min_stem
  EXTRACT putative_minimum_palindrome FROM ref VIA start, end

  IF putative_minimum_palindrome is palindromic:
    CALL longest_palidromic_sequence BY extanding putative_minimum_palindrome edges
    ADD putative_minimum_palindrome_with_spacer

  ELSE:
    FOR each putative_spacer_length in max_spacer:
      IF putative_minimum_palindrome_with_spacer is palindromic:
        CALL longest_palidromic_sequence BY extanding putative_minimum_palindrome_with_spacer edges
        ADD putative_minimum_palindrome_with_spacer
  ENDIF
```

Note that, based on the literature, a pre-filtering of extracted palindromes si done. The main criteria are:

* A minimum length for the seed palindrome fixed to 8 nucleotides.
* An identity of at leaste 85% between both stems.
* A maximum spacer length fixed to 20 nucleotides.

These values can be modified by user but they have a biological relevance ([*Wang and al.*]((https://www.sciencedirect.com/science/article/pii/S0014579306000986))).

The output BED file contians the following fields:

* Chromosome name.
* Start in reference genome (0 based).
* Stop in reference genome.
* Palindrome type: Perfect or Near (less than 100% identity or Spacer).
* Palindrome length with spacer.
* Spacer sequence length.
* Mismatch rate in %.
* AT percentage in palindrome.


### Filtering step

From the paper [*long inverted repeats in eukaryotic genomes: Recombinogenic motifs determine genomic plasticity*](https://www.sciencedirect.com/science/article/pii/S0014579306000986), the following formule was used to select putative unstable palindromes based on their structural characteristics:

```math
stem\_length/spacer\_length >= mismatch\_ratio
```

A higher stem length relative to the spacer length suggests a disposition to intrastrand base pairing initially and a stable secondary structure subsequently, especially if the mismatch ratio is low.

New fields are added to the BED file and will consitute the main databases to queary for palindromic fragiles sites:

* Palindrome type: Perfect, Near and Spacer (if there is a spacer sequence).
* Recombinogenicity score which corresponds to ```stem\_length/spacer\_length```.
* Size groups to query by size group the database: 0-50 bp/51-99 bp/100-200 bp/>200 bp.

INSinPAL, as default settings, uses palindromic fragile sites within ```100-200 bp and >200 bp``` as longer palindromes are more prone to form secondary structures.

## SV callers

### Benchmarking

### INSurVeyor

### Manta

### Basil

## Annotations

### MEIs annotation

### Source of inserted sequence

### Length of inserted sequence

### Functionnal annotations


## Installation


## Execution
