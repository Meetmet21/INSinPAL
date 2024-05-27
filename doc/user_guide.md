INSinPAL User Guide
===================

## Table of content

[//]: #

* [Introduction](#Introduction)
* [Palindromic fragile sites](#Palindromic-fragile-sites)
  * [Palindrome mining algorithm](#Palindrome-mining-algorithm)
  * [Filtering step](#Filtering-step)
  * [Interesting sites](#Interesting-sites)
* [SV callers](#SV-callers)
  * [INSurVeyor](#INSurVeyor)
  * [Manta](#Manta)
  * [Basil](#Basil)
  * [Benchmarking](#Benchmarking)
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

In a genomic context, a palindrome consists of two inverted repeated DNA sequences, where the reverse complement of one sequence matches the sequence of the other. These regions are prone to forming DNA secondary structures through intrastrand base pairing. The resulting structures, such as hairpins on single-strand DNA or cruciforms on double-strand DNA, create hotspots for chromosomal rearrangements. This occurs because they can lead to double-strand breaks (DSBs) through replication fork stalling ([Mikleni and al.](https://www.mdpi.com/1422-0067/22/6/2840)). If DSBs are not repaired or are misrepaired, it can lead to deletions. translocations, inversions and large insertions. These alterations can distrupt the genetic functions of nearby genes, potentially leading to various genetic disorders and rare diseases ([J. Bissler](https://www.imrpress.com/journal/FBL/3/4/10.2741/A284)).


### Palindrome mining algorithm

### Filtering step

### Interesting sites


## SV callers

### INSurVeyor

### Manta

### Basil

### Benchmarking

## Annotations

### MEIs annotation

### Source of inserted sequence

### Length of inserted sequence

### Functionnal annotations


## Installation


## Execution
